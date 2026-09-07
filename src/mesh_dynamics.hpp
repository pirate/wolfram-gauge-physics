#pragma once

#include "equivariant_pair_search.hpp"
#include "shared_edge.hpp"

namespace wgphysics::research {

// Fixed oriented two-complex, shared link degrees of freedom, and complete
// incident-face accounting. Unlike CellChain, faces do not own exclusive links.
class MeshDynamics {
public:
    enum class Lift { ExclusiveClosingLinks, SharedEdge };
    struct DirectedLink { std::size_t edge; bool reversed; };
    using Path = std::vector<DirectedLink>;
    struct Patch {
        CellPairPatch specification;
        Path first, second, connector;
        std::vector<std::size_t> reads, writes, affected_faces;
        Lift lift=Lift::ExclusiveClosingLinks;
    };
    struct Event {
        std::size_t patch, depth;
        std::vector<std::size_t> parents;
    };

    MeshDynamics(const OrientedCellComplex& complex, const FiberBundleConnection& connection,
                 const PairTable& table)
        : complex_(complex), fiber_(connection.fiber()), group_(connection.local_gauge_group()),
          forward_(table) {
        if (complex.base().vertices()!=connection.base().vertices() ||
            complex.base().edges()!=connection.base().edges())
            throw std::invalid_argument("mesh and connection disagree");
        validate_pair_table(group_,table);
        inverse_.resize(table.size());
        for (std::size_t i=0;i<table.size();++i) inverse_[table[i]]=i;
        identity_=group_.index_of(Permutation::identity(fiber_.vertex_count()));
        for (const auto edge : complex.base().edges()) {
            edge_ids_[edge]=edges_.size(); edges_.push_back(edge);
            values_.push_back(group_.index_of(connection.edge_transport(edge.first,edge.second)));
        }
        incidence_.resize(edges_.size()); last_writer_.resize(edges_.size());
        for (const auto& face : complex.faces()) {
            faces_.push_back(compile_path(face));
            const auto id=faces_.size()-1;
            for (const auto link : faces_.back()) incidence_[link.edge].push_back(id);
            face_values_.push_back(transport(faces_.back()));
        }
        sectors_.resize(group_.order());
        for (uint16_t a=0;a<group_.order();++a) {
            sectors_[a]=a;
            for (uint16_t g=0;g<group_.order();++g) sectors_[a]=std::min(sectors_[a],conjugate(g,a));
        }
    }

    std::size_t add_patch(const CellPairPatch& specification) {
        return append_patch(specification,interaction_support(complex_.base(),specification),Lift::ExclusiveClosingLinks);
    }

    std::size_t add_patch(const SharedEdgePatch& specification) {
        const auto support=shared_edge_support(complex_,specification);
        return append_patch({specification.first_loop,specification.second_loop,{specification.first_loop.front()}},
                            support,Lift::SharedEdge);
    }

private:
    std::size_t append_patch(const CellPairPatch& specification,const InteractionSupport& support,Lift lift) {
        // Validate actual faces, not arbitrarily chosen cycles of the one-skeleton.
        for (const auto* loop : {&specification.first_loop,&specification.second_loop}) {
            const auto canonical=canonical_oriented_face(*loop);
            if (std::find(complex_.faces().begin(),complex_.faces().end(),canonical)==complex_.faces().end())
                throw std::invalid_argument("update loop is not an oriented mesh face");
        }
        Patch patch{specification,compile_path(specification.first_loop),
                    compile_path(specification.second_loop),compile_path(specification.connector),{},{},{},lift};
        for (const auto edge : support.reads) patch.reads.push_back(edge_ids_.at(edge));
        std::set<std::size_t> affected;
        for (const auto edge : support.writes) {
            const auto id=edge_ids_.at(edge); patch.writes.push_back(id);
            affected.insert(incidence_[id].begin(),incidence_[id].end());
        }
        patch.affected_faces.assign(affected.begin(),affected.end());
        patches_.push_back(std::move(patch));
        return patches_.size()-1;
    }

public:
    void update(std::size_t id, HurwitzDirection direction=HurwitzDirection::Forward, bool record=false) {
        const auto& p=patches_.at(id);
        if (record && unrecorded_) throw std::logic_error("cannot append provenance after unrecorded events");
        if (record) {
            std::set<std::size_t> parents;
            for (const auto edge : p.reads) if (last_writer_[edge]) parents.insert(*last_writer_[edge]);
            Event event{id,0,{parents.begin(),parents.end()}};
            for (const auto parent : parents) event.depth=std::max(event.depth,events_[parent].depth+1);
            for (const auto edge : p.writes) last_writer_[edge]=events_.size();
            events_.push_back(std::move(event));
        } else unrecorded_=true;
        const auto a=transport(p.first),b=transport(p.second),t=transport(p.connector);
        const auto based_b=conjugate(t,b);
        const auto& table=direction==HurwitzDirection::Forward ? forward_ : inverse_;
        const auto target=table[a*group_.order()+based_b];
        const uint16_t next_a=target/group_.order(), next_b=conjugate(group_.inverse(t),target%group_.order());
        if(p.lift==Lift::SharedEdge) {
            const auto shared=p.first.front();
            // S' = S A^-1 X; independent oracle instead computes Q^-1 X.
            set_link(shared,group_.multiply(link_value(shared),group_.multiply(group_.inverse(a),next_a)));
        } else {
            const auto first=p.first.back(),second=p.second.back();
            const auto u=group_.multiply(group_.multiply(next_a,group_.inverse(a)),link_value(first));
            const auto v=group_.multiply(group_.multiply(next_b,group_.inverse(b)),link_value(second));
            set_link(first,u); set_link(second,v);
        }
        for (const auto face : p.affected_faces) face_values_[face]=transport(faces_[face]);
    }

    bool independent(std::size_t a,std::size_t b) const {
        const auto& p=patches_.at(a); const auto& q=patches_.at(b);
        for (const auto edge : p.writes)
            if (std::binary_search(q.reads.begin(),q.reads.end(),edge)) return false;
        for (const auto edge : q.writes)
            if (std::binary_search(p.reads.begin(),p.reads.end(),edge)) return false;
        return true;
    }

    // Color actual link read/write conflicts, not just face adjacency. These are
    // commuting physical link updates. Concurrent implementations must separately
    // refresh the union of affected face caches after all link writes complete.
    std::vector<std::vector<std::size_t>> commuting_layers(const std::vector<std::size_t>& order) const {
        std::vector<std::vector<std::size_t>> layers;
        std::vector<std::set<std::size_t>> reads,writes;
        std::set<std::size_t> seen;
        for (const auto id : order) {
            const auto& p=patches_.at(id);
            if (!seen.insert(id).second) throw std::invalid_argument("duplicate scheduled patch");
            std::size_t layer=0;
            for (;layer<layers.size();++layer) {
                bool conflict=false;
                for (const auto e : p.writes) if (reads[layer].contains(e)) conflict=true;
                for (const auto e : p.reads) if (writes[layer].contains(e)) conflict=true;
                if (!conflict) break;
            }
            if (layer==layers.size()) { layers.emplace_back(); reads.emplace_back(); writes.emplace_back(); }
            layers[layer].push_back(id);
            reads[layer].insert(p.reads.begin(),p.reads.end()); writes[layer].insert(p.writes.begin(),p.writes.end());
        }
        return layers;
    }

    FiberBundleConnection snapshot() const {
        FiberBundleConnection result(complex_.base(),fiber_);
        for (std::size_t i=0;i<edges_.size();++i)
            result.set_transport(edges_[i].first,edges_[i].second,group_.elements()[values_[i]]);
        return result;
    }
    const auto& complex() const { return complex_; }
    const auto& group() const { return group_; }
    const auto& values() const { return values_; }
    const auto& face_values() const { return face_values_; }
    const auto& patches() const { return patches_; }
    const auto& events() const { return events_; }
    const auto& incidence() const { return incidence_; }
    std::pair<uint16_t,uint16_t> based_pair(std::size_t id) const {
        const auto& p=patches_.at(id);
        return {transport(p.first),conjugate(transport(p.connector),transport(p.second))};
    }
    uint16_t sector(uint16_t value) const { return sectors_.at(value); }
    uint16_t identity() const { return identity_; }

private:
    Path compile_path(const std::vector<Vertex>& vertices) const {
        if (vertices.empty()) throw std::invalid_argument("empty mesh path");
        if (!complex_.base().has_vertex(vertices.front())) throw std::out_of_range("path vertex absent");
        Path result;
        for (std::size_t i=1;i<vertices.size();++i) {
            const auto u=vertices[i-1],v=vertices[i];
            result.push_back({edge_ids_.at(infragauge::canonical_edge(u,v)),u>v});
        }
        return result;
    }
    uint16_t link_value(DirectedLink link) const {
        return link.reversed ? group_.inverse(values_[link.edge]) : values_[link.edge];
    }
    void set_link(DirectedLink link,uint16_t value) {
        values_[link.edge]=link.reversed ? group_.inverse(value) : value;
    }
    uint16_t transport(const Path& path) const {
        auto value=identity_;
        for (const auto link : path) value=group_.multiply(link_value(link),value);
        return value;
    }
    uint16_t conjugate(uint16_t g,uint16_t a) const {
        return group_.multiply(g,group_.multiply(a,group_.inverse(g)));
    }
    OrientedCellComplex complex_;
    FiberGraph fiber_;
    AutomorphismTables group_;
    PairTable forward_,inverse_;
    uint16_t identity_{};
    std::vector<Edge> edges_;
    std::map<Edge,std::size_t> edge_ids_;
    std::vector<uint16_t> values_,face_values_,sectors_;
    std::vector<Path> faces_;
    std::vector<std::vector<std::size_t>> incidence_;
    std::vector<Patch> patches_;
    std::vector<std::optional<std::size_t>> last_writer_;
    std::vector<Event> events_;
    bool unrecorded_{};
};

inline OrientedCellComplex triangulated_torus(std::size_t side) {
    if (side<3 || side>128) throw std::invalid_argument("torus side must lie in [3,128]");
    std::vector<Vertex> vertices(side*side); std::iota(vertices.begin(),vertices.end(),0);
    std::set<Edge> edges; std::vector<std::vector<Vertex>> faces;
    const auto vertex=[&](std::size_t x,std::size_t y) { return static_cast<Vertex>((y%side)*side+x%side); };
    for (std::size_t y=0;y<side;++y) for (std::size_t x=0;x<side;++x) {
        const auto a=vertex(x,y),b=vertex(x+1,y),c=vertex(x+1,y+1),d=vertex(x,y+1);
        for (auto face : {std::vector<Vertex>{a,b,c,a},std::vector<Vertex>{a,c,d,a}}) {
            for (std::size_t i=1;i<face.size();++i) edges.insert(infragauge::canonical_edge(face[i-1],face[i]));
            faces.push_back(std::move(face));
        }
    }
    return OrientedCellComplex(BaseGraph(vertices,{edges.begin(),edges.end()}),faces);
}

// Every ordered adjacent triangle pair, every exclusive closing-link choice,
// and its unique connector in the union with both written edges removed.
inline std::vector<CellPairPatch> triangular_mesh_patches(const OrientedCellComplex& complex) {
    std::map<Edge,std::vector<std::size_t>> incidence;
    for (std::size_t f=0;f<complex.faces().size();++f) {
        const auto& face=complex.faces()[f];
        if (face.size()!=4) throw std::invalid_argument("patch enumeration currently requires triangles");
        for (std::size_t k=0;k<3;++k) incidence[infragauge::canonical_edge(face[k],face[k+1])].push_back(f);
    }
    std::vector<CellPairPatch> result;
    for (const auto& [shared,faces] : incidence) {
        if (faces.size()>2) throw std::invalid_argument("nonmanifold triangle edge");
        if (faces.size()!=2) continue;
        for (const auto first : faces) for (const auto second : faces) if (first!=second) {
            const auto& a=complex.faces()[first]; const auto& b=complex.faces()[second];
            for (std::size_t i=0;i<3;++i) for (std::size_t j=0;j<3;++j) {
                std::vector<Vertex> p{a[i],a[(i+1)%3],a[(i+2)%3],a[i]};
                std::vector<Vertex> q{b[j],b[(j+1)%3],b[(j+2)%3],b[j]};
                const auto u=infragauge::canonical_edge(p[2],p[3]),v=infragauge::canonical_edge(q[2],q[3]);
                if (u==shared || v==shared) continue;
                std::map<Vertex,std::set<Vertex>> neighbors;
                for (const auto* loop : {&p,&q}) for (std::size_t k=0;k<3;++k) {
                    const auto edge=infragauge::canonical_edge((*loop)[k],(*loop)[k+1]);
                    if (edge==u || edge==v) continue;
                    neighbors[edge.first].insert(edge.second); neighbors[edge.second].insert(edge.first);
                }
                std::map<Vertex,Vertex> parents{{q[0],q[0]}};
                std::vector<Vertex> queue{q[0]};
                for (std::size_t k=0;k<queue.size();++k) for (const auto next : neighbors[queue[k]])
                    if (!parents.contains(next)) { parents[next]=queue[k]; queue.push_back(next); }
                if (!parents.contains(p[0])) continue;
                std::vector<Vertex> path{p[0]};
                while (path.back()!=q[0]) path.push_back(parents.at(path.back()));
                std::reverse(path.begin(),path.end());
                result.push_back({std::move(p),std::move(q),std::move(path)});
            }
        }
    }
    return result;
}

} // namespace wgphysics::research
