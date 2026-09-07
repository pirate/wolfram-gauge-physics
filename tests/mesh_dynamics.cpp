#include <iostream>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;
void require(bool ok,const char* text) { if (!ok) throw std::runtime_error(text); }

bool same(const FiberBundleConnection& a,const FiberBundleConnection& b) {
    if (a.base().edges()!=b.base().edges()) return false;
    for (const auto [u,v] : a.base().edges()) if (a.edge_transport(u,v)!=b.edge_transport(u,v)) return false;
    return true;
}

int main() {
    try {
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        PairTable rule;
        for (const auto& a : group.elements()) for (const auto& b : group.elements()) {
            const auto [c,d]=pair_targets(a,b,PairInteraction::BoundaryShear);
            rule.push_back(group.index_of(c)*group.order()+group.index_of(d));
        }
        const auto complex=triangulated_torus(3);
        require(complex.base().edges().size()==27 && complex.faces().size()==18,"torus incidence counts");
        std::map<Edge,int> oriented_incidence;
        for(const auto& face : complex.faces()) for(std::size_t k=1;k<face.size();++k)
            oriented_incidence[wgphysics::infragauge::canonical_edge(face[k-1],face[k])]+=face[k-1]<face[k] ? 1 : -1;
        for(const auto& [edge,sum] : oriented_incidence) { (void)edge; require(sum==0,"torus orientations disagree"); }
        FiberBundleConnection initial(complex.base(),fiber);
        std::size_t index=0;
        for (const auto [u,v] : complex.base().edges()) initial.set_transport(u,v,group.elements()[(index++*3+1)%8]);
        MeshDynamics mesh(complex,initial,rule);
        const auto specifications=triangular_mesh_patches(complex);
        require(specifications.size()==216,"complete ordered patch count");
        for (const auto& p : specifications) mesh.add_patch(p);
        for (const auto& faces : mesh.incidence()) require(faces.size()==2,"closed torus edge incidence");
        auto oracle=initial;
        for (std::size_t i=0;i<specifications.size();++i) {
            oracle=apply_pair_table_interaction(oracle,specifications[i],rule);
            mesh.update(i,HurwitzDirection::Forward,true);
            require(same(mesh.snapshot(),oracle),"shared-link kernel differs from independent link oracle");
            for (std::size_t f=0;f<complex.faces().size();++f)
                require(group.elements()[mesh.face_values()[f]]==oracle.holonomy(complex.faces()[f]),
                        "incident-face cache missed backreaction");
        }
        // Every read's latest writer is a parent, including spectator/connector links.
        std::vector<std::optional<std::size_t>> last(complex.base().edges().size());
        for (std::size_t i=0;i<mesh.events().size();++i) {
            const auto& event=mesh.events()[i]; const auto& p=mesh.patches()[event.patch];
            std::set<std::size_t> parents;
            for (const auto edge : p.reads) if (last[edge]) parents.insert(*last[edge]);
            require(event.parents==std::vector<std::size_t>(parents.begin(),parents.end()),"missing causal read dependency");
            for (const auto edge : p.writes) last[edge]=i;
        }
        for (std::size_t i=specifications.size();i-->0;) mesh.update(i,HurwitzDirection::Inverse);
        require(same(mesh.snapshot(),initial),"whole-mesh inverse replay failed");
        bool rejected=false;
        try { mesh.update(0,HurwitzDirection::Forward,true); } catch (const std::logic_error&) { rejected=true; }
        require(rejected,"provenance appended after unrecorded history");
        // Independent local frames at every vertex, on every generated patch.
        std::map<Vertex,Permutation> frames;
        for (const auto vertex : complex.base().vertices()) frames.emplace(vertex,group.elements()[(5*vertex+3)%8]);
        MeshDynamics plain(complex,initial,rule),gauged(complex,initial.gauge_transform(frames),rule);
        for (const auto& p : specifications) { plain.add_patch(p); gauged.add_patch(p); }
        for (std::size_t i=0;i<specifications.size();++i) {
            plain.update(i); gauged.update(i);
            require(same(plain.snapshot().gauge_transform(frames),gauged.snapshot()),"local gauge covariance failed");
        }
        std::vector<std::size_t> order(specifications.size()); std::iota(order.begin(),order.end(),0);
        const auto layers=plain.commuting_layers(order);
        std::size_t scheduled=0;
        for (const auto& layer : layers) {
            for (const auto a : layer) for (const auto b : layer) if (a!=b)
                require(plain.independent(a,b),"conflicting patches in claimed parallel layer");
            MeshDynamics forward=plain,backward=plain;
            for (const auto p : layer) forward.update(p);
            for (auto it=layer.rbegin();it!=layer.rend();++it) backward.update(*it);
            require(forward.values()==backward.values(),"parallel layer changes under reordering");
            scheduled+=layer.size();
        }
        require(scheduled==specifications.size(),"scheduler dropped patches");
        MeshDynamics flat(complex,FiberBundleConnection(complex.base(),fiber).gauge_transform(frames),rule);
        for (const auto& p : specifications) flat.add_patch(p);
        const auto flat_before=flat.values();
        for (const auto id : order) flat.update(id);
        require(flat.values()==flat_before,"pure-gauge flat vacuum changed");
        rejected=false;
        try { mesh.add_patch({{0,1,4,3,0},{0,4,3,0},{0}}); } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"non-face loop accepted");
        std::cout<<specifications.size()<<" shared-face oracle updates and gauge checks; "<<layers.size()
                 <<" conflict-free layers; cache, inverse, flatness and causal dependencies passed\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
