#pragma once

#include "gauge_interactions.hpp"

namespace wgphysics::research {

// A different link realization from CellPairPatch: the first loop STARTS with
// u->v, the second ENDS with v->u, and both are based at u. Their ordered product
// is the transport around the exterior boundary, with the common edge canceled.
struct SharedEdgePatch {
    std::vector<Vertex> first_loop, second_loop;
};

inline InteractionSupport shared_edge_support(const OrientedCellComplex& complex,
                                              const SharedEdgePatch& patch) {
    const auto& a=patch.first_loop; const auto& b=patch.second_loop;
    for(const auto* loop : {&a,&b}) {
        const auto canonical=canonical_oriented_face(*loop);
        if(std::find(complex.faces().begin(),complex.faces().end(),canonical)==complex.faces().end())
            throw std::invalid_argument("shared-edge loop is not an oriented mesh face");
    }
    if(a.size()<4 || b.size()<4 || a.front()!=b.front() || a[1]!=b[b.size()-2])
        throw std::invalid_argument("shared-edge loops need opposite traversals at a common root");
    const auto shared=infragauge::canonical_edge(a[0],a[1]);
    std::set<Vertex> common;
    for(const auto v : a) if(std::find(b.begin(),b.end(),v)!=b.end()) common.insert(v);
    if(common!=std::set<Vertex>{a[0],a[1]})
        throw std::invalid_argument("shared-edge patch must be a two-face disk");
    std::size_t incident=0;
    for(const auto& face : complex.faces()) if(loop_contains_edge(face,shared)) ++incident;
    if(incident!=2) throw std::invalid_argument("shared-edge write requires exactly two incident faces");
    InteractionSupport result;
    for(const auto* loop : {&a,&b}) for(std::size_t i=1;i<loop->size();++i)
        result.reads.insert(infragauge::canonical_edge((*loop)[i-1],(*loop)[i]));
    result.writes.insert(shared);
    return result;
}

template<class PairMap>
inline FiberBundleConnection realize_shared_edge_interaction(
    const OrientedCellComplex& complex,const FiberBundleConnection& connection,
    const SharedEdgePatch& patch,PairMap targets) {
    (void)shared_edge_support(complex,patch);
    if(complex.base().edges()!=connection.base().edges() ||
       complex.base().vertices()!=connection.base().vertices())
        throw std::invalid_argument("shared-edge mesh and connection disagree");
    const auto mul=Permutation::compose;
    const auto a=connection.holonomy(patch.first_loop),b=connection.holonomy(patch.second_loop);
    const auto [x,y]=targets(a,b);
    if(mul(x,y)!=mul(a,b)) throw std::invalid_argument("shared-edge target changes boundary transport");
    // A = Q S and B = S^-1 P. Q and P are exterior paths, hence unchanged.
    const std::vector<Vertex> exterior(patch.first_loop.begin()+1,patch.first_loop.end());
    const auto q=connection.parallel_transport(exterior);
    auto result=connection;
    result.set_transport(patch.first_loop[0],patch.first_loop[1],mul(q.inverse(),x));
    if(result.holonomy(patch.first_loop)!=x || result.holonomy(patch.second_loop)!=y)
        throw std::logic_error("shared-edge realization failed its holonomy targets");
    return result;
}

// Both oriented rooted versions of each interior edge. They are distinct
// candidate operators, not asserted equivalent by a basepoint or braid argument.
inline std::vector<SharedEdgePatch> shared_edge_patches(const OrientedCellComplex& complex) {
    std::map<Edge,std::vector<std::size_t>> incidence;
    for(std::size_t f=0;f<complex.faces().size();++f) {
        const auto& face=complex.faces()[f];
        for(std::size_t i=1;i<face.size();++i)
            incidence[infragauge::canonical_edge(face[i-1],face[i])].push_back(f);
    }
    const auto rotate=[](const std::vector<Vertex>& face,Vertex root) {
        const auto begin=std::find(face.begin(),face.end()-1,root)-face.begin();
        const auto length=face.size()-1;
        std::vector<Vertex> loop;
        for(std::size_t k=0;k<=length;++k) loop.push_back(face[(begin+k)%length]);
        return loop;
    };
    std::vector<SharedEdgePatch> result;
    for(const auto& [edge,faces] : incidence) {
        if(faces.size()>2) throw std::invalid_argument("nonmanifold shared-edge incidence");
        if(faces.size()!=2) continue;
        for(std::size_t order=0;order<2;++order) {
            const auto& first=complex.faces()[faces[order]];
            for(std::size_t i=1;i<first.size();++i)
                if(infragauge::canonical_edge(first[i-1],first[i])==edge) {
                    SharedEdgePatch patch{rotate(first,first[i-1]),rotate(complex.faces()[faces[1-order]],first[i-1])};
                    (void)shared_edge_support(complex,patch);
                    result.push_back(std::move(patch));
                }
        }
    }
    return result;
}

} // namespace wgphysics::research
