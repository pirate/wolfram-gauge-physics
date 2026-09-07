#pragma once

#include "triple_rules.hpp"

namespace wgphysics::research {

// Oriented disk fan: (u,v0,v1,u), (u,v1,v2,u), (u,v2,v3,u).
// The tuple supplied to the rule is (H3,H2,H1), whose product is the boundary.
struct ThreeFacePatch { std::array<std::vector<Vertex>,3> faces; };

inline InteractionSupport three_face_support(const OrientedCellComplex& complex,const ThreeFacePatch& patch) {
    const auto u=patch.faces[0].empty() ? Vertex{} : patch.faces[0][0];
    std::vector<Vertex> rim;
    InteractionSupport support;
    for(std::size_t i=0;i<3;++i) {
        const auto& f=patch.faces[i];
        if(f.size()!=4 || f[0]!=u || f[3]!=u || (i && rim.back()!=f[1]))
            throw std::invalid_argument("three-face patch is not an oriented triangular fan");
        if(std::find(complex.faces().begin(),complex.faces().end(),canonical_oriented_face(f))==complex.faces().end())
            throw std::invalid_argument("three-face loop is not an actual oriented cell");
        if(!i) rim.push_back(f[1]);
        rim.push_back(f[2]);
        for(std::size_t k=1;k<f.size();++k) support.reads.insert(infragauge::canonical_edge(f[k-1],f[k]));
    }
    const std::set<Vertex> distinct(rim.begin(),rim.end());
    if(distinct.size()!=4 || distinct.contains(u)) throw std::invalid_argument("three-face fan is not a five-vertex disk");
    for(std::size_t k=1;k<3;++k) {
        const auto e=infragauge::canonical_edge(u,rim[k]);
        std::size_t incident=0;
        for(const auto& face : complex.faces()) if(loop_contains_edge(face,e)) ++incident;
        if(incident!=2) throw std::invalid_argument("three-face write has an unaccounted incident face");
        support.writes.insert(e);
    }
    return support;
}

inline std::vector<ThreeFacePatch> three_face_patches(const OrientedCellComplex& complex) {
    std::map<Vertex,std::map<Vertex,Vertex>> links;
    for(const auto& face : complex.faces()) {
        if(face.size()!=4) throw std::invalid_argument("three-face enumeration requires triangles");
        for(std::size_t i=0;i<3;++i) {
            auto& next=links[face[i]];
            if(!next.emplace(face[(i+1)%3],face[(i+2)%3]).second)
                throw std::invalid_argument("inconsistent oriented fan incidence");
        }
    }
    std::vector<ThreeFacePatch> result;
    for(const auto& [u,next] : links) for(const auto& [start,second] : next) {
        (void)second;
        std::vector<Vertex> rim{start};
        while(rim.size()<4 && next.contains(rim.back())) rim.push_back(next.at(rim.back()));
        if(rim.size()!=4 || std::set<Vertex>(rim.begin(),rim.end()).size()!=4) continue;
        ThreeFacePatch p{{std::vector<Vertex>{u,rim[0],rim[1],u},
                          std::vector<Vertex>{u,rim[1],rim[2],u},
                          std::vector<Vertex>{u,rim[2],rim[3],u}}};
        (void)three_face_support(complex,p); result.push_back(std::move(p));
    }
    return result;
}

inline FiberBundleConnection realize_three_face_interaction(
    const OrientedCellComplex& complex,const FiberBundleConnection& input,const ThreeFacePatch& patch,const TripleTable& table) {
    (void)three_face_support(complex,patch);
    if(complex.base().edges()!=input.base().edges() || complex.base().vertices()!=input.base().vertices())
        throw std::invalid_argument("three-face mesh and connection disagree");
    const AutomorphismTables group(input.local_gauge_group());
    validate_triple_table(group,table);
    Triple state;
    for(std::size_t i=0;i<3;++i) state[2-i]=group.index_of(input.holonomy(patch.faces[i]));
    const auto targets=decode_triple(table.entries[encode_triple(state,group.order())],group.order());
    const auto mul=Permutation::compose;
    auto output=input;
    // Solve from each unchanged outer-prefix path and the requested cumulative
    // holonomy, rather than sharing the compiled runner's spoke recurrence.
    std::vector<Vertex> exterior{patch.faces[0][0],patch.faces[0][1]};
    auto cumulative=Permutation::identity(input.fiber().vertex_count());
    for(std::size_t i=0;i<2;++i) {
        exterior.push_back(patch.faces[i][2]);
        cumulative=mul(group.elements()[targets[2-i]],cumulative);
        output.set_transport(exterior.front(),exterior.back(),mul(input.parallel_transport(exterior),cumulative.inverse()));
    }
    for(std::size_t i=0;i<3;++i)
        if(output.holonomy(patch.faces[i])!=group.elements()[targets[2-i]])
            throw std::logic_error("three-face realization missed a holonomy target");
    return output;
}

} // namespace wgphysics::research
