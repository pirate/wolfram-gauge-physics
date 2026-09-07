#include <iostream>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;
void require(bool condition,const char* message) { if(!condition) throw std::runtime_error(message); }
bool same(const FiberBundleConnection& a,const FiberBundleConnection& b) {
    for(const auto [u,v] : a.base().edges()) if(a.edge_transport(u,v)!=b.edge_transport(u,v)) return false;
    return true;
}
TripleTable selected_rule(bool transport) {
    TripleTable result{std::vector<std::size_t>(512)}; std::iota(result.entries.begin(),result.entries.end(),0);
    for(const auto [a,b] : {std::pair{13,76},{37,289},{100,352},{265,328}}) std::swap(result.entries[a],result.entries[b]);
    if(transport) for(std::size_t g=1;g<8;++g) std::swap(result.entries[g],result.entries[64*g]);
    return result;
}

int main() {
    try {
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const OrientedCellComplex disk(BaseGraph({0,1,2,3,4,5},
            {{0,1},{1,2},{0,2},{2,3},{0,3},{3,4},{0,4},{0,5},{1,5}}),
            {{0,1,2,0},{0,2,3,0},{0,3,4,0},{1,0,5,1}});
        const ThreeFacePatch patch{{std::vector<Vertex>{0,1,2,0},std::vector<Vertex>{0,2,3,0},std::vector<Vertex>{0,3,4,0}}};
        const auto support=three_face_support(disk,patch);
        const auto initial=[&](std::size_t code) {
            FiberBundleConnection c(disk.base(),fiber); const auto t=decode_triple(code,8);
            c.set_transport(0,2,group.elements()[group.inverse(t[2])]);
            c.set_transport(0,3,group.elements()[group.multiply(group.inverse(t[2]),group.inverse(t[1]))]);
            c.set_transport(0,4,group.elements()[group.inverse(triple_product(group,code))]);
            c.set_transport(0,5,group.elements()[3]); c.set_transport(1,5,group.elements()[1]);
            return c;
        };
        std::size_t checks=0;
        for(const bool combined : {false,true}) {
            const auto rule=selected_rule(combined); validate_triple_table(group,rule);
            for(std::size_t code=0;code<512;++code) {
                const auto raw=initial(code);
                const auto raw_after=realize_three_face_interaction(disk,raw,patch,rule);
                for(std::size_t frame_case=0;frame_case<3;++frame_case) {
                    std::map<Vertex,Permutation> frames;
                    for(const auto v : disk.base().vertices())
                        frames.emplace(v,group.elements()[frame_case ? (code+v*(2*frame_case+1)+frame_case)%8 : 0]);
                    const auto input=raw.gauge_transform(frames);
                    MeshDynamics mesh(disk,input,rule); mesh.add_patch(patch);
                    require(mesh.patches()[0].writes.size()==2 && mesh.patches()[0].affected_faces.size()==3,
                            "three-face support omits or adds an affected face");
                    mesh.update(0);
                    const auto actual=mesh.snapshot();
                    require(same(actual,raw_after.gauge_transform(frames)),"three-face local frame covariance failed");
                    require(same(actual,realize_three_face_interaction(disk,input,patch,rule)),"prefix oracle disagrees with recurrence kernel");
                    for(const auto [u,v] : disk.base().edges()) if(!support.writes.contains({u,v}))
                        require(actual.edge_transport(u,v)==input.edge_transport(u,v),"an exterior or spectator link changed");
                    for(std::size_t f=0;f<disk.faces().size();++f)
                        require(group.elements()[mesh.face_values()[f]]==actual.holonomy(disk.faces()[f]),"three-face cache failed");
                    require(actual.holonomy(disk.faces()[3])==input.holonomy(disk.faces()[3]),"spectator face changed");
                    mesh.update(0,HurwitzDirection::Inverse); require(same(mesh.snapshot(),input),"raw link inverse failed");
                    ++checks;
                }
            }
        }
        const auto left=initial(76),right=initial(97);
        require(left.gauge_invariant_signature()!=right.gauge_invariant_signature(),"feedback witness is gauge redundancy");
        for(const auto [u,v] : disk.base().edges()) if(!support.writes.contains({u,v}))
            require(left.edge_transport(u,v)==right.edge_transport(u,v),"feedback witnesses have different exteriors");
        const auto torus=triangulated_torus(3); const auto rule=selected_rule(true);
        FiberBundleConnection input(torus.base(),fiber);
        std::size_t k=0;
        for(const auto [u,v] : torus.base().edges()) input.set_transport(u,v,group.elements()[(3*(k++)+1)%8]);
        MeshDynamics mesh(torus,input,rule); const auto specifications=three_face_patches(torus);
        require(specifications.size()==54,"torus three-face wedge count failed");
        for(const auto& p : specifications) mesh.add_patch(p);
        std::vector<std::size_t> order(specifications.size()); std::iota(order.begin(),order.end(),0);
        auto oracle=input;
        for(const auto p : order) {
            oracle=realize_three_face_interaction(torus,oracle,specifications[p],rule);
            mesh.update(p,HurwitzDirection::Forward,true); require(same(mesh.snapshot(),oracle),"overlapping fan oracle mismatch");
        }
        std::vector<std::optional<std::size_t>> last(torus.base().edges().size());
        for(std::size_t i=0;i<mesh.events().size();++i) {
            const auto& event=mesh.events()[i]; const auto& p=mesh.patches()[event.patch];
            std::set<std::size_t> parents;
            for(const auto e : p.reads) if(last[e]) parents.insert(*last[e]);
            require(event.parents==std::vector<std::size_t>(parents.begin(),parents.end()),"fan causal parents wrong");
            for(const auto e : p.writes) last[e]=i;
        }
        for(const auto& layer : mesh.commuting_layers(order)) {
            auto a=mesh,b=mesh;
            for(const auto p : layer) a.update(p);
            for(auto it=layer.rbegin();it!=layer.rend();++it) b.update(*it);
            require(a.values()==b.values(),"fan commuting layer depends on order");
        }
        for(auto it=order.rbegin();it!=order.rend();++it) mesh.update(*it,HurwitzDirection::Inverse);
        require(same(mesh.snapshot(),input),"overlapping fan inverse replay failed");
        bool rejected=false;
        try { mesh.add_patch(SharedEdgePatch{{0,1,4,0},{0,4,3,0}}); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"pair patch accepted by triple dynamics");
        rejected=false;
        try { three_face_support(disk,ThreeFacePatch{{patch.faces[0],patch.faces[1],patch.faces[0]}}); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"non-disk three-face patch accepted");
        rejected=false;
        const OrientedCellComplex nonmanifold(BaseGraph({0,1,2,3,4,5},
            {{0,1},{1,2},{0,2},{2,3},{0,3},{3,4},{0,4},{0,5},{2,5}}),
            {{0,1,2,0},{0,2,3,0},{0,3,4,0},{0,2,5,0}});
        try { three_face_support(nonmanifold,patch); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"extra face on an internal spoke accepted");
        std::cout<<checks<<" fan realization/frame/inverse cases; witness, spectator, polygon boundary and torus controls passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
