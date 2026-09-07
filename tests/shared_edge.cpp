#include <iostream>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
bool same(const FiberBundleConnection& a,const FiberBundleConnection& b) {
    for(const auto [u,v] : a.base().edges()) if(a.edge_transport(u,v)!=b.edge_transport(u,v)) return false;
    return true;
}

int main() {
    try {
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto n=group.order();
        PairTable table;
        for(const auto& a : group.elements()) for(const auto& b : group.elements()) {
            const auto [x,y]=pair_targets(a,b,PairInteraction::BoundaryShear);
            table.push_back(group.index_of(x)*n+group.index_of(y));
        }
        const OrientedCellComplex disk(BaseGraph({0,1,2,3},{{0,1},{0,2},{0,3},{1,2},{1,3}}),
                                       {{0,1,2,0},{0,3,1,0}});
        const auto patches=shared_edge_patches(disk);
        require(patches.size()==2,"wrong rooted disk patch count");
        const std::vector<Vertex> boundary{0,3,1,2,0};
        // Every assignment to all five links: includes all local gauge frames,
        // not merely gauge-fixed holonomies. Test a non-involutive word rule.
        for(std::size_t code=0;code<32768;++code) {
            FiberBundleConnection input(disk.base(),fiber);
            auto remaining=code;
            for(const auto [u,v] : disk.base().edges()) {
                input.set_transport(u,v,group.elements()[remaining%n]); remaining/=n;
            }
            for(const auto& patch : patches) {
                MeshDynamics mesh(disk,input,table); mesh.add_patch(patch);
                require(mesh.patches()[0].writes.size()==1 && mesh.patches()[0].affected_faces.size()==2,
                        "shared-edge primitive has hidden spectators");
                const auto before=mesh.values();
                const auto oracle=realize_shared_edge_interaction(disk,input,patch,[&](const auto& a,const auto& b) {
                    return pair_targets(a,b,PairInteraction::BoundaryShear);
                });
                mesh.update(0);
                const auto after=mesh.snapshot();
                require(same(after,oracle),"compiled shared-edge kernel differs from exterior-path oracle");
                require(after.holonomy(boundary)==input.holonomy(boundary),"boundary holonomy changed");
                for(const auto [u,v] : disk.base().edges()) if(wgphysics::infragauge::canonical_edge(u,v)!=
                    wgphysics::infragauge::canonical_edge(patch.first_loop[0],patch.first_loop[1]))
                    require(after.edge_transport(u,v)==input.edge_transport(u,v),"exterior link changed");
                for(std::size_t f=0;f<disk.faces().size();++f)
                    require(group.elements()[mesh.face_values()[f]]==after.holonomy(disk.faces()[f]),"face cache wrong");
                mesh.update(0,HurwitzDirection::Inverse);
                require(mesh.values()==before,"shared-edge inverse failed on raw links");
            }
        }
        const auto torus=triangulated_torus(3);
        FiberBundleConnection input(torus.base(),fiber);
        std::size_t cursor=0;
        for(const auto [u,v] : torus.base().edges()) input.set_transport(u,v,group.elements()[(3*(cursor++)+1)%n]);
        std::map<Vertex,Permutation> frames;
        for(const auto v : torus.base().vertices()) frames.emplace(v,group.elements()[(5*v+3)%n]);
        MeshDynamics plain(torus,input,table),gauged(torus,input.gauge_transform(frames),table);
        const auto torus_patches=shared_edge_patches(torus);
        require(torus_patches.size()==54,"closed mesh must expose two rooted patches per edge");
        for(const auto& p : torus_patches) { plain.add_patch(p); gauged.add_patch(p); }
        std::vector<std::size_t> order(torus_patches.size()); std::iota(order.begin(),order.end(),0);
        for(const auto p : order) {
            plain.update(p,HurwitzDirection::Forward,true); gauged.update(p);
            require(same(plain.snapshot().gauge_transform(frames),gauged.snapshot()),"local frame covariance failed");
        }
        std::vector<std::optional<std::size_t>> last(torus.base().edges().size());
        for(std::size_t i=0;i<plain.events().size();++i) {
            const auto& event=plain.events()[i]; const auto& patch=plain.patches()[event.patch];
            std::set<std::size_t> parents;
            for(const auto e : patch.reads) if(last[e]) parents.insert(*last[e]);
            require(event.parents==std::vector<std::size_t>(parents.begin(),parents.end()),"missing causal read dependency");
            for(const auto e : patch.writes) last[e]=i;
        }
        for(const auto& layer : plain.commuting_layers(order)) {
            auto a=plain,b=plain;
            for(const auto p : layer) a.update(p);
            for(auto it=layer.rbegin();it!=layer.rend();++it) b.update(*it);
            require(a.values()==b.values(),"claimed commuting layer depends on ordering");
        }
        bool rejected=false;
        const OrientedCellComplex nonmanifold(BaseGraph({0,1,2,3,4},
            {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4}}),{{0,1,2,0},{0,3,1,0},{0,1,4,0}});
        try { (void)shared_edge_support(nonmanifold,patches[0]); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"third incident face silently ignored");
        rejected=false;
        try { (void)shared_edge_patches(OrientedCellComplex(disk.base(),{{0,1,2,0},{0,1,3,0}})); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"incompatible face orientations accepted");
        const OrientedCellComplex polygons(BaseGraph({0,1,2,3,4,5},
            {{0,1},{1,2},{2,3},{0,3},{0,4},{4,5},{1,5}}),{{0,1,2,3,0},{0,4,5,1,0}});
        FiberBundleConnection polygon_input(polygons.base(),fiber);
        cursor=0;
        for(const auto [u,v] : polygons.base().edges()) polygon_input.set_transport(u,v,group.elements()[(cursor++)%n]);
        for(const auto& patch : shared_edge_patches(polygons)) {
            MeshDynamics mesh(polygons,polygon_input,table); mesh.add_patch(patch); mesh.update(0);
            const auto oracle=realize_shared_edge_interaction(polygons,polygon_input,patch,[&](const auto& a,const auto& b) {
                return pair_targets(a,b,PairInteraction::BoundaryShear);
            });
            require(same(mesh.snapshot(),oracle),"polygonal shared-edge realization failed");
        }
        std::cout<<"65536 raw-link disk updates checked; torus frames, causality, commuting layers and incidence guards passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
