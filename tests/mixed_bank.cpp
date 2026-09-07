#include <iostream>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;
void require(bool b,const char* message) { if(!b) throw std::runtime_error(message); }
bool same(const FiberBundleConnection& a,const FiberBundleConnection& b) {
    for(const auto [u,v] : a.base().edges()) if(a.edge_transport(u,v)!=b.edge_transport(u,v)) return false;
    return true;
}
int main() {
    try {
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto census=search_triple_involutions(group);
        const auto complex=triangulated_torus(3);
        FiberBundleConnection input(complex.base(),fiber);
        std::size_t k=0;
        for(const auto [u,v] : complex.base().edges()) input.set_transport(u,v,group.elements()[(k++*3+1)%8]);
        PairTable transport(64); std::iota(transport.begin(),transport.end(),0);
        for(std::size_t a=1;a<8;++a) std::swap(transport[a],transport[8*a]);
        MeshDynamics mesh(complex,input,transport);
        const auto r1=mesh.register_rule(census.at(48)),r2=mesh.register_rule(census.at(55));
        require(r1!=r2 && mesh.register_rule(census.at(48))==r1,"rule registry identity failed");
        require(mesh.register_rule(transport)==0,"pair rule deduplication failed");
        const auto fans=three_face_patches(complex); const auto pairs=shared_edge_patches(complex);
        for(const auto& p : fans) mesh.add_patch(p,r1);
        for(const auto& p : pairs) mesh.add_patch(p,0);
        auto default_dispatch=mesh,explicit_dispatch=mesh;
        default_dispatch.update(0); explicit_dispatch.update_with_rule(0,r1);
        require(default_dispatch.values()==explicit_dispatch.values(),"nonzero default rule dispatch failed");
        const auto original=mesh.values();
        const auto valid_events=mesh.events().size();
        bool rejected=false;
        try { mesh.update_with_rule(0,0,HurwitzDirection::Forward,true); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected && mesh.values()==original && mesh.events().size()==valid_events,"arity error mutated state or provenance");
        rejected=false;
        try { mesh.add_patch(pairs[0],r2); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"pair support accepted a triple default rule");
        rejected=false;
        try { mesh.update_with_rule(0,999); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"unknown rule accepted");
        auto bad=census[48]; bad.entries[0]=1;
        rejected=false;
        try { mesh.register_rule(bad); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected && mesh.rules().size()==3,"bad rule changed the registry");
        std::map<Vertex,Permutation> frames;
        for(const auto v : complex.base().vertices()) frames.emplace(v,group.elements()[(5*v+3)%8]);
        MeshDynamics gauged(complex,input.gauge_transform(frames),transport);
        gauged.register_rule(census[48]); gauged.register_rule(census[55]);
        for(const auto& p : fans) gauged.add_patch(p,r1);
        for(const auto& p : pairs) gauged.add_patch(p,0);
        auto oracle=input;
        std::vector<std::pair<std::size_t,std::size_t>> history;
        for(std::size_t step=0;step<324;++step) {
            const auto index=(step*17+step/3)%fans.size();
            const auto rule=step%3==0 ? 0 : step%3==1 ? r1 : r2;
            const auto p=rule ? index : fans.size()+index;
            if(rule) oracle=realize_three_face_interaction(complex,oracle,fans[index],census[rule==r1?48:55]);
            else oracle=realize_shared_edge_interaction(complex,oracle,pairs[index],[&](auto a,auto b) {
                const auto target=transport[group.index_of(a)*8+group.index_of(b)];
                return std::pair{group.elements()[target/8],group.elements()[target%8]};
            });
            mesh.update_with_rule(p,rule,HurwitzDirection::Forward,true);
            gauged.update_with_rule(p,rule);
            require(same(mesh.snapshot(),oracle),"mixed bank diverges from independent connection lifts");
            require(same(mesh.snapshot().gauge_transform(frames),gauged.snapshot()),"mixed-bank local frame covariance failed");
            require(mesh.events().back().rule==rule,"event omitted the selected rule");
            for(std::size_t f=0;f<complex.faces().size();++f)
                require(group.elements()[mesh.face_values()[f]]==oracle.holonomy(complex.faces()[f]),"mixed face cache failed");
            history.emplace_back(p,rule);
        }
        std::vector<std::optional<std::size_t>> last(complex.base().edges().size());
        for(std::size_t i=0;i<mesh.events().size();++i) {
            const auto& event=mesh.events()[i]; const auto& p=mesh.patches()[event.patch];
            std::set<std::size_t> parents;
            for(const auto edge : p.reads) if(last[edge]) parents.insert(*last[edge]);
            require(event.parents==std::vector<std::size_t>(parents.begin(),parents.end()),"mixed-rule causal parent mismatch");
            for(const auto edge : p.writes) last[edge]=i;
        }
        for(auto it=history.rbegin();it!=history.rend();++it) mesh.update_with_rule(it->first,it->second,HurwitzDirection::Inverse);
        require(mesh.values()==original,"mixed-rule raw-link reversal failed");
        std::vector<std::size_t> order(mesh.patches().size()); std::iota(order.begin(),order.end(),0);
        for(const auto& layer : mesh.commuting_layers(order)) {
            auto a=mesh,b=mesh;
            const auto selected=[&](std::size_t p) { return p>=fans.size() ? 0 : p%2 ? r1 : r2; };
            for(const auto p : layer) a.update_with_rule(p,selected(p));
            for(auto it=layer.rbegin();it!=layer.rend();++it) b.update_with_rule(*it,selected(*it));
            require(a.values()==b.values(),"mixed-rule commuting layer changed under reversal");
        }
        std::cout<<"Mixed arities/rules: 324 oracle/frame/cache checks; exact inverse, registry, guards and causal rule provenance passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
