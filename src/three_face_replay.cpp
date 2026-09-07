#include <iostream>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;
template<class Values> void array(const Values& values) {
    std::cout<<'['; bool first=true;
    for(const auto value : values) { if(!first) std::cout<<','; first=false; std::cout<<value; }
    std::cout<<']';
}

int main() {
    try {
        std::size_t side;
        if(!(std::cin>>side) || side<3 || side>24) throw std::invalid_argument("side must be 3..24");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        TripleTable rule{std::vector<std::size_t>(512)};
        for(auto& entry : rule.entries) if(!(std::cin>>entry)) throw std::invalid_argument("missing triple table");
        validate_triple_table(group,rule);
        const auto complex=triangulated_torus(side);
        FiberBundleConnection input(complex.base(),fiber);
        for(const auto [u,v] : complex.base().edges()) {
            std::size_t value;
            if(!(std::cin>>value) || value>=group.order()) throw std::invalid_argument("invalid link value");
            input.set_transport(u,v,group.elements()[value]);
        }
        std::size_t count;
        if(!(std::cin>>count) || count>1000000) throw std::invalid_argument("event count exceeds replay bound");
        MeshDynamics mesh(complex,input,rule);
        for(const auto& spec : three_face_patches(complex)) mesh.add_patch(spec);
        std::vector<std::size_t> path(count);
        for(auto& p : path)
            if(!(std::cin>>p) || p>=mesh.patches().size()) throw std::invalid_argument("invalid replay patch");
        std::string extra;
        if(std::cin>>extra) throw std::invalid_argument("unexpected trailing input");
        const auto original=mesh.values();
        std::cout<<"{\"steps\":[";
        for(std::size_t i=0;i<path.size();++i) {
            const auto p=path[i],before=encode_triple(mesh.triple_values(p),8);
            mesh.update(p);
            if(i) std::cout<<',';
            std::cout<<'['<<p<<','<<before<<','<<encode_triple(mesh.triple_values(p),8)<<']';
        }
        std::cout<<"],\"final_links\":"; array(mesh.values());
        std::cout<<",\"final_face_holonomies\":"; array(mesh.face_values());
        for(auto it=path.rbegin();it!=path.rend();++it) mesh.update(*it,HurwitzDirection::Inverse);
        if(mesh.values()!=original) throw std::logic_error("replay inverse failed");
        std::cout<<",\"exact_link_inverse\":true}\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
