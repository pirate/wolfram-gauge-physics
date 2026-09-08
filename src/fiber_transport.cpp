#include <array>
#include <iostream>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;

template<class Values> void print_array(const Values& values) {
    std::cout << '[';
    bool first=true;
    for(const auto value : values) { if(!first) std::cout << ','; first=false; std::cout << value; }
    std::cout << ']';
}

int main() {
    try {
        std::size_t side,vertices,edge_count,attempts;
        if(!(std::cin>>side>>vertices>>edge_count>>attempts) || side<3 || side>24 ||
           vertices<1 || vertices>6 || edge_count>vertices*(vertices-1)/2 || attempts>100000)
            throw std::invalid_argument("invalid bounded fiber transport dimensions");
        std::set<Edge> unique;
        for(std::size_t i=0;i<edge_count;++i) {
            Vertex u,v;
            if(!(std::cin>>u>>v) || u>=vertices || v>=vertices || u==v ||
               !unique.insert(wgphysics::infragauge::canonical_edge(u,v)).second)
                throw std::invalid_argument("invalid simple fiber edge");
        }
        const FiberGraph fiber(static_cast<uint32_t>(vertices),{unique.begin(),unique.end()});
        const auto automorphisms=fiber.automorphisms();
        if(automorphisms.size()>24) throw std::invalid_argument("bounded transport supports at most 24 derived automorphisms");
        const AutomorphismTables group(automorphisms);
        const auto n=group.order(),identity=group.index_of(Permutation::identity(vertices));
        PairTable vacancy(n*n);
        for(std::size_t a=0;a<n;++a) for(std::size_t b=0;b<n;++b)
            vacancy[a*n+b]=((a==identity)!=(b==identity)) ? b*n+a : a*n+b;
        validate_pair_table(group,vacancy);
        const auto complex=triangulated_torus(side);
        FiberBundleConnection input(complex.base(),fiber);
        for(const auto [u,v] : complex.base().edges()) {
            std::size_t value;
            if(!(std::cin>>value) || value>=n) throw std::invalid_argument("invalid derived-group link");
            input.set_transport(u,v,group.elements()[value]);
        }
        const auto patches=shared_edge_patches(complex);
        std::vector<std::size_t> schedule(attempts);
        for(auto& patch : schedule)
            if(!(std::cin>>patch) || patch>=patches.size()) throw std::invalid_argument("invalid rooted transport patch");
        std::string extra;
        if(std::cin>>extra) throw std::invalid_argument("trailing fiber transport input");

        MeshDynamics mesh(complex,input,vacancy);
        for(const auto& patch : patches) mesh.add_patch(patch);
        const auto initial=mesh.values(),initial_faces=mesh.face_values();
        const auto initial_signature=input.gauge_invariant_signature();
        std::vector<std::array<std::size_t,4>> events;
        for(std::size_t tick=0;tick<attempts;++tick) {
            const auto patch=schedule[tick];
            const auto [a,b]=mesh.based_pair(patch);
            const std::size_t code=a*n+b,target=vacancy[code];
            mesh.update(patch);
            if(code!=target) events.push_back({tick+1,patch,code,target});
        }
        const auto final_links=mesh.values(),final_faces=mesh.face_values();
        const auto final=mesh.snapshot();
        const auto final_signature=final.gauge_invariant_signature();
        for(std::size_t f=0;f<complex.faces().size();++f)
            if(group.elements()[final_faces[f]]!=final.holonomy(complex.faces()[f]))
                throw std::logic_error("generic transport face cache differs from raw connection");
        for(auto i=schedule.rbegin();i!=schedule.rend();++i)
            mesh.update(*i,HurwitzDirection::Inverse);
        if(mesh.values()!=initial || mesh.face_values()!=initial_faces)
            throw std::logic_error("generic transport inverse failed");

        std::cout << "{\"automorphisms\":[";
        for(std::size_t i=0;i<n;++i) { if(i) std::cout << ','; print_array(group.elements()[i].image()); }
        std::cout << "],\"events\":[";
        for(std::size_t i=0;i<events.size();++i) { if(i) std::cout << ','; print_array(events[i]); }
        std::cout << "],\"initial_faces\":"; print_array(initial_faces);
        std::cout << ",\"final_faces\":"; print_array(final_faces);
        std::cout << ",\"final_links\":"; print_array(final_links);
        std::cout << ",\"initial_signature\":"; print_array(initial_signature);
        std::cout << ",\"final_signature\":"; print_array(final_signature);
        std::cout << ",\"exact_link_inverse\":true,\"face_cache_matches_raw_connection\":true}\n";
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
