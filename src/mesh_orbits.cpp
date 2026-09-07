#include <iostream>
#include <random>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;

template<class T> void json_array(const T& values) {
    std::cout<<'['; bool first=true;
    for(const auto v : values) { if(!first) std::cout<<','; first=false; std::cout<<v; }
    std::cout<<']';
}

void json_patches(const std::vector<CellPairPatch>& patches) {
    std::cout<<'[';
    for(std::size_t i=0;i<patches.size();++i) {
        if(i) std::cout<<',';
        std::cout<<"{\"first\":"; json_array(patches[i].first_loop);
        std::cout<<",\"second\":"; json_array(patches[i].second_loop);
        std::cout<<",\"connector\":"; json_array(patches[i].connector); std::cout<<'}';
    }
    std::cout<<']';
}

std::vector<std::size_t> features(const MeshDynamics& mesh,const std::vector<std::size_t>& pair_class,std::size_t count) {
    std::vector<std::size_t> result(count);
    for(std::size_t p=0;p<mesh.patches().size();++p) {
        const auto [a,b]=mesh.based_pair(p);
        ++result[pair_class[a*mesh.group().order()+b]];
    }
    return result;
}

OrientedCellComplex octahedron() {
    std::set<Edge> edges; std::vector<std::vector<Vertex>> faces;
    for(Vertex a=0;a<2;++a) for(Vertex b=2;b<4;++b) for(Vertex c=4;c<6;++c) {
        std::vector<Vertex> face=(a+b+c)%2 ? std::vector<Vertex>{a,c,b,a} : std::vector<Vertex>{a,b,c,a};
        for(std::size_t i=1;i<face.size();++i) edges.insert(wgphysics::infragauge::canonical_edge(face[i-1],face[i]));
        faces.push_back(face);
    }
    return OrientedCellComplex(BaseGraph({0,1,2,3,4,5},{edges.begin(),edges.end()}),faces);
}

int main() {
    try {
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto n=group.order();
        PairTable table(n*n);
        for(auto& value : table) if(!(std::cin>>value)) throw std::invalid_argument("missing pair table");
        validate_pair_table(group,table);
        std::map<std::size_t,std::size_t> class_ids;
        std::vector<std::size_t> pair_class(n*n),class_representatives;
        for(std::size_t pair=0;pair<table.size();++pair) {
            auto best=pair;
            for(uint16_t g=0;g<n;++g) best=std::min(best,conjugated_pair(group,pair,g));
            if(!class_ids.contains(best)) { class_ids[best]=class_ids.size(); class_representatives.push_back(best); }
            pair_class[pair]=class_ids.at(best);
        }
        const OrientedCellComplex tetra(BaseGraph({0,1,2,3},{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}),
                                        {{0,1,2,0},{0,2,3,0},{0,3,1,0},{1,3,2,1}});
        const auto patches=triangular_mesh_patches(tetra);
        std::map<std::vector<Vertex>,std::size_t> state_ids;
        std::vector<std::vector<uint16_t>> representatives;
        for(uint16_t a=0;a<n;++a) for(uint16_t b=0;b<n;++b) for(uint16_t c=0;c<n;++c) {
            FiberBundleConnection connection(tetra.base(),fiber);
            connection.set_transport(1,2,group.elements()[a]);
            connection.set_transport(1,3,group.elements()[b]);
            connection.set_transport(2,3,group.elements()[c]);
            const auto key=connection.gauge_invariant_signature();
            if(!state_ids.contains(key)) { state_ids[key]=representatives.size(); representatives.push_back({a,b,c}); }
        }
        std::cout<<"{\"schema\":1,\"gauge_quotient\":\"local frames on a fixed labeled base; no vertex-isomorphism quotient\","
                    "\"pair_class_representatives\":"; json_array(class_representatives);
        std::cout<<",\"pair_class_by_raw_pair\":"; json_array(pair_class);
        std::cout<<",\"tetrahedron\":{\"rooted_states\":"<<n*n*n<<",\"patches\":"<<patches.size()<<",\"patch_specs\":";
        json_patches(patches); std::cout<<",\"states\":[\n";
        for(std::size_t id=0;id<representatives.size();++id) {
            const auto& r=representatives[id];
            FiberBundleConnection connection(tetra.base(),fiber);
            connection.set_transport(1,2,group.elements()[r[0]]);
            connection.set_transport(1,3,group.elements()[r[1]]);
            connection.set_transport(2,3,group.elements()[r[2]]);
            MeshDynamics mesh(tetra,connection,table);
            for(const auto& p : patches) mesh.add_patch(p);
            const auto before=mesh.values();
            std::vector<std::size_t> next,stabilizer;
            for(uint16_t g=0;g<n;++g)
                if(std::all_of(r.begin(),r.end(),[&](uint16_t a) { return group.multiply(g,a)==group.multiply(a,g); }))
                    stabilizer.push_back(g);
            const auto f=features(mesh,pair_class,class_representatives.size());
            std::vector<uint16_t> sectors;
            for(const auto h : mesh.face_values()) sectors.push_back(mesh.sector(h));
            for(std::size_t p=0;p<patches.size();++p) {
                mesh.update(p);
                const auto key=mesh.snapshot().gauge_invariant_signature();
                if(!state_ids.contains(key)) throw std::logic_error("tetrahedron transition escaped exhaustive quotient");
                next.push_back(state_ids.at(key));
                mesh.update(p,HurwitzDirection::Inverse);
                if(mesh.values()!=before) throw std::logic_error("tetrahedron inverse failed");
            }
            if(id) std::cout<<",\n";
            std::cout<<"{\"id\":"<<id<<",\"rooted_chords\":"; json_array(r);
            std::cout<<",\"stabilizer\":"; json_array(stabilizer);
            std::cout<<",\"face_sectors\":"; json_array(sectors);
            std::cout<<",\"pair_features\":"; json_array(f);
            std::cout<<",\"successors\":"; json_array(next); std::cout<<'}';
        }
        std::cout<<"]},\"octahedron\":{\"seed\":993121,\"samples\":128,\"steps_per_sample\":16,\"edges\":[";
        const auto oct=octahedron(); const auto oct_patches=triangular_mesh_patches(oct);
        bool first_edge=true;
        for(const auto [u,v] : oct.base().edges()) { if(!first_edge) std::cout<<','; first_edge=false; std::cout<<'['<<u<<','<<v<<']'; }
        std::cout<<"],\"patches\":"; json_patches(oct_patches);
        std::cout<<",\"checks\":[\n";
        std::mt19937_64 random(993121);
        for(std::size_t sample=0;sample<128;++sample) {
            FiberBundleConnection connection(oct.base(),fiber);
            for(const auto [u,v] : oct.base().edges()) connection.set_transport(u,v,group.elements()[random()%n]);
            MeshDynamics mesh(oct,connection,table);
            for(const auto& p : oct_patches) mesh.add_patch(p);
            auto f=features(mesh,pair_class,class_representatives.size());
            for(std::size_t step=0;step<16;++step) {
                const auto p=random()%oct_patches.size(); const auto links=mesh.values();
                mesh.update(p); const auto after=features(mesh,pair_class,class_representatives.size());
                if(sample || step) std::cout<<",\n";
                std::cout<<"{\"sample\":"<<sample<<",\"step\":"<<step<<",\"patch\":"<<p<<",\"before_links\":"; json_array(links);
                std::cout<<",\"before_features\":"; json_array(f);
                std::cout<<",\"after_features\":"; json_array(after); std::cout<<'}';
                f=after;
            }
        }
        std::cout<<"]}}\n";
        std::cerr<<representatives.size()<<" exact tetrahedron gauge states; "<<patches.size()
                 <<" generators; 2048 octahedron checks exported\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
