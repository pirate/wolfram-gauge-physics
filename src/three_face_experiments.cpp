#include <iostream>
#include <random>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;
template<class Values> void array(const Values& values) {
    std::cout<<'['; bool first=true;
    for(const auto value : values) { if(!first) std::cout<<','; first=false; std::cout<<value; }
    std::cout<<']';
}
std::vector<std::size_t> histogram(const MeshDynamics& mesh) {
    std::vector<std::size_t> result(mesh.group().order());
    for(const auto a : mesh.face_values()) ++result[mesh.sector(a)];
    return result;
}

int main(int argc,char** argv) {
    try {
        std::size_t side=12,layers=256,trials=8; uint64_t seed=6290831;
        for(int i=1;i<argc;i+=2) {
            if(i+1>=argc) throw std::invalid_argument("missing option value");
            const std::string key=argv[i]; const auto value=std::stoull(argv[i+1]);
            if(key=="--side") side=value; else if(key=="--layers") layers=value;
            else if(key=="--trials") trials=value; else if(key=="--seed") seed=value;
            else throw std::invalid_argument("unknown option");
        }
        if(!layers || layers>8192 || !trials || trials>64) throw std::invalid_argument("bounded positive layers/trials required");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto n=group.order();
        std::array<TripleTable,3> rules;
        for(auto& rule : rules) {
            rule.entries.resize(n*n*n);
            for(auto& value : rule.entries) if(!(std::cin>>value)) throw std::invalid_argument("missing triple tables");
            validate_triple_table(group,rule);
        }
        for(std::size_t x=0;x<rules[0].entries.size();++x) {
            if(rules[0].entries[x]!=x && rules[1].entries[x]!=x) throw std::invalid_argument("collision and transport supports overlap");
            if(rules[2].entries[x]!=rules[1].entries[rules[0].entries[x]]) throw std::invalid_argument("combined table is not the stated composition");
        }
        const auto complex=triangulated_torus(side); const auto patches=three_face_patches(complex);
        MeshDynamics topology(complex,FiberBundleConnection(complex.base(),fiber),rules[2]);
        for(const auto& p : patches) topology.add_patch(p);
        const Vertex center=(side/2)*side+side/2;
        const auto prefix=static_cast<std::size_t>(std::find_if(patches.begin(),patches.end(),[&](const auto& p) {
            return p.faces[0][0]==center;
        })-patches.begin());
        if(prefix==patches.size()) throw std::logic_error("missing central fan");
        std::vector<std::set<std::size_t>> adjacent(complex.faces().size());
        for(const auto& fs : topology.incidence()) {
            if(fs.size()!=2) throw std::logic_error("expected closed manifold");
            adjacent[fs[0]].insert(fs[1]); adjacent[fs[1]].insert(fs[0]);
        }
        std::vector<std::vector<std::vector<std::size_t>>> schedules;
        for(std::size_t trial=0;trial<trials;++trial) {
            std::vector<std::size_t> order(patches.size()); std::iota(order.begin(),order.end(),0);
            std::mt19937_64 random(seed+trial*100003);
            for(std::size_t k=order.size();k>1;--k) std::swap(order[k-1],order[random()%k]);
            schedules.push_back(topology.commuting_layers(order));
        }
        std::cout<<"{\"schema\":1,\"side\":"<<side<<",\"layers\":"<<layers<<",\"trials\":"<<trials
                 <<",\"seed\":"<<seed<<",\"faces\":"<<complex.faces().size()<<",\"edges\":"<<complex.base().edges().size()
                 <<",\"patches\":"<<patches.size()<<",\"prefix_patch\":"<<prefix<<",\"schedules\":[";
        for(std::size_t t=0;t<trials;++t) {
            if(t) std::cout<<','; std::cout<<'[';
            for(std::size_t i=0;i<schedules[t].size();++i) { if(i) std::cout<<','; array(schedules[t][i]); }
            std::cout<<']';
        }
        std::cout<<"],\"trajectory_columns\":[\"layer\",\"collision_events\",\"transport_events\",\"nonflat_faces\",\"max_seed_distance\","
                    "\"sector0\",\"sector1\",\"sector2\",\"sector3\",\"sector4\",\"sector5\",\"sector6\",\"sector7\"],\"runs\":[";
        bool first=true;
        const std::array<std::string,3> modes{"collision","transport","combined"};
        const Vertex other=(side/2)*side+(side/2+std::max<std::size_t>(1,side/3))%side;
        const Vertex other_end=(other/side)*side+(other+1)%side;
        for(std::size_t trial=0;trial<trials;++trial) {
            for(const std::string condition : {"flat","reflection_link","rotation_link","reflection_pair","opposite_reflection_pair","witness_left","witness_right","random_links"}) {
                FiberBundleConnection initial(complex.base(),fiber);
                if(condition=="reflection_link" || condition=="reflection_pair" || condition=="opposite_reflection_pair")
                    initial.set_transport(center,center+1,group.elements()[1]);
                if(condition=="rotation_link") initial.set_transport(center,center+1,group.elements()[3]);
                if(condition=="reflection_pair") initial.set_transport(other,other_end,group.elements()[1]);
                if(condition=="opposite_reflection_pair") initial.set_transport(other,other_end,group.elements()[4]);
                if(condition=="witness_left" || condition=="witness_right") {
                    const Triple state=condition=="witness_left" ? Triple{1,1,4} : Triple{1,4,1};
                    const auto& p=patches[prefix]; const auto u=p.faces[0][0];
                    initial.set_transport(u,p.faces[0][2],group.elements()[group.inverse(state[2])]);
                    initial.set_transport(u,p.faces[1][2],group.elements()[group.multiply(group.inverse(state[2]),group.inverse(state[1]))]);
                    initial.set_transport(u,p.faces[2][2],group.elements()[group.inverse(triple_product(group,encode_triple(state,n)))]);
                }
                if(condition=="random_links") {
                    std::mt19937_64 random(seed+700001+trial);
                    for(const auto [u,v] : complex.base().edges()) initial.set_transport(u,v,group.elements()[random()%n]);
                }
                for(std::size_t mode=0;mode<modes.size();++mode) {
                    MeshDynamics mesh(complex,initial,rules[mode]);
                    for(const auto& p : patches) mesh.add_patch(p);
                    const auto original=mesh.values();
                    const auto missing=std::numeric_limits<std::size_t>::max();
                    std::vector<std::size_t> distance(complex.faces().size(),missing),queue;
                    for(std::size_t f=0;f<distance.size();++f) if(mesh.face_values()[f]!=mesh.identity()) { distance[f]=0; queue.push_back(f); }
                    for(std::size_t i=0;i<queue.size();++i) for(const auto next : adjacent[queue[i]])
                        if(distance[next]==missing) { distance[next]=distance[queue[i]]+1; queue.push_back(next); }
                    if(!first) std::cout<<','; first=false;
                    std::cout<<"{\"trial\":"<<trial<<",\"condition\":\""<<condition<<"\",\"mode\":\""<<modes[mode]<<"\",\"initial_links\":";
                    array(original); std::cout<<",\"trajectory\":[";
                    const auto layer_at=[&](std::size_t tick) {
                        return tick==1 ? std::vector<std::size_t>{prefix} : schedules[trial][(tick-2)%schedules[trial].size()];
                    };
                    std::vector<std::vector<std::size_t>> history;
                    std::size_t updates=0;
                    for(std::size_t tick=0;tick<=layers;++tick) {
                        std::size_t collisions=0,transports=0;
                        if(tick) for(const auto p : layer_at(tick)) {
                            const auto state=encode_triple(mesh.triple_values(p),n),next=rules[mode].entries[state];
                            if(state!=next) {
                                if(rules[0].entries[state]!=state) ++collisions; else ++transports;
                            }
                            mesh.update(p); ++updates;
                        }
                        const auto hist=histogram(mesh); history.push_back(hist);
                        std::size_t radius=0;
                        for(std::size_t f=0;f<distance.size();++f) if(mesh.face_values()[f]!=mesh.identity()) {
                            if(distance[f]==missing || distance[f]>2*tick) throw std::logic_error("fan propagation bound failed");
                            radius=std::max(radius,distance[f]);
                        }
                        if(tick%16==0 || tick==layers) {
                            const auto oracle=mesh.snapshot();
                            for(std::size_t f=0;f<complex.faces().size();++f)
                                if(oracle.holonomy(complex.faces()[f])!=group.elements()[mesh.face_values()[f]])
                                    throw std::logic_error("fan face cache disagrees with raw links");
                        }
                        if(tick) std::cout<<',';
                        std::cout<<'['<<tick<<','<<collisions<<','<<transports<<','<<complex.faces().size()-hist[mesh.identity()]<<','<<radius;
                        for(const auto h : hist) std::cout<<','<<h; std::cout<<']';
                    }
                    std::cout<<"],\"final_links\":"; array(mesh.values());
                    for(std::size_t tick=layers;tick>0;--tick) {
                        const auto layer=layer_at(tick);
                        for(auto it=layer.rbegin();it!=layer.rend();++it) mesh.update(*it,HurwitzDirection::Inverse);
                        if(histogram(mesh)!=history[tick-1]) throw std::logic_error("fan reverse histogram echo failed");
                    }
                    if(mesh.values()!=original) throw std::logic_error("fan reverse raw-link echo failed");
                    std::cout<<",\"updates\":"<<updates<<",\"exact_inverse_replay\":true,\"reverse_histogram_echo\":true}";
                }
            }
            std::cerr<<"trial "<<trial<<": eight conditions, three rule controls, "<<schedules[trial].size()<<" schedule layers; inverses verified\n";
        }
        std::cout<<"]}\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
