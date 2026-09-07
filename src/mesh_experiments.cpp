#include <iostream>
#include <random>
#include "mesh_dynamics.hpp"

using namespace wgphysics::research;

template<class T> void array(std::ostream& out,const T& values) {
    out<<'['; bool first=true;
    for (const auto value : values) { if (!first) out<<','; first=false; out<<value; }
    out<<']';
}

void verify_snapshot(const MeshDynamics& mesh) {
    const auto oracle=mesh.snapshot();
    for (std::size_t f=0;f<mesh.face_values().size();++f)
        if (oracle.holonomy(mesh.complex().faces()[f])!=mesh.group().elements()[mesh.face_values()[f]])
            throw std::logic_error("shared-face cache differs from complete recomputation");
}

void census(const FiberGraph& fiber,const AutomorphismTables& group,const PairTable& table,bool shared=false) {
    const OrientedCellComplex complex(BaseGraph({0,1,2,3,4,5},
        {{0,1},{0,2},{0,3},{0,4},{1,2},{1,4},{2,3},{2,5},{3,5}}),
        {{0,1,2,0},{0,2,3,0},{0,4,1,0},{2,5,3,2}});
    const CellPairPatch patch{{1,2,0,1},{3,0,2,3},{3,0,2,1}};
    const SharedEdgePatch edge_patch{{0,2,3,0},{0,1,2,0}};
    const std::vector<Vertex> outer{0,4,1,2,5,3,0};
    const auto n=group.order();
    const auto make_connection=[&](uint16_t a,uint16_t b,uint16_t c,uint16_t d) {
        FiberBundleConnection connection(complex.base(),fiber);
        connection.set_transport(1,2,group.elements()[a]);
        connection.set_transport(2,3,group.elements()[b]);
        connection.set_transport(1,4,group.elements()[group.inverse(c)]);
        connection.set_transport(3,5,group.elements()[group.inverse(group.multiply(b,d))]);
        return connection;
    };
    std::cout<<"\"local_census\":{\"columns\":[\"A\",\"B\",\"C\",\"D\",\"A_after\",\"B_after\",\"C_after\",\"D_after\"],"
                "\"faces\":[[0,1,2,0],[0,2,3,0],[0,4,1,0],[2,5,3,2]],"
                "\"patch\":";
    if(shared) std::cout<<"{\"first\":[0,2,3,0],\"second\":[0,1,2,0],\"connector\":[0]}";
    else std::cout<<"{\"first\":[1,2,0,1],\"second\":[3,0,2,3],\"connector\":[3,0,2,1]}";
    std::cout<<",\"outer_boundary\":[0,4,1,2,5,3,0],\"transitions\":[\n";
    bool first=true;
    for (uint16_t a=0;a<n;++a) for (uint16_t b=0;b<n;++b)
        for (uint16_t c=0;c<n;++c) for (uint16_t d=0;d<n;++d) {
            const auto connection=make_connection(a,b,c,d);
            MeshDynamics mesh(complex,connection,table);
            if (mesh.face_values()!=std::vector<uint16_t>{a,b,c,d}) throw std::logic_error("census gauge parametrization failed");
            if(shared) mesh.add_patch(edge_patch); else mesh.add_patch(patch);
            mesh.update(0);
            const auto targets=[&](const auto& x,const auto& y) {
                const auto result=table[group.index_of(x)*n+group.index_of(y)];
                return std::pair{group.elements()[result/n],group.elements()[result%n]};
            };
            const auto oracle=shared ? realize_shared_edge_interaction(complex,connection,edge_patch,targets)
                                     : realize_pair_interaction(connection,patch,targets);
            const auto actual=mesh.snapshot();
            for (const auto [u,v] : complex.base().edges())
                if (actual.edge_transport(u,v)!=oracle.edge_transport(u,v)) throw std::logic_error("local link oracle disagrees");
            if (actual.holonomy(outer)!=connection.holonomy(outer)) throw std::logic_error("outer boundary changed");
            verify_snapshot(mesh);
            if (!first) std::cout<<",\n"; first=false;
            std::cout<<'['<<a<<','<<b<<','<<c<<','<<d;
            for (const auto value : mesh.face_values()) std::cout<<','<<value;
            std::cout<<']';
            mesh.update(0,HurwitzDirection::Inverse);
            if (mesh.face_values()!=std::vector<uint16_t>{a,b,c,d}) throw std::logic_error("local inverse failed");
            const auto reversed=mesh.snapshot();
            for (const auto [u,v] : complex.base().edges())
                if (reversed.edge_transport(u,v)!=connection.edge_transport(u,v)) throw std::logic_error("local inverse changed a link");
        }
    if(shared) {
        std::cout<<"]}";
        std::cerr<<"4096 shared-edge states: target faces, unchanged spectators, boundary and inverse checks passed\n";
        return;
    }
    std::cout<<"],\"fixed_boundary_witness\":{";
    const auto r=group.index_of(Permutation({1,2,3,0}));
    const auto t=group.index_of(Permutation({0,3,2,1}));
    const auto id=group.index_of(Permutation::identity(4));
    const auto opposite=group.multiply(group.multiply(r,r),t);
    const auto first_connection=make_connection(id,t,t,t);
    const auto other_connection=make_connection(id,t,t,opposite);
    // Match every exterior link, not merely the boundary's conjugacy class.
    // g_v = U_first(u,v) g_u U_other(u,v)^-1 along the boundary.
    std::map<Vertex,Permutation> frames;
    bool found=false;
    for (const auto& g : group.elements())
        if (Permutation::compose(g,Permutation::compose(other_connection.holonomy(outer),g.inverse()))==first_connection.holonomy(outer)) {
            frames.emplace(outer.front(),g); found=true; break;
        }
    if (!found) throw std::logic_error("witness boundary holonomies not conjugate");
    for (std::size_t i=1;i+1<outer.size();++i) {
        const auto u=outer[i-1],v=outer[i];
        frames.emplace(v,Permutation::compose(first_connection.edge_transport(u,v),
            Permutation::compose(frames.at(u),other_connection.edge_transport(u,v).inverse())));
    }
    const auto second_connection=other_connection.gauge_transform(frames);
    for (std::size_t i=1;i<outer.size();++i)
        if (first_connection.edge_transport(outer[i-1],outer[i])!=second_connection.edge_transport(outer[i-1],outer[i]))
            throw std::logic_error("witness exterior links differ");
    if (first_connection.gauge_invariant_signature()==second_connection.gauge_invariant_signature())
        throw std::logic_error("witness is only a gauge-frame difference");
    std::cout<<"\"edges\":[";
    bool first_edge=true;
    for(const auto [u,v] : complex.base().edges()) {
        if(!first_edge) std::cout<<','; first_edge=false; std::cout<<'['<<u<<','<<v<<']';
    }
    std::cout<<"],\"states\":[";
    std::vector<uint16_t> initial_sectors,first_out;
    bool first_state=true;
    for(const auto* input : {&first_connection,&second_connection}) {
        MeshDynamics state(complex,*input,table); state.add_patch(patch);
        std::vector<uint16_t> before;
        for(const auto a : state.face_values()) before.push_back(state.sector(a));
        if (!first_state && initial_sectors!=before) throw std::logic_error("witness incoming face sectors differ");
        if(first_state) initial_sectors=before;
        if(!first_state) std::cout<<',';
        std::cout<<"{\"initial_links\":"; array(std::cout,state.values());
        std::cout<<",\"initial_face_sectors\":"; array(std::cout,before);
        state.update(0);
        std::vector<uint16_t> after;
        for(const auto a : state.face_values()) after.push_back(state.sector(a));
        if(first_state) first_out=after;
        else if(first_out==after) throw std::logic_error("selected law has no fixed-boundary witness response");
        std::cout<<",\"final_links\":"; array(std::cout,state.values());
        std::cout<<",\"final_face_sectors\":"; array(std::cout,after); std::cout<<'}';
        first_state=false;
    }
    std::cout<<"],\"identical_exterior_links\":true,\"initial_states_gauge_equivalent\":false}}";
    std::cerr<<"4096 shared-incidence states: microscopic oracle, boundary and inverse checks passed\n";
}

int main(int argc,char** argv) {
    try {
        std::size_t side=12,layers=256,trials=8;
        uint64_t seed=819031;
        bool shared=false;
        for (int i=1;i<argc;i+=2) {
            if (i+1>=argc) throw std::invalid_argument("missing option value");
            const std::string key=argv[i];
            if(key=="--lift") {
                const std::string value=argv[i+1];
                if(value!="shared-edge" && value!="exclusive-closing") throw std::invalid_argument("unknown mesh lift");
                shared=value=="shared-edge"; continue;
            }
            const auto value=std::stoull(argv[i+1]);
            if (key=="--side") side=value;
            else if (key=="--layers") layers=value;
            else if (key=="--trials") trials=value;
            else if (key=="--seed") seed=value;
            else throw std::invalid_argument("unknown option");
        }
        if (layers>8192 || trials>64) throw std::invalid_argument("bounded runner requires layers<=8192, trials<=64");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        PairTable table(group.order()*group.order());
        for (auto& value : table) if (!(std::cin>>value)) throw std::invalid_argument("missing pair table");
        validate_pair_table(group,table);
        const auto complex=triangulated_torus(side);
        const auto specifications=shared ? std::vector<CellPairPatch>{} : triangular_mesh_patches(complex);
        const auto edge_specifications=shared ? shared_edge_patches(complex) : std::vector<SharedEdgePatch>{};
        const auto patch_count=specifications.size()+edge_specifications.size();
        const auto add_patches=[&](MeshDynamics& mesh) {
            for(const auto& p : specifications) mesh.add_patch(p);
            for(const auto& p : edge_specifications) mesh.add_patch(p);
        };
        const Permutation reflection({0,3,2,1}),rotation({1,2,3,0});
        std::cout<<"{\"schema\":1,\"side\":"<<side<<",\"layers\":"<<layers<<",\"trials\":"<<trials
                 <<",\"seed\":"<<seed<<",\"vertices\":"<<complex.base().vertices().size()
                 <<",\"edges\":"<<complex.base().edges().size()<<",\"faces\":"<<complex.faces().size()
                 <<",\"patches\":"<<patch_count<<",\"seed_group_elements\":{\"reflection\":"
                 <<group.index_of(reflection)<<",\"rotation\":"<<group.index_of(rotation)<<"},\"group_permutations\":[";
        for (std::size_t i=0;i<group.elements().size();++i) { if(i) std::cout<<','; array(std::cout,group.elements()[i].image()); }
        std::cout<<"],";
        if(shared) std::cout<<"\"lift\":\"shared-edge\",";
        census(fiber,group,table,shared);
        std::cout<<",\"trajectory_columns\":[\"layer\",\"active_faces\",\"max_seed_dual_distance\",\"components\",\"largest_component\","
                    "\"sector0\",\"sector1\",\"sector2\",\"sector3\",\"sector4\",\"sector5\",\"sector6\",\"sector7\"],\"runs\":[\n";
        bool first_run=true;
        const Vertex center=(side/2)*side+side/2;
        const Vertex other=(side/2)*side+(side/2+side/4)%side;
        if (Permutation::compose(reflection,rotation)==Permutation::compose(rotation,reflection))
            throw std::logic_error("pair seeds commute");
        for (std::size_t trial=0;trial<trials;++trial) {
            MeshDynamics topology(complex,FiberBundleConnection(complex.base(),fiber),table);
            add_patches(topology);
            std::vector<std::size_t> order(patch_count); std::iota(order.begin(),order.end(),0);
            std::mt19937_64 scheduler(seed+trial*100003);
            for (std::size_t i=order.size();i>1;--i) std::swap(order[i-1],order[scheduler()%i]);
            const auto schedule=topology.commuting_layers(order);
            std::vector<std::set<std::size_t>> adjacent(complex.faces().size());
            for (const auto& fs : topology.incidence()) {
                if (fs.size()!=2) throw std::logic_error("torus has boundary/nonmanifold edge");
                adjacent[fs[0]].insert(fs[1]); adjacent[fs[1]].insert(fs[0]);
            }
            for (const std::string condition : {"flat","reflection","rotation","noncommuting_pair","random_links"}) {
                FiberBundleConnection initial(complex.base(),fiber);
                if (condition=="reflection" || condition=="noncommuting_pair") initial.set_transport(center,center+1,reflection);
                if (condition=="rotation" || condition=="noncommuting_pair")
                    initial.set_transport(other,(other+side)%(side*side),rotation);
                if (condition=="random_links") {
                    std::mt19937_64 random(seed+700001+trial);
                    for (const auto [u,v] : complex.base().edges()) initial.set_transport(u,v,group.elements()[random()%group.order()]);
                }
                MeshDynamics mesh(complex,initial,table);
                add_patches(mesh);
                const auto initial_values=mesh.values();
                const auto missing=std::numeric_limits<std::size_t>::max();
                std::vector<std::size_t> distance(complex.faces().size(),missing),frontier;
                for (std::size_t f=0;f<distance.size();++f) if (mesh.face_values()[f]!=mesh.identity()) {
                    distance[f]=0; frontier.push_back(f);
                }
                for (std::size_t i=0;i<frontier.size();++i) for (const auto next : adjacent[frontier[i]])
                    if (distance[next]==missing) { distance[next]=distance[frontier[i]]+1; frontier.push_back(next); }
                if (!first_run) std::cout<<",\n"; first_run=false;
                std::cout<<"{\"trial\":"<<trial<<",\"condition\":\""<<condition<<"\",\"schedule_layers\":"<<schedule.size()
                         <<",\"seed_edges\":[["<<center<<','<<center+1<<"],["<<other<<','<<(other+side)%(side*side)
                         <<"]]";
                if(shared) {
                    std::cout<<",\"initial_links\":"; array(std::cout,initial_values);
                    std::cout<<",\"schedule\":[";
                    for(std::size_t i=0;i<schedule.size();++i) { if(i) std::cout<<','; array(std::cout,schedule[i]); }
                    std::cout<<']';
                }
                std::cout<<",\"trajectory\":[";
                std::size_t updates=0;
                std::vector<std::vector<std::size_t>> forward_histograms;
                for (std::size_t tick=0;tick<=layers;++tick) {
                    if (tick) for (const auto p : schedule[(tick-1)%schedule.size()]) { mesh.update(p); ++updates; }
                    std::vector<std::size_t> histogram(group.order());
                    std::vector<bool> active(distance.size()),visited(distance.size());
                    std::size_t count=0,radius=0,components=0,largest=0;
                    for (std::size_t f=0;f<distance.size();++f) {
                        ++histogram[mesh.sector(mesh.face_values()[f])];
                        if (mesh.face_values()[f]!=mesh.identity()) {
                            active[f]=true; ++count;
                            if (distance[f]==missing || distance[f]>(shared ? tick : 2*tick))
                                throw std::logic_error("curvature escaped the local propagation bound");
                            radius=std::max(radius,distance[f]);
                        }
                    }
                    for (std::size_t f=0;f<distance.size();++f) if (active[f] && !visited[f]) {
                        std::vector<std::size_t> component{f}; visited[f]=true; ++components;
                        for (std::size_t i=0;i<component.size();++i) for (const auto next : adjacent[component[i]])
                            if (active[next] && !visited[next]) { visited[next]=true; component.push_back(next); }
                        largest=std::max(largest,component.size());
                    }
                    if (condition=="flat" && count) throw std::logic_error("flat mesh gained curvature");
                    if (tick%16==0 || tick==layers) verify_snapshot(mesh);
                    forward_histograms.push_back(histogram);
                    if(tick) std::cout<<',';
                    std::cout<<'['<<tick<<','<<count<<','<<radius<<','<<components<<','<<largest;
                    for(const auto value : histogram) std::cout<<','<<value;
                    std::cout<<']';
                }
                std::cout<<']';
                if(shared) { std::cout<<",\"final_links\":"; array(std::cout,mesh.values()); }
                for (std::size_t tick=layers;tick>0;--tick) {
                    const auto& layer=schedule[(tick-1)%schedule.size()];
                    for(auto it=layer.rbegin();it!=layer.rend();++it) mesh.update(*it,HurwitzDirection::Inverse);
                    std::vector<std::size_t> histogram(group.order());
                    for (const auto value : mesh.face_values()) ++histogram[mesh.sector(value)];
                    if (histogram!=forward_histograms[tick-1]) throw std::logic_error("reverse echo differs at intermediate layer");
                }
                if (mesh.values()!=initial_values) throw std::logic_error("mesh experiment inverse replay failed");
                std::cout<<",\"updates\":"<<updates<<",\"exact_inverse_replay\":true,\"reverse_histogram_echo\":true}";
            }
            std::cerr<<"trial "<<trial<<": five mesh conditions, "<<schedule.size()<<" conflict-free schedule layers, inverse replay passed\n";
        }
        std::cout<<"]}\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
