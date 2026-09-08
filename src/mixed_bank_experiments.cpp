#include <iostream>
#include "mesh_dynamics.hpp"
#include "cycle_fiber_cli.hpp"

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
        const bool raw_events=argc>1 && std::string_view(argv[argc-1])=="--raw-events";
        if(raw_events) --argc;
        const auto fiber=cycle_fiber_from_arguments(argc,argv);
        std::size_t emitted_link_values=0;
        std::size_t side,rule_count,condition_count,attempts,stride;
        if(!(std::cin>>side>>rule_count>>condition_count>>attempts>>stride) || side<3 || side>48 ||
           !rule_count || rule_count>144 || !condition_count || condition_count>16 ||
           !attempts || attempts>1000000 || !stride || attempts%stride)
            throw std::invalid_argument("invalid bounded bank experiment dimensions");
        const AutomorphismTables group(fiber.automorphisms());
        const auto n=group.order();
        PairTable pair(n*n);
        for(auto& x : pair) if(!(std::cin>>x)) throw std::invalid_argument("missing pair table");
        validate_pair_table(group,pair);
        std::vector<TripleTable> triples(rule_count,TripleTable{std::vector<std::size_t>(n*n*n)});
        for(auto& rule : triples) {
            for(auto& x : rule.entries) if(!(std::cin>>x)) throw std::invalid_argument("missing triple table");
            validate_triple_table(group,rule);
        }
        const auto complex=triangulated_torus(side);
        const auto fans=three_face_patches(complex);
        const auto pairs=shared_edge_patches(complex);
        if(fans.size()!=pairs.size()) throw std::logic_error("runner expects equal rooted support counts");
        std::vector<FiberBundleConnection> inputs;
        for(std::size_t c=0;c<condition_count;++c) {
            FiberBundleConnection input(complex.base(),fiber);
            for(const auto [u,v] : complex.base().edges()) {
                std::size_t a;
                if(!(std::cin>>a) || a>=n) throw std::invalid_argument("invalid initial link");
                input.set_transport(u,v,group.elements()[a]);
            }
            inputs.push_back(std::move(input));
        }
        std::vector<std::size_t> schedule(attempts);
        for(auto& event : schedule)
            if(!(std::cin>>event) || event>=(rule_count+1)*fans.size()) throw std::invalid_argument("invalid bank schedule");
        std::string extra;
        if(std::cin>>extra) throw std::invalid_argument("trailing bank experiment input");
        std::cout<<"{\"runs\":["; bool first_run=true;
        for(std::size_t c=0;c<inputs.size();++c) for(const bool combined : {false,true}) {
            MeshDynamics mesh(complex,inputs[c],pair);
            std::vector<std::size_t> rule_ids{0};
            for(const auto& rule : triples) rule_ids.push_back(mesh.register_rule(rule));
            for(const auto& p : fans) mesh.add_patch(p,rule_ids[1]);
            for(const auto& p : pairs) mesh.add_patch(p,0);
            const auto initial=mesh.values();
            std::vector<std::vector<std::size_t>> history{histogram(mesh)};
            if(!first_run) std::cout<<','; first_run=false;
            std::cout<<"{\"condition\":"<<c<<",\"mode\":\""<<(combined?"combined":"transport")<<"\",\"events\":[";
            bool first_event=true;
            const auto decode=[&](std::size_t encoded) {
                const auto rule=encoded/fans.size(),p=encoded%fans.size();
                return std::pair{rule,rule ? p : fans.size()+p};
            };
            for(std::size_t tick=0;tick<attempts;++tick) {
                const auto [rule,p]=decode(schedule[tick]);
                if(combined || !rule) {
                    std::size_t code;
                    if(rule) code=encode_triple(mesh.triple_values(p),n);
                    else { const auto [a,b]=mesh.based_pair(p); code=a*n+b; }
                    const auto target=mesh.rules()[rule_ids[rule]].forward[code];
                    mesh.update_with_rule(p,rule_ids[rule]);
                    if(target!=code) {
                        if(!first_event) std::cout<<','; first_event=false;
                        std::cout<<'['<<tick+1<<','<<rule<<','<<p%fans.size()<<','<<code<<','<<target;
                        if(raw_events) {
                            // Stream bounded diagnostic snapshots without storing trajectories
                            // in the engine. Default output remains byte-for-byte unchanged.
                            if(mesh.values().size()>8000000-emitted_link_values)
                                throw std::invalid_argument("raw event output exceeds eight million link values");
                            emitted_link_values+=mesh.values().size();
                            std::cout<<','; array(mesh.values());
                        }
                        std::cout<<']';
                    }
                }
                if((tick+1)%stride==0) history.push_back(histogram(mesh));
            }
            std::cout<<"],\"histograms\":[";
            for(std::size_t i=0;i<history.size();++i) { if(i) std::cout<<','; array(history[i]); }
            std::cout<<"],\"final_links\":"; array(mesh.values());
            const auto final=mesh.snapshot();
            for(std::size_t f=0;f<complex.faces().size();++f)
                if(group.elements()[mesh.face_values()[f]]!=final.holonomy(complex.faces()[f]))
                    throw std::logic_error("bank face cache differs from raw connection");
            for(std::size_t tick=attempts;tick>0;--tick) {
                const auto [rule,p]=decode(schedule[tick-1]);
                if(combined || !rule) mesh.update_with_rule(p,rule_ids[rule],HurwitzDirection::Inverse);
                if((tick-1)%stride==0 && histogram(mesh)!=history[(tick-1)/stride])
                    throw std::logic_error("bank intermediate reverse histogram echo failed");
            }
            if(mesh.values()!=initial) throw std::logic_error("bank raw-link inverse failed");
            std::cout<<",\"exact_link_inverse\":true,\"reverse_histogram_echo\":true}";
        }
        std::cout<<"]}\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
