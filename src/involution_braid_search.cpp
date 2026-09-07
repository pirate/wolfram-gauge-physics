#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "equivariant_pair_search.hpp"

using namespace wgphysics::research;

int main(int argc,char** argv) {
    try {
        std::string path="out/involution-braid-search.json";
        if (argc==3 && std::string(argv[1])=="--output") path=argv[2];
        else if (argc!=1) throw std::invalid_argument("usage: wgphysics_involution_braid_search [--output FILE]");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto n=group.order();
        std::vector<std::vector<uint16_t>> action(n,std::vector<uint16_t>(n));
        for (uint16_t g=0;g<n;++g) for (uint16_t a=0;a<n;++a)
            action[g][a]=group.multiply(g,group.multiply(a,group.inverse(g)));
        std::vector<std::set<uint16_t>> subgroups;
        for (std::size_t x=0;x<static_cast<std::size_t>(n)*n;++x)
            subgroups.push_back(pair_generated_subgroup(group,x));
        struct Solution { PairTable table; bool strict; std::size_t growth; };
        std::vector<Solution> solutions;
        std::vector<uint16_t> spectator_left,spectator_right,spectator_input;
        std::size_t spectator_solution{};
        std::size_t growth_candidates=0;
        const auto count=visit_pair_involutions(group,[&](const PairTable& table) {
            std::size_t grown=0;
            for (std::size_t x=0;x<table.size();++x)
                if (subgroups[table[x]].size()>subgroups[x].size()) ++grown;
            if (grown) ++growth_candidates;
            bool strict=true,physical=true;
            for (uint16_t a=0;a<n && physical;++a) for (uint16_t b=0;b<n && physical;++b)
                for (uint16_t c=0;c<n && physical;++c) {
                    std::array<uint16_t,3> left{a,b,c},right=left;
                    const auto apply=[&](auto& values,std::size_t i) {
                        const auto result=table[values[i]*n+values[i+1]];
                        values[i]=result/n; values[i+1]=result%n;
                    };
                    apply(left,0); apply(left,1); apply(left,0);
                    apply(right,1); apply(right,0); apply(right,1);
                    if (left==right) continue;
                    strict=false; bool equivalent=false;
                    for (uint16_t g=0;g<n && !equivalent;++g)
                        equivalent=action[g][left[0]]==right[0] && action[g][left[1]]==right[1]
                                   && action[g][left[2]]==right[2];
                    physical=equivalent;
                }
            if (physical) {
                validate_pair_table(group,table); solutions.push_back({table,strict,grown});
                if (!strict && spectator_left.empty()) {
                    // A generating spectator pair removes the freedom to conjugate
                    // an isolated triple without affecting the surrounding state.
                    const auto r=group.index_of(Permutation({1,2,3,0}));
                    const auto s=group.index_of(Permutation({0,3,2,1}));
                    for (uint16_t a=0;a<n && spectator_left.empty();++a)
                        for (uint16_t b=0;b<n && spectator_left.empty();++b)
                            for (uint16_t c=0;c<n && spectator_left.empty();++c) {
                                std::vector<uint16_t> left{a,b,c,r,s},right=left;
                                const auto act=[&](auto& values,std::size_t i) {
                                    const auto y=table[values[i]*n+values[i+1]];
                                    values[i]=y/n; values[i+1]=y%n;
                                };
                                act(left,0);act(left,1);act(left,0);
                                act(right,1);act(right,0);act(right,1);
                                if (left==right) continue;
                                bool equivalent=false;
                                for (uint16_t g=0;g<n;++g) {
                                    bool matches=true;
                                    for (std::size_t i=0;i<left.size();++i)
                                        matches=matches && action[g][left[i]]==right[i];
                                    equivalent=equivalent || matches;
                                }
                                if (equivalent) throw std::logic_error("generating spectators failed to resolve gauge mismatch");
                                spectator_left=left;spectator_right=right;
                                spectator_input={a,b,c,r,s};spectator_solution=solutions.size()-1;
                            }
                }
            }
        });
        const auto parent=std::filesystem::path(path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot open output");
        out << "{\"schema\":1,\"fiber\":\"C4\",\"exhausted\":true,\"candidate_count\":" << count
            << ",\"subgroup_growing_candidates\":" << growth_candidates
            << ",\"isolated_triple_braid_solutions\":" << solutions.size() << ",\"solutions\":[\n";
        std::size_t strict_count=0,physical_growing=0,strict_growing=0;
        for (std::size_t i=0;i<solutions.size();++i) {
            const auto& solution=solutions[i];
            if (solution.strict) ++strict_count;
            if (solution.growth) ++physical_growing;
            if (solution.strict && solution.growth) ++strict_growing;
            if (i) out << ",\n";
            out << "{\"id\":" << i << ",\"strict_braid\":" << (solution.strict ? "true" : "false")
                << ",\"subgroup_grown_pairs\":" << solution.growth << ",\"table\":[";
            for (std::size_t j=0;j<solution.table.size();++j) { if (j) out << ','; out << solution.table[j]; }
            out << "]}";
        }
        out << "\n],\"strict_braid_solutions\":" << strict_count
            << ",\"subgroup_growing_isolated_triple_solutions\":" << physical_growing
            << ",\"subgroup_growing_strict_solutions\":" << strict_growing
            << ",\"spectator_counterexample\":{\"solution_id\":" << spectator_solution << ",\"input\":[";
        for (std::size_t i=0;i<spectator_input.size();++i) { if (i) out << ','; out << spectator_input[i]; }
        out << "],\"left\":[";
        for (std::size_t i=0;i<spectator_left.size();++i) { if (i) out << ','; out << spectator_left[i]; }
        out << "],\"right\":[";
        for (std::size_t i=0;i<spectator_right.size();++i) { if (i) out << ','; out << spectator_right[i]; }
        out << "]}}\n";
        if (!out) throw std::runtime_error("failed writing output");
        std::cout << count << " candidates; " << strict_count << " strict / " << solutions.size()
                  << " gauge-quotiented braid solutions; " << physical_growing << " grow subgroups\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
