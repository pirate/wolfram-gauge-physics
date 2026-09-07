#include <filesystem>
#include <fstream>
#include <iostream>
#include "cell_chain.hpp"

using namespace wgphysics::research;

int main(int argc,char** argv) {
    try {
        std::string path="out/pair-involutions.json";
        if (argc==3 && std::string(argv[1])=="--output") path=argv[2];
        else if (argc!=1) throw std::invalid_argument("usage: wgphysics_involution_census [--output FILE]");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto rules=search_pair_involutions(group);
        const auto n=group.order(), id=group.index_of(Permutation::identity(4));
        std::vector<std::set<uint16_t>> subgroups;
        for (std::size_t x=0;x<static_cast<std::size_t>(n)*n;++x)
            subgroups.push_back(pair_generated_subgroup(group,x));
        const auto parent=std::filesystem::path(path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot open output");
        out << "{\"schema\":1,\"fiber\":\"C4\",\"group_order\":" << n
            << ",\"construction\":\"minimal transpositions closed under conjugation and orientation reversal\","
               "\"claim\":\"exhaustive within this construction only\",\"rules\":[\n";
        std::size_t escaping_rules=0;
        for (std::size_t index=0;index<rules.size();++index) {
            const auto& rule=rules[index];
            std::size_t moved=0,grown=0,changed=0,braid_failures=0;
            for (std::size_t x=0;x<rule.size();++x) {
                if (rule[x]!=x) ++moved;
                if (subgroups[rule[x]]!=subgroups[x]) ++changed;
                if (subgroups[rule[x]].size()>subgroups[x].size()) ++grown;
            }
            if (grown) ++escaping_rules;
            for (uint16_t a=0;a<n;++a) for (uint16_t b=0;b<n;++b) for (uint16_t c=0;c<n;++c) {
                std::vector<uint16_t> left{a,b,c},right=left;
                const auto act=[&](auto& values,std::size_t i) {
                    const auto result=rule[values[i]*n+values[i+1]];
                    values[i]=result/n; values[i+1]=result%n;
                };
                act(left,0); act(left,1); act(left,0);
                act(right,1); act(right,0); act(right,1);
                if (left!=right) ++braid_failures;
            }
            if (index) out << ",\n";
            out << "{\"id\":" << index << ",\"moved_pairs\":" << moved
                << ",\"subgroup_changed_pairs\":" << changed << ",\"subgroup_grown_pairs\":" << grown
                << ",\"braid_failures_of_512\":" << braid_failures << ",\"table\":[";
            for (std::size_t x=0;x<rule.size();++x) { if (x) out << ','; out << rule[x]; }
            out << "],\"single_seed_runs\":[";
            bool first=true;
            for (uint16_t seed=0;seed<n;++seed) if (seed!=id) for (std::size_t phase=0;phase<2;++phase) {
                std::vector<uint16_t> initial(128,id); initial[64]=seed;
                CellChain chain(fiber,initial,std::vector<uint16_t>(127,id),rule);
                const auto total=chain.total();
                if (!first) out << ',';
                first=false;
                out << "{\"seed_element\":" << seed << ",\"phase\":" << phase << ",\"samples\":[";
                for (std::size_t tick=0;tick<=32;++tick) {
                    if (tick) for (std::size_t i=(tick-1+phase)%2;i+1<128;i+=2)
                        chain.update(i,HurwitzDirection::Forward,false);
                    if (chain.total()!=total) throw std::logic_error("ordered product changed in chain");
                    if (tick%4) continue;
                    std::size_t active=0; long left=0,right=0;
                    for (std::size_t i=0;i<128;++i) if (chain.values()[i]!=id) {
                        if (!active) left=static_cast<long>(i)-64;
                        right=static_cast<long>(i)-64; ++active;
                    }
                    if (tick) out << ',';
                    out << "{\"layer\":" << tick << ",\"active\":" << active << ",\"left\":" << left
                        << ",\"right\":" << right << ",\"subgroup_order\":" << chain.generated_subgroup().size() << '}';
                }
                out << "]}";
            }
            out << "]}";
        }
        out << "\n],\"rule_count\":" << rules.size() << ",\"subgroup_growing_rules\":" << escaping_rules << "}\n";
        if (!out) throw std::runtime_error("failed to write output");
        std::cout << rules.size() << " rules, " << escaping_rules << " grow generated subgroups; wrote " << path << '\n';
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
