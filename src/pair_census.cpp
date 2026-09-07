#include <filesystem>
#include <fstream>
#include <iostream>
#include "gauge_interactions.hpp"

using namespace wgphysics::research;
using Code = std::pair<std::size_t, std::size_t>;

int main(int argc, char** argv) {
    try {
        std::string path = "out/pair-interactions.json";
        if (argc == 3 && std::string(argv[1]) == "--output") path = argv[2];
        else if (argc != 1) throw std::invalid_argument("usage: wgphysics_pair_census [--output FILE]");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const auto group = fiber.automorphisms();
        const AutomorphismTables tables(group);
        const auto mul = Permutation::compose;
        const auto key = [&](const Permutation& a, const Permutation& b) {
            Code best{group.size(),group.size()};
            for (const auto& g : group)
                best = std::min(best, Code{tables.index_of(mul(g,mul(a,g.inverse()))),
                                         tables.index_of(mul(g,mul(b,g.inverse())))});
            return best;
        };
        const auto sector_multiset = [&](const Permutation& a, const Permutation& b) {
            std::vector<std::vector<Vertex>> sectors{
                canonical_conjugacy_image(group,a), canonical_conjugacy_image(group,b)};
            std::sort(sectors.begin(),sectors.end());
            return sectors;
        };
        const auto parent = std::filesystem::path(path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("could not open output");
        out << "{\n  \"schema\":1,\n  \"fiber\":\"C4\",\n  \"group_order\":" << group.size()
            << ",\n  \"quotient\":\"ordered holonomy pairs under simultaneous Aut(F) conjugation\","
               "\n  \"rules\":[\n";
        bool first = true;
        for (const auto rule : {PairInteraction::Hurwitz,PairInteraction::BoundaryShear}) {
            std::map<Code,Code> quotient;
            std::size_t changed_multisets = 0, braid_failures = 0;
            for (const auto& a : group) for (const auto& b : group) {
                const auto next = pair_targets(a,b,rule);
                if (mul(next.first,next.second) != mul(a,b))
                    throw std::logic_error("boundary holonomy changed");
                const auto before = key(a,b), after = key(next.first,next.second);
                const auto [it, inserted] = quotient.emplace(before,after);
                if (!inserted && it->second != after)
                    throw std::logic_error("update is not well-defined on gauge quotient");
                if (sector_multiset(a,b) != sector_multiset(next.first,next.second))
                    ++changed_multisets;
                for (const auto& c : group) {
                    std::vector<Permutation> left{a,b,c},right=left;
                    const auto act = [&](auto& values,std::size_t i) {
                        auto result = pair_targets(values[i],values[i+1],rule);
                        values[i]=result.first; values[i+1]=result.second;
                    };
                    act(left,0); act(left,1); act(left,0);
                    act(right,1); act(right,0); act(right,1);
                    if (left != right) ++braid_failures;
                }
            }
            std::set<Code> image, visited;
            for (const auto& [from,to] : quotient) image.insert(to);
            if (image.size() != quotient.size()) throw std::logic_error("quotient not reversible");
            std::map<std::size_t,std::size_t> cycles;
            for (const auto& [start,unused] : quotient) {
                if (visited.contains(start)) continue;
                auto at=start; std::size_t length=0;
                do { visited.insert(at); ++length; at=quotient.at(at); } while (at != start);
                ++cycles[length];
            }
            if (!first) out << ",\n";
            first=false;
            out << "    {\"name\":\"" << (rule == PairInteraction::Hurwitz ? "hurwitz" : "boundary_shear")
                << "\",\"raw_pairs\":" << group.size()*group.size()
                << ",\"physical_pairs\":" << quotient.size()
                << ",\"changed_sector_multisets\":" << changed_multisets
                << ",\"braid_triples\":" << group.size()*group.size()*group.size()
                << ",\"braid_relation_failures\":" << braid_failures
                << ",\"quotient_cycle_counts_by_length\":{";
            bool first_cycle=true;
            for (const auto& [length,count] : cycles) {
                if (!first_cycle) out << ',';
                first_cycle=false;
                out << '"' << length << "\":" << count;
            }
            out << "}}";
        }
        out << "\n  ]\n}\n";
        if (!out) throw std::runtime_error("could not finish output");
        std::cout << "wrote " << path << '\n';
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
