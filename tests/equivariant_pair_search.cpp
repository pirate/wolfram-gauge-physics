#include <iostream>
#include "cell_chain.hpp"

using namespace wgphysics::research;
void require(bool ok,const char* msg) { if (!ok) throw std::runtime_error(msg); }

int main() {
    try {
        const FiberGraph edge(2,{{0,1}});
        const AutomorphismTables small(edge.automorphisms());
        PairTable permutation{0,1,2,3};
        std::set<PairTable> brute;
        do {
            try { validate_pair_table(small,permutation); brute.insert(permutation); }
            catch (const std::invalid_argument&) {}
        } while (std::next_permutation(permutation.begin(),permutation.end()));
        const auto minimal=search_pair_involutions(small);
        require(brute==std::set<PairTable>(minimal.begin(),minimal.end()),"C2 search disagrees with brute force");
        std::set<PairTable> unions;
        const auto visited=visit_pair_involutions(small,[&](const auto& table) { unions.insert(table); });
        require(unions==brute && visited==brute.size(),"full involution visitor repeats or omits C2 laws");
        bool capped=false;
        try { visit_pair_involutions(small,[](const auto&) {},1); }
        catch (const std::runtime_error&) { capped=true; }
        require(capped,"capped enumeration silently claimed exhaustion");
        const FiberGraph square(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(square.automorphisms());
        const auto rules=search_pair_involutions(group);
        const BaseGraph base({0,1,2,3,4,5},{{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{3,0}});
        const CellPairPatch patch{{0,1,2,0},{3,4,5,3},{3,0}};
        std::size_t grew=0,checked=0;
        for (const auto& rule : rules) {
            validate_pair_table(group,rule);
            for (std::size_t x=0;x<rule.size();++x) {
                require(rule[rule[x]]==x,"constructed rule not involutive");
                for (uint16_t g=0;g<group.order();++g)
                    require((conjugated_pair(group,x,g)==x)==(conjugated_pair(group,rule[x],g)==rule[x]),
                            "equivariant bijection changed the point stabilizer");
                if (pair_generated_subgroup(group,rule[x]).size()>pair_generated_subgroup(group,x).size()) ++grew;
                FiberBundleConnection state(base,square);
                state.set_transport(2,0,group.elements()[x/group.order()]);
                state.set_transport(5,3,group.elements()[x%group.order()]);
                state.set_transport(3,0,group.elements()[(x+3)%group.order()]);
                const auto result=apply_pair_table_interaction(state,patch,rule);
                const auto restored=apply_pair_table_interaction(result,patch,rule);
                CellChain chain(square,{static_cast<uint16_t>(x/group.order()),static_cast<uint16_t>(x%group.order())},
                                {static_cast<uint16_t>((x+3)%group.order())},rule);
                chain.update(0);
                require(result.holonomy(patch.first_loop)==group.elements()[chain.values()[0]] &&
                        result.holonomy(patch.second_loop)==group.elements()[chain.values()[1]],"compiled/link mismatch");
                const std::map<Vertex,Permutation> frames{{0,group.elements()[(x+1)%group.order()]},
                    {3,group.elements()[(x+2)%group.order()]},{2,group.elements()[(x+4)%group.order()]}};
                const auto framed=apply_pair_table_interaction(state.gauge_transform(frames),patch,rule);
                const auto expected=result.gauge_transform(frames);
                for (const auto [u,v] : base.edges()) {
                    require(restored.edge_transport(u,v)==state.edge_transport(u,v),"link inverse failed");
                    require(framed.edge_transport(u,v)==expected.edge_transport(u,v),"local frame covariance failed");
                }
                ++checked;
            }
        }
        require(grew>0,"search found no escape from word-map subgroup conservation");
        std::cout << rules.size() << " D4 rules, " << checked << " link/oracle cases, " << grew
                  << " subgroup-growing transitions; C2 brute-force check passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
