#include <iostream>
#include "gauge_interactions.hpp"

using namespace wgphysics::research;

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

bool same(const FiberBundleConnection& a, const FiberBundleConnection& b) {
    for (const auto [u, v] : a.base().edges())
        if (a.edge_transport(u, v) != b.edge_transport(u, v)) return false;
    return true;
}

int main() {
    try {
        const FiberGraph square(4, {{0,1},{1,2},{2,3},{3,0}});
        const auto group = square.automorphisms();
        const BaseGraph base({0,1,2,3,4,5,6,7},
                             {{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{3,6},{6,0},{6,7}});
        const CellPairPatch patch{{0,1,2,0},{3,4,5,3},{3,6,0}};
        std::size_t covariance_checks = 0, braid_checks = 0;
        for (const auto rule : {PairInteraction::Hurwitz, PairInteraction::BoundaryShear}) {
            for (const auto& a : group) for (const auto& b : group) for (const auto& t : group) {
                FiberBundleConnection state(base, square);
                state.set_transport(2,0,a);
                state.set_transport(5,3,b);
                state.set_transport(3,6,t);
                const auto next = apply_pair_interaction(state, patch, rule);
                require(same(state, apply_pair_interaction(next, patch, rule,
                             HurwitzDirection::Inverse)), "link-level inverse failed");
                const auto support = interaction_support(state, patch);
                for (const auto [u,v] : base.edges()) if (!support.writes.contains({u,v}))
                    require(next.edge_transport(u,v) == state.edge_transport(u,v),
                            "changed an off-support link");
                for (const auto& g : group) for (const auto& h : group) {
                    const std::map<Vertex, Permutation> frames{
                        {0,g},{3,h},{1,t},{2,b},{4,a},{5,t},{6,b},{7,h}};
                    require(same(apply_pair_interaction(state.gauge_transform(frames), patch, rule),
                                 next.gauge_transform(frames)), "local frame covariance failed");
                    ++covariance_checks;
                }
                // Orientation reversal exchanges cells as well as inverting both holonomies.
                const auto reversed = pair_targets(b.inverse(), a.inverse(), rule);
                const auto inverse = pair_targets(a,b,rule,HurwitzDirection::Inverse);
                require(reversed.second.inverse() == inverse.first &&
                        reversed.first.inverse() == inverse.second,
                        "orientation/time reversal identity failed");
            }
        }
        // The existing primitive is a braid action on based holonomies.
        for (const auto& a : group) for (const auto& b : group) for (const auto& c : group) {
            auto left = std::vector<Permutation>{a,b,c}, right = left;
            const auto act = [](auto& values, std::size_t i) {
                const auto pair = pair_targets(values[i],values[i+1],PairInteraction::Hurwitz);
                values[i] = pair.first; values[i+1] = pair.second;
            };
            act(left,0); act(left,1); act(left,0);
            act(right,1); act(right,0); act(right,1);
            require(left == right, "Hurwitz braid relation failed");
            ++braid_checks;
        }
        const Permutation rotation({1,2,3,0}), identity = Permutation::identity(4);
        const BaseGraph three({0,1,2,3,4,5,6,7,8},
            {{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{6,7},{7,8},{8,6},{3,0},{6,3}});
        const CellPairPatch adjacent_first{{0,1,2,0},{3,4,5,3},{3,0}};
        const CellPairPatch adjacent_second{{3,4,5,3},{6,7,8,6},{6,3}};
        std::size_t link_braid_checks = 0;
        for (const auto rule : {PairInteraction::Hurwitz,PairInteraction::BoundaryShear})
            for (const auto& a : group) for (const auto& b : group) for (const auto& c : group) {
                FiberBundleConnection start(three,square);
                start.set_transport(2,0,a); start.set_transport(5,3,b); start.set_transport(8,6,c);
                start.set_transport(3,0,rotation); start.set_transport(6,3,group.back());
                auto left=start,right=start;
                for (const auto* spec : {&adjacent_first,&adjacent_second,&adjacent_first})
                    left=apply_pair_interaction(left,*spec,rule);
                for (const auto* spec : {&adjacent_second,&adjacent_first,&adjacent_second})
                    right=apply_pair_interaction(right,*spec,rule);
                require(same(left,right),"transported link-level braid identity failed");
                ++link_braid_checks;
            }
        const auto shear = pair_targets(rotation, identity, PairInteraction::BoundaryShear);
        require(shear.first == Permutation::compose(rotation, rotation) &&
                shear.second == rotation.inverse(), "shear did not change curvature sectors");
        require(canonical_conjugacy_image(group, shear.first) !=
                canonical_conjugacy_image(group, rotation) && shear.first != identity,
                "claimed sector change is only a frame change");
        // Disjoint write/read supports commute even with a shared read-only connector.
        const BaseGraph four({0,1,2,3,4,5,6,7,8,9},
            {{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{3,0},
             {0,6},{6,7},{7,0},{3,8},{8,9},{9,3}});
        const CellPairPatch p{{0,1,2,0},{3,4,5,3},{3,0}};
        const CellPairPatch q{{0,6,7,0},{3,8,9,3},{3,0}};
        FiberBundleConnection initial(four,square);
        initial.set_transport(2,0,rotation); initial.set_transport(9,3,group.back());
        require(independent_supports(interaction_support(initial,p),interaction_support(initial,q)),
                "shared reads wrongly conflict");
        require(!independent_supports(interaction_support(initial,p),interaction_support(initial,p)),
                "overlapping writes wrongly independent");
        const auto apply = [&](const auto& state, const auto& spec) {
            return apply_pair_interaction(state,spec,PairInteraction::BoundaryShear);
        };
        require(same(apply(apply(initial,p),q),apply(apply(initial,q),p)),
                "support-independent updates do not commute");
        // Caller-selected paths and exclusive closing links are real preconditions.
        for (const auto& bad : std::vector<CellPairPatch>{
                {{0,1,2,0},{3,4,5,3},{0,6,3}},
                {{0,1,2,0},{3,4,5,3},{3,6,0,2,0}},
                {{0,1,2,0,1,2,0},{3,4,5,3},{3,6,0}}}) {
            bool rejected = false;
            try { (void)apply_pair_interaction(FiberBundleConnection(base,square),bad,
                                              PairInteraction::BoundaryShear); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected,"invalid patch accepted");
        }
        std::cout << covariance_checks << " link-level gauge covariance checks; "
                  << braid_checks << " holonomy braid identities; " << link_braid_checks
                  << " transported link braid identities; inverse, support, sector tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
