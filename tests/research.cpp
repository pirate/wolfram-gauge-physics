#include "research.hpp"

#include <cmath>
#include <iostream>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    try {
        using namespace wgphysics::infragauge;
        using namespace wgphysics::research;

        // #1: a wheel's center has a C4 link, so its local symmetry is inferred
        // as D4 without supplying a group or calling the link a gauge fiber.
        const BaseGraph wheel(
            {0, 1, 2, 3, 4},
            {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 3}, {3, 4}, {4, 1}});
        const auto inferred = infer_link_fibers(wheel);
        const auto c4_class = std::find_if(inferred.begin(), inferred.end(), [](const auto& item) {
            return item.representative.vertex_count() == 4 && item.automorphism_order == 8;
        });
        require(c4_class != inferred.end(), "failed to infer the C4 link fiber");
        require(c4_class->centers == std::vector<Vertex>{0}, "C4 link was assigned to the wrong center");

        const FiberGraph square(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        const auto group = square.automorphisms();
        const Permutation quarter_turn({1, 2, 3, 0});
        const Permutation reflection({0, 3, 2, 1});
        const BaseGraph triangle({0, 1, 2}, {{0, 1}, {1, 2}, {2, 0}});
        FiberBundleConnection curved(triangle, square);
        curved.set_transport(2, 0, quarter_turn);

        // #2: connection state participates in exact multiway identity. All
        // eight fresh-frame extensions register as one physical child.
        GaugeSubdivisionEvolution evolution(curved);
        const auto outputs = evolution.subdivide(0, 0, 1, 3);
        require(outputs.size() == 1, "gauge copies became distinct rewrite children");
        require(evolution.state_count() == 2, "unexpected physical state count after subdivision");
        require(evolution.events().size() == 1 && evolution.events().front().raw_extensions == 8,
                "rewrite event lost its raw-to-physical orbit count");
        const auto physical_identity =
            canonical_connection_state(evolution.state(outputs.front()));
        for (const auto& raw_extension : curved.subdivide_edge_extensions(0, 1, 3)) {
            require(canonical_connection_state(raw_extension) == physical_identity,
                    "physical-only subdivision disagrees with an exhaustive orbit member");
        }
        const auto framed = curved.gauge_transform({{0, reflection}, {1, quarter_turn}});
        require(canonical_connection_state(curved) == canonical_connection_state(framed),
                "joint state identity depends on local fiber frames");
        FiberBundleConnection relabeled(
            BaseGraph({10, 11, 12}, {{10, 11}, {11, 12}, {12, 10}}), square);
        relabeled.set_transport(12, 10, quarter_turn);
        require(canonical_connection_state(curved) == canonical_connection_state(relabeled),
                "joint state identity depends on base labels");

        // #3: non-isomorphic fibers and partial horizontal lifts. A triangle
        // embeds only two selected adjacent vertices into a square; an absent
        // image is an explicit failed lift rather than an invented transport.
        const BaseGraph path({0, 1, 2}, {{0, 1}, {1, 2}});
        const FiberGraph triangle_fiber(3, {{0, 1}, {1, 2}, {2, 0}});
        const FiberGraph path_fiber(3, {{0, 1}, {1, 2}});
        VariableFiberBundle variable(path);
        variable.set_fiber(0, triangle_fiber);
        variable.set_fiber(1, square);
        variable.set_fiber(2, path_fiber);
        variable.set_partial_transport(
            0, 1, PartialFiberMap(3, 4, {Vertex{0}, Vertex{1}, std::nullopt}));
        variable.set_partial_transport(
            1, 2, PartialFiberMap(4, 3, {Vertex{0}, Vertex{1}, Vertex{2}, std::nullopt}));
        require(variable.lift({0, 1, 2}, 0) == std::optional<Vertex>{0},
                "valid partial horizontal lift failed");
        require(!variable.lift({0, 1, 2}, 2), "undefined partial lift was silently completed");

        // #4: enumerate the complete finite family of reversible local maps
        // that preserve identity, inversion, and conjugation equivariance.
        const AutomorphismTables tables(group);
        const auto dynamics = search_local_gauge_dynamics(tables);
        require(!dynamics.empty(), "no gauge-equivariant reversible dynamics found");
        const auto identity_rule = std::find_if(dynamics.begin(), dynamics.end(), [&](const auto& rule) {
            for (uint16_t value = 0; value < tables.order(); ++value) {
                if (rule.image[value] != value) return false;
            }
            return true;
        });
        require(identity_rule != dynamics.end(), "dynamics census omitted the identity update");
        require(std::any_of(dynamics.begin(), dynamics.end(), [&](const auto& rule) {
                    return rule.fixed_elements < tables.order();
                }),
                "dynamics census found no nontrivial reversible update");

        // #5: score curvature changes against causal edges. Event 1 is reached
        // from declared source 0; event 2 changes independently and is exposed.
        FiberBundleConnection flat(triangle, square);
        FiberBundleConnection reflected(triangle, square);
        reflected.set_transport(2, 0, reflection);
        const std::vector<std::vector<Vertex>> loops{{0, 1, 2, 0}};
        const auto flat_sector = curvature_sector(flat, loops);
        const auto curved_sector = curvature_sector(curved, loops);
        const auto reflected_sector = curvature_sector(reflected, loops);
        require(flat_sector != curved_sector && curved_sector != reflected_sector,
                "curvature sectors collapsed before causal analysis");
        const auto report = analyze_causal_curvature(
            {{0, flat_sector, curved_sector},
             {1, curved_sector, reflected_sector},
             {2, flat_sector, reflected_sector}},
            {{0, 1}},
            {0});
        require(report.changed_events == 3 && report.source_events == 1,
                "causal curvature event counts are wrong");
        require(report.causally_reached_changes == 1 && report.off_causal_changes == 1,
                "causal curvature reachability is wrong");
        require(report.maximum_causal_depth == 1 && std::abs(report.causal_alignment - 0.5) < 1e-12,
                "causal curvature score is wrong");

        // #6, deliberately: diagnose before compressing. The subdivision
        // isometry has a flat rank-|D4| Schmidt spectrum, so rank four must
        // discard exactly half the norm. There is no lossless low-rank shortcut.
        const auto orbit = SubdivisionWavefunction::from_total_transport(
            group, Permutation::identity(square.vertex_count()));
        const auto spectrum = orbit.schmidt_spectrum();
        require(spectrum.size() == group.size(), "subdivision Schmidt rank is wrong");
        require(std::all_of(spectrum.begin(), spectrum.end(), [&](const auto value) {
                    return std::abs(value - 1.0 / std::sqrt(8.0)) < 1e-12;
                }),
                "subdivision Schmidt spectrum is not flat");
        require(std::abs(orbit.best_rank_discarded_norm(4) - 0.5) < 1e-12,
                "rank-four discarded norm is wrong");

        // #8: all unlabeled simple graphs through four vertices are included
        // exactly once (1 + 2 + 4 + 11 = 18), with finite derived observables.
        const auto census = enumerate_fiber_census(4);
        require(census.size() == 18, "small-fiber census has the wrong graph count");
        require(std::any_of(census.begin(), census.end(), [](const auto& entry) {
                    return entry.vertices == 4 && entry.automorphism_order == 8 &&
                           entry.nonabelian && entry.reversible_equivariant_dynamics == 8;
                }),
                "small-fiber census omitted the D4 case");
        require(std::all_of(census.begin(), census.end(), [](const auto& entry) {
                    return entry.subdivision_physical_children == 1 &&
                           entry.exact_subdivision_schmidt_rank == entry.automorphism_order;
                }),
                "census rewrite or compression invariant is inconsistent");

        std::cout << "Research milestones 1-5: PASS\n"
                  << "inferred fiber classes: " << inferred.size() << '\n'
                  << "D4 reversible equivariant dynamics: " << dynamics.size() << '\n'
                  << "causal alignment fixture: " << report.causal_alignment << '\n'
                  << "subdivision Schmidt rank: " << spectrum.size() << '\n'
                  << "unlabeled fibers through n=4: " << census.size() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Research milestones 1-5: FAIL: " << error.what() << '\n';
        return 1;
    }
}
