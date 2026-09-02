#include "product_evolution.hpp"

#include <hypergraph/parallel_evolution.hpp>
#include <hypergraph/pattern.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void evolve_subdivisions(hypergraph::Hypergraph& graph, uint32_t steps, uint32_t threads = 1) {
    graph.set_state_canonicalization_mode(hypergraph::StateCanonicalizationMode::None);
    hypergraph::ParallelEvolutionEngine engine(&graph, threads);
    engine.set_explore_from_canonical_states_only(false);
    engine.add_rule(hypergraph::make_rule(0).lhs({0, 1}).rhs({0, 2}).rhs({2, 1}).build());
    engine.evolve({{0, 1}, {1, 2}, {2, 0}}, steps);
}

}  // namespace

int main() {
    try {
        using namespace wgphysics::infragauge;
        using namespace wgphysics::product;
        using wgphysics::research::canonical_connection_state;

        hypergraph::Hypergraph graph;
        evolve_subdivisions(graph, 2);
        const FiberGraph fiber(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        const BaseGraph base({0, 1, 2}, {{0, 1}, {1, 2}, {2, 0}});
        FiberBundleConnection initial(base, fiber);
        initial.set_transport(2, 0, Permutation({1, 2, 3, 0}));

        const auto result = evolve_product(graph, 0, initial);
        require(result.events.size() == graph.num_published_events(),
                "not every real engine event received a product transition");
        require(result.physical_state_by_raw_state.size() == graph.num_published_states(),
                "not every raw engine state received a product identity slot");
        require(std::all_of(
                    result.physical_state_by_raw_state.begin(),
                    result.physical_state_by_raw_state.end(),
                    [](const auto id) { return id != std::numeric_limits<std::size_t>::max(); }),
                "a reachable raw engine state has no product identity");
        require(result.physical_states.size() == 3,
                "triangle subdivisions should quotient to one physical state at each depth");
        require(!result.causal_edges.empty(), "actual engine run produced no causal edges");
        require(result.causal_curvature.changed_events == 0 &&
                    result.causal_curvature.causal_alignment == 1.0,
                "transport-preserving subdivision changed its curvature sector");
        require(std::all_of(result.events.begin(), result.events.end(), [](const auto& event) {
                    return !event.curvature_sector_changed && event.gauge_orbit_size == 8 &&
                           std::abs(event.gauge_orbit_size *
                                        event.amplitude_per_labeled_factorization *
                                        event.amplitude_per_labeled_factorization -
                                    1.0) < 1e-12;
                }),
                "an engine event has a spurious curvature change or non-isometric gauge orbit");

        const auto census = census_subdivision_rule(graph, 0, fiber);
        require(census.initial_curvature_sectors == 5,
                "D4 product census has the wrong curvature-sector count");
        require(census.physical_product_states_across_sectors == 15,
                "D4 product census has the wrong joint closure size");
        require(census.curvature_violations == 0 &&
                    census.maximum_gauge_orbit_norm_error < 1e-12,
                "D4 product census violated curvature preservation or orbit norm");

        const auto framed = initial.gauge_transform({
            {0, Permutation({0, 3, 2, 1})},
            {1, Permutation({1, 2, 3, 0})}});
        const auto framed_result = evolve_product(graph, 0, framed);
        require(framed_result.physical_states.size() == result.physical_states.size(),
                "local frame choice changed the product closure size");
        for (std::size_t id = 0; id < result.physical_states.size(); ++id) {
            require(canonical_connection_state(result.physical_states[id]) ==
                        canonical_connection_state(framed_result.physical_states[id]),
                    "local frame choice changed a product-state identity");
        }

        hypergraph::Hypergraph parallel_graph;
        evolve_subdivisions(parallel_graph, 2, 4);
        const auto parallel_result = evolve_product(parallel_graph, 0, initial);
        require(parallel_graph.num_published_states() == graph.num_published_states() &&
                    parallel_graph.num_published_events() == graph.num_published_events(),
                "raw engine closure depends on worker count");
        require(parallel_result.physical_states.size() == result.physical_states.size() &&
                    parallel_result.events.size() == result.events.size() &&
                    parallel_result.causal_edges.size() == result.causal_edges.size(),
                "product closure depends on worker count");
        std::set<std::vector<Vertex>> serial_identities;
        std::set<std::vector<Vertex>> parallel_identities;
        for (const auto& state : result.physical_states) {
            serial_identities.insert(canonical_connection_state(state));
        }
        for (const auto& state : parallel_result.physical_states) {
            parallel_identities.insert(canonical_connection_state(state));
        }
        require(serial_identities == parallel_identities,
                "product state identities depend on worker scheduling");

        hypergraph::Hypergraph unsupported;
        unsupported.set_state_canonicalization_mode(hypergraph::StateCanonicalizationMode::None);
        hypergraph::ParallelEvolutionEngine unsupported_engine(&unsupported, 1);
        unsupported_engine.add_rule(
            hypergraph::make_rule(0)
                .lhs({0, 1})
                .rhs({0, 2})
                .rhs({2, 3})
                .rhs({3, 1})
                .build());
        unsupported_engine.evolve({{0, 1}, {1, 2}, {2, 0}}, 1);
        bool rejected = false;
        try {
            (void)evolve_product(unsupported, 0, initial);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "unsupported rule morphology was silently assigned a connection update");

        std::cout << "Official-engine gauge product evolution: PASS\n"
                  << "raw engine states: " << graph.num_published_states() << '\n'
                  << "physical product states: " << result.physical_states.size() << '\n'
                  << "engine events with sectors: " << result.events.size() << '\n'
                  << "engine causal edges: " << result.causal_edges.size() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Official-engine gauge product evolution: FAIL: " << error.what() << '\n';
        return 1;
    }
}
