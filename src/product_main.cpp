#include <filesystem>
#include <fstream>
#include <hypergraph/parallel_evolution.hpp>
#include <hypergraph/pattern.hpp>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "product_evolution.hpp"

namespace {

template <typename T>
void write_values(std::ostream& output, const std::vector<T>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) output << ',';
        output << values[index];
    }
    output << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        uint32_t steps = 2;
        std::string output_path = "out/product-evolution.json";
        std::optional<std::size_t> dynamics_index;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--help") {
                std::cout << "wgphysics_product_evolve [--steps N] [--output FILE] "
                             "[--cell-dynamics-index N]\n";
                return 0;
            }
            if (index + 1 >= argc) throw std::invalid_argument("missing value after " + argument);
            const std::string value = argv[++index];
            if (argument == "--steps")
                steps = static_cast<uint32_t>(std::stoul(value));
            else if (argument == "--output")
                output_path = value;
            else if (argument == "--cell-dynamics-index") {
                dynamics_index = static_cast<std::size_t>(std::stoull(value));
            } else
                throw std::invalid_argument("unknown argument: " + argument);
        }

        hypergraph::Hypergraph graph;
        graph.set_state_canonicalization_mode(hypergraph::StateCanonicalizationMode::None);
        hypergraph::ParallelEvolutionEngine engine(&graph, 1);
        engine.set_explore_from_canonical_states_only(false);
        engine.add_rule(hypergraph::make_rule(0).lhs({0, 1}).rhs({0, 2}).rhs({2, 1}).build());
        engine.evolve({{0, 1}, {1, 2}, {2, 0}}, steps);

        using namespace wgphysics::infragauge;
        const FiberGraph fiber(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        FiberBundleConnection initial(BaseGraph({0, 1, 2}, {{0, 1}, {1, 2}, {2, 0}}), fiber);
        initial.set_transport(2, 0, Permutation({1, 2, 3, 0}));
        const wgphysics::infragauge::AutomorphismTables tables(fiber.automorphisms());
        const auto dynamics = wgphysics::research::search_local_gauge_dynamics(tables);
        const wgphysics::research::LocalGaugeDynamics* selected_dynamics = nullptr;
        if (dynamics_index) {
            if (*dynamics_index >= dynamics.size()) {
                throw std::out_of_range(
                    "cell dynamics index is outside the exact candidate census");
            }
            selected_dynamics = &dynamics[*dynamics_index];
        }
        const auto product =
            wgphysics::product::evolve_product(graph, 0, initial, selected_dynamics);

        const auto parent = std::filesystem::path(output_path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("could not open output: " + output_path);
        output << std::setprecision(17);
        output << "{\n  \"schema\":3,\n"
               << "  \"semantics\":\"raw engine evolution plus exact joint "
                  "base/gauge quotient\",\n"
               << "  \"steps\":" << steps << ",\n"
               << "  \"cell_dynamics\":";
        if (selected_dynamics) {
            output << "{\"index\":" << *dynamics_index << ",\"element_map\":";
            write_values(output, selected_dynamics->image);
            output << ",\"curvature_sector_map\":";
            write_values(output, wgphysics::research::induced_conjugacy_sector_map(
                                     tables, *selected_dynamics));
            output << '}';
        } else {
            output << "null";
        }
        output << ",\n"
               << "  \"counts\":{\"raw_engine_states\":" << graph.num_published_states()
               << ",\"physical_product_states\":" << product.physical_states.size()
               << ",\"events\":" << product.events.size()
               << ",\"causal_edges\":" << product.causal_edges.size()
               << ",\"curvature_changes\":" << product.causal_curvature.changed_events
               << ",\"source_events\":" << product.causal_curvature.source_events
               << ",\"causally_reached_changes\":"
               << product.causal_curvature.causally_reached_changes
               << ",\"off_causal_changes\":" << product.causal_curvature.off_causal_changes
               << ",\"maximum_causal_depth\":" << product.causal_curvature.maximum_causal_depth
               << ",\"causal_alignment\":" << product.causal_curvature.causal_alignment
               << "},\n  \"physical_states\":[\n";
        for (std::size_t index = 0; index < product.physical_states.size(); ++index) {
            const auto& state = product.physical_states[index];
            output << "    {\"id\":" << index
                   << ",\"base_vertices\":" << state.base().vertices().size()
                   << ",\"base_edges\":" << state.base().edges().size() << ",\"gauge_signature\":";
            write_values(output, state.gauge_invariant_signature());
            output << '}';
            if (index + 1 != product.physical_states.size()) output << ',';
            output << '\n';
        }
        output << "  ],\n  \"raw_state_products\":";
        write_values(output, product.physical_state_by_raw_state);
        output << ",\n  \"event_sectors\":[\n";
        for (std::size_t index = 0; index < product.events.size(); ++index) {
            const auto& event = product.events[index];
            output << "    {\"engine_event\":" << event.engine_event
                   << ",\"raw_input\":" << event.raw_input_state
                   << ",\"raw_output\":" << event.raw_output_state
                   << ",\"physical_input\":" << event.physical_input_state
                   << ",\"physical_output\":" << event.physical_output_state << ",\"subdivision\":["
                   << event.subdivided_from << ',' << event.fresh_vertex << ','
                   << event.subdivided_to << ']'
                   << ",\"gauge_orbit_size\":" << event.gauge_orbit_size
                   << ",\"amplitude_per_labeled_factorization\":"
                   << event.amplitude_per_labeled_factorization << ",\"curvature_changed\":"
                   << (event.curvature_sector_changed ? "true" : "false") << '}';
            if (index + 1 != product.events.size()) output << ',';
            output << '\n';
        }
        output << "  ],\n  \"causal_edges\":[";
        for (std::size_t index = 0; index < product.causal_edges.size(); ++index) {
            if (index) output << ',';
            output << '[' << product.causal_edges[index].first << ','
                   << product.causal_edges[index].second << ']';
        }
        output << "]\n}\n";
        std::cout << "wrote " << output_path << " with " << graph.num_published_states()
                  << " raw states and " << product.physical_states.size()
                  << " physical product states\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
