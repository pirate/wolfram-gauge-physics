#include "product_evolution.hpp"

#include <hypergraph/parallel_evolution.hpp>
#include <hypergraph/pattern.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void write_code(std::ostream& output, const std::vector<wgphysics::infragauge::Vertex>& code) {
    output << '[';
    for (std::size_t index = 0; index < code.size(); ++index) {
        if (index) output << ',';
        output << code[index];
    }
    output << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::string output_path = "data/rule-fiber-census-subdivision.json";
        if (argc == 3 && std::string(argv[1]) == "--output") output_path = argv[2];
        else if (argc != 1) throw std::invalid_argument("usage: wgphysics_product_census [--output FILE]");

        hypergraph::Hypergraph graph;
        graph.set_state_canonicalization_mode(hypergraph::StateCanonicalizationMode::None);
        hypergraph::ParallelEvolutionEngine engine(&graph, 1);
        engine.set_explore_from_canonical_states_only(false);
        engine.add_rule(
            hypergraph::make_rule(0).lhs({0, 1}).rhs({0, 2}).rhs({2, 1}).build());
        engine.evolve({{0, 1}, {1, 2}, {2, 0}}, 2);

        const auto fibers = wgphysics::research::enumerate_fiber_census(4);
        std::vector<wgphysics::product::RuleFiberCensusEntry> entries;
        entries.reserve(fibers.size());
        for (const auto& fiber : fibers) {
            entries.push_back(wgphysics::product::census_subdivision_rule(
                graph, 0, wgphysics::research::fiber_from_canonical_code(fiber.graph_code)));
        }

        const auto parent = std::filesystem::path(output_path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("could not open output: " + output_path);
        output << "{\n  \"schema\":1,\n"
               << "  \"rule\":\"{x,y}->{x,w},{w,y}\",\n"
               << "  \"initial_base\":\"triangle\",\n"
               << "  \"steps\":2,\n"
               << "  \"engine_state_semantics\":\"full raw provenance\",\n"
               << "  \"product_identity\":\"base isomorphism plus local gauge quotient\",\n"
               << "  \"entries\":[\n";
        for (std::size_t index = 0; index < entries.size(); ++index) {
            const auto& entry = entries[index];
            output << "    {\"fiber_graph_code\":";
            write_code(output, entry.fiber_graph_code);
            output << ",\"automorphism_order\":" << entry.automorphism_order
                   << ",\"initial_curvature_sectors\":" << entry.initial_curvature_sectors
                   << ",\"raw_engine_states_per_sector\":"
                   << entry.raw_engine_states_per_sector
                   << ",\"physical_product_states_across_sectors\":"
                   << entry.physical_product_states_across_sectors
                   << ",\"engine_events_across_sectors\":"
                   << entry.engine_events_across_sectors
                   << ",\"causal_edges_across_sectors\":"
                   << entry.causal_edges_across_sectors
                   << ",\"curvature_violations\":" << entry.curvature_violations
                   << ",\"gauge_orbit_norm_verified\":"
                   << (entry.maximum_gauge_orbit_norm_error < 1e-12 ? "true" : "false") << '}';
            if (index + 1 != entries.size()) output << ',';
            output << '\n';
        }
        output << "  ]\n}\n";
        std::cout << "wrote " << output_path << " with " << entries.size()
                  << " rule-by-fiber product rows\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
