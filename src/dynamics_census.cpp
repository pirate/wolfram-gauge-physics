#include <filesystem>
#include <fstream>
#include <hypergraph/parallel_evolution.hpp>
#include <hypergraph/pattern.hpp>
#include <iomanip>
#include <iostream>
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
        std::string output_path = "data/d4-cell-dynamics-census.json";
        if (argc == 3 && std::string(argv[1]) == "--output")
            output_path = argv[2];
        else if (argc != 1) {
            throw std::invalid_argument("usage: wgphysics_dynamics_census [--output FILE]");
        }

        hypergraph::Hypergraph graph;
        graph.set_state_canonicalization_mode(hypergraph::StateCanonicalizationMode::None);
        hypergraph::ParallelEvolutionEngine engine(&graph, 1);
        engine.set_explore_from_canonical_states_only(false);
        engine.add_rule(hypergraph::make_rule(0).lhs({0, 1}).rhs({0, 2}).rhs({2, 1}).build());
        engine.evolve({{0, 1}, {1, 2}, {2, 0}}, 2);

        const wgphysics::infragauge::FiberGraph fiber(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        const auto entries = wgphysics::product::census_cell_dynamics(graph, 0, fiber);
        const auto sector_count =
            wgphysics::research::conjugacy_representatives(fiber.automorphisms()).size();

        const auto parent = std::filesystem::path(output_path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("could not open output: " + output_path);
        output << std::setprecision(17);
        output << "{\n  \"schema\":1,\n"
               << "  \"fiber\":\"C4\",\n"
               << "  \"derived_gauge_group\":\"Aut(C4)=D4\",\n"
               << "  \"rule\":\"{x,y}->{x,w},{w,y}\",\n"
               << "  \"initial_base\":\"triangle\",\n"
               << "  \"steps\":2,\n"
               << "  \"curvature_sectors\":" << sector_count << ",\n"
               << "  \"candidate_constraints\":[\"bijection\",\"identity_fixed\","
                  "\"inversion_equivariant\",\"conjugation_equivariant\"],\n"
               << "  \"entries\":[\n";
        for (std::size_t index = 0; index < entries.size(); ++index) {
            const auto& entry = entries[index];
            output << "    {\"rule_index\":" << entry.rule_index << ",\"element_map\":";
            write_values(output, entry.element_map);
            output << ",\"curvature_sector_map\":";
            write_values(output, entry.curvature_sector_map);
            output << ",\"fixed_elements\":" << entry.fixed_elements
                   << ",\"fixed_curvature_sectors\":" << entry.fixed_curvature_sectors
                   << ",\"physical_product_states\":" << entry.physical_product_states
                   << ",\"changed_events\":" << entry.changed_events
                   << ",\"source_events\":" << entry.source_events
                   << ",\"causally_reached_changes\":" << entry.causally_reached_changes
                   << ",\"off_causal_changes\":" << entry.off_causal_changes
                   << ",\"maximum_causal_depth\":" << entry.maximum_causal_depth
                   << ",\"causal_alignment\":" << entry.causal_alignment << '}';
            if (index + 1 != entries.size()) output << ',';
            output << '\n';
        }
        output << "  ]\n}\n";
        std::cout << "wrote " << output_path << " with " << entries.size()
                  << " exact local dynamics candidates\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
