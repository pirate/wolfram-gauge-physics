#include "research.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct Config {
    std::string output = "data/fiber-census-n4.json";
    wgphysics::infragauge::Vertex maximum_vertices = 4;
};

Config parse_args(int argc, char** argv) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help") {
            std::cout << "wgphysics_census [--max-vertices N] [--output FILE]\n";
            std::exit(0);
        }
        if (i + 1 >= argc) throw std::invalid_argument("missing value after " + argument);
        const std::string value = argv[++i];
        if (argument == "--max-vertices") {
            config.maximum_vertices = static_cast<wgphysics::infragauge::Vertex>(std::stoul(value));
        } else if (argument == "--output") {
            config.output = value;
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    return config;
}

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
        const auto config = parse_args(argc, argv);
        const auto census = wgphysics::research::enumerate_fiber_census(config.maximum_vertices);
        const auto parent = std::filesystem::path(config.output).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream output(config.output);
        if (!output) throw std::runtime_error("could not open output: " + config.output);
        output << "{\n  \"schema\":1,\n"
               << "  \"scope\":\"all unlabeled simple graphs with 1 through "
               << config.maximum_vertices << " vertices\",\n"
               << "  \"dynamics_search_max_group_order\":8,\n"
               << "  \"entries\":[\n";
        for (std::size_t index = 0; index < census.size(); ++index) {
            const auto& entry = census[index];
            output << "    {\"vertices\":" << entry.vertices
                   << ",\"edges\":" << entry.edges << ",\"graph_code\":";
            write_code(output, entry.graph_code);
            output << ",\"automorphism_order\":" << entry.automorphism_order
                   << ",\"nonabelian\":" << (entry.nonabelian ? "true" : "false")
                   << ",\"triangle_curvature_sectors\":" << entry.conjugacy_classes
                   << ",\"reversible_equivariant_dynamics\":";
            if (entry.reversible_equivariant_dynamics) {
                output << *entry.reversible_equivariant_dynamics;
            } else {
                output << "null";
            }
            output << ",\"subdivision_raw_extensions\":" << entry.subdivision_raw_extensions
                   << ",\"subdivision_physical_children\":"
                   << entry.subdivision_physical_children
                   << ",\"exact_subdivision_schmidt_rank\":"
                   << entry.exact_subdivision_schmidt_rank
                   << ",\"half_rank_discarded_norm\":"
                   << entry.half_rank_discarded_norm << '}';
            if (index + 1 != census.size()) output << ',';
            output << '\n';
        }
        output << "  ]\n}\n";
        std::cout << "wrote " << config.output << " with " << census.size()
                  << " unlabeled fiber graphs\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
