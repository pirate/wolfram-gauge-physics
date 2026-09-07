#include <hypergraph/hypergraph.hpp>
#include <hypergraph/parallel_evolution.hpp>
#include <hypergraph/pattern.hpp>

#include "observables.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace hypergraph;

namespace {

using EdgeSpec = std::vector<uint32_t>;
using EdgeList = std::vector<EdgeSpec>;

struct Config {
    std::string rule = "0,1;0,2->0,2;0,3;1,3;2,3";
    std::string init = "1,2;1,3";
    std::string output = "out/evolution.json";
    uint32_t steps = 3;
    uint32_t threads = 1;
    uint32_t probe_radius = 8;
    uint32_t diffusion_steps = 16;
    uint32_t probe_sources = 64;
};

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    std::stringstream stream(text);
    std::string part;
    while (std::getline(stream, part, delimiter)) {
        if (!part.empty()) parts.push_back(part);
    }
    return parts;
}

EdgeList parse_edges(const std::string& text) {
    EdgeList edges;
    for (const auto& edge_text : split(text, ';')) {
        EdgeSpec edge;
        for (const auto& vertex : split(edge_text, ',')) {
            size_t consumed = 0;
            const auto value = std::stoul(vertex, &consumed);
            if (consumed != vertex.size()) {
                throw std::invalid_argument("invalid vertex: " + vertex);
            }
            edge.push_back(static_cast<uint32_t>(value));
        }
        if (edge.empty()) throw std::invalid_argument("empty hyperedge");
        edges.push_back(std::move(edge));
    }
    if (edges.empty()) throw std::invalid_argument("empty edge list");
    return edges;
}

std::pair<EdgeList, EdgeList> parse_rule(const std::string& text) {
    const auto arrow = text.find("->");
    if (arrow == std::string::npos || text.find("->", arrow + 2) != std::string::npos) {
        throw std::invalid_argument("rule must contain exactly one ->");
    }
    return {parse_edges(text.substr(0, arrow)), parse_edges(text.substr(arrow + 2))};
}

uint64_t mix(uint64_t hash, uint64_t value) {
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
    return hash;
}

struct Snapshot {
    StateId id{};
    uint32_t step{};
    EdgeList edges;
    std::vector<uint32_t> vertices;
    std::unordered_map<uint32_t, std::vector<uint32_t>> adjacency;
    std::unordered_map<uint32_t, uint32_t> incident_edges;
};

Snapshot snapshot(const Hypergraph& graph, StateId state_id) {
    const auto& state = graph.get_state(state_id);
    Snapshot out{};
    out.id = state_id;
    out.step = state.step;
    std::set<uint32_t> vertices;
    state.edges.for_each([&](EdgeId edge_id) {
        const auto& edge = graph.get_edge(edge_id);
        EdgeSpec values(edge.vertices, edge.vertices + edge.arity);
        out.edges.push_back(values);
        for (const auto vertex : values) {
            vertices.insert(vertex);
            ++out.incident_edges[vertex];
            out.adjacency[vertex];
        }
        for (size_t i = 0; i < values.size(); ++i) {
            for (size_t j = i + 1; j < values.size(); ++j) {
                if (values[i] == values[j]) continue;
                out.adjacency[values[i]].push_back(values[j]);
                out.adjacency[values[j]].push_back(values[i]);
            }
        }
    });
    out.vertices.assign(vertices.begin(), vertices.end());
    for (auto& [_, neighbors] : out.adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
    std::sort(out.edges.begin(), out.edges.end());
    return out;
}

uint32_t component_count(const Snapshot& state) {
    std::set<uint32_t> unseen(state.vertices.begin(), state.vertices.end());
    uint32_t components = 0;
    while (!unseen.empty()) {
        ++components;
        std::queue<uint32_t> frontier;
        frontier.push(*unseen.begin());
        unseen.erase(unseen.begin());
        while (!frontier.empty()) {
            const auto vertex = frontier.front();
            frontier.pop();
            const auto found = state.adjacency.find(vertex);
            if (found == state.adjacency.end()) continue;
            for (const auto neighbor : found->second) {
                if (unseen.erase(neighbor)) frontier.push(neighbor);
            }
        }
    }
    return components;
}

uint32_t diameter(const Snapshot& state) {
    uint32_t maximum = 0;
    for (const auto source : state.vertices) {
        std::unordered_map<uint32_t, uint32_t> distance{{source, 0}};
        std::queue<uint32_t> frontier;
        frontier.push(source);
        while (!frontier.empty()) {
            const auto vertex = frontier.front();
            frontier.pop();
            for (const auto neighbor : state.adjacency.at(vertex)) {
                if (distance.contains(neighbor)) continue;
                distance[neighbor] = distance[vertex] + 1;
                maximum = std::max(maximum, distance[neighbor]);
                frontier.push(neighbor);
            }
        }
    }
    return maximum;
}

std::map<uint64_t, uint32_t> local_fingerprints(const Snapshot& state) {
    std::unordered_map<uint32_t, uint64_t> colors;
    for (const auto vertex : state.vertices) {
        uint64_t hash = 0xcbf29ce484222325ULL;
        hash = mix(hash, state.adjacency.at(vertex).size());
        hash = mix(hash, state.incident_edges.at(vertex));
        colors[vertex] = hash;
    }
    for (int round = 0; round < 2; ++round) {
        auto next = colors;
        for (const auto vertex : state.vertices) {
            std::vector<uint64_t> neighbors;
            for (const auto neighbor : state.adjacency.at(vertex)) {
                neighbors.push_back(colors.at(neighbor));
            }
            std::sort(neighbors.begin(), neighbors.end());
            uint64_t hash = mix(colors.at(vertex), round + 1);
            for (const auto color : neighbors) hash = mix(hash, color);
            next[vertex] = hash;
        }
        colors = std::move(next);
    }
    std::map<uint64_t, uint32_t> histogram;
    for (const auto& [_, color] : colors) ++histogram[color];
    return histogram;
}

wgphysics::infragauge::BaseGraph projected_base_graph(const Snapshot& state) {
    using wgphysics::infragauge::BaseGraph;
    using wgphysics::infragauge::Edge;
    using wgphysics::infragauge::Vertex;
    std::set<Edge> edges;
    for (const auto& [vertex, neighbors] : state.adjacency) {
        for (const auto neighbor : neighbors) {
            if (vertex != neighbor) edges.insert(wgphysics::infragauge::canonical_edge(vertex, neighbor));
        }
    }
    return BaseGraph(std::vector<Vertex>(state.vertices.begin(), state.vertices.end()),
                     std::vector<Edge>(edges.begin(), edges.end()));
}

void json_doubles(std::ostream& out, const std::vector<double>& values) {
    out << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) out << ',';
        out << std::setprecision(9) << values[index];
    }
    out << ']';
}

void json_optional_doubles(std::ostream& out,
                           const std::vector<std::optional<double>>& values) {
    out << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) out << ',';
        if (values[index]) out << std::setprecision(9) << *values[index];
        else out << "null";
    }
    out << ']';
}

void json_plateau(std::ostream& out,
                  const std::optional<wgphysics::observables::ScalePlateau>& plateau) {
    if (!plateau) {
        out << "null";
        return;
    }
    out << "{\"first_scale\":" << plateau->first_scale
        << ",\"last_scale\":" << plateau->last_scale
        << ",\"mean\":" << std::setprecision(9) << plateau->mean
        << ",\"relative_span\":" << plateau->relative_span << '}';
}

void json_edges(std::ostream& out, const EdgeList& edges) {
    out << '[';
    for (size_t i = 0; i < edges.size(); ++i) {
        if (i) out << ',';
        out << '[';
        for (size_t j = 0; j < edges[i].size(); ++j) {
            if (j) out << ',';
            out << edges[i][j];
        }
        out << ']';
    }
    out << ']';
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            std::cout << "wgphysics_evolve [--rule LHS->RHS] [--init EDGES] [--steps N] "
                         "[--threads N] [--output FILE] [--probe-radius N] "
                         "[--diffusion-steps N] [--probe-sources N]\n\n"
                         "Edges use comma-separated vertices and semicolon-separated hyperedges.\n"
                         "Example: --rule '0,1;0,2->0,2;0,3;1,3;2,3' --init '1,2;1,3'\n";
            std::exit(0);
        }
        if (i + 1 >= argc) throw std::invalid_argument("missing value after " + arg);
        const std::string value = argv[++i];
        if (arg == "--rule") config.rule = value;
        else if (arg == "--init") config.init = value;
        else if (arg == "--steps") config.steps = std::stoul(value);
        else if (arg == "--threads") config.threads = std::stoul(value);
        else if (arg == "--output") config.output = value;
        else if (arg == "--probe-radius") config.probe_radius = std::stoul(value);
        else if (arg == "--diffusion-steps") config.diffusion_steps = std::stoul(value);
        else if (arg == "--probe-sources") config.probe_sources = std::stoul(value);
        else throw std::invalid_argument("unknown argument: " + arg);
    }
    if (config.threads == 0) throw std::invalid_argument("threads must be positive");
    if (config.probe_sources == 0) throw std::invalid_argument("probe sources must be positive");
    return config;
}

void write_result(const Config& config, const Hypergraph& graph, std::ostream& out) {
    const auto causal = graph.causal_graph().get_causal_edges();
    const auto branchial = graph.causal_graph().get_branchial_edges();
    std::set<StateId> canonical_states;
    for (StateId id = 0; id < graph.num_published_states(); ++id) {
        const auto& state = graph.get_state(id);
        canonical_states.insert(state.canonical_id == INVALID_ID ? id : state.canonical_id);
    }
    out << "{\n  \"schema\":2,\n  \"semantics\":\"unlabeled ordered hypergraph rewriting with intrinsic probes\",\n";
    out << "  \"rule\":\"" << config.rule << "\",\n  \"initial_state\":\"" << config.init << "\",\n";
    out << "  \"steps\":" << config.steps << ",\n";
    out << "  \"counts\":{\"raw_states\":" << graph.num_published_states()
        << ",\"canonical_states\":" << canonical_states.size()
        << ",\"events\":" << graph.num_published_events() << ",\"causal_edges\":" << causal.size()
        << ",\"branchial_edges\":" << branchial.size() << "},\n  \"states\":[\n";
    bool first_state = true;
    for (StateId id = 0; id < graph.num_published_states(); ++id) {
        const auto state = snapshot(graph, id);
        if (!first_state) out << ",\n";
        first_state = false;
        const auto fingerprints = local_fingerprints(state);
        const auto intrinsic = wgphysics::observables::probe_intrinsic_geometry(
            projected_base_graph(state),
            {config.probe_radius, config.diffusion_steps, config.probe_sources});
        double mean_degree = 0.0;
        for (const auto vertex : state.vertices) mean_degree += state.adjacency.at(vertex).size();
        if (!state.vertices.empty()) mean_degree /= state.vertices.size();
        const auto& raw_state = graph.get_state(id);
        const auto canonical_id = raw_state.canonical_id == INVALID_ID ? id : raw_state.canonical_id;
        out << "    {\"id\":" << id << ",\"canonical_id\":" << canonical_id
            << ",\"step\":" << state.step << ",\"edges\":";
        json_edges(out, state.edges);
        out << ",\"observables\":{\"vertices\":" << state.vertices.size()
            << ",\"edges\":" << state.edges.size() << ",\"components\":" << component_count(state)
            << ",\"diameter\":" << diameter(state) << ",\"mean_degree\":"
            << std::fixed << std::setprecision(6) << mean_degree << "},\"local_fingerprints\":{";
        bool first = true;
        for (const auto& [hash, count] : fingerprints) {
            if (!first) out << ',';
            first = false;
            std::ostringstream key;
            key << std::hex << hash;
            out << '\"' << key.str() << "\":" << count;
        }
        out << "},\"intrinsic_geometry\":{\"projection\":\"hypergraph_2_section\""
            << ",\"sampled_sources\":" << intrinsic.sampled_sources
            << ",\"mean_shell_volume\":";
        json_doubles(out, intrinsic.mean_shell_volume);
        out << ",\"mean_ball_volume\":";
        json_doubles(out, intrinsic.mean_ball_volume);
        out << ",\"ball_volume_cv\":";
        json_doubles(out, intrinsic.ball_volume_coefficient_of_variation);
        out << ",\"volume_dimension\":";
        json_optional_doubles(out, intrinsic.volume_dimension);
        out << ",\"lazy_return_probability\":";
        json_doubles(out, intrinsic.lazy_return_probability);
        out << ",\"spectral_dimension\":";
        json_optional_doubles(out, intrinsic.spectral_dimension);
        out << ",\"volume_dimension_plateau\":";
        json_plateau(out, intrinsic.volume_dimension_plateau);
        out << ",\"spectral_dimension_plateau\":";
        json_plateau(out, intrinsic.spectral_dimension_plateau);
        out << "}}";
    }
    out << "\n  ],\n  \"events\":[";
    bool first = true;
    for (EventId id = 0; id < graph.num_published_events(); ++id) {
        const auto& event = graph.get_event(id);
        if (event.id == INVALID_ID) continue;
        if (!first) out << ',';
        first = false;
        out << "{\"id\":" << id << ",\"rule\":" << event.rule_index
            << ",\"input\":" << event.input_state << ",\"output\":" << event.output_state << '}';
    }
    out << "],\n  \"causal_edges\":[";
    for (size_t i = 0; i < causal.size(); ++i) {
        if (i) out << ',';
        out << '[' << causal[i].producer << ',' << causal[i].consumer << ']';
    }
    out << "],\n  \"branchial_edges\":[";
    for (size_t i = 0; i < branchial.size(); ++i) {
        if (i) out << ',';
        out << '[' << branchial[i].event1 << ',' << branchial[i].event2 << ']';
    }
    out << "]\n}\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto config = parse_args(argc, argv);
        const auto [lhs, rhs] = parse_rule(config.rule);
        const auto initial = parse_edges(config.init);

        auto builder = make_rule(0);
        for (const auto& edge : lhs) builder.lhs(edge);
        for (const auto& edge : rhs) builder.rhs(edge);

        Hypergraph graph;
        graph.set_state_canonicalization_mode(StateCanonicalizationMode::Full);
        ParallelEvolutionEngine engine(&graph, config.threads);
        engine.set_explore_from_canonical_states_only(true);
        engine.add_rule(builder.build());
        engine.evolve(initial, config.steps);

        const auto parent = std::filesystem::path(config.output).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream output(config.output);
        if (!output) throw std::runtime_error("could not open output: " + config.output);
        write_result(config, graph, output);
        std::cout << "wrote " << config.output << " with " << graph.num_published_states()
                  << " raw states and " << graph.num_published_events() << " events\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
