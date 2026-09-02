#pragma once

#include "research.hpp"

#include <hypergraph/hypergraph.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wgphysics::product {

using infragauge::BaseGraph;
using infragauge::Edge;
using infragauge::FiberBundleConnection;
using infragauge::FiberGraph;
using infragauge::Vertex;
using research::CausalCurvatureReport;
using research::CurvatureEvent;

struct ProductEvent {
    hypergraph::EventId engine_event{};
    hypergraph::StateId raw_input_state{};
    hypergraph::StateId raw_output_state{};
    std::size_t physical_input_state{};
    std::size_t physical_output_state{};
    Vertex subdivided_from{};
    Vertex subdivided_to{};
    Vertex fresh_vertex{};
    std::size_t gauge_orbit_size{};
    double amplitude_per_labeled_factorization{};
    bool curvature_sector_changed{};
};

struct ProductEvolutionResult {
    std::vector<FiberBundleConnection> physical_states;
    std::vector<std::size_t> physical_state_by_raw_state;
    std::vector<ProductEvent> events;
    std::vector<std::pair<std::size_t, std::size_t>> causal_edges;
    CausalCurvatureReport causal_curvature;
};

inline BaseGraph binary_base_graph(
    const hypergraph::Hypergraph& graph, hypergraph::StateId state_id) {
    const auto& state = graph.get_state(state_id);
    std::set<Vertex> vertices;
    std::set<Edge> edges;
    std::size_t raw_edge_count = 0;
    state.edges.for_each([&](hypergraph::EdgeId edge_id) {
        const auto& edge = graph.get_edge(edge_id);
        if (edge.arity != 2 || edge.vertices[0] == edge.vertices[1]) {
            throw std::invalid_argument(
                "gauge product state currently requires loop-free binary hyperedges");
        }
        const auto u = static_cast<Vertex>(edge.vertices[0]);
        const auto v = static_cast<Vertex>(edge.vertices[1]);
        vertices.insert(u);
        vertices.insert(v);
        edges.insert(infragauge::canonical_edge(u, v));
        ++raw_edge_count;
    });
    if (edges.size() != raw_edge_count) {
        throw std::invalid_argument(
            "ordered or duplicate binary hyperedges cannot be collapsed into an undirected connection");
    }
    return BaseGraph(
        {vertices.begin(), vertices.end()}, {edges.begin(), edges.end()});
}

inline bool same_base(const BaseGraph& left, const BaseGraph& right) {
    return left.vertices() == right.vertices() && left.edges() == right.edges();
}

struct SubdivisionMorphism {
    Vertex from{};
    Vertex to{};
    Vertex midpoint{};
};

inline SubdivisionMorphism subdivision_morphism(
    const hypergraph::Hypergraph& graph,
    const hypergraph::Event& event,
    const BaseGraph& input_base) {
    if (event.num_consumed != 1 || event.num_produced != 2) {
        throw std::invalid_argument(
            "unsupported gauge rewrite: expected one consumed and two produced edges");
    }
    const auto& consumed = graph.get_edge(event.consumed_edges[0]);
    if (consumed.arity != 2 || consumed.vertices[0] == consumed.vertices[1]) {
        throw std::invalid_argument("unsupported gauge rewrite: consumed edge is not simple binary");
    }
    const auto from = static_cast<Vertex>(consumed.vertices[0]);
    const auto to = static_cast<Vertex>(consumed.vertices[1]);
    if (!input_base.has_edge(from, to)) {
        throw std::invalid_argument("rewrite provenance consumed an edge absent from its input base");
    }

    std::set<Vertex> produced_vertices;
    std::set<Edge> produced_edges;
    for (std::size_t index = 0; index < event.num_produced; ++index) {
        const auto& produced = graph.get_edge(event.produced_edges[index]);
        if (produced.arity != 2 || produced.vertices[0] == produced.vertices[1]) {
            throw std::invalid_argument("unsupported gauge rewrite: produced edge is not simple binary");
        }
        const auto u = static_cast<Vertex>(produced.vertices[0]);
        const auto v = static_cast<Vertex>(produced.vertices[1]);
        produced_vertices.insert(u);
        produced_vertices.insert(v);
        produced_edges.insert(infragauge::canonical_edge(u, v));
    }
    std::vector<Vertex> fresh;
    std::set_difference(
        produced_vertices.begin(), produced_vertices.end(),
        input_base.vertices().begin(), input_base.vertices().end(),
        std::back_inserter(fresh));
    if (fresh.size() != 1) {
        throw std::invalid_argument("unsupported gauge rewrite: subdivision needs one fresh vertex");
    }
    const auto midpoint = fresh.front();
    const std::set<Edge> expected{
        infragauge::canonical_edge(from, midpoint),
        infragauge::canonical_edge(midpoint, to)};
    if (produced_edges != expected) {
        throw std::invalid_argument(
            "unsupported gauge rewrite: produced edges do not factor the consumed edge");
    }
    return {from, to, midpoint};
}

// Exact sidecar product evolution for the supported subdivision morphism.
// The engine must retain raw states/events: quotienting the bare base before
// connection data is attached can erase transitions that a connection makes
// inequivalent. We therefore consume every raw event and perform the joint
// base-isomorphism/local-gauge quotient here.
inline ProductEvolutionResult evolve_product(
    const hypergraph::Hypergraph& graph,
    hypergraph::StateId initial_state,
    FiberBundleConnection initial_connection) {
    if (initial_state >= graph.num_published_states()) {
        throw std::out_of_range("initial product state is absent from the engine run");
    }
    if (graph.state_canonicalization_mode() != hypergraph::StateCanonicalizationMode::None) {
        throw std::invalid_argument(
            "exact product evolution requires unquotiented raw engine states");
    }
    if (!same_base(initial_connection.base(), binary_base_graph(graph, initial_state))) {
        throw std::invalid_argument("initial connection base does not match the engine state");
    }

    ProductEvolutionResult result;
    result.physical_state_by_raw_state.assign(
        graph.num_published_states(), std::numeric_limits<std::size_t>::max());
    std::vector<std::optional<FiberBundleConnection>> connection_by_raw_state(
        graph.num_published_states());
    std::map<std::vector<Vertex>, std::size_t> physical_ids;
    auto register_physical = [&](const FiberBundleConnection& connection) {
        const auto identity = research::canonical_connection_state(connection);
        const auto found = physical_ids.find(identity);
        if (found != physical_ids.end()) return found->second;
        const auto id = result.physical_states.size();
        physical_ids.emplace(identity, id);
        result.physical_states.push_back(connection);
        return id;
    };

    connection_by_raw_state[initial_state] = initial_connection;
    result.physical_state_by_raw_state[initial_state] = register_physical(initial_connection);

    std::set<hypergraph::EventId> pending;
    for (hypergraph::EventId id = 0; id < graph.num_published_events(); ++id) pending.insert(id);
    std::vector<CurvatureEvent> curvature_events;
    while (!pending.empty()) {
        bool progressed = false;
        for (auto iterator = pending.begin(); iterator != pending.end();) {
            const auto event_id = *iterator;
            const auto& event = graph.get_event(event_id);
            if (event.id == hypergraph::INVALID_ID ||
                event.input_state >= connection_by_raw_state.size() ||
                !connection_by_raw_state[event.input_state]) {
                ++iterator;
                continue;
            }
            if (event.output_state >= connection_by_raw_state.size()) {
                throw std::out_of_range("event output state is absent from the engine run");
            }
            const auto& input_connection = *connection_by_raw_state[event.input_state];
            if (!same_base(input_connection.base(), binary_base_graph(graph, event.input_state))) {
                throw std::logic_error("product connection drifted from its raw input state");
            }
            const auto morphism = subdivision_morphism(graph, event, input_connection.base());
            auto output_connection = input_connection.subdivide_edge_representative(
                morphism.from, morphism.to, morphism.midpoint);
            if (!same_base(output_connection.base(), binary_base_graph(graph, event.output_state))) {
                throw std::logic_error("induced connection base disagrees with the raw output state");
            }
            if (connection_by_raw_state[event.output_state]) {
                throw std::logic_error("raw engine state has multiple product-state parents");
            }
            const auto input_sector = input_connection.gauge_invariant_signature();
            const auto output_sector = output_connection.gauge_invariant_signature();
            const auto physical_output = register_physical(output_connection);
            const auto orbit_size = input_connection.local_gauge_group().size();
            const auto labeling_amplitude = 1.0 / std::sqrt(static_cast<double>(orbit_size));
            connection_by_raw_state[event.output_state] = output_connection;
            result.physical_state_by_raw_state[event.output_state] = physical_output;
            result.events.push_back({
                event_id,
                event.input_state,
                event.output_state,
                result.physical_state_by_raw_state[event.input_state],
                physical_output,
                morphism.from,
                morphism.to,
                morphism.midpoint,
                orbit_size,
                labeling_amplitude,
                input_sector != output_sector});
            curvature_events.push_back({event_id, input_sector, output_sector});
            iterator = pending.erase(iterator);
            progressed = true;
        }
        if (!progressed) {
            throw std::logic_error(
                "engine events are disconnected from the initial raw product state");
        }
    }

    const auto engine_causal = graph.causal_graph().get_causal_edges();
    result.causal_edges.reserve(engine_causal.size());
    for (const auto& edge : engine_causal) {
        result.causal_edges.emplace_back(edge.producer, edge.consumer);
    }
    result.causal_curvature = research::analyze_causal_curvature(
        curvature_events, result.causal_edges, {});
    return result;
}

struct RuleFiberCensusEntry {
    std::vector<Vertex> fiber_graph_code;
    std::size_t automorphism_order{};
    std::size_t initial_curvature_sectors{};
    std::size_t raw_engine_states_per_sector{};
    std::size_t physical_product_states_across_sectors{};
    std::size_t engine_events_across_sectors{};
    std::size_t causal_edges_across_sectors{};
    std::size_t curvature_violations{};
    double maximum_gauge_orbit_norm_error{};
};

inline RuleFiberCensusEntry census_subdivision_rule(
    const hypergraph::Hypergraph& graph,
    hypergraph::StateId initial_state,
    const FiberGraph& fiber) {
    const auto base = binary_base_graph(graph, initial_state);
    if (base.edges().empty()) throw std::invalid_argument("census initial base has no edge");
    const auto group = fiber.automorphisms();
    const auto sectors = research::conjugacy_representatives(group);
    std::set<std::vector<Vertex>> physical_identities;
    RuleFiberCensusEntry entry;
    entry.fiber_graph_code = research::canonical_graph_code(fiber);
    entry.automorphism_order = group.size();
    entry.initial_curvature_sectors = sectors.size();
    entry.raw_engine_states_per_sector = graph.num_published_states();
    const auto seeded_edge = *base.edges().begin();
    for (const auto& sector : sectors) {
        FiberBundleConnection initial(base, fiber);
        initial.set_transport(seeded_edge.first, seeded_edge.second, sector);
        const auto evolution = evolve_product(graph, initial_state, std::move(initial));
        for (const auto& state : evolution.physical_states) {
            physical_identities.insert(research::canonical_connection_state(state));
        }
        entry.engine_events_across_sectors += evolution.events.size();
        entry.causal_edges_across_sectors += evolution.causal_edges.size();
        entry.curvature_violations += evolution.causal_curvature.changed_events;
        for (const auto& event : evolution.events) {
            const auto orbit_norm = static_cast<double>(event.gauge_orbit_size) *
                event.amplitude_per_labeled_factorization *
                event.amplitude_per_labeled_factorization;
            entry.maximum_gauge_orbit_norm_error = std::max(
                entry.maximum_gauge_orbit_norm_error, std::abs(orbit_norm - 1.0));
        }
    }
    entry.physical_product_states_across_sectors = physical_identities.size();
    return entry;
}

}  // namespace wgphysics::product
