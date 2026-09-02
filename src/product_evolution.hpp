#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <hypergraph/hypergraph.hpp>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "research.hpp"

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

template <typename State>
struct BasicProductEvolutionResult {
    std::vector<State> physical_states;
    std::vector<std::size_t> physical_state_by_raw_state;
    std::vector<ProductEvent> events;
    std::vector<std::pair<std::size_t, std::size_t>> causal_edges;
    CausalCurvatureReport causal_curvature;
};

using ProductEvolutionResult = BasicProductEvolutionResult<FiberBundleConnection>;

struct CellProductState {
    research::OrientedCellComplex complex;
    FiberBundleConnection connection;
};

using CellProductEvolutionResult = BasicProductEvolutionResult<CellProductState>;

inline BaseGraph binary_base_graph(const hypergraph::Hypergraph& graph,
                                   hypergraph::StateId state_id) {
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
            "ordered or duplicate binary hyperedges cannot "
            "be collapsed into an undirected connection");
    }
    return BaseGraph({vertices.begin(), vertices.end()}, {edges.begin(), edges.end()});
}

inline bool same_base(const BaseGraph& left, const BaseGraph& right) {
    return left.vertices() == right.vertices() && left.edges() == right.edges();
}

struct SubdivisionMorphism {
    Vertex from{};
    Vertex to{};
    Vertex midpoint{};
};

inline SubdivisionMorphism subdivision_morphism(const hypergraph::Hypergraph& graph,
                                                const hypergraph::Event& event,
                                                const BaseGraph& input_base) {
    if (event.num_consumed != 1 || event.num_produced != 2) {
        throw std::invalid_argument(
            "unsupported gauge rewrite: expected one "
            "consumed and two produced edges");
    }
    const auto& consumed = graph.get_edge(event.consumed_edges[0]);
    if (consumed.arity != 2 || consumed.vertices[0] == consumed.vertices[1]) {
        throw std::invalid_argument(
            "unsupported gauge rewrite: consumed edge is not simple binary");
    }
    const auto from = static_cast<Vertex>(consumed.vertices[0]);
    const auto to = static_cast<Vertex>(consumed.vertices[1]);
    if (!input_base.has_edge(from, to)) {
        throw std::invalid_argument(
            "rewrite provenance consumed an edge absent from its input base");
    }

    std::set<Vertex> produced_vertices;
    std::set<Edge> produced_edges;
    for (std::size_t index = 0; index < event.num_produced; ++index) {
        const auto& produced = graph.get_edge(event.produced_edges[index]);
        if (produced.arity != 2 || produced.vertices[0] == produced.vertices[1]) {
            throw std::invalid_argument(
                "unsupported gauge rewrite: produced edge is not simple binary");
        }
        const auto u = static_cast<Vertex>(produced.vertices[0]);
        const auto v = static_cast<Vertex>(produced.vertices[1]);
        produced_vertices.insert(u);
        produced_vertices.insert(v);
        produced_edges.insert(infragauge::canonical_edge(u, v));
    }
    std::vector<Vertex> fresh;
    std::set_difference(produced_vertices.begin(), produced_vertices.end(),
                        input_base.vertices().begin(), input_base.vertices().end(),
                        std::back_inserter(fresh));
    if (fresh.size() != 1) {
        throw std::invalid_argument(
            "unsupported gauge rewrite: subdivision needs one fresh vertex");
    }
    const auto midpoint = fresh.front();
    const std::set<Edge> expected{infragauge::canonical_edge(from, midpoint),
                                  infragauge::canonical_edge(midpoint, to)};
    if (produced_edges != expected) {
        throw std::invalid_argument(
            "unsupported gauge rewrite: produced edges do "
            "not factor the consumed edge");
    }
    return {from, to, midpoint};
}

// Shared exact event-provenance runner. State operations are supplied by the
// connection-only and explicit-cell wrappers below, so the scheduling, causal
// attachment, and safety checks have one implementation.
template <typename State, typename ConnectionOf, typename IdentityOf, typename SectorOf,
          typename Transition>
inline BasicProductEvolutionResult<State> evolve_product_impl(
    const hypergraph::Hypergraph& graph, hypergraph::StateId initial_state, State initial,
    ConnectionOf connection_of, IdentityOf identity_of, SectorOf sector_of,
    Transition transition) {
    if (initial_state >= graph.num_published_states()) {
        throw std::out_of_range("initial product state is absent from the engine run");
    }
    if (graph.state_canonicalization_mode() != hypergraph::StateCanonicalizationMode::None) {
        throw std::invalid_argument(
            "exact product evolution requires unquotiented raw engine states");
    }
    if (!same_base(connection_of(initial).base(), binary_base_graph(graph, initial_state))) {
        throw std::invalid_argument("initial connection base does not match the engine state");
    }

    BasicProductEvolutionResult<State> result;
    result.physical_state_by_raw_state.assign(graph.num_published_states(),
                                              std::numeric_limits<std::size_t>::max());
    std::vector<std::optional<State>> state_by_raw_state(graph.num_published_states());
    std::map<std::vector<Vertex>, std::size_t> physical_ids;
    auto register_physical = [&](const State& state) {
        const auto identity = identity_of(state);
        const auto found = physical_ids.find(identity);
        if (found != physical_ids.end()) return found->second;
        const auto id = result.physical_states.size();
        physical_ids.emplace(identity, id);
        result.physical_states.push_back(state);
        return id;
    };

    state_by_raw_state[initial_state] = initial;
    result.physical_state_by_raw_state[initial_state] = register_physical(initial);

    std::set<hypergraph::EventId> pending;
    for (hypergraph::EventId id = 0; id < graph.num_published_events(); ++id) pending.insert(id);
    std::vector<CurvatureEvent> curvature_events;
    while (!pending.empty()) {
        bool progressed = false;
        for (auto iterator = pending.begin(); iterator != pending.end();) {
            const auto event_id = *iterator;
            const auto& event = graph.get_event(event_id);
            if (event.id == hypergraph::INVALID_ID ||
                event.input_state >= state_by_raw_state.size() ||
                !state_by_raw_state[event.input_state]) {
                ++iterator;
                continue;
            }
            if (event.output_state >= state_by_raw_state.size()) {
                throw std::out_of_range("event output state is absent from the engine run");
            }
            const auto& input = *state_by_raw_state[event.input_state];
            const auto& input_connection = connection_of(input);
            if (!same_base(input_connection.base(), binary_base_graph(graph, event.input_state))) {
                throw std::logic_error("product connection drifted from its raw input state");
            }
            const auto morphism = subdivision_morphism(graph, event, input_connection.base());
            auto output = transition(input, morphism);
            const auto& output_connection = connection_of(output);
            if (!same_base(output_connection.base(),
                           binary_base_graph(graph, event.output_state))) {
                throw std::logic_error(
                    "induced connection base disagrees with the raw output state");
            }
            if (state_by_raw_state[event.output_state]) {
                throw std::logic_error("raw engine state has multiple product-state parents");
            }
            const auto input_sector = sector_of(input);
            const auto output_sector = sector_of(output);
            const auto physical_output = register_physical(output);
            const auto orbit_size = input_connection.local_gauge_group().size();
            const auto labeling_amplitude = 1.0 / std::sqrt(static_cast<double>(orbit_size));
            state_by_raw_state[event.output_state] = output;
            result.physical_state_by_raw_state[event.output_state] = physical_output;
            result.events.push_back({event_id, event.input_state, event.output_state,
                                     result.physical_state_by_raw_state[event.input_state],
                                     physical_output, morphism.from, morphism.to, morphism.midpoint,
                                     orbit_size, labeling_amplitude,
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
    std::set<std::size_t> source_events;
    const auto initial_step = graph.get_state(initial_state).step;
    for (const auto& event : curvature_events) {
        if (event.before != event.after &&
            graph.get_state(graph.get_event(event.id).input_state).step == initial_step) {
            source_events.insert(event.id);
        }
    }
    result.causal_curvature =
        research::analyze_causal_curvature(curvature_events, result.causal_edges, source_events);
    return result;
}

// Exact sidecar product evolution for the supported subdivision morphism.
// The engine must retain raw states/events: quotienting the bare base before
// connection data is attached can erase transitions that a connection makes
// inequivalent. We therefore consume every raw event and perform the joint
// base-isomorphism/local-gauge quotient here.
inline ProductEvolutionResult evolve_product(
    const hypergraph::Hypergraph& graph, hypergraph::StateId initial_state,
    FiberBundleConnection initial_connection,
    const research::LocalGaugeDynamics* cell_dynamics = nullptr) {
    return evolve_product_impl(
        graph, initial_state, std::move(initial_connection),
        [](const FiberBundleConnection& state) -> const FiberBundleConnection& { return state; },
        [](const FiberBundleConnection& state) {
            return research::canonical_connection_state(state);
        },
        [](const FiberBundleConnection& state) { return state.gauge_invariant_signature(); },
        [cell_dynamics](const FiberBundleConnection& input,
                        const SubdivisionMorphism& morphism) {
            auto output = input.subdivide_edge_representative(
                morphism.from, morphism.to, morphism.midpoint);
            if (cell_dynamics) {
                output = research::apply_holonomy_dynamics(
                    output, research::simple_cycle_loop(output.base()), *cell_dynamics);
            }
            return output;
        });
}

// The same official-engine closure with explicit oriented faces carried
// through each event. Cell incidence participates in physical identity, and
// curvature is measured on those faces rather than an inferred cycle basis.
inline CellProductEvolutionResult evolve_cell_product(
    const hypergraph::Hypergraph& graph, hypergraph::StateId initial_state,
    CellProductState initial) {
    if (!same_base(initial.complex.base(), initial.connection.base())) {
        throw std::invalid_argument("initial cell complex and connection bases differ");
    }
    return evolve_product_impl(
        graph, initial_state, std::move(initial),
        [](const CellProductState& state) -> const FiberBundleConnection& {
            return state.connection;
        },
        [](const CellProductState& state) {
            return research::canonical_cell_connection_state(state.complex, state.connection);
        },
        [](const CellProductState& state) {
            return research::curvature_sector(state.connection, state.complex.faces());
        },
        [](const CellProductState& input, const SubdivisionMorphism& morphism) {
            return CellProductState{
                input.complex.subdivide_edge(
                    morphism.from, morphism.to, morphism.midpoint),
                input.connection.subdivide_edge_representative(
                    morphism.from, morphism.to, morphism.midpoint)};
        });
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

struct DynamicsCensusEntry {
    std::size_t rule_index{};
    std::vector<uint16_t> element_map;
    std::vector<uint16_t> curvature_sector_map;
    std::size_t fixed_elements{};
    std::size_t fixed_curvature_sectors{};
    std::size_t physical_product_states{};
    std::size_t changed_events{};
    std::size_t source_events{};
    std::size_t causally_reached_changes{};
    std::size_t off_causal_changes{};
    std::size_t maximum_causal_depth{};
    double causal_alignment{};
};

inline std::vector<DynamicsCensusEntry> census_cell_dynamics(const hypergraph::Hypergraph& graph,
                                                             hypergraph::StateId initial_state,
                                                             const FiberGraph& fiber) {
    const auto base = binary_base_graph(graph, initial_state);
    if (base.edges().empty()) throw std::invalid_argument("dynamics census base has no edge");
    const infragauge::AutomorphismTables tables(fiber.automorphisms());
    const auto rules = research::search_local_gauge_dynamics(tables);
    const auto sectors = research::conjugacy_representatives(tables.elements());
    const auto seeded_edge = *base.edges().begin();
    std::vector<DynamicsCensusEntry> entries;
    entries.reserve(rules.size());
    for (std::size_t rule_index = 0; rule_index < rules.size(); ++rule_index) {
        const auto& rule = rules[rule_index];
        DynamicsCensusEntry entry;
        entry.rule_index = rule_index;
        entry.element_map = rule.image;
        entry.curvature_sector_map = research::induced_conjugacy_sector_map(tables, rule);
        entry.fixed_elements = rule.fixed_elements;
        for (std::size_t sector = 0; sector < entry.curvature_sector_map.size(); ++sector) {
            if (entry.curvature_sector_map[sector] == sector) ++entry.fixed_curvature_sectors;
        }
        std::set<std::vector<Vertex>> physical_identities;
        for (const auto& representative : sectors) {
            FiberBundleConnection initial(base, fiber);
            initial.set_transport(seeded_edge.first, seeded_edge.second, representative);
            const auto evolution = evolve_product(graph, initial_state, std::move(initial), &rule);
            for (const auto& state : evolution.physical_states) {
                physical_identities.insert(research::canonical_connection_state(state));
            }
            entry.changed_events += evolution.causal_curvature.changed_events;
            entry.source_events += evolution.causal_curvature.source_events;
            entry.causally_reached_changes += evolution.causal_curvature.causally_reached_changes;
            entry.off_causal_changes += evolution.causal_curvature.off_causal_changes;
            entry.maximum_causal_depth = std::max(entry.maximum_causal_depth,
                                                  evolution.causal_curvature.maximum_causal_depth);
        }
        entry.physical_product_states = physical_identities.size();
        const auto non_source_changes = entry.changed_events - entry.source_events;
        entry.causal_alignment =
            non_source_changes == 0
                ? 1.0
                : static_cast<double>(entry.causally_reached_changes) / non_source_changes;
        entries.push_back(std::move(entry));
    }
    return entries;
}

inline RuleFiberCensusEntry census_subdivision_rule(const hypergraph::Hypergraph& graph,
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
            entry.maximum_gauge_orbit_norm_error =
                std::max(entry.maximum_gauge_orbit_norm_error, std::abs(orbit_norm - 1.0));
        }
    }
    entry.physical_product_states_across_sectors = physical_identities.size();
    return entry;
}

}  // namespace wgphysics::product
