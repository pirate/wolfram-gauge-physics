#pragma once

#include "infragauge.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wgphysics::research {

using infragauge::AutomorphismTables;
using infragauge::BaseGraph;
using infragauge::Edge;
using infragauge::FiberBundleConnection;
using infragauge::FiberGraph;
using infragauge::Permutation;
using infragauge::Vertex;

inline std::vector<Vertex> canonical_graph_code(const FiberGraph& graph) {
    std::vector<Vertex> labels(graph.vertex_count());
    std::iota(labels.begin(), labels.end(), 0);
    std::vector<Vertex> best;
    bool first = true;
    do {
        std::vector<Vertex> candidate{graph.vertex_count()};
        for (Vertex u = 0; u < graph.vertex_count(); ++u) {
            for (Vertex v = u + 1; v < graph.vertex_count(); ++v) {
                candidate.push_back(graph.adjacent(labels[u], labels[v]) ? 1U : 0U);
            }
        }
        if (first || candidate < best) {
            best = std::move(candidate);
            first = false;
        }
    } while (std::next_permutation(labels.begin(), labels.end()));
    return best;
}

inline FiberGraph fiber_from_canonical_code(const std::vector<Vertex>& code) {
    if (code.empty()) throw std::invalid_argument("fiber graph code is empty");
    const auto count = code.front();
    if (code.size() != 1 + static_cast<std::size_t>(count) * (count - 1) / 2) {
        throw std::invalid_argument("fiber graph code has the wrong length");
    }
    std::vector<Edge> edges;
    std::size_t position = 1;
    for (Vertex u = 0; u < count; ++u) {
        for (Vertex v = u + 1; v < count; ++v) {
            if (code[position] > 1) throw std::invalid_argument("fiber graph code is not binary");
            if (code[position++]) edges.emplace_back(u, v);
        }
    }
    return FiberGraph(count, std::move(edges));
}

inline std::vector<Permutation> conjugacy_representatives(
    const std::vector<Permutation>& group) {
    if (group.empty()) throw std::invalid_argument("group cannot be empty");
    std::map<std::vector<Vertex>, Permutation> representatives;
    for (const auto& value : group) {
        std::vector<Vertex> canonical;
        bool first = true;
        for (const auto& frame : group) {
            const auto conjugated = Permutation::compose(
                frame, Permutation::compose(value, frame.inverse()));
            if (first || conjugated.image() < canonical) {
                canonical = conjugated.image();
                first = false;
            }
        }
        representatives.try_emplace(std::move(canonical), value);
    }
    std::vector<Permutation> result;
    result.reserve(representatives.size());
    for (const auto& [_, value] : representatives) result.push_back(value);
    return result;
}

struct InferredFiberClass {
    FiberGraph representative;
    std::vector<Vertex> centers;
    std::vector<Vertex> canonical_code;
    std::size_t automorphism_order{};
};

// Infer a candidate fiber from the graph induced by each vertex's open link.
// Exact graph canonicalization groups repeated neighborhoods without using
// coordinates, supplied group names, or a learned embedding.
inline std::vector<InferredFiberClass> infer_link_fibers(
    const BaseGraph& base, std::size_t minimum_multiplicity = 1) {
    std::vector<InferredFiberClass> classes;
    for (const auto center : base.vertices()) {
        const auto neighbors = base.neighbors(center);
        if (neighbors.empty()) continue;
        std::map<Vertex, Vertex> local_label;
        for (Vertex i = 0; i < neighbors.size(); ++i) local_label.emplace(neighbors[i], i);
        std::vector<Edge> link_edges;
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
            for (std::size_t j = i + 1; j < neighbors.size(); ++j) {
                if (base.has_edge(neighbors[i], neighbors[j])) {
                    link_edges.emplace_back(static_cast<Vertex>(i), static_cast<Vertex>(j));
                }
            }
        }
        FiberGraph link(static_cast<Vertex>(neighbors.size()), std::move(link_edges));
        auto code = canonical_graph_code(link);
        const auto found = std::find_if(classes.begin(), classes.end(), [&](const auto& value) {
            return value.canonical_code == code;
        });
        if (found == classes.end()) {
            classes.push_back(
                {link, {center}, std::move(code), link.automorphisms().size()});
        } else {
            found->centers.push_back(center);
        }
    }
    classes.erase(
        std::remove_if(classes.begin(), classes.end(), [&](const auto& value) {
            return value.centers.size() < minimum_multiplicity;
        }),
        classes.end());
    std::sort(classes.begin(), classes.end(), [](const auto& left, const auto& right) {
        return left.canonical_code < right.canonical_code;
    });
    return classes;
}

// Exact identity for small uniform-fiber states. It quotients arbitrary base
// relabelings and all local fiber-frame choices. A canonical fiber template is
// part of the identity, so different microscopic fibers cannot collide.
inline std::vector<Vertex> canonical_connection_state(
    const FiberBundleConnection& connection) {
    const std::vector<Vertex> old_vertices(
        connection.base().vertices().begin(), connection.base().vertices().end());
    std::vector<Vertex> new_labels(old_vertices.size());
    std::iota(new_labels.begin(), new_labels.end(), 0);

    // Lexicographic minimization factors exactly: first retain only labelings
    // with the least bare-base adjacency code, then minimize the connection
    // quotient over ties. This avoids reconstructing a connection for all n!
    // labelings. On an n-cycle only its 2n canonical dihedral labelings survive.
    std::vector<Vertex> best_base_code;
    std::vector<std::vector<Vertex>> canonical_labelings;
    bool first_base = true;
    do {
        std::vector<Vertex> old_at_new(old_vertices.size());
        for (std::size_t i = 0; i < old_vertices.size(); ++i) {
            old_at_new[new_labels[i]] = old_vertices[i];
        }
        std::vector<Vertex> base_code;
        for (Vertex u = 0; u < old_vertices.size(); ++u) {
            for (Vertex v = u + 1; v < old_vertices.size(); ++v) {
                base_code.push_back(
                    connection.base().has_edge(old_at_new[u], old_at_new[v]) ? 1U : 0U);
            }
        }
        if (first_base || base_code < best_base_code) {
            best_base_code = std::move(base_code);
            canonical_labelings.clear();
            canonical_labelings.push_back(new_labels);
            first_base = false;
        } else if (base_code == best_base_code) {
            canonical_labelings.push_back(new_labels);
        }
    } while (std::next_permutation(new_labels.begin(), new_labels.end()));

    std::vector<Vertex> best_gauge_code;
    bool first_gauge = true;
    for (const auto& labeling : canonical_labelings) {
        std::map<Vertex, Vertex> relabel;
        for (std::size_t i = 0; i < old_vertices.size(); ++i) {
            relabel.emplace(old_vertices[i], labeling[i]);
        }
        std::vector<Edge> relabeled_edges;
        for (const auto [u, v] : connection.base().edges()) {
            relabeled_edges.emplace_back(relabel.at(u), relabel.at(v));
        }
        std::vector<Vertex> vertices(old_vertices.size());
        std::iota(vertices.begin(), vertices.end(), 0);
        FiberBundleConnection relabeled(
            BaseGraph(vertices, relabeled_edges), connection.fiber());
        for (const auto [u, v] : connection.base().edges()) {
            relabeled.set_transport(
                relabel.at(u), relabel.at(v), connection.edge_transport(u, v));
        }
        const auto gauge_code = relabeled.gauge_invariant_signature();
        if (first_gauge || gauge_code < best_gauge_code) {
            best_gauge_code = gauge_code;
            first_gauge = false;
        }
    }

    std::vector<Vertex> result{static_cast<Vertex>(old_vertices.size())};
    const auto fiber_code = canonical_graph_code(connection.fiber());
    result.insert(result.end(), fiber_code.begin(), fiber_code.end());
    result.push_back(std::numeric_limits<Vertex>::max());
    result.insert(result.end(), best_base_code.begin(), best_base_code.end());
    result.push_back(std::numeric_limits<Vertex>::max());
    result.insert(result.end(), best_gauge_code.begin(), best_gauge_code.end());
    return result;
}

struct GaugeRewriteEvent {
    std::size_t input_state{};
    std::size_t output_state{};
    std::size_t raw_extensions{};
};

class GaugeSubdivisionEvolution {
public:
    explicit GaugeSubdivisionEvolution(FiberBundleConnection initial) {
        register_state(std::move(initial));
    }

    std::size_t state_count() const { return states_.size(); }
    const FiberBundleConnection& state(std::size_t id) const { return states_.at(id); }
    const std::vector<GaugeRewriteEvent>& events() const { return events_; }

    std::vector<std::size_t> subdivide(
        std::size_t input, Vertex from, Vertex to, Vertex midpoint) {
        const auto raw_orbit_size = states_.at(input).local_gauge_group().size();
        auto representative = states_.at(input).subdivide_edge_representative(from, to, midpoint);
        const auto output = register_state(std::move(representative));
        events_.push_back({input, output, raw_orbit_size});
        return {output};
    }

private:
    std::vector<FiberBundleConnection> states_;
    std::map<std::vector<Vertex>, std::size_t> identities_;
    std::vector<GaugeRewriteEvent> events_;

    std::size_t register_state(FiberBundleConnection state) {
        auto identity = canonical_connection_state(state);
        const auto found = identities_.find(identity);
        if (found != identities_.end()) return found->second;
        const auto id = states_.size();
        identities_.emplace(std::move(identity), id);
        states_.push_back(std::move(state));
        return id;
    }
};

class PartialFiberMap {
public:
    PartialFiberMap(
        Vertex source_size, Vertex target_size, std::vector<std::optional<Vertex>> image)
        : source_size_(source_size), target_size_(target_size), image_(std::move(image)) {
        if (image_.size() != source_size_) {
            throw std::invalid_argument("partial map has the wrong source size");
        }
        std::set<Vertex> used;
        for (const auto value : image_) {
            if (!value) continue;
            if (*value >= target_size_) throw std::out_of_range("partial map target is absent");
            if (!used.insert(*value).second) throw std::invalid_argument("partial map is not injective");
        }
    }

    Vertex source_size() const { return source_size_; }
    Vertex target_size() const { return target_size_; }
    std::optional<Vertex> operator()(Vertex source) const {
        if (source >= source_size_) throw std::out_of_range("partial map source is absent");
        return image_[source];
    }
    std::size_t domain_size() const {
        return static_cast<std::size_t>(std::count_if(
            image_.begin(), image_.end(), [](const auto& value) { return value.has_value(); }));
    }
    bool preserves_induced_adjacency(
        const FiberGraph& source, const FiberGraph& target) const {
        if (source.vertex_count() != source_size_ || target.vertex_count() != target_size_) {
            return false;
        }
        for (Vertex u = 0; u < source_size_; ++u) {
            if (!image_[u]) continue;
            for (Vertex v = u + 1; v < source_size_; ++v) {
                if (!image_[v]) continue;
                if (source.adjacent(u, v) != target.adjacent(*image_[u], *image_[v])) return false;
            }
        }
        return true;
    }
    PartialFiberMap inverse() const {
        std::vector<std::optional<Vertex>> result(target_size_);
        for (Vertex source = 0; source < source_size_; ++source) {
            if (image_[source]) result[*image_[source]] = source;
        }
        return PartialFiberMap(target_size_, source_size_, std::move(result));
    }
    static PartialFiberMap compose(
        const PartialFiberMap& after, const PartialFiberMap& before) {
        if (before.target_size_ != after.source_size_) {
            throw std::invalid_argument("partial-map intermediate fibers differ");
        }
        std::vector<std::optional<Vertex>> result(before.source_size_);
        for (Vertex source = 0; source < before.source_size_; ++source) {
            const auto middle = before(source);
            if (middle) result[source] = after(*middle);
        }
        return PartialFiberMap(before.source_size_, after.target_size_, std::move(result));
    }

private:
    Vertex source_size_;
    Vertex target_size_;
    std::vector<std::optional<Vertex>> image_;
};

class VariableFiberBundle {
public:
    explicit VariableFiberBundle(BaseGraph base) : base_(std::move(base)) {}

    void set_fiber(Vertex base_vertex, FiberGraph fiber) {
        if (!base_.has_vertex(base_vertex)) throw std::out_of_range("fiber base vertex is absent");
        fibers_.insert_or_assign(base_vertex, std::move(fiber));
    }
    void set_partial_transport(Vertex from, Vertex to, const PartialFiberMap& transport) {
        if (!base_.has_edge(from, to)) throw std::out_of_range("transport edge is absent");
        const auto& source = fibers_.at(from);
        const auto& target = fibers_.at(to);
        if (!transport.preserves_induced_adjacency(source, target)) {
            throw std::invalid_argument("partial lift is not an induced fiber embedding");
        }
        transports_.insert_or_assign({from, to}, transport);
        transports_.insert_or_assign({to, from}, transport.inverse());
    }
    std::optional<Vertex> lift(const std::vector<Vertex>& path, Vertex start) const {
        if (path.empty()) throw std::invalid_argument("lift path cannot be empty");
        if (!fibers_.contains(path.front()) || start >= fibers_.at(path.front()).vertex_count()) {
            throw std::out_of_range("lift starts outside its fiber");
        }
        std::optional<Vertex> current = start;
        for (std::size_t i = 1; i < path.size() && current; ++i) {
            const auto found = transports_.find({path[i - 1], path[i]});
            if (found == transports_.end()) throw std::out_of_range("path has no partial lift");
            current = found->second(*current);
        }
        return current;
    }

private:
    BaseGraph base_;
    std::map<Vertex, FiberGraph> fibers_;
    std::map<std::pair<Vertex, Vertex>, PartialFiberMap> transports_;
};

struct LocalGaugeDynamics {
    std::vector<uint16_t> image;
    std::size_t fixed_elements{};
};

inline std::vector<LocalGaugeDynamics> search_local_gauge_dynamics(
    const AutomorphismTables& group) {
    std::vector<uint16_t> candidate(group.order());
    std::iota(candidate.begin(), candidate.end(), 0);
    const auto identity = group.index_of(Permutation::identity(group.fiber_size()));
    std::vector<LocalGaugeDynamics> result;
    do {
        if (candidate[identity] != identity) continue;
        bool valid = true;
        for (uint16_t value = 0; value < group.order() && valid; ++value) {
            if (candidate[group.inverse(value)] != group.inverse(candidate[value])) valid = false;
            for (uint16_t frame = 0; frame < group.order() && valid; ++frame) {
                const auto conjugated = group.multiply(
                    frame, group.multiply(value, group.inverse(frame)));
                const auto transformed_conjugate = group.multiply(
                    frame, group.multiply(candidate[value], group.inverse(frame)));
                if (candidate[conjugated] != transformed_conjugate) valid = false;
            }
        }
        if (valid) {
            std::size_t actual_fixed = 0;
            for (uint16_t value = 0; value < group.order(); ++value) {
                if (candidate[value] == value) ++actual_fixed;
            }
            result.push_back({candidate, actual_fixed});
        }
    } while (std::next_permutation(candidate.begin(), candidate.end()));
    return result;
}

struct FiberCensusEntry {
    Vertex vertices{};
    std::size_t edges{};
    std::vector<Vertex> graph_code;
    std::size_t automorphism_order{};
    bool nonabelian{};
    std::size_t conjugacy_classes{};
    std::optional<std::size_t> reversible_equivariant_dynamics;
    std::size_t subdivision_raw_extensions{};
    std::size_t subdivision_physical_children{};
    std::size_t exact_subdivision_schmidt_rank{};
    double half_rank_discarded_norm{};
};

inline std::vector<FiberCensusEntry> enumerate_fiber_census(
    Vertex maximum_vertices = 4, std::size_t maximum_dynamics_order = 8) {
    if (maximum_vertices == 0 || maximum_vertices > 5) {
        throw std::invalid_argument("exact fiber census supports one through five vertices");
    }
    std::map<std::vector<Vertex>, FiberGraph> representatives;
    for (Vertex count = 1; count <= maximum_vertices; ++count) {
        std::vector<Edge> possible_edges;
        for (Vertex u = 0; u < count; ++u) {
            for (Vertex v = u + 1; v < count; ++v) possible_edges.emplace_back(u, v);
        }
        const auto graph_count = std::uint64_t{1} << possible_edges.size();
        for (std::uint64_t mask = 0; mask < graph_count; ++mask) {
            std::vector<Edge> edges;
            for (std::size_t bit = 0; bit < possible_edges.size(); ++bit) {
                if ((mask >> bit) & 1U) edges.push_back(possible_edges[bit]);
            }
            FiberGraph graph(count, std::move(edges));
            representatives.try_emplace(canonical_graph_code(graph), graph);
        }
    }

    std::vector<FiberCensusEntry> result;
    result.reserve(representatives.size());
    for (const auto& [code, graph] : representatives) {
        const auto automorphisms = graph.automorphisms();
        bool nonabelian = false;
        for (const auto& left : automorphisms) {
            for (const auto& right : automorphisms) {
                if (Permutation::compose(left, right) != Permutation::compose(right, left)) {
                    nonabelian = true;
                }
            }
        }
        const auto conjugacy_classes = conjugacy_representatives(automorphisms);
        std::optional<std::size_t> dynamics;
        if (automorphisms.size() <= maximum_dynamics_order) {
            dynamics = search_local_gauge_dynamics(AutomorphismTables(automorphisms)).size();
        }
        const auto retained_rank = automorphisms.size() / 2;
        const auto discarded = 1.0 -
            static_cast<double>(retained_rank) / static_cast<double>(automorphisms.size());
        result.push_back({
            graph.vertex_count(), graph.edges().size(), code, automorphisms.size(), nonabelian,
            conjugacy_classes.size(), dynamics, automorphisms.size(), 1, automorphisms.size(),
            discarded});
    }
    return result;
}

inline std::vector<Vertex> curvature_sector(
    const FiberBundleConnection& connection,
    const std::vector<std::vector<Vertex>>& loops) {
    std::vector<Vertex> result;
    const auto& group = connection.local_gauge_group();
    for (const auto& loop : loops) {
        const auto value = connection.holonomy(loop);
        std::vector<Vertex> canonical;
        bool first = true;
        for (const auto& frame : group) {
            const auto conjugated = Permutation::compose(
                frame, Permutation::compose(value, frame.inverse()));
            if (first || conjugated.image() < canonical) {
                canonical = conjugated.image();
                first = false;
            }
        }
        result.push_back(static_cast<Vertex>(canonical.size()));
        result.insert(result.end(), canonical.begin(), canonical.end());
    }
    return result;
}

struct CurvatureEvent {
    std::size_t id{};
    std::vector<Vertex> before;
    std::vector<Vertex> after;
};

struct CausalCurvatureReport {
    std::size_t changed_events{};
    std::size_t source_events{};
    std::size_t causally_reached_changes{};
    std::size_t off_causal_changes{};
    std::size_t maximum_causal_depth{};
    double causal_alignment{};
};

inline CausalCurvatureReport analyze_causal_curvature(
    const std::vector<CurvatureEvent>& events,
    const std::vector<std::pair<std::size_t, std::size_t>>& causal_edges,
    const std::set<std::size_t>& declared_sources) {
    std::map<std::size_t, CurvatureEvent> by_id;
    std::set<std::size_t> changed;
    for (const auto& event : events) {
        if (!by_id.emplace(event.id, event).second) throw std::invalid_argument("duplicate event id");
        if (event.before != event.after) changed.insert(event.id);
    }
    for (const auto source : declared_sources) {
        if (!changed.contains(source)) throw std::invalid_argument("source event has no curvature change");
    }
    std::map<std::size_t, std::vector<std::size_t>> outgoing;
    for (const auto [producer, consumer] : causal_edges) {
        if (!by_id.contains(producer) || !by_id.contains(consumer)) {
            throw std::out_of_range("causal edge references an absent event");
        }
        outgoing[producer].push_back(consumer);
    }
    std::map<std::size_t, std::size_t> distance;
    std::queue<std::size_t> frontier;
    for (const auto source : declared_sources) {
        distance[source] = 0;
        frontier.push(source);
    }
    while (!frontier.empty()) {
        const auto producer = frontier.front();
        frontier.pop();
        for (const auto consumer : outgoing[producer]) {
            if (!distance.contains(consumer)) {
                distance[consumer] = distance[producer] + 1;
                frontier.push(consumer);
            }
        }
    }
    CausalCurvatureReport report;
    report.changed_events = changed.size();
    report.source_events = declared_sources.size();
    for (const auto event : changed) {
        if (declared_sources.contains(event)) continue;
        const auto found = distance.find(event);
        if (found == distance.end()) {
            ++report.off_causal_changes;
        } else {
            ++report.causally_reached_changes;
            report.maximum_causal_depth = std::max(report.maximum_causal_depth, found->second);
        }
    }
    const auto non_source = report.changed_events - report.source_events;
    report.causal_alignment = non_source == 0
        ? 1.0
        : static_cast<double>(report.causally_reached_changes) / non_source;
    return report;
}

}  // namespace wgphysics::research
