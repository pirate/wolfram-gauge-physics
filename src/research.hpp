#pragma once

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

#include "infragauge.hpp"

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

inline std::vector<Permutation> conjugacy_representatives(const std::vector<Permutation>& group) {
    if (group.empty()) throw std::invalid_argument("group cannot be empty");
    std::map<std::vector<Vertex>, Permutation> representatives;
    for (const auto& value : group) {
        std::vector<Vertex> canonical;
        bool first = true;
        for (const auto& frame : group) {
            const auto conjugated =
                Permutation::compose(frame, Permutation::compose(value, frame.inverse()));
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
inline std::vector<InferredFiberClass> infer_link_fibers(const BaseGraph& base,
                                                         std::size_t minimum_multiplicity = 1) {
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
            classes.push_back({link, {center}, std::move(code), link.automorphisms().size()});
        } else {
            found->centers.push_back(center);
        }
    }
    classes.erase(std::remove_if(classes.begin(), classes.end(),
                                 [&](const auto& value) {
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
inline std::vector<Vertex> canonical_connection_state(const FiberBundleConnection& connection) {
    const std::vector<Vertex> old_vertices(connection.base().vertices().begin(),
                                           connection.base().vertices().end());
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
                base_code.push_back(connection.base().has_edge(old_at_new[u], old_at_new[v]) ? 1U
                                                                                             : 0U);
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
        FiberBundleConnection relabeled(BaseGraph(vertices, relabeled_edges), connection.fiber());
        for (const auto [u, v] : connection.base().edges()) {
            relabeled.set_transport(relabel.at(u), relabel.at(v), connection.edge_transport(u, v));
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

    std::vector<std::size_t> subdivide(std::size_t input, Vertex from, Vertex to, Vertex midpoint) {
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

inline std::vector<Vertex> canonical_oriented_face(std::vector<Vertex> boundary) {
    if (boundary.size() < 4 || boundary.front() != boundary.back()) {
        throw std::invalid_argument("an oriented face must be a closed nondegenerate boundary");
    }
    boundary.pop_back();
    std::vector<Vertex> best;
    bool first = true;
    for (std::size_t shift = 0; shift < boundary.size(); ++shift) {
        std::vector<Vertex> candidate;
        candidate.reserve(boundary.size() + 1);
        for (std::size_t index = 0; index < boundary.size(); ++index) {
            candidate.push_back(boundary[(shift + index) % boundary.size()]);
        }
        candidate.push_back(candidate.front());
        if (first || candidate < best) {
            best = std::move(candidate);
            first = false;
        }
    }
    return best;
}

// An oriented combinatorial two-complex. Faces are first-class cyclic boundary
// words; they are never guessed from cycles in the one-skeleton.
class OrientedCellComplex {
   public:
    OrientedCellComplex(BaseGraph base, std::vector<std::vector<Vertex>> faces)
        : base_(std::move(base)) {
        faces_.reserve(faces.size());
        std::set<std::vector<Vertex>> unique_faces;
        for (auto& face : faces) {
            auto canonical = canonical_oriented_face(std::move(face));
            std::set<Vertex> distinct(canonical.begin(), canonical.end() - 1);
            if (distinct.size() + 1 != canonical.size()) {
                throw std::invalid_argument("cell boundary repeats a vertex");
            }
            for (std::size_t index = 1; index < canonical.size(); ++index) {
                if (!base_.has_edge(canonical[index - 1], canonical[index])) {
                    throw std::invalid_argument("cell boundary traverses an absent base edge");
                }
            }
            if (!unique_faces.insert(canonical).second) {
                throw std::invalid_argument("duplicate oriented face");
            }
            faces_.push_back(std::move(canonical));
        }
    }

    const BaseGraph& base() const { return base_; }
    const std::vector<std::vector<Vertex>>& faces() const { return faces_; }

    OrientedCellComplex subdivide_edge(Vertex from, Vertex to, Vertex midpoint) const {
        const auto changed_edge = infragauge::canonical_edge(from, to);
        std::vector<std::vector<Vertex>> subdivided_faces;
        subdivided_faces.reserve(faces_.size());
        for (const auto& face : faces_) {
            std::vector<Vertex> boundary{face.front()};
            for (std::size_t index = 1; index < face.size(); ++index) {
                if (infragauge::canonical_edge(face[index - 1], face[index]) == changed_edge) {
                    boundary.push_back(midpoint);
                }
                boundary.push_back(face[index]);
            }
            subdivided_faces.push_back(std::move(boundary));
        }
        return OrientedCellComplex(base_.subdivide(from, to, midpoint),
                                   std::move(subdivided_faces));
    }

   private:
    BaseGraph base_;
    std::vector<std::vector<Vertex>> faces_;
};

inline std::vector<Vertex> canonical_cell_connection_state(
    const OrientedCellComplex& complex, const FiberBundleConnection& connection) {
    if (complex.base().vertices() != connection.base().vertices() ||
        complex.base().edges() != connection.base().edges()) {
        throw std::invalid_argument("cell complex and connection have different base graphs");
    }
    const std::vector<Vertex> old_vertices(complex.base().vertices().begin(),
                                           complex.base().vertices().end());
    std::vector<Vertex> labels(old_vertices.size());
    std::iota(labels.begin(), labels.end(), 0);
    std::vector<Vertex> best;
    bool first = true;
    do {
        std::map<Vertex, Vertex> relabel;
        for (std::size_t index = 0; index < old_vertices.size(); ++index) {
            relabel.emplace(old_vertices[index], labels[index]);
        }
        std::vector<Edge> edges;
        for (const auto [from, to] : complex.base().edges()) {
            edges.emplace_back(relabel.at(from), relabel.at(to));
        }
        std::vector<Vertex> vertices(old_vertices.size());
        std::iota(vertices.begin(), vertices.end(), 0);
        FiberBundleConnection relabeled_connection(BaseGraph(vertices, edges), connection.fiber());
        for (const auto [from, to] : connection.base().edges()) {
            relabeled_connection.set_transport(relabel.at(from), relabel.at(to),
                                               connection.edge_transport(from, to));
        }
        std::vector<std::vector<Vertex>> faces;
        for (const auto& face : complex.faces()) {
            std::vector<Vertex> relabeled_face;
            relabeled_face.reserve(face.size());
            for (const auto vertex : face) relabeled_face.push_back(relabel.at(vertex));
            faces.push_back(canonical_oriented_face(std::move(relabeled_face)));
        }
        std::sort(faces.begin(), faces.end());

        std::vector<Vertex> candidate{static_cast<Vertex>(vertices.size())};
        const auto fiber_code = canonical_graph_code(connection.fiber());
        candidate.insert(candidate.end(), fiber_code.begin(), fiber_code.end());
        candidate.push_back(std::numeric_limits<Vertex>::max());
        for (Vertex from = 0; from < vertices.size(); ++from) {
            for (Vertex to = from + 1; to < vertices.size(); ++to) {
                candidate.push_back(relabeled_connection.base().has_edge(from, to) ? 1U : 0U);
            }
        }
        candidate.push_back(static_cast<Vertex>(faces.size()));
        for (const auto& face : faces) {
            candidate.push_back(static_cast<Vertex>(face.size() - 1));
            candidate.insert(candidate.end(), face.begin(), face.end() - 1);
        }
        candidate.push_back(std::numeric_limits<Vertex>::max());
        const auto gauge_code = relabeled_connection.gauge_invariant_signature();
        candidate.insert(candidate.end(), gauge_code.begin(), gauge_code.end());
        if (first || candidate < best) {
            best = std::move(candidate);
            first = false;
        }
    } while (std::next_permutation(labels.begin(), labels.end()));
    return best;
}

class PartialFiberMap {
   public:
    PartialFiberMap(Vertex source_size, Vertex target_size,
                    std::vector<std::optional<Vertex>> image)
        : source_size_(source_size), target_size_(target_size), image_(std::move(image)) {
        if (image_.size() != source_size_) {
            throw std::invalid_argument("partial map has the wrong source size");
        }
        std::set<Vertex> used;
        for (const auto value : image_) {
            if (!value) continue;
            if (*value >= target_size_) throw std::out_of_range("partial map target is absent");
            if (!used.insert(*value).second)
                throw std::invalid_argument("partial map is not injective");
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
    bool preserves_induced_adjacency(const FiberGraph& source, const FiberGraph& target) const {
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
    static PartialFiberMap compose(const PartialFiberMap& after, const PartialFiberMap& before) {
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

inline void validate_local_gauge_dynamics(const AutomorphismTables& group,
                                          const LocalGaugeDynamics& dynamics) {
    if (dynamics.image.size() != group.order()) {
        throw std::invalid_argument("holonomy dynamics uses a different gauge group");
    }
    std::vector<uint16_t> sorted = dynamics.image;
    std::sort(sorted.begin(), sorted.end());
    for (uint16_t value = 0; value < group.order(); ++value) {
        if (sorted[value] != value) {
            throw std::invalid_argument("holonomy dynamics is not reversible");
        }
    }
    const auto identity = group.index_of(Permutation::identity(group.fiber_size()));
    if (dynamics.image[identity] != identity) {
        throw std::invalid_argument("holonomy dynamics does not preserve flatness");
    }
    for (uint16_t value = 0; value < group.order(); ++value) {
        if (dynamics.image[group.inverse(value)] != group.inverse(dynamics.image[value])) {
            throw std::invalid_argument("holonomy dynamics depends on loop orientation");
        }
        for (uint16_t frame = 0; frame < group.order(); ++frame) {
            const auto conjugated =
                group.multiply(frame, group.multiply(value, group.inverse(frame)));
            const auto transformed_image =
                group.multiply(frame, group.multiply(dynamics.image[value], group.inverse(frame)));
            if (dynamics.image[conjugated] != transformed_image) {
                throw std::invalid_argument("holonomy dynamics is not gauge equivariant");
            }
        }
    }
}

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
                const auto conjugated =
                    group.multiply(frame, group.multiply(value, group.inverse(frame)));
                const auto transformed_conjugate =
                    group.multiply(frame, group.multiply(candidate[value], group.inverse(frame)));
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

inline std::vector<Vertex> canonical_conjugacy_image(const std::vector<Permutation>& group,
                                                     const Permutation& value) {
    std::vector<Vertex> canonical;
    bool first = true;
    for (const auto& frame : group) {
        const auto conjugated =
            Permutation::compose(frame, Permutation::compose(value, frame.inverse()));
        if (first || conjugated.image() < canonical) {
            canonical = conjugated.image();
            first = false;
        }
    }
    return canonical;
}

inline std::vector<uint16_t> induced_conjugacy_sector_map(const AutomorphismTables& group,
                                                          const LocalGaugeDynamics& dynamics) {
    validate_local_gauge_dynamics(group, dynamics);
    const auto representatives = conjugacy_representatives(group.elements());
    std::map<std::vector<Vertex>, uint16_t> sector_by_code;
    for (uint16_t index = 0; index < representatives.size(); ++index) {
        sector_by_code.emplace(canonical_conjugacy_image(group.elements(), representatives[index]),
                               index);
    }
    std::vector<uint16_t> result;
    result.reserve(representatives.size());
    for (const auto& representative : representatives) {
        const auto element = group.index_of(representative);
        const auto code =
            canonical_conjugacy_image(group.elements(), group.elements()[dynamics.image[element]]);
        result.push_back(sector_by_code.at(code));
    }
    return result;
}

inline std::vector<Vertex> simple_cycle_loop(const BaseGraph& base) {
    if (base.vertices().size() < 3 || base.edges().size() != base.vertices().size()) {
        throw std::invalid_argument("cell update requires one simple cycle");
    }
    for (const auto vertex : base.vertices()) {
        if (base.neighbors(vertex).size() != 2) {
            throw std::invalid_argument("cell update requires degree two at every base vertex");
        }
    }
    const auto start = *base.vertices().begin();
    const auto start_neighbors = base.neighbors(start);
    Vertex previous = start;
    Vertex current = start_neighbors.front();
    std::vector<Vertex> loop{start, current};
    while (current != start) {
        const auto neighbors = base.neighbors(current);
        const auto next = neighbors.front() == previous ? neighbors.back() : neighbors.front();
        previous = current;
        current = next;
        loop.push_back(current);
        if (loop.size() > base.vertices().size() + 1) {
            throw std::logic_error("cycle traversal did not close");
        }
    }
    if (loop.size() != base.vertices().size() + 1) {
        throw std::invalid_argument("base contains more than one cycle component");
    }
    return loop;
}

// Apply a conjugation-equivariant holonomy rule to one cell. If H is the
// current based holonomy and phi(H) the target, multiplying the closing edge
// by phi(H) H^-1 changes exactly that loop transport. The correction transforms
// at the loop basepoint, so the edge update remains gauge covariant.
inline FiberBundleConnection apply_holonomy_dynamics(const FiberBundleConnection& connection,
                                                     const std::vector<Vertex>& loop,
                                                     const LocalGaugeDynamics& dynamics) {
    const AutomorphismTables tables(connection.local_gauge_group());
    validate_local_gauge_dynamics(tables, dynamics);
    const auto current = connection.holonomy(loop);
    const auto current_index = tables.index_of(current);
    const auto target_index = dynamics.image[current_index];
    if (target_index >= tables.order())
        throw std::out_of_range("holonomy dynamics target is absent");
    const auto& target = tables.elements()[target_index];
    const auto correction = Permutation::compose(target, current.inverse());
    const auto closing_from = loop[loop.size() - 2];
    const auto closing_to = loop.back();
    FiberBundleConnection result = connection;
    result.set_transport(
        closing_from, closing_to,
        Permutation::compose(correction, connection.edge_transport(closing_from, closing_to)));
    if (result.holonomy(loop) != target) {
        throw std::logic_error("cell update failed to realize its target holonomy");
    }
    return result;
}

enum class HurwitzDirection { Forward, Inverse };

inline bool loop_contains_edge(const std::vector<Vertex>& loop, Edge edge) {
    for (std::size_t index = 1; index < loop.size(); ++index) {
        if (infragauge::canonical_edge(loop[index - 1], loop[index]) == edge) return true;
    }
    return false;
}

// The elementary braid action on two adjacent, commonly based cells:
//   (A, B) -> (A B A^-1, A).
// It is reversible, commutes with simultaneous conjugation at the basepoint,
// and preserves the ordered total holonomy A B. Unlike unary equivariant maps,
// it can transport a nontrivial conjugacy sector from one cell to the other.
inline FiberBundleConnection apply_hurwitz_cell_exchange(
    const FiberBundleConnection& connection, const std::vector<Vertex>& first_loop,
    const std::vector<Vertex>& second_loop,
    HurwitzDirection direction = HurwitzDirection::Forward) {
    if (first_loop.size() < 4 || second_loop.size() < 4 ||
        first_loop.front() != first_loop.back() || second_loop.front() != second_loop.back() ||
        first_loop.front() != second_loop.front()) {
        throw std::invalid_argument(
            "Hurwitz exchange requires two nondegenerate loops at one basepoint");
    }
    const auto first_closing =
        infragauge::canonical_edge(first_loop[first_loop.size() - 2], first_loop.back());
    const auto second_closing =
        infragauge::canonical_edge(second_loop[second_loop.size() - 2], second_loop.back());
    if (first_closing == second_closing || loop_contains_edge(second_loop, first_closing) ||
        loop_contains_edge(first_loop, second_closing)) {
        throw std::invalid_argument(
            "Hurwitz realization requires one exclusive closing edge per cell");
    }

    const auto first = connection.holonomy(first_loop);
    const auto second = connection.holonomy(second_loop);
    const auto targets =
        direction == HurwitzDirection::Forward
            ? std::pair{Permutation::compose(first, Permutation::compose(second, first.inverse())),
                        first}
            : std::pair{second, Permutation::compose(second.inverse(),
                                                     Permutation::compose(first, second))};
    const auto& [target_first, target_second] = targets;

    FiberBundleConnection result = connection;
    const auto realize = [&](const std::vector<Vertex>& loop, const Permutation& current,
                             const Permutation& target) {
        const auto from = loop[loop.size() - 2];
        const auto to = loop.back();
        const auto correction = Permutation::compose(target, current.inverse());
        result.set_transport(from, to,
                             Permutation::compose(correction, connection.edge_transport(from, to)));
    };
    realize(first_loop, first, target_first);
    realize(second_loop, second, target_second);
    if (result.holonomy(first_loop) != target_first ||
        result.holonomy(second_loop) != target_second) {
        throw std::logic_error("Hurwitz exchange failed to realize its target holonomies");
    }
    const auto total_before = Permutation::compose(first, second);
    const auto total_after = Permutation::compose(target_first, target_second);
    if (total_before != total_after) {
        throw std::logic_error("Hurwitz exchange violated total holonomy conservation");
    }
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
        const auto discarded =
            1.0 - static_cast<double>(retained_rank) / static_cast<double>(automorphisms.size());
        result.push_back({graph.vertex_count(), graph.edges().size(), code, automorphisms.size(),
                          nonabelian, conjugacy_classes.size(), dynamics, automorphisms.size(), 1,
                          automorphisms.size(), discarded});
    }
    return result;
}

inline std::vector<Vertex> curvature_sector(const FiberBundleConnection& connection,
                                            const std::vector<std::vector<Vertex>>& loops) {
    std::vector<Vertex> result;
    const auto& group = connection.local_gauge_group();
    for (const auto& loop : loops) {
        const auto value = connection.holonomy(loop);
        const auto canonical = canonical_conjugacy_image(group, value);
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
        if (!by_id.emplace(event.id, event).second)
            throw std::invalid_argument("duplicate event id");
        if (event.before != event.after) changed.insert(event.id);
    }
    for (const auto source : declared_sources) {
        if (!changed.contains(source))
            throw std::invalid_argument("source event has no curvature change");
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
    report.causal_alignment =
        non_source == 0 ? 1.0 : static_cast<double>(report.causally_reached_changes) / non_source;
    return report;
}

}  // namespace wgphysics::research
