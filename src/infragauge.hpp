#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wgphysics::infragauge {

using Vertex = uint32_t;
using Edge = std::pair<Vertex, Vertex>;

inline Edge canonical_edge(Vertex u, Vertex v) {
    if (u == v) throw std::invalid_argument("self-edges are not supported");
    return u < v ? Edge{u, v} : Edge{v, u};
}

class Permutation {
public:
    explicit Permutation(std::vector<Vertex> image) : image_(std::move(image)) {
        auto sorted = image_;
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t i = 0; i < sorted.size(); ++i) {
            if (sorted[i] != i) throw std::invalid_argument("invalid permutation");
        }
    }

    static Permutation identity(std::size_t size) {
        std::vector<Vertex> image(size);
        std::iota(image.begin(), image.end(), 0);
        return Permutation(std::move(image));
    }

    std::size_t size() const { return image_.size(); }

    Vertex operator()(Vertex value) const {
        if (value >= image_.size()) throw std::out_of_range("permutation input outside fiber");
        return image_[value];
    }

    const std::vector<Vertex>& image() const { return image_; }

    Permutation inverse() const {
        std::vector<Vertex> result(size());
        for (Vertex i = 0; i < size(); ++i) result[image_[i]] = i;
        return Permutation(std::move(result));
    }

    // compose(after, before)(x) = after(before(x)).
    static Permutation compose(const Permutation& after, const Permutation& before) {
        require_same_size(after, before);
        std::vector<Vertex> result(after.size());
        for (Vertex i = 0; i < after.size(); ++i) result[i] = after(before(i));
        return Permutation(std::move(result));
    }

    std::vector<uint32_t> cycle_signature() const {
        std::vector<bool> visited(size(), false);
        std::vector<uint32_t> result;
        for (Vertex start = 0; start < size(); ++start) {
            if (visited[start]) continue;
            uint32_t length = 0;
            auto current = start;
            do {
                visited[current] = true;
                current = (*this)(current);
                ++length;
            } while (!visited[current]);
            result.push_back(length);
        }
        std::sort(result.begin(), result.end());
        return result;
    }

    bool operator==(const Permutation&) const = default;

private:
    std::vector<Vertex> image_;

    static void require_same_size(const Permutation& left, const Permutation& right) {
        if (left.size() != right.size()) throw std::invalid_argument("permutation sizes differ");
    }
};

class FiberGraph {
public:
    FiberGraph(uint32_t vertex_count, std::vector<Edge> edges)
        : vertex_count_(vertex_count) {
        if (vertex_count_ == 0) throw std::invalid_argument("fiber must contain a vertex");
        for (const auto [u, v] : edges) {
            if (u >= vertex_count_ || v >= vertex_count_) {
                throw std::out_of_range("fiber edge endpoint outside fiber");
            }
            edges_.insert(canonical_edge(u, v));
        }
    }

    uint32_t vertex_count() const { return vertex_count_; }
    const std::set<Edge>& edges() const { return edges_; }
    bool adjacent(Vertex u, Vertex v) const {
        if (u >= vertex_count_ || v >= vertex_count_ || u == v) return false;
        return edges_.contains(canonical_edge(u, v));
    }

    bool is_automorphism(const Permutation& permutation) const {
        if (permutation.size() != vertex_count_) return false;
        for (Vertex u = 0; u < vertex_count_; ++u) {
            for (Vertex v = u + 1; v < vertex_count_; ++v) {
                if (adjacent(u, v) != adjacent(permutation(u), permutation(v))) return false;
            }
        }
        return true;
    }

    // Exact enumeration is intentional. These fibers are microscopic and small;
    // exhaustive automorphisms make the emergent local group falsifiable.
    std::vector<Permutation> automorphisms() const {
        std::vector<Vertex> candidate(vertex_count_);
        std::iota(candidate.begin(), candidate.end(), 0);
        std::vector<Permutation> result;
        do {
            Permutation permutation(candidate);
            if (is_automorphism(permutation)) result.push_back(std::move(permutation));
        } while (std::next_permutation(candidate.begin(), candidate.end()));
        return result;
    }

private:
    uint32_t vertex_count_;
    std::set<Edge> edges_;
};

// Dense table representation of Aut(F). This is the device-facing form: a
// connection value is a small integer, and composition, inverse, and action are
// allocation-free table lookups rather than permutation manipulation.
class AutomorphismTables {
public:
    explicit AutomorphismTables(std::vector<Permutation> group) : group_(std::move(group)) {
        if (group_.empty()) throw std::invalid_argument("automorphism group cannot be empty");
        if (group_.size() > std::numeric_limits<uint16_t>::max()) {
            throw std::invalid_argument("automorphism group exceeds compact table index");
        }
        const auto fiber_size = group_.front().size();
        multiplication_.resize(group_.size() * group_.size());
        inverses_.resize(group_.size());
        action_.resize(group_.size() * fiber_size);
        for (std::size_t left = 0; left < group_.size(); ++left) {
            if (group_[left].size() != fiber_size) {
                throw std::invalid_argument("automorphisms act on different fibers");
            }
            inverses_[left] = index_of(group_[left].inverse());
            for (std::size_t value = 0; value < fiber_size; ++value) {
                action_[left * fiber_size + value] =
                    static_cast<uint16_t>(group_[left](static_cast<Vertex>(value)));
            }
            for (std::size_t right = 0; right < group_.size(); ++right) {
                multiplication_[left * group_.size() + right] =
                    index_of(Permutation::compose(group_[left], group_[right]));
            }
        }
    }

    uint16_t order() const { return static_cast<uint16_t>(group_.size()); }
    uint16_t fiber_size() const { return static_cast<uint16_t>(group_.front().size()); }
    const std::vector<Permutation>& elements() const { return group_; }
    const std::vector<uint16_t>& multiplication_data() const { return multiplication_; }
    const std::vector<uint16_t>& inverse_data() const { return inverses_; }
    const std::vector<uint16_t>& action_data() const { return action_; }
    uint16_t multiply(uint16_t left, uint16_t right) const {
        require_group_index(left);
        require_group_index(right);
        return multiplication_[static_cast<std::size_t>(left) * group_.size() + right];
    }
    uint16_t inverse(uint16_t value) const {
        require_group_index(value);
        return inverses_[value];
    }
    uint16_t act(uint16_t group_element, uint16_t fiber_vertex) const {
        require_group_index(group_element);
        if (fiber_vertex >= fiber_size()) throw std::out_of_range("fiber action input is absent");
        return action_[static_cast<std::size_t>(group_element) * fiber_size() + fiber_vertex];
    }
    uint16_t index_of(const Permutation& value) const {
        const auto found = std::find(group_.begin(), group_.end(), value);
        if (found == group_.end()) throw std::invalid_argument("permutation is outside group");
        return static_cast<uint16_t>(std::distance(group_.begin(), found));
    }

private:
    std::vector<Permutation> group_;
    std::vector<uint16_t> multiplication_;
    std::vector<uint16_t> inverses_;
    std::vector<uint16_t> action_;

    void require_group_index(uint16_t value) const {
        if (value >= group_.size()) throw std::out_of_range("group-table index is absent");
    }
};

class BaseGraph {
public:
    BaseGraph(std::vector<Vertex> vertices, std::vector<Edge> edges)
        : vertices_(vertices.begin(), vertices.end()) {
        if (vertices_.size() != vertices.size()) throw std::invalid_argument("duplicate base vertex");
        for (const auto [u, v] : edges) {
            if (!vertices_.contains(u) || !vertices_.contains(v)) {
                throw std::out_of_range("base edge endpoint is not a base vertex");
            }
            edges_.insert(canonical_edge(u, v));
        }
    }

    bool has_vertex(Vertex vertex) const { return vertices_.contains(vertex); }
    bool has_edge(Vertex u, Vertex v) const { return edges_.contains(canonical_edge(u, v)); }
    const std::set<Vertex>& vertices() const { return vertices_; }
    const std::set<Edge>& edges() const { return edges_; }

    std::vector<Vertex> neighbors(Vertex vertex) const {
        if (!has_vertex(vertex)) throw std::out_of_range("base vertex is absent");
        std::vector<Vertex> result;
        for (const auto [u, v] : edges_) {
            if (u == vertex) result.push_back(v);
            if (v == vertex) result.push_back(u);
        }
        return result;
    }

    BaseGraph subdivide(Vertex from, Vertex to, Vertex midpoint) const {
        if (!has_edge(from, to)) throw std::out_of_range("cannot subdivide a missing edge");
        if (has_vertex(midpoint)) throw std::invalid_argument("midpoint must be a fresh vertex");
        auto vertices = vertices_;
        auto edges = edges_;
        vertices.insert(midpoint);
        edges.erase(canonical_edge(from, to));
        edges.insert(canonical_edge(from, midpoint));
        edges.insert(canonical_edge(midpoint, to));
        return BaseGraph(
            std::vector<Vertex>(vertices.begin(), vertices.end()),
            std::vector<Edge>(edges.begin(), edges.end()));
    }

private:
    std::set<Vertex> vertices_;
    std::set<Edge> edges_;
};

class FiberBundleConnection {
public:
    FiberBundleConnection(BaseGraph base, FiberGraph fiber)
        : base_(std::move(base)), fiber_(std::move(fiber)),
          gauge_group_(fiber_.automorphisms()) {
        const auto identity = Permutation::identity(fiber_.vertex_count());
        for (const auto edge : base_.edges()) transports_.emplace(edge, identity);
    }

    const BaseGraph& base() const { return base_; }
    const FiberGraph& fiber() const { return fiber_; }
    const std::vector<Permutation>& local_gauge_group() const { return gauge_group_; }

    void set_transport(Vertex from, Vertex to, const Permutation& transport) {
        if (!base_.has_edge(from, to)) throw std::out_of_range("transport edge is absent");
        require_automorphism(transport);
        transports_.insert_or_assign(
            canonical_edge(from, to), from < to ? transport : transport.inverse());
    }

    Permutation edge_transport(Vertex from, Vertex to) const {
        if (!base_.has_edge(from, to)) throw std::out_of_range("transport edge is absent");
        const auto& stored = transports_.at(canonical_edge(from, to));
        return from < to ? stored : stored.inverse();
    }

    Permutation parallel_transport(const std::vector<Vertex>& path) const {
        if (path.empty()) throw std::invalid_argument("path cannot be empty");
        auto result = Permutation::identity(fiber_.vertex_count());
        for (std::size_t i = 1; i < path.size(); ++i) {
            result = Permutation::compose(edge_transport(path[i - 1], path[i]), result);
        }
        return result;
    }

    Permutation holonomy(const std::vector<Vertex>& loop) const {
        if (loop.size() < 2 || loop.front() != loop.back()) {
            throw std::invalid_argument("holonomy path must be a closed loop");
        }
        return parallel_transport(loop);
    }

    bool flat_on(const std::vector<std::vector<Vertex>>& loops) const {
        const auto identity = Permutation::identity(fiber_.vertex_count());
        return std::all_of(loops.begin(), loops.end(), [&](const auto& loop) {
            return holonomy(loop) == identity;
        });
    }

    // A change of local fiber frame at each base vertex. Connection transport
    // U_uv transforms as g_v U_uv g_u^-1; closed-loop holonomy is conjugated.
    FiberBundleConnection gauge_transform(
        const std::map<Vertex, Permutation>& local_frames) const {
        auto result = *this;
        const auto identity = Permutation::identity(fiber_.vertex_count());
        for (const auto& [base_vertex, frame] : local_frames) {
            if (!base_.has_vertex(base_vertex)) throw std::out_of_range("frame base vertex absent");
            require_automorphism(frame);
        }
        for (const auto [u, v] : base_.edges()) {
            const auto found_u = local_frames.find(u);
            const auto found_v = local_frames.find(v);
            const auto& g_u = found_u == local_frames.end() ? identity : found_u->second;
            const auto& g_v = found_v == local_frames.end() ? identity : found_v->second;
            const auto transformed = Permutation::compose(
                g_v, Permutation::compose(edge_transport(u, v), g_u.inverse()));
            result.set_transport(u, v, transformed);
        }
        return result;
    }

    bool covariantly_constant_section(const std::map<Vertex, Vertex>& section) const {
        if (section.size() != base_.vertices().size()) return false;
        for (const auto vertex : base_.vertices()) {
            const auto found = section.find(vertex);
            if (found == section.end() || found->second >= fiber_.vertex_count()) return false;
        }
        for (const auto [u, v] : base_.edges()) {
            if (edge_transport(u, v)(section.at(u)) != section.at(v)) return false;
        }
        return true;
    }

    std::vector<Vertex> holonomy_fixed_points(const std::vector<Vertex>& loop) const {
        const auto transport = holonomy(loop);
        std::vector<Vertex> result;
        for (Vertex value = 0; value < fiber_.vertex_count(); ++value) {
            if (transport(value) == value) result.push_back(value);
        }
        return result;
    }

    // Exact quotient by all local changes of fiber frame on a fixed base graph.
    // A deterministic spanning forest removes tree transports. Each chord gives
    // one root-based holonomy; simultaneous conjugation at a component root is
    // removed by choosing the lexicographically least representative.
    std::vector<Vertex> gauge_invariant_signature() const {
        std::set<Vertex> visited;
        std::vector<Vertex> signature;
        for (const auto root : base_.vertices()) {
            if (visited.contains(root)) continue;
            std::map<Vertex, Vertex> parent{{root, root}};
            std::vector<Vertex> queue{root};
            std::set<Edge> tree_edges;
            visited.insert(root);
            for (std::size_t position = 0; position < queue.size(); ++position) {
                const auto current = queue[position];
                for (const auto neighbor : base_.neighbors(current)) {
                    if (visited.insert(neighbor).second) {
                        parent.emplace(neighbor, current);
                        tree_edges.insert(canonical_edge(current, neighbor));
                        queue.push_back(neighbor);
                    }
                }
            }

            std::vector<Permutation> fundamental_holonomies;
            for (const auto [u, v] : base_.edges()) {
                if (!parent.contains(u) || !parent.contains(v) ||
                    tree_edges.contains(canonical_edge(u, v))) {
                    continue;
                }
                auto loop = path_from_root(parent, root, u);
                loop.push_back(v);
                const auto root_to_v = path_from_root(parent, root, v);
                for (std::size_t i = root_to_v.size() - 1; i > 0; --i) {
                    loop.push_back(root_to_v[i - 1]);
                }
                fundamental_holonomies.push_back(holonomy(loop));
            }

            std::vector<Vertex> component_signature;
            bool first_candidate = true;
            for (const auto& global_frame : gauge_group_) {
                std::vector<Vertex> candidate;
                for (const auto& value : fundamental_holonomies) {
                    const auto conjugated = Permutation::compose(
                        global_frame,
                        Permutation::compose(value, global_frame.inverse()));
                    candidate.insert(
                        candidate.end(), conjugated.image().begin(), conjugated.image().end());
                }
                if (first_candidate || candidate < component_signature) {
                    component_signature = std::move(candidate);
                    first_candidate = false;
                }
            }
            // Delimit components and retain how many independent cycles each has.
            signature.push_back(static_cast<Vertex>(fundamental_holonomies.size()));
            signature.insert(signature.end(), component_signature.begin(), component_signature.end());
        }
        return signature;
    }

    // Extend the base rewrite (u--v) -> (u--w--v). Every factorization
    // U_uv = U_wv U_uw is enumerated from Aut(F). The resulting connections
    // are one orbit under changes of frame at the fresh fiber w.
    std::vector<FiberBundleConnection> subdivide_edge_extensions(
        Vertex from, Vertex to, Vertex midpoint) const {
        const auto total_transport = edge_transport(from, to);
        const auto new_base = base_.subdivide(from, to, midpoint);
        std::vector<FiberBundleConnection> result;
        result.reserve(gauge_group_.size());
        for (const auto& first : gauge_group_) {
            FiberBundleConnection extension(new_base, fiber_);
            for (const auto [u, v] : base_.edges()) {
                if (canonical_edge(u, v) != canonical_edge(from, to)) {
                    extension.set_transport(u, v, edge_transport(u, v));
                }
            }
            const auto second = Permutation::compose(total_transport, first.inverse());
            extension.set_transport(from, midpoint, first);
            extension.set_transport(midpoint, to, second);
            result.push_back(std::move(extension));
        }
        return result;
    }

    // One representative of the physical subdivision orbit. The fresh frame
    // is fixed to identity, so callers that quotient gauge immediately never
    // need to materialize |Aut(F)| equivalent connection copies.
    FiberBundleConnection subdivide_edge_representative(
        Vertex from, Vertex to, Vertex midpoint) const {
        const auto total_transport = edge_transport(from, to);
        FiberBundleConnection result(base_.subdivide(from, to, midpoint), fiber_);
        for (const auto [u, v] : base_.edges()) {
            if (canonical_edge(u, v) != canonical_edge(from, to)) {
                result.set_transport(u, v, edge_transport(u, v));
            }
        }
        result.set_transport(from, midpoint, Permutation::identity(fiber_.vertex_count()));
        result.set_transport(midpoint, to, total_transport);
        return result;
    }

private:
    BaseGraph base_;
    FiberGraph fiber_;
    std::vector<Permutation> gauge_group_;
    std::map<Edge, Permutation> transports_;

    void require_automorphism(const Permutation& permutation) const {
        if (!fiber_.is_automorphism(permutation)) {
            throw std::invalid_argument("connection/frame is not a fiber automorphism");
        }
    }

    static std::vector<Vertex> path_from_root(
        const std::map<Vertex, Vertex>& parent, Vertex root, Vertex target) {
        std::vector<Vertex> reversed{target};
        while (target != root) {
            target = parent.at(target);
            reversed.push_back(target);
        }
        return std::vector<Vertex>(reversed.rbegin(), reversed.rend());
    }
};

// Quantum amplitude over connection factorizations created by an edge
// subdivision. The basis is Aut(F) x Aut(F), derived from the fiber graph.
class SubdivisionWavefunction {
public:
    using Amplitude = std::complex<double>;

    static SubdivisionWavefunction from_total_transport(
        std::vector<Permutation> group, const Permutation& total_transport) {
        if (group.empty()) throw std::invalid_argument("local gauge group cannot be empty");
        SubdivisionWavefunction result(std::move(group));
        const auto magnitude = 1.0 / std::sqrt(static_cast<double>(result.group_.size()));
        for (std::size_t first = 0; first < result.group_.size(); ++first) {
            const auto second_value = Permutation::compose(
                total_transport, result.group_[first].inverse());
            const auto second = result.group_index(second_value);
            result.amplitudes_[result.index(first, second)] = magnitude;
        }
        return result;
    }

    std::size_t group_size() const { return group_.size(); }

    double norm_squared() const {
        double result = 0.0;
        for (const auto amplitude : amplitudes_) result += std::norm(amplitude);
        return result;
    }

    std::size_t support_size(double tolerance = 1e-12) const {
        return static_cast<std::size_t>(std::count_if(
            amplitudes_.begin(), amplitudes_.end(),
            [=](const auto amplitude) { return std::abs(amplitude) > tolerance; }));
    }

    // The coefficient matrix for a subdivision orbit is monomial: each first
    // transport determines exactly one second transport. Its nonzero entry
    // magnitudes are therefore the exact Schmidt singular values. Computing
    // them from the stored amplitudes keeps this diagnostic honest if the
    // construction changes later.
    std::vector<double> schmidt_spectrum(double tolerance = 1e-12) const {
        std::vector<double> spectrum;
        std::vector<std::size_t> column_counts(group_.size(), 0);
        for (std::size_t row = 0; row < group_.size(); ++row) {
            std::size_t row_count = 0;
            for (std::size_t column = 0; column < group_.size(); ++column) {
                const auto magnitude = std::abs(amplitudes_[index(row, column)]);
                if (magnitude <= tolerance) continue;
                ++row_count;
                ++column_counts[column];
                spectrum.push_back(magnitude);
            }
            if (row_count != 1) {
                throw std::logic_error("subdivision state is no longer a monomial orbit matrix");
            }
        }
        if (std::any_of(column_counts.begin(), column_counts.end(), [](const auto count) {
                return count != 1;
            })) {
            throw std::logic_error("subdivision orbit matrix is not bijective");
        }
        std::sort(spectrum.begin(), spectrum.end(), std::greater<>());
        return spectrum;
    }

    double best_rank_discarded_norm(std::size_t retained_rank) const {
        const auto spectrum = schmidt_spectrum();
        if (retained_rank > spectrum.size()) {
            throw std::out_of_range("retained Schmidt rank exceeds exact rank");
        }
        double discarded = 0.0;
        for (std::size_t i = retained_rank; i < spectrum.size(); ++i) {
            discarded += spectrum[i] * spectrum[i];
        }
        return discarded;
    }

    SubdivisionWavefunction gauge_transform_midpoint(const Permutation& frame) const {
        const auto frame_index = group_index(frame);
        (void)frame_index;
        SubdivisionWavefunction result(group_);
        for (std::size_t first = 0; first < group_.size(); ++first) {
            for (std::size_t second = 0; second < group_.size(); ++second) {
                const auto transformed_first = group_index(
                    Permutation::compose(frame, group_[first]));
                const auto transformed_second = group_index(
                    Permutation::compose(group_[second], frame.inverse()));
                result.amplitudes_[result.index(transformed_first, transformed_second)] +=
                    amplitudes_[index(first, second)];
            }
        }
        return result;
    }

    Amplitude inner_product(const SubdivisionWavefunction& other) const {
        require_same_group(other);
        Amplitude result = 0.0;
        for (std::size_t i = 0; i < amplitudes_.size(); ++i) {
            result += std::conj(amplitudes_[i]) * other.amplitudes_[i];
        }
        return result;
    }

    bool approximately_equal(const SubdivisionWavefunction& other, double tolerance = 1e-12) const {
        if (!same_group(other)) return false;
        for (std::size_t i = 0; i < amplitudes_.size(); ++i) {
            if (std::abs(amplitudes_[i] - other.amplitudes_[i]) > tolerance) return false;
        }
        return true;
    }

private:
    explicit SubdivisionWavefunction(std::vector<Permutation> group)
        : group_(std::move(group)), amplitudes_(group_.size() * group_.size(), 0.0) {}

    std::vector<Permutation> group_;
    std::vector<Amplitude> amplitudes_;

    std::size_t index(std::size_t first, std::size_t second) const {
        return first * group_.size() + second;
    }

    std::size_t group_index(const Permutation& value) const {
        const auto found = std::find(group_.begin(), group_.end(), value);
        if (found == group_.end()) throw std::invalid_argument("permutation is outside gauge group");
        return static_cast<std::size_t>(std::distance(group_.begin(), found));
    }

    bool same_group(const SubdivisionWavefunction& other) const {
        return group_ == other.group_;
    }

    void require_same_group(const SubdivisionWavefunction& other) const {
        if (!same_group(other)) throw std::invalid_argument("wavefunctions use different groups");
    }
};

}  // namespace wgphysics::infragauge
