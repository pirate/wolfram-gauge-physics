#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "infragauge.hpp"

namespace wgphysics::observables {

using infragauge::BaseGraph;
using infragauge::FiberBundleConnection;
using infragauge::Permutation;
using infragauge::Vertex;

struct ProbeConfig {
    std::size_t maximum_radius{8};
    std::size_t diffusion_steps{16};
    std::size_t maximum_sources{64};
    std::size_t minimum_plateau_scales{3};
    double maximum_plateau_relative_span{0.35};
};

struct ScalePlateau {
    std::size_t first_scale{};
    std::size_t last_scale{};
    double mean{};
    double relative_span{};
};

struct IntrinsicGeometryProfile {
    std::size_t sampled_sources{};
    std::vector<double> mean_shell_volume;
    std::vector<double> mean_ball_volume;
    std::vector<double> ball_volume_coefficient_of_variation;
    std::vector<std::optional<double>> volume_dimension;
    std::vector<double> lazy_return_probability;
    std::vector<std::optional<double>> spectral_dimension;
    std::optional<ScalePlateau> volume_dimension_plateau;
    std::optional<ScalePlateau> spectral_dimension_plateau;
};

namespace detail {

struct IndexedGraph {
    std::vector<Vertex> vertices;
    std::map<Vertex, std::size_t> index;
    std::vector<std::vector<std::size_t>> neighbors;
};

inline IndexedGraph index_graph(const BaseGraph& graph) {
    IndexedGraph result;
    result.vertices.assign(graph.vertices().begin(), graph.vertices().end());
    result.neighbors.resize(result.vertices.size());
    for (std::size_t i = 0; i < result.vertices.size(); ++i) {
        result.index.emplace(result.vertices[i], i);
    }
    for (const auto [left, right] : graph.edges()) {
        const auto u = result.index.at(left);
        const auto v = result.index.at(right);
        result.neighbors[u].push_back(v);
        result.neighbors[v].push_back(u);
    }
    return result;
}

inline std::vector<std::size_t> deterministic_sources(std::size_t vertex_count,
                                                       std::size_t maximum_sources) {
    if (maximum_sources == 0) throw std::invalid_argument("maximum probe sources must be positive");
    const auto count = std::min(vertex_count, maximum_sources);
    std::vector<std::size_t> result;
    result.reserve(count);
    for (std::size_t sample = 0; sample < count; ++sample) {
        result.push_back((sample * vertex_count) / count);
    }
    return result;
}

inline std::vector<std::size_t> distances_from(const IndexedGraph& graph, std::size_t source) {
    const auto unreachable = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> distance(graph.vertices.size(), unreachable);
    std::queue<std::size_t> frontier;
    distance[source] = 0;
    frontier.push(source);
    while (!frontier.empty()) {
        const auto current = frontier.front();
        frontier.pop();
        for (const auto neighbor : graph.neighbors[current]) {
            if (distance[neighbor] != unreachable) continue;
            distance[neighbor] = distance[current] + 1;
            frontier.push(neighbor);
        }
    }
    return distance;
}

inline std::optional<ScalePlateau> detect_plateau(
    const std::vector<std::optional<double>>& values, std::size_t minimum_scales,
    double maximum_relative_span) {
    if (minimum_scales == 0 || maximum_relative_span < 0.0) {
        throw std::invalid_argument("invalid plateau detector configuration");
    }
    std::optional<ScalePlateau> best;
    for (std::size_t first = 0; first < values.size(); ++first) {
        if (!values[first] || !std::isfinite(*values[first])) continue;
        double sum = 0.0;
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        for (std::size_t last = first; last < values.size(); ++last) {
            if (!values[last] || !std::isfinite(*values[last])) break;
            sum += *values[last];
            minimum = std::min(minimum, *values[last]);
            maximum = std::max(maximum, *values[last]);
            const auto length = last - first + 1;
            if (length < minimum_scales) continue;
            const auto mean = sum / static_cast<double>(length);
            if (std::abs(mean) < 1e-12) continue;
            const auto relative_span = (maximum - minimum) / std::abs(mean);
            if (relative_span > maximum_relative_span) continue;
            if (!best || length > best->last_scale - best->first_scale + 1 ||
                (length == best->last_scale - best->first_scale + 1 &&
                 relative_span < best->relative_span)) {
                best = ScalePlateau{first, last, mean, relative_span};
            }
        }
    }
    return best;
}

}  // namespace detail

// All geometry here is intrinsic. Ball volumes use shortest-path distance in
// the graph, and diffusion uses a lazy random walk. No embedding coordinates,
// target dimension, or expected continuum shape enter the calculation.
inline IntrinsicGeometryProfile probe_intrinsic_geometry(const BaseGraph& graph,
                                                          const ProbeConfig& config = {}) {
    const auto indexed = detail::index_graph(graph);
    IntrinsicGeometryProfile result;
    result.mean_shell_volume.assign(config.maximum_radius + 1, 0.0);
    result.mean_ball_volume.assign(config.maximum_radius + 1, 0.0);
    result.ball_volume_coefficient_of_variation.assign(config.maximum_radius + 1, 0.0);
    result.volume_dimension.assign(config.maximum_radius + 1, std::nullopt);
    result.lazy_return_probability.assign(config.diffusion_steps + 1, 0.0);
    result.spectral_dimension.assign(config.diffusion_steps + 1, std::nullopt);
    if (indexed.vertices.empty()) return result;

    const auto sources = detail::deterministic_sources(indexed.vertices.size(),
                                                        config.maximum_sources);
    result.sampled_sources = sources.size();
    std::vector<std::vector<double>> ball_samples(config.maximum_radius + 1);
    for (auto& samples : ball_samples) samples.reserve(sources.size());

    for (const auto source : sources) {
        const auto distances = detail::distances_from(indexed, source);
        std::vector<double> shell(config.maximum_radius + 1, 0.0);
        for (const auto distance : distances) {
            if (distance <= config.maximum_radius) shell[distance] += 1.0;
        }
        double ball = 0.0;
        for (std::size_t radius = 0; radius <= config.maximum_radius; ++radius) {
            ball += shell[radius];
            result.mean_shell_volume[radius] += shell[radius];
            result.mean_ball_volume[radius] += ball;
            ball_samples[radius].push_back(ball);
        }

        std::vector<double> probability(indexed.vertices.size(), 0.0);
        probability[source] = 1.0;
        result.lazy_return_probability[0] += 1.0;
        for (std::size_t step = 1; step <= config.diffusion_steps; ++step) {
            std::vector<double> next(indexed.vertices.size(), 0.0);
            for (std::size_t vertex = 0; vertex < indexed.vertices.size(); ++vertex) {
                if (indexed.neighbors[vertex].empty()) {
                    next[vertex] += probability[vertex];
                    continue;
                }
                next[vertex] += 0.5 * probability[vertex];
                const auto share = 0.5 * probability[vertex] /
                                   static_cast<double>(indexed.neighbors[vertex].size());
                for (const auto neighbor : indexed.neighbors[vertex]) next[neighbor] += share;
            }
            probability = std::move(next);
            result.lazy_return_probability[step] += probability[source];
        }
    }

    const auto denominator = static_cast<double>(sources.size());
    for (std::size_t radius = 0; radius <= config.maximum_radius; ++radius) {
        result.mean_shell_volume[radius] /= denominator;
        result.mean_ball_volume[radius] /= denominator;
        double squared_error = 0.0;
        for (const auto value : ball_samples[radius]) {
            const auto delta = value - result.mean_ball_volume[radius];
            squared_error += delta * delta;
        }
        const auto standard_deviation = std::sqrt(squared_error / denominator);
        if (result.mean_ball_volume[radius] > 0.0) {
            result.ball_volume_coefficient_of_variation[radius] =
                standard_deviation / result.mean_ball_volume[radius];
        }
        if (radius >= 2 && result.mean_ball_volume[radius - 1] > 0.0 &&
            result.mean_ball_volume[radius] > result.mean_ball_volume[radius - 1]) {
            result.volume_dimension[radius] =
                std::log(result.mean_ball_volume[radius] /
                         result.mean_ball_volume[radius - 1]) /
                std::log(static_cast<double>(radius) / static_cast<double>(radius - 1));
        }
    }

    for (auto& probability : result.lazy_return_probability) probability /= denominator;
    for (std::size_t step = 2; step + 1 <= config.diffusion_steps; ++step) {
        const auto before = result.lazy_return_probability[step - 1];
        const auto after = result.lazy_return_probability[step + 1];
        if (before > 0.0 && after > 0.0) {
            result.spectral_dimension[step] =
                -2.0 * std::log(after / before) /
                std::log(static_cast<double>(step + 1) / static_cast<double>(step - 1));
        }
    }

    result.volume_dimension_plateau = detail::detect_plateau(
        result.volume_dimension, config.minimum_plateau_scales,
        config.maximum_plateau_relative_span);
    result.spectral_dimension_plateau = detail::detect_plateau(
        result.spectral_dimension, config.minimum_plateau_scales,
        config.maximum_plateau_relative_span);
    return result;
}

struct GaugeLoopObservable {
    std::size_t edge_length{};
    std::vector<std::uint32_t> conjugacy_signature;
    double normalized_permutation_character{};
    double curvature_weight{};
    bool trivial{};
};

// The normalized trace of the permutation representation is a class function,
// hence invariant under local frame changes that conjugate closed holonomy.
inline GaugeLoopObservable probe_gauge_loop(const FiberBundleConnection& connection,
                                            const std::vector<Vertex>& loop) {
    if (loop.size() < 2 || loop.front() != loop.back()) {
        throw std::invalid_argument("gauge probe requires a closed loop");
    }
    const auto holonomy = connection.holonomy(loop);
    std::size_t fixed_points = 0;
    for (Vertex value = 0; value < holonomy.size(); ++value) {
        if (holonomy(value) == value) ++fixed_points;
    }
    const auto character = static_cast<double>(fixed_points) /
                           static_cast<double>(connection.fiber().vertex_count());
    // Cycle type classifies conjugacy in the full symmetric group, which can
    // merge distinct sectors of the smaller derived group Aut(F).
    auto sector = holonomy.image();
    for (const auto& frame : connection.local_gauge_group()) {
        const auto conjugated = Permutation::compose(
            frame, Permutation::compose(holonomy, frame.inverse()));
        sector = std::min(sector, conjugated.image());
    }
    return {loop.size() - 1, sector, character, 1.0 - character,
            holonomy == Permutation::identity(connection.fiber().vertex_count())};
}

struct GaugeFieldProfile {
    std::vector<GaugeLoopObservable> loops;
    std::size_t nontrivial_loops{};
    double curvature_density{};
    double mean_normalized_character{};
    double effective_curvature_support{};
    double normalized_localization{};
};

inline GaugeFieldProfile probe_gauge_field(
    const FiberBundleConnection& connection,
    const std::vector<std::vector<Vertex>>& loops) {
    GaugeFieldProfile result;
    result.loops.reserve(loops.size());
    double weight_sum = 0.0;
    double weight_squared_sum = 0.0;
    for (const auto& loop : loops) {
        auto observable = probe_gauge_loop(connection, loop);
        if (!observable.trivial) ++result.nontrivial_loops;
        result.mean_normalized_character += observable.normalized_permutation_character;
        weight_sum += observable.curvature_weight;
        weight_squared_sum += observable.curvature_weight * observable.curvature_weight;
        result.loops.push_back(std::move(observable));
    }
    if (loops.empty()) return result;
    result.curvature_density = static_cast<double>(result.nontrivial_loops) /
                               static_cast<double>(loops.size());
    result.mean_normalized_character /= static_cast<double>(loops.size());
    if (weight_squared_sum > 0.0) {
        result.effective_curvature_support = weight_sum * weight_sum / weight_squared_sum;
        if (loops.size() > 1) {
            result.normalized_localization =
                (static_cast<double>(loops.size()) / result.effective_curvature_support - 1.0) /
                (static_cast<double>(loops.size()) - 1.0);
        } else {
            result.normalized_localization = 1.0;
        }
    }
    return result;
}

struct LoopLengthBin {
    std::size_t edge_length{};
    std::size_t samples{};
    double mean_character{};
    double character_variance{};
};

inline std::vector<LoopLengthBin> aggregate_loop_characters(
    const std::vector<GaugeLoopObservable>& observations) {
    std::map<std::size_t, std::vector<double>> by_length;
    for (const auto& observation : observations) {
        by_length[observation.edge_length].push_back(
            observation.normalized_permutation_character);
    }
    std::vector<LoopLengthBin> result;
    result.reserve(by_length.size());
    for (const auto& [length, values] : by_length) {
        const auto mean = std::accumulate(values.begin(), values.end(), 0.0) /
                          static_cast<double>(values.size());
        double variance = 0.0;
        for (const auto value : values) {
            const auto delta = value - mean;
            variance += delta * delta;
        }
        variance /= static_cast<double>(values.size());
        result.push_back({length, values.size(), mean, variance});
    }
    return result;
}

}  // namespace wgphysics::observables
