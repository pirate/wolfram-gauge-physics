#pragma once

#include "research.hpp"

namespace wgphysics::research {

// These are explicit candidate word automorphisms, not inferred equations of motion.
enum class PairInteraction { Hurwitz, BoundaryShear };

inline std::pair<Permutation, Permutation> pair_targets(
    const Permutation& a, const Permutation& b, PairInteraction rule,
    HurwitzDirection direction = HurwitzDirection::Forward) {
    const auto mul = Permutation::compose;
    switch (rule) {
        case PairInteraction::Hurwitz:
            if (direction == HurwitzDirection::Forward)
                return {mul(a, mul(b, a.inverse())), a};
            return {b, mul(b.inverse(), mul(a, b))};
        case PairInteraction::BoundaryShear:
            // Preserve P = A B by taking A' = P A and B' = A^-1.
            // Recover A = B'^-1 and B = B' A' B'.
            if (direction == HurwitzDirection::Forward)
                return {mul(mul(a, b), a), a.inverse()};
            return {b.inverse(), mul(b, mul(a, b))};
    }
    throw std::invalid_argument("unknown pair interaction");
}

struct InteractionSupport {
    std::set<Edge> reads;
    std::set<Edge> writes;
};

inline bool independent_supports(const InteractionSupport& a, const InteractionSupport& b) {
    for (const auto& edge : a.writes)
        if (b.reads.contains(edge) || b.writes.contains(edge)) return false;
    for (const auto& edge : b.writes)
        if (a.reads.contains(edge) || a.writes.contains(edge)) return false;
    return true;
}

struct CellPairPatch {
    std::vector<Vertex> first_loop;
    std::vector<Vertex> second_loop;
    // Directed from the SECOND basepoint to the FIRST; part of the operator.
    std::vector<Vertex> connector;
};

inline InteractionSupport interaction_support(const FiberBundleConnection& connection,
                                               const CellPairPatch& patch) {
    const auto validate_loop = [&](const std::vector<Vertex>& loop) {
        if (loop.size() < 4 || loop.front() != loop.back())
            throw std::invalid_argument("interaction requires simple closed cell loops");
        const std::set<Vertex> distinct(loop.begin(), loop.end() - 1);
        if (distinct.size() != loop.size() - 1)
            throw std::invalid_argument("cell loop repeats a vertex");
    };
    validate_loop(patch.first_loop);
    validate_loop(patch.second_loop);
    if (patch.connector.empty() || patch.connector.front() != patch.second_loop.front() ||
        patch.connector.back() != patch.first_loop.front())
        throw std::invalid_argument("connector must run from second basepoint to first");
    InteractionSupport support;
    const auto add_path = [&](const std::vector<Vertex>& path) {
        for (std::size_t i = 1; i < path.size(); ++i) {
            if (!connection.base().has_edge(path[i - 1], path[i]))
                throw std::invalid_argument("interaction path leaves the base graph");
            support.reads.insert(infragauge::canonical_edge(path[i - 1], path[i]));
        }
    };
    add_path(patch.first_loop);
    add_path(patch.second_loop);
    add_path(patch.connector);
    const auto closing = [](const std::vector<Vertex>& loop) {
        return infragauge::canonical_edge(loop[loop.size() - 2], loop.back());
    };
    const auto first = closing(patch.first_loop), second = closing(patch.second_loop);
    if (loop_contains_edge(patch.second_loop, first) ||
        loop_contains_edge(patch.first_loop, second) ||
        loop_contains_edge(patch.connector, first) ||
        loop_contains_edge(patch.connector, second))
        throw std::invalid_argument("closing links must be exclusive and off the connector");
    support.writes = {first, second};
    return support;
}

inline FiberBundleConnection apply_pair_interaction(
    const FiberBundleConnection& connection, const CellPairPatch& patch,
    PairInteraction rule, HurwitzDirection direction = HurwitzDirection::Forward) {
    (void)interaction_support(connection, patch);
    const auto mul = Permutation::compose;
    const auto t = connection.parallel_transport(patch.connector);
    const auto a = connection.holonomy(patch.first_loop);
    const auto b = connection.holonomy(patch.second_loop);
    const auto based_b = mul(t, mul(b, t.inverse()));
    const auto [target_a, target_based_b] = pair_targets(a, based_b, rule, direction);
    const auto target_b = mul(t.inverse(), mul(target_based_b, t));
    FiberBundleConnection result = connection;
    const auto realize = [&](const std::vector<Vertex>& loop, const Permutation& before,
                             const Permutation& after) {
        const auto from = loop[loop.size() - 2], to = loop.back();
        result.set_transport(from, to, mul(mul(after, before.inverse()),
                                          connection.edge_transport(from, to)));
    };
    realize(patch.first_loop, a, target_a);
    realize(patch.second_loop, b, target_b);
    if (result.holonomy(patch.first_loop) != target_a ||
        result.holonomy(patch.second_loop) != target_b ||
        result.parallel_transport(patch.connector) != t ||
        mul(target_a, target_based_b) != mul(a, based_b))
        throw std::logic_error("pair interaction violated its exact postconditions");
    return result;
}

}  // namespace wgphysics::research
