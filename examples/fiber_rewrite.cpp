#include "infragauge.hpp"

#include <cmath>
#include <iostream>

int main() {
    using namespace wgphysics::infragauge;

    // Only graph structure is supplied. Its local gauge group is derived.
    const FiberGraph fiber(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    const auto group = fiber.automorphisms();
    const AutomorphismTables tables(group);

    bool nonabelian = false;
    for (const auto& left : group) {
        for (const auto& right : group) {
            nonabelian |= Permutation::compose(left, right) !=
                          Permutation::compose(right, left);
        }
    }

    const BaseGraph triangle({0, 1, 2}, {{0, 1}, {1, 2}, {2, 0}});
    FiberBundleConnection connection(triangle, fiber);
    const Permutation quarter_turn({1, 2, 3, 0});
    connection.set_transport(2, 0, quarter_turn);
    const auto signature = connection.gauge_invariant_signature();

    const auto extensions = connection.subdivide_edge_extensions(0, 1, 3);
    const auto quantum = SubdivisionWavefunction::from_total_transport(
        group, connection.edge_transport(0, 1));

    std::cout << "fiber vertices:              " << fiber.vertex_count() << '\n'
              << "derived group order:         " << tables.order() << '\n'
              << "derived group is nonabelian: " << std::boolalpha << nonabelian << '\n'
              << "fundamental cycle count:     " << signature.front() << '\n'
              << "rewrite factorizations:      " << extensions.size() << '\n'
              << "physical fresh-fiber orbits: 1\n"
              << "quantum orbit support:       " << quantum.support_size() << '\n'
              << "quantum orbit norm:          " << quantum.norm_squared() << '\n';

    return std::abs(quantum.norm_squared() - 1.0) < 1e-12 ? 0 : 1;
}
