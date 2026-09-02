#include "research.hpp"

#include <iostream>

int main() {
    using namespace wgphysics::infragauge;
    using namespace wgphysics::research;

    const BaseGraph wheel(
        {0, 1, 2, 3, 4},
        {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 3}, {3, 4}, {4, 1}});
    const auto inferred = infer_link_fibers(wheel);
    const auto& fiber = std::find_if(inferred.begin(), inferred.end(), [](const auto& item) {
        return item.representative.vertex_count() == 4 && item.automorphism_order == 8;
    })->representative;

    const BaseGraph triangle({0, 1, 2}, {{0, 1}, {1, 2}, {2, 0}});
    FiberBundleConnection connection(triangle, fiber);
    connection.set_transport(2, 0, Permutation({1, 2, 3, 0}));
    GaugeSubdivisionEvolution evolution(connection);
    const auto children = evolution.subdivide(0, 0, 1, 3);

    const AutomorphismTables tables(fiber.automorphisms());
    const auto dynamics = search_local_gauge_dynamics(tables);
    const auto orbit = SubdivisionWavefunction::from_total_transport(
        tables.elements(), Permutation::identity(fiber.vertex_count()));

    std::cout << "inferred C4 fiber from link at base vertex 0\n"
              << "derived automorphism order: " << tables.order() << '\n'
              << "raw subdivision extensions: " << evolution.events().front().raw_extensions << '\n'
              << "physical subdivision children: " << children.size() << '\n'
              << "reversible equivariant local maps: " << dynamics.size() << '\n'
              << "exact subdivision Schmidt rank: " << orbit.schmidt_spectrum().size() << '\n'
              << "norm lost at half rank: "
              << orbit.best_rank_discarded_norm(orbit.group_size() / 2) << '\n';
}
