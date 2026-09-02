#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

#include "research.hpp"

namespace {

void print_permutation(const wgphysics::infragauge::Permutation& value) {
    std::cout << '[';
    for (std::size_t index = 0; index < value.image().size(); ++index) {
        if (index) std::cout << ',';
        std::cout << value.image()[index];
    }
    std::cout << ']';
}

}  // namespace

int main() {
    using namespace wgphysics::infragauge;
    using namespace wgphysics::research;

    const FiberGraph fiber(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    const BaseGraph base({0, 1, 2, 3}, {{0, 1}, {1, 2}, {2, 0}, {1, 3}, {3, 0}});
    const std::vector<Vertex> first_cell{0, 1, 2, 0};
    const std::vector<Vertex> second_cell{0, 1, 3, 0};
    const OrientedCellComplex complex(base, {first_cell, second_cell});
    FiberBundleConnection initial(base, fiber);
    initial.set_transport(2, 0, Permutation({1, 2, 3, 0}));

    const auto evolved = apply_hurwitz_cell_exchange(initial, first_cell, second_cell);
    const auto restored =
        apply_hurwitz_cell_exchange(evolved, first_cell, second_cell, HurwitzDirection::Inverse);
    bool inverse_exact = true;
    for (const auto [from, to] : base.edges()) {
        inverse_exact &= restored.edge_transport(from, to) == initial.edge_transport(from, to);
    }
    const auto total_before =
        Permutation::compose(initial.holonomy(first_cell), initial.holonomy(second_cell));
    const auto total_after =
        Permutation::compose(evolved.holonomy(first_cell), evolved.holonomy(second_cell));
    const auto subdivided = complex.subdivide_edge(0, 1, 4);

    std::cout << "derived group order: " << fiber.automorphisms().size() << '\n';
    std::cout << "cell holonomies before: ";
    print_permutation(initial.holonomy(first_cell));
    std::cout << ' ';
    print_permutation(initial.holonomy(second_cell));
    std::cout << "\ncell holonomies after:  ";
    print_permutation(evolved.holonomy(first_cell));
    std::cout << ' ';
    print_permutation(evolved.holonomy(second_cell));
    std::cout << "\ntotal holonomy conserved: " << std::boolalpha << (total_before == total_after)
              << '\n';
    std::cout << "inverse restores links:    " << inverse_exact << '\n';
    std::cout << "incident faces updated:   "
              << static_cast<std::size_t>(
                     std::count_if(subdivided.faces().begin(), subdivided.faces().end(),
                                   [](const auto& face) {
                                       return std::find(face.begin(), face.end(), 4) != face.end();
                                   }))
              << '\n';
    return 0;
}
