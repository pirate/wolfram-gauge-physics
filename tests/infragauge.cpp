#include "infragauge.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {

using wgphysics::infragauge::Permutation;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Permutation permutation(std::initializer_list<uint32_t> image) {
    return Permutation(std::vector<uint32_t>(image));
}

}  // namespace

int main() {
    try {
        using namespace wgphysics::infragauge;

        // No gauge group is declared. A square fiber is supplied as bare graph
        // structure; its eight-element, nonabelian D4 automorphism group follows.
        const FiberGraph fiber(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        const auto group = fiber.automorphisms();
        require(group.size() == 8, "square fiber did not derive its full automorphism group");
        bool found_noncommuting_pair = false;
        for (const auto& left : group) {
            for (const auto& right : group) {
                if (Permutation::compose(left, right) != Permutation::compose(right, left)) {
                    found_noncommuting_pair = true;
                }
            }
        }
        require(found_noncommuting_pair, "derived square-fiber group was unexpectedly abelian");
        const AutomorphismTables tables(group);
        require(tables.order() == 8 && tables.fiber_size() == 4,
                "compact automorphism tables have the wrong dimensions");

        const BaseGraph triangle({0, 1, 2}, {{0, 1}, {1, 2}, {2, 0}});
        FiberBundleConnection connection(triangle, fiber);
        const auto quarter_turn = permutation({1, 2, 3, 0});
        const auto reflection = permutation({0, 3, 2, 1});
        connection.set_transport(2, 0, quarter_turn);

        const std::vector<uint32_t> loop{0, 1, 2, 0};
        const auto holonomy = connection.holonomy(loop);
        require(holonomy == quarter_turn, "parallel transport composed in the wrong order");
        require(holonomy.cycle_signature() == std::vector<uint32_t>{4},
                "quarter-turn holonomy had the wrong conjugacy invariant");
        require(!connection.flat_on({loop}), "curved connection was reported flat");
        require(connection.holonomy_fixed_points(loop).empty(),
                "quarter-turn holonomy unexpectedly admitted a horizontal leaf");

        // Local fiber-frame changes conjugate closed holonomy. Its raw matrix
        // changes, while its cycle structure—the observable—does not.
        const auto transformed = connection.gauge_transform({
            {0, reflection}, {1, quarter_turn}, {2, group.back()}});
        require(transformed.holonomy(loop).cycle_signature() == holonomy.cycle_signature(),
                "holonomy conjugacy class changed under local fiber automorphisms");
        require(transformed.gauge_invariant_signature() == connection.gauge_invariant_signature(),
                "spanning-tree gauge quotient changed under local fiber automorphisms");
        require(tables.multiply(tables.index_of(reflection), tables.index_of(quarter_turn)) ==
                    tables.index_of(Permutation::compose(reflection, quarter_turn)),
                "compact group multiplication table disagreed with permutation composition");
        require(tables.act(tables.index_of(quarter_turn), 0) == 1,
                "compact group action table disagreed with fiber action");

        FiberBundleConnection flat(triangle, fiber);
        require(flat.flat_on({loop}), "identity connection was reported curved");
        require(flat.gauge_invariant_signature() != connection.gauge_invariant_signature(),
                "gauge quotient conflated flat and curved connections");
        require(flat.covariantly_constant_section({{0, 2}, {1, 2}, {2, 2}}),
                "flat product bundle rejected a constant section");

        // Coupling to a Wolfram-style base rewrite: factor the old horizontal
        // lift through a fresh fiber and verify the complete gauge orbit.
        const auto extensions = connection.subdivide_edge_extensions(0, 1, 3);
        require(extensions.size() == group.size(),
                "edge rewrite did not enumerate Aut(F) factorizations");
        for (std::size_t i = 0; i < extensions.size(); ++i) {
            const auto& extension = extensions[i];
            require(extension.parallel_transport({0, 3, 1}) == connection.edge_transport(0, 1),
                    "edge rewrite changed boundary parallel transport");
            const auto gauge_copy = extensions.front().gauge_transform({{3, group[i]}});
            require(gauge_copy.edge_transport(0, 3) == extension.edge_transport(0, 3) &&
                        gauge_copy.edge_transport(3, 1) == extension.edge_transport(3, 1),
                    "rewrite factorizations did not form one fresh-fiber gauge orbit");
        }

        // Quantize the rewrite orbit without choosing a named continuous group.
        const auto identity = Permutation::identity(4);
        const auto state = SubdivisionWavefunction::from_total_transport(group, identity);
        require(std::abs(state.norm_squared() - 1.0) < 1e-12,
                "automorphism-orbit state was not normalized");
        require(state.support_size() == group.size(),
                "automorphism-orbit state had the wrong support");
        for (const auto& frame : group) {
            require(state.gauge_transform_midpoint(frame).approximately_equal(state),
                    "quantum rewrite state depended on the fresh fiber frame");
        }
        const auto curved_sector =
            SubdivisionWavefunction::from_total_transport(group, quarter_turn);
        require(std::abs(state.inner_product(curved_sector)) < 1e-12,
                "distinct boundary holonomy sectors were not orthogonal");

        std::cout << "InfraGauge-derived combinatorial gauge invariants: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "InfraGauge-derived combinatorial gauge invariants: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
