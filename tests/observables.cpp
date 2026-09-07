#include "observables.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

wgphysics::infragauge::BaseGraph cycle_graph(std::size_t size, uint32_t label_offset = 0) {
    using namespace wgphysics::infragauge;
    std::vector<Vertex> vertices(size);
    std::iota(vertices.begin(), vertices.end(), label_offset);
    std::vector<Edge> edges;
    for (std::size_t index = 0; index < size; ++index) {
        edges.emplace_back(vertices[index], vertices[(index + 1) % size]);
    }
    return BaseGraph(std::move(vertices), std::move(edges));
}

wgphysics::infragauge::BaseGraph torus_graph(std::size_t width, std::size_t height) {
    using namespace wgphysics::infragauge;
    std::vector<Vertex> vertices(width * height);
    std::iota(vertices.begin(), vertices.end(), 0);
    std::vector<Edge> edges;
    auto vertex = [=](std::size_t x, std::size_t y) {
        return static_cast<Vertex>((y % height) * width + (x % width));
    };
    for (std::size_t y = 0; y < height; ++y) {
        for (std::size_t x = 0; x < width; ++x) {
            edges.emplace_back(vertex(x, y), vertex(x + 1, y));
            edges.emplace_back(vertex(x, y), vertex(x, y + 1));
        }
    }
    return BaseGraph(std::move(vertices), std::move(edges));
}

double mean_window(const std::vector<std::optional<double>>& values,
                   std::size_t first, std::size_t last) {
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t index = first; index <= last; ++index) {
        if (values.at(index)) {
            sum += *values[index];
            ++count;
        }
    }
    if (count == 0) throw std::runtime_error("dimension window contains no estimates");
    return sum / static_cast<double>(count);
}

void require_close(double left, double right, double tolerance, const char* message) {
    require(std::abs(left - right) <= tolerance, message);
}

}  // namespace

int main() {
    try {
        using namespace wgphysics::infragauge;
        using namespace wgphysics::observables;

        const ProbeConfig config{16, 32, 256, 3, 0.35};
        const auto cycle = probe_intrinsic_geometry(cycle_graph(128), config);
        const auto torus = probe_intrinsic_geometry(torus_graph(16, 16), config);

        // The detector is not told either graph's dimension. It recovers the
        // expected intrinsic regimes from graph-distance volume growth.
        const auto cycle_volume_dimension = mean_window(cycle.volume_dimension, 6, 12);
        const auto torus_volume_dimension = mean_window(torus.volume_dimension, 3, 6);
        require(cycle_volume_dimension > 0.8 && cycle_volume_dimension < 1.1,
                "cycle volume growth did not recover one dimension");
        require(torus_volume_dimension > 1.5 && torus_volume_dimension < 2.2,
                "torus volume growth did not recover two dimensions");
        require(torus_volume_dimension > cycle_volume_dimension + 0.6,
                "intrinsic volume probe failed to distinguish graph dimensions");

        // Lazy diffusion removes the even/odd artifact of bipartite graphs.
        // Intermediate-time return scaling also distinguishes 1D from 2D.
        const auto cycle_spectral_dimension = mean_window(cycle.spectral_dimension, 8, 20);
        const auto torus_spectral_dimension = mean_window(torus.spectral_dimension, 4, 10);
        require(cycle_spectral_dimension > 0.7 && cycle_spectral_dimension < 1.3,
                "cycle diffusion did not recover one spectral dimension");
        require(torus_spectral_dimension > 1.2 && torus_spectral_dimension < 2.5,
                "torus diffusion did not recover a higher spectral dimension");
        require(torus_spectral_dimension > cycle_spectral_dimension + 0.3,
                "spectral probe failed to distinguish cycle and torus");

        require(std::all_of(cycle.ball_volume_coefficient_of_variation.begin(),
                            cycle.ball_volume_coefficient_of_variation.end(),
                            [](double value) { return value < 1e-12; }),
                "vertex-transitive cycle was reported as inhomogeneous");
        require(std::all_of(torus.ball_volume_coefficient_of_variation.begin(),
                            torus.ball_volume_coefficient_of_variation.end(),
                            [](double value) { return value < 1e-12; }),
                "vertex-transitive torus was reported as inhomogeneous");

        const auto relabeled = probe_intrinsic_geometry(cycle_graph(128, 1000), config);
        require(cycle.mean_ball_volume == relabeled.mean_ball_volume &&
                    cycle.lazy_return_probability == relabeled.lazy_return_probability,
                "intrinsic geometry depends on vertex labels");

        // Gauge probes use only conjugacy data and normalized characters.
        // One curved face among two should be maximally localized.
        const FiberGraph square(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        const Permutation quarter_turn({1, 2, 3, 0});
        const Permutation reflection({0, 3, 2, 1});
        const BaseGraph two_faces({0, 1, 2, 3},
                                  {{0, 1}, {1, 2}, {2, 0}, {1, 3}, {3, 0}});
        const std::vector<std::vector<Vertex>> loops{{0, 1, 2, 0}, {0, 1, 3, 0}};
        FiberBundleConnection localized(two_faces, square);
        localized.set_transport(2, 0, quarter_turn);
        const auto gauge_profile = probe_gauge_field(localized, loops);
        require(gauge_profile.nontrivial_loops == 1 &&
                    std::abs(gauge_profile.curvature_density - 0.5) < 1e-12,
                "gauge probe failed to count the localized curved face");
        require_close(gauge_profile.effective_curvature_support, 1.0, 1e-12,
                      "curvature support is not localized on one face");
        require_close(gauge_profile.normalized_localization, 1.0, 1e-12,
                      "curvature localization score is wrong");

        const auto framed = localized.gauge_transform({{0, reflection}, {1, quarter_turn}});
        const auto framed_profile = probe_gauge_field(framed, loops);
        for (std::size_t index = 0; index < loops.size(); ++index) {
            require(gauge_profile.loops[index].conjugacy_signature ==
                            framed_profile.loops[index].conjugacy_signature &&
                        gauge_profile.loops[index].normalized_permutation_character ==
                            framed_profile.loops[index].normalized_permutation_character,
                    "loop observable depends on a local gauge frame");
        }

        auto half_turn_state = localized, edge_reflection_state = localized;
        const Permutation half_turn({2,3,0,1}), edge_reflection({1,0,3,2});
        require(half_turn.cycle_signature() == edge_reflection.cycle_signature(),
                "control must have identical symmetric-group cycle type");
        half_turn_state.set_transport(2,0,half_turn);
        edge_reflection_state.set_transport(2,0,edge_reflection);
        require(probe_gauge_loop(half_turn_state,loops[0]).conjugacy_signature !=
                    probe_gauge_loop(edge_reflection_state,loops[0]).conjugacy_signature,
                "probe merged distinct Aut(C4) curvature sectors");

        const auto bins = aggregate_loop_characters(gauge_profile.loops);
        require(bins.size() == 1 && bins.front().edge_length == 3 &&
                    bins.front().samples == 2 &&
                    std::abs(bins.front().mean_character - 0.5) < 1e-12,
                "loop-character aggregation is wrong");

        std::cout << "Intrinsic geometry and gauge observables: PASS\n"
                  << "cycle volume dimension: " << cycle_volume_dimension << '\n'
                  << "torus volume dimension: " << torus_volume_dimension << '\n'
                  << "cycle spectral dimension: " << cycle_spectral_dimension << '\n'
                  << "torus spectral dimension: " << torus_spectral_dimension << '\n'
                  << "localized curvature support: "
                  << gauge_profile.effective_curvature_support << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Intrinsic geometry and gauge observables: FAIL: " << error.what() << '\n';
        return 1;
    }
}
