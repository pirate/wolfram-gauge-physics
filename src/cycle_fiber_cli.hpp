#pragma once

#include <charconv>
#include <string_view>
#include "research.hpp"

namespace wgphysics::research {

// Keep the historical square-fiber protocol by default. The optional cycle
// changes adjacency only; every runner still derives its automorphism group.
inline FiberGraph cycle_fiber_from_arguments(int argc,char** argv) {
    uint32_t vertices=4;
    if(argc!=1) {
        if(argc!=3 || std::string_view(argv[1])!="--cycle")
            throw std::invalid_argument("usage: runner [--cycle 3..6]");
        const std::string_view argument(argv[2]);
        const auto [end,error]=std::from_chars(argument.data(),argument.data()+argument.size(),vertices);
        if(error!=std::errc{} || end!=argument.data()+argument.size() || vertices<3 || vertices>6)
            throw std::invalid_argument("cycle fiber requires 3 through 6 vertices");
    }
    std::vector<Edge> edges;
    for(uint32_t i=0;i<vertices;++i) edges.emplace_back(i,(i+1)%vertices);
    return FiberGraph(vertices,std::move(edges));
}

} // namespace wgphysics::research
