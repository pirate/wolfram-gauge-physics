#include <iostream>
#include "d4_loop_observer.hpp"

using namespace wgphysics::infragauge;
using namespace wgphysics::observables;
template<class T> void array(const T& values) {
    std::cout<<'['; bool first=true;
    for(const auto x : values) { if(!first) std::cout<<','; first=false; std::cout<<x; }
    std::cout<<']';
}
template<class T> void matrix(const T& values) {
    std::cout<<'['; bool first=true;
    for(const auto& x : values) { if(!first) std::cout<<','; first=false; array(x); }
    std::cout<<']';
}
int main() {
    try {
        std::size_t nv,ne,ns;
        if(!(std::cin>>nv>>ne>>ns) || nv>1000 || ne>3000 || !ns || ns>128)
            throw std::invalid_argument("invalid bounded observer dimensions");
        std::vector<Vertex> vertices(nv); std::vector<Edge> edges(ne);
        for(auto& v : vertices) if(!(std::cin>>v)) throw std::invalid_argument("missing vertex");
        for(auto& [u,v] : edges) if(!(std::cin>>u>>v) || u>=v) throw std::invalid_argument("missing canonical edge");
        const BaseGraph base(vertices,edges);
        if(std::vector<Edge>(base.edges().begin(),base.edges().end())!=edges)
            throw std::invalid_argument("edges must be distinct and sorted");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const D4LoopForest forest(base,group);
        std::cout<<"{\"states\":[";
        for(std::size_t s=0;s<ns;++s) {
            std::vector<uint16_t> links(ne);
            for(auto& x : links) { unsigned int value; if(!(std::cin>>value) || value>=8) throw std::invalid_argument("invalid link"); x=value; }
            const auto state=forest.probe(links);
            const auto representative=forest.representative(state.normalized);
            FiberBundleConnection raw(base,fiber),canonical(base,fiber);
            for(std::size_t e=0;e<ne;++e) {
                const auto [u,v]=edges[e];
                raw.set_transport(u,v,group.elements()[links[e]]);
                canonical.set_transport(u,v,group.elements()[representative[e]]);
            }
            if(raw.gauge_invariant_signature()!=canonical.gauge_invariant_signature())
                throw std::logic_error("linear observer differs from legacy gauge quotient");
            if(s) std::cout<<',';
            std::cout<<"{\"based_loops\":"; matrix(state.based_loops);
            std::cout<<",\"normalized\":"; matrix(state.normalized);
            std::cout<<",\"ranks\":"; array(state.ranks);
            std::cout<<",\"representative_links\":"; array(representative);
            std::cout<<",\"legacy_quotient_agrees\":true}";
        }
        std::string extra; if(std::cin>>extra) throw std::invalid_argument("trailing observer input");
        std::cout<<"]}\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
