#include <chrono>
#include <atomic>
#include <iostream>
#include "d4_loop_observer.hpp"

using namespace wgphysics::infragauge;
using namespace wgphysics::observables;
using Clock=std::chrono::steady_clock;
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double,std::milli>(Clock::now()-start).count();
}
BaseGraph torus(std::size_t side) {
    std::vector<Vertex> vertices(side*side); std::iota(vertices.begin(),vertices.end(),0);
    std::vector<Edge> edges; edges.reserve(3*side*side);
    const auto v=[&](std::size_t x,std::size_t y){return static_cast<Vertex>((y%side)*side+x%side);};
    for(std::size_t y=0;y<side;++y) for(std::size_t x=0;x<side;++x) {
        edges.push_back(canonical_edge(v(x,y),v(x+1,y)));
        edges.push_back(canonical_edge(v(x,y),v(x,y+1)));
        edges.push_back(canonical_edge(v(x,y),v(x+1,y+1)));
    }
    return BaseGraph(vertices,edges);
}
int main(int argc,char** argv) {
    try {
        const bool smoke=argc==2 && std::string(argv[1])=="--smoke";
        if(argc!=1 && !smoke) throw std::invalid_argument("usage: wgphysics_loop_bench [--smoke]");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        std::cout<<"{\"schema\":1,\"clock\":\"steady_clock milliseconds; CPU only\",\"runs\":[";
        const std::vector<std::size_t> sizes=smoke ? std::vector<std::size_t>{16,32} : std::vector<std::size_t>{16,32,128,600};
        bool first=true;
        for(const auto side : sizes) {
            const auto base=torus(side);
            const auto start=Clock::now();
            const D4LoopForest forest(base,group);
            const auto compile_ms=milliseconds(start);
            std::vector<uint16_t> links; links.reserve(base.edges().size());
            uint32_t rng=619973;
            for(std::size_t e=0;e<base.edges().size();++e) { rng=1664525*rng+1013904223; links.push_back((rng>>24)%8); }
            const auto expected=forest.probe(links).normalized;
            std::vector<double> samples;
            for(int repeat=0;repeat<5;++repeat) {
                const auto begin=Clock::now();
                std::atomic_signal_fence(std::memory_order_seq_cst);
                const auto state=forest.probe(links);
                std::atomic_signal_fence(std::memory_order_seq_cst);
                samples.push_back(milliseconds(begin));
                if(state.normalized!=expected) throw std::logic_error("repeated observer disagrees");
            }
            std::vector<uint16_t> changed; changed.reserve(links.size()); std::size_t i=0;
            for(const auto [u,v] : base.edges()) {
                const auto frame=[](Vertex x){return static_cast<uint16_t>((x*3+x/7+1)%8);};
                changed.push_back(group.multiply(frame(v),group.multiply(links[i++],group.inverse(frame(u)))));
            }
            if(forest.probe(changed).normalized!=expected || forest.probe(forest.representative(expected)).normalized!=expected)
                throw std::logic_error("large observer failed local-frame or reconstruction test");
            std::vector<double> legacy_samples;
            if(side<=32) {
                FiberBundleConnection raw(base,fiber),canonical(base,fiber);
                const auto representative=forest.representative(expected); i=0;
                for(const auto [u,v] : base.edges()) {
                    raw.set_transport(u,v,group.elements()[links[i]]);
                    canonical.set_transport(u,v,group.elements()[representative[i++]]);
                }
                const auto reference=canonical.gauge_invariant_signature();
                for(int repeat=0;repeat<3;++repeat) {
                    const auto begin=Clock::now();
                    std::atomic_signal_fence(std::memory_order_seq_cst);
                    const auto value=raw.gauge_invariant_signature();
                    std::atomic_signal_fence(std::memory_order_seq_cst);
                    legacy_samples.push_back(milliseconds(begin));
                    if(value!=reference) throw std::logic_error("legacy gauge quotient disagrees");
                }
            }
            if(!first) std::cout<<','; first=false;
            std::cout<<"{\"side\":"<<side<<",\"vertices\":"<<base.vertices().size()<<",\"edges\":"<<base.edges().size()
                     <<",\"independent_loops\":"<<expected[0].size()<<",\"forest_compile_ms\":"<<compile_ms<<",\"probe_ms\":[";
            for(std::size_t j=0;j<samples.size();++j) { if(j) std::cout<<','; std::cout<<samples[j]; }
            std::cout<<"],\"legacy_ms\":[";
            for(std::size_t j=0;j<legacy_samples.size();++j) { if(j) std::cout<<','; std::cout<<legacy_samples[j]; }
            std::cout<<"],\"local_frame_check\":true,\"representative_check\":true,\"legacy_check\":"<<(legacy_samples.empty()?"null":"true")<<'}';
        }
        std::cout<<"]}\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
