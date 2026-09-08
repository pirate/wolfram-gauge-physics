#include <iostream>
#include "triple_rules.hpp"
#include "cycle_fiber_cli.hpp"

using namespace wgphysics::research;
int main(int argc,char** argv) {
    try {
        const auto fiber=cycle_fiber_from_arguments(argc,argv);
        const AutomorphismTables group(fiber.automorphisms());
        const auto rules=search_triple_involutions(group);
        std::cout<<"{\"schema\":1,\"group_order\":"<<group.order()<<",\"minimal_rules\":[";
        for(std::size_t id=0;id<rules.size();++id) {
            if(id) std::cout<<',';
            std::cout<<"{\"id\":"<<id<<",\"transpositions\":["; bool first=true;
            for(std::size_t x=0;x<rules[id].entries.size();++x) if(x<rules[id].entries[x]) {
                if(!first) std::cout<<','; first=false; std::cout<<'['<<x<<','<<rules[id].entries[x]<<']';
            }
            std::cout<<"]}";
        }
        std::cout<<"]}\n";
        std::cerr<<rules.size()<<" minimal triple involutions verified; not a census of all disjoint unions\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
