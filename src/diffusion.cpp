#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include "cell_chain.hpp"

using namespace wgphysics::research;

int main(int argc,char** argv) {
    try {
        std::size_t trials=4096,cells=1040,layers=256;
        uint64_t seed=731291;
        for (int i=1;i<argc;i+=2) {
            if (i+1>=argc) throw std::invalid_argument("missing option value");
            const std::string key=argv[i]; const auto value=std::stoull(argv[i+1]);
            if (key=="--trials") trials=value;
            else if (key=="--cells") cells=value;
            else if (key=="--layers") layers=value;
            else if (key=="--seed") seed=value;
            else throw std::invalid_argument("unknown option");
        }
        if (trials<32 || trials%32 || cells<80 || cells%4 || layers<16 || layers%2 || layers>(cells-16)/4)
            throw std::invalid_argument("need trials divisible by 32, cells divisible by 4, even layers>=16, and cells>=4*layers+16");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables group(fiber.automorphisms());
        const auto n=group.order(),id=group.index_of(Permutation::identity(4));
        PairTable table(n*n);
        std::vector<uint16_t> symbols(n);
        for (auto& x : table) if (!(std::cin>>x)) throw std::invalid_argument("missing pair table input");
        for (auto& x : symbols) if (!(std::cin>>x) || x>3) throw std::invalid_argument("invalid charge symbol input");
        validate_pair_table(group,table);
        std::array<std::vector<uint16_t>,4> representatives;
        for (uint16_t a=0;a<n;++a) representatives[symbols[a]].push_back(a);
        for (const auto& entries : representatives) if (entries.empty()) throw std::invalid_argument("charge map not surjective");
        // Prove the effective oracle used below is this raw rule's factor map.
        for (uint16_t a=0;a<n;++a) for (uint16_t b=0;b<n;++b) {
            const auto result=table[a*n+b];
            auto left=symbols[a],right=symbols[b];
            if (!left || !right) std::swap(left,right);
            if (symbols[result/n]!=left || symbols[result%n]!=right)
                throw std::invalid_argument("selected table does not factor to zero-vacancy exclusion");
        }
        std::vector<std::size_t> times{0};
        for (std::size_t t=4;t<layers;t*=2) times.push_back(t);
        times.push_back(layers);
        std::cout << std::setprecision(15) << "{\"schema\":1,\"cells\":" << cells
                  << ",\"layers\":" << layers << ",\"seed\":" << seed
                  << ",\"clock\":\"one disjoint matching layer\",\"ensembles\":[\n";
        std::size_t ensemble_index=0;
        for (const double density : {0.0,0.25,0.5,0.75,1.0}) {
            const auto samples=(density==0 || density==1) ? std::min<std::size_t>(trials,256) : trials;
            const auto block_size=samples/32;
            std::vector<std::vector<std::array<double,4>>> blocks(32,
                std::vector<std::array<double,4>>(times.size(),{0,0,0,0}));
            std::vector<std::vector<std::size_t>> hist(times.size(),std::vector<std::size_t>(2*layers+1));
            std::mt19937_64 random(seed+ensemble_index*100003);
            const auto uniform=[&]() { return (random()>>11)*0x1.0p-53; };
            std::size_t initial_occupied=0;
            for (std::size_t trial=0;trial<samples;++trial) {
                std::vector<uint16_t> values(cells),coarse(cells);
                for (std::size_t i=0;i<cells;++i) {
                    coarse[i]=uniform()<density ? 2 : 0;
                    const auto& choices=representatives[coarse[i]];
                    values[i]=choices[random()%choices.size()];
                    if (i!=cells/2 && coarse[i]) ++initial_occupied;
                }
                coarse[cells/2]=1; values[cells/2]=representatives[1][trial%representatives[1].size()];
                CellChain microscopic(fiber,values,std::vector<uint16_t>(cells-1,id),table);
                const auto initial_total=microscopic.total();
                std::size_t tag=cells/2,checkpoint=0;
                for (std::size_t tick=0;tick<=layers;++tick) {
                    if (tick) {
                        for (std::size_t i=(tick-1+trial%2)%2;i+1<cells;i+=2) {
                            microscopic.update(i,HurwitzDirection::Forward,false);
                            if (!coarse[i] || !coarse[i+1]) std::swap(coarse[i],coarse[i+1]);
                            if (symbols[microscopic.values()[i]]!=coarse[i] ||
                                symbols[microscopic.values()[i+1]]!=coarse[i+1])
                                throw std::logic_error("microscopic/quotient disagreement");
                            if (coarse[i]==1) tag=i;
                            if (coarse[i+1]==1) tag=i+1;
                        }
                    }
                    if (tick!=times[checkpoint]) continue;
                    const auto x=static_cast<long>(tag)-static_cast<long>(cells/2);
                    auto& moments=blocks[trial/block_size][checkpoint];
                    double power=x;
                    for (auto& value : moments) { value+=power; power*=x; }
                    ++hist[checkpoint][x+static_cast<long>(layers)];
                    ++checkpoint;
                }
                if (microscopic.total()!=initial_total) throw std::logic_error("microscopic total product changed");
            }
            if (ensemble_index++) std::cout << ",\n";
            std::cout << "{\"density\":" << density << ",\"trials\":" << samples
                      << ",\"measured_initial_density\":" << double(initial_occupied)/(samples*(cells-1))
                      << ",\"block_size\":" << block_size << ",\"times\":[";
            for (std::size_t i=0;i<times.size();++i) { if(i) std::cout<<','; std::cout<<times[i]; }
            std::cout << "],\"block_raw_moment_sums\":[";
            for (std::size_t b=0;b<blocks.size();++b) {
                if(b) std::cout<<','; std::cout<<'[';
                for (std::size_t t=0;t<times.size();++t) {
                    if(t) std::cout<<','; std::cout<<'[';
                    for (std::size_t k=0;k<4;++k) { if(k) std::cout<<','; std::cout<<blocks[b][t][k]; }
                    std::cout<<']';
                }
                std::cout<<']';
            }
            std::cout << "],\"displacement_histograms\":[";
            for (std::size_t t=0;t<times.size();++t) {
                if(t) std::cout<<','; std::cout<<'[';
                for (std::size_t x=0;x<hist[t].size();++x) { if(x) std::cout<<','; std::cout<<hist[t][x]; }
                std::cout<<']';
            }
            std::cout << "]}" << std::flush;
            std::cerr << "density " << density << ": " << samples << " microscopic trials verified\n";
        }
        std::cout << "\n]}\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
