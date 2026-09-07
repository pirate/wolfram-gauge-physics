#include <filesystem>
#include <fstream>
#include <iostream>
#include "cell_chain.hpp"

using namespace wgphysics::research;

int main(int argc,char** argv) {
    try {
        std::string path="out/chain-experiments.json";
        if (argc==3 && std::string(argv[1])=="--output") path=argv[2];
        else if (argc!=1) throw std::invalid_argument("usage: wgphysics_chain_experiments [--output FILE]");
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables tables(fiber.automorphisms());
        const auto id=tables.index_of(Permutation::identity(4));
        const auto rotation=tables.index_of(Permutation({1,2,3,0}));
        const auto reflection=tables.index_of(Permutation({0,3,2,1}));
        const auto parent=std::filesystem::path(path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot open experiment output");
        out << "{\"schema\":1,\"geometry\":\"supplied open chain of triangular cells\","
               "\"clock\":\"alternating disjoint pair layers\",\"fiber\":\"C4\","
               "\"runs\":[\n";
        bool first_run=true;
        std::size_t runs=0,updates=0;
        for (const auto rule : {PairInteraction::Hurwitz,PairInteraction::BoundaryShear})
            for (const std::size_t count : {64,128,256,1024})
                for (const auto seed : {"flat","rotation","reflection","neutral_pair","noncommuting_pair"})
                    for (std::size_t phase=0;phase<2;++phase) {
                        std::vector<uint16_t> initial(count,id);
                        const std::size_t center=count/2;
                        const std::string seed_name=seed;
                        if (seed_name=="rotation") initial[center]=rotation;
                        if (seed_name=="reflection") initial[center]=reflection;
                        if (seed_name=="neutral_pair") {
                            initial[center]=rotation; initial[center+1]=tables.inverse(rotation);
                        }
                        if (seed_name=="noncommuting_pair") {
                            initial[center]=rotation; initial[center+1]=reflection;
                        }
                        CellChain chain(fiber,initial,std::vector<uint16_t>(count-1,id),rule);
                        const auto initial_total=chain.total();
                        const auto initial_subgroup=chain.generated_subgroup();
                        if (!first_run) out << ",\n";
                        first_run=false;
                        out << "{\"rule\":\"" << (rule==PairInteraction::Hurwitz ? "hurwitz" : "boundary_shear")
                            << "\",\"cells\":" << count << ",\"seed\":\"" << seed
                            << "\",\"schedule_phase\":" << phase << ",\"samples\":[";
                        for (std::size_t tick=0;tick<=count/4;++tick) {
                            if (tick) for (std::size_t i=(tick-1+phase)%2;i+1<count;i+=2) {
                                chain.update(i,HurwitzDirection::Forward,false); ++updates;
                            }
                            if (chain.total()!=initial_total) throw std::logic_error("lost total holonomy");
                            if (tick%4 && tick!=count/4) continue;
                            std::size_t active=0; long left=0,right=0; double sum=0,sum2=0;
                            std::map<uint16_t,std::size_t> histogram;
                            for (std::size_t i=0;i<count;++i) {
                                ++histogram[chain.sector(chain.values()[i])];
                                if (chain.values()[i]==id) continue;
                                const auto position=static_cast<long>(i)-static_cast<long>(center);
                                if (!active) left=position;
                                right=position; ++active; sum+=position; sum2+=position*position;
                                if (std::abs(position)>static_cast<long>(tick)+1)
                                    throw std::logic_error("changed outside the imposed dependency cone");
                            }
                            if (tick) out << ',';
                            out << "{\"layer\":" << tick << ",\"active_cells\":" << active
                                << ",\"left\":" << left << ",\"right\":" << right
                                << ",\"mean_position\":" << (active ? sum/active : 0)
                                << ",\"mean_squared_position\":" << (active ? sum2/active : 0)
                                << ",\"sector_counts\":{";
                            bool first_sector=true;
                            for (const auto [sector,n] : histogram) {
                                if (!first_sector) out << ',';
                                first_sector=false;
                                out << '"' << sector << "\":" << n;
                            }
                            out << '}';
                            if (count==128) {
                                out << ",\"sector_profile\":[";
                                for (std::size_t i=0;i<count;++i) {
                                    if (i) out << ',';
                                    out << chain.sector(chain.values()[i]);
                                }
                                out << ']';
                            }
                            out << '}';
                        }
                        if (chain.generated_subgroup()!=initial_subgroup)
                            throw std::logic_error("generated subgroup changed");
                        out << "],\"total_holonomy_conserved\":true,\"generated_subgroup_order\":"
                            << initial_subgroup.size() << ",\"generated_subgroup_conserved\":true}";
                        ++runs;
                    }
        out << "\n]}\n";
        if (!out) throw std::runtime_error("failed writing experiment output");
        std::cout << runs << " runs; " << updates << " exact local updates; wrote " << path << '\n';
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
