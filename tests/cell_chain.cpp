#include <iostream>
#include "cell_chain.hpp"

using namespace wgphysics::research;
void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }

int main() {
    try {
        const FiberGraph fiber(4,{{0,1},{1,2},{2,3},{3,0}});
        const AutomorphismTables tables(fiber.automorphisms());
        const auto identity=tables.index_of(Permutation::identity(4));
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        constexpr std::size_t count=8;
        for (Vertex i=0;i<count;++i) {
            for (Vertex k=0;k<3;++k) vertices.push_back(3*i+k);
            edges.insert(edges.end(),{{3*i,3*i+1},{3*i+1,3*i+2},{3*i+2,3*i}});
            if (i) edges.emplace_back(3*i,3*(i-1));
        }
        const BaseGraph base(vertices,edges);
        std::size_t checked=0;
        for (const auto rule : {PairInteraction::Hurwitz,PairInteraction::BoundaryShear}) {
            FiberBundleConnection reference(base,fiber);
            std::vector<uint16_t> initial(count), connectors(count-1);
            for (Vertex i=0;i<count;++i) {
                initial[i]=(3*i+1)%tables.order();
                reference.set_transport(3*i+2,3*i,tables.elements()[initial[i]]);
                if (i) {
                    connectors[i-1]=(5*i+2)%tables.order();
                    reference.set_transport(3*i,3*(i-1),tables.elements()[connectors[i-1]]);
                }
            }
            CellChain chain(fiber,initial,connectors,rule);
            const auto total=chain.total();
            const auto subgroup=chain.generated_subgroup();
            for (std::size_t tick=0;tick<16;++tick) {
                for (std::size_t i=tick%2;i+1<count;i+=2) {
                    const Vertex a=3*i,b=3*(i+1);
                    reference=apply_pair_interaction(reference,{{a,a+1,a+2,a},{b,b+1,b+2,b},{b,a}},rule);
                    chain.update(i);
                    for (Vertex j=0;j<count;++j)
                        require(tables.elements()[chain.values()[j]]==reference.holonomy({3*j,3*j+1,3*j+2,3*j}),
                                "compiled chain differs from link oracle");
                    require(chain.total()==total,"global transported product not conserved");
                    require(chain.generated_subgroup()==subgroup,"generated holonomy subgroup changed");
                    ++checked;
                }
            }
            const auto events=chain.events();
            // Replay a different topological ordering of the actual dependency DAG.
            CellChain reordered(fiber,initial,connectors,rule);
            std::vector<bool> done(events.size(),false);
            std::size_t reordered_positions=0;
            for (std::size_t position=0;position<events.size();++position) {
                std::optional<std::size_t> selected;
                for (std::size_t id=0;id<events.size();++id) if (!done[id]) {
                    bool ready=true;
                    for (const auto p : events[id].parents) ready=ready && done[p];
                    if (ready) selected=id; // highest ready ID reverses independent work
                }
                require(selected.has_value(),"dependency graph contains a cycle");
                const auto id=*selected, left=events[id].left;
                require(std::pair{reordered.values()[left],reordered.values()[left+1]}==events[id].before,
                        "topological replay read a different input version");
                reordered.update(left,HurwitzDirection::Forward,false);
                done[id]=true;
                if (id!=position) ++reordered_positions;
            }
            require(reordered_positions>0 && reordered.values()==chain.values(),
                    "independent schedule replay did not preserve the final connection");
            for (std::size_t i=events.size();i-->0;)
                chain.update(events[i].left,HurwitzDirection::Inverse,false);
            require(chain.values()==initial,"reversing recorded schedule did not recover links");
            CellChain dependency(fiber,std::vector<uint16_t>(4,identity),
                                  std::vector<uint16_t>(3,identity),rule);
            dependency.update(0); dependency.update(2); dependency.update(1);
            require(dependency.events()[0].parents.empty() && dependency.events()[1].parents.empty(),
                    "independent updates acquired dependencies");
            require(dependency.events()[2].parents==std::vector<std::size_t>{0,1},
                    "read/write dependency provenance incorrect");
            require(dependency.events()[2].depth==1,"causal depth incorrect");
        }
        // Independent check on a cyclic subgroup: shear reduces to a linear map mod 4.
        const Permutation rotation({1,2,3,0});
        std::vector<uint16_t> powers{identity};
        for (int i=1;i<4;++i) powers.push_back(tables.multiply(powers.back(),tables.index_of(rotation)));
        std::vector<unsigned> exponent(128); exponent[64]=1;
        std::vector<uint16_t> values(128,identity); values[64]=powers[1];
        CellChain linear(fiber,values,std::vector<uint16_t>(127,identity),PairInteraction::BoundaryShear);
        for (std::size_t tick=0;tick<32;++tick) {
            for (std::size_t i=tick%2;i+1<128;i+=2) {
                const auto a=exponent[i],b=exponent[i+1];
                exponent[i]=(2*a+b)%4; exponent[i+1]=(4-a)%4;
                linear.update(i,HurwitzDirection::Forward,false);
            }
            for (std::size_t i=0;i<128;++i)
                require(linear.values()[i]==powers[exponent[i]],"cyclic-subgroup reduction failed");
        }
        // Full nonabelian two-layer reduction: the even sublattice streams
        // exactly; the odd sublattice is driven by the streaming variables.
        const FiberGraph triangle(3,{{0,1},{1,2},{2,0}});
        const AutomorphismTables nonabelian(triangle.automorphisms());
        std::vector<uint16_t> input(128);
        for (std::size_t i=0;i<input.size();++i) input[i]=(5*i+i/3+1)%nonabelian.order();
        const auto trivial=nonabelian.index_of(Permutation::identity(3));
        CellChain split(triangle,input,std::vector<uint16_t>(127,trivial),PairInteraction::BoundaryShear);
        for (std::size_t i=0;i+1<input.size();i+=2) split.update(i,HurwitzDirection::Forward,false);
        for (std::size_t i=1;i+1<input.size();i+=2) split.update(i,HurwitzDirection::Forward,false);
        for (std::size_t j=1;j+1<input.size()/2;++j) {
            require(split.values()[2*j]==input[2*(j-1)],"exact streaming sublattice reduction failed");
            const auto a=nonabelian.inverse(input[2*j]),c=input[2*(j+1)],b=input[2*(j+1)+1];
            const auto target=nonabelian.multiply(a,nonabelian.multiply(c,
                nonabelian.multiply(b,nonabelian.multiply(c,a))));
            require(split.values()[2*j+1]==target,"nonabelian driven sublattice reduction failed");
        }
        std::cout << checked << " compiled/link oracle steps, reverse replay, dependency DAG, "
                     "subgroup conservation, cyclic and nonabelian sublattice reductions passed\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
