#pragma once

#include <array>
#include <span>
#include "infragauge.hpp"

namespace wgphysics::observables {

// Derived and exhaustively verified from Aut(F), never an unchecked group label.
class D4CentralCoordinates {
public:
    explicit D4CentralCoordinates(const infragauge::AutomorphismTables& group) : group_(group) {
        if(group_.order()!=8) throw std::invalid_argument("central coordinates require order eight");
        identity_=group_.index_of(infragauge::Permutation::identity(group_.fiber_size()));
        std::vector<uint16_t> center;
        for(uint16_t a=0;a<8;++a) {
            bool central=true;
            for(uint16_t b=0;b<8;++b) central &= group_.multiply(a,b)==group_.multiply(b,a);
            if(central) center.push_back(a);
        }
        if(center.size()!=2) throw std::invalid_argument("central coordinates require a two-element center");
        z_=center[0]==identity_ ? center[1] : center[0];
        bool found=false;
        for(uint16_t r=0;r<8 && !found;++r) for(uint16_t s=0;s<8 && !found;++s)
            if(group_.multiply(r,r)==identity_ && group_.multiply(s,s)==identity_ &&
               group_.multiply(r,s)!=group_.multiply(s,r)) { r_=r; s_=s; found=true; }
        if(!found) throw std::invalid_argument("no noncommuting involution generators");
        std::set<uint16_t> covered;
        for(uint16_t code=0;code<8;++code) {
            const auto a=code&1,b=(code>>1)&1,c=(code>>2)&1;
            const auto value=group_.multiply(c ? z_ : identity_,group_.multiply(a ? r_ : identity_,b ? s_ : identity_));
            element_[code]=value; bits_[value]=code; covered.insert(value);
        }
        if(covered.size()!=8) throw std::invalid_argument("central normal form is incomplete");
        for(uint16_t x=0;x<8;++x) for(uint16_t y=0;y<8;++y) {
            const auto a=bits_[x],b=bits_[y];
            const auto product=(a^b)^((((a>>1)&1)&(b&1))<<2);
            if(element_[product]!=group_.multiply(x,y)) throw std::logic_error("central multiplication law failed");
            const auto conjugate=group_.multiply(x,group_.multiply(y,group_.inverse(x)));
            const auto changed=b^((((a&1)&((b>>1)&1))^(((a>>1)&1)&(b&1)))<<2);
            if(bits_[conjugate]!=changed) throw std::logic_error("central conjugation law failed");
        }
    }

    std::pair<std::vector<uint16_t>,std::size_t> normalize(std::span<const uint16_t> loops) const {
        std::array<std::size_t,2> pivots{}; std::size_t rank=0;
        for(std::size_t i=0;i<loops.size();++i) {
            if(loops[i]>=8) throw std::out_of_range("invalid based holonomy");
            const auto kind=bits_[loops[i]]&3;
            if(kind && rank<2 && (rank==0 || kind!=(bits_[loops[pivots[0]]]&3))) pivots[rank++]=i;
        }
        uint16_t frame=0;
        for(;frame<4;++frame) {
            bool zero=true;
            for(std::size_t j=0;j<rank;++j) zero &= (change(bits_[loops[pivots[j]]],frame)&4)==0;
            if(zero) break;
        }
        if(frame==4) throw std::logic_error("independent frame pivots were not solvable");
        std::vector<uint16_t> result; result.reserve(loops.size());
        for(const auto value : loops) result.push_back(element_[change(bits_[value],frame)]);
        return {std::move(result),rank};
    }

    uint16_t identity() const { return identity_; }
    const infragauge::AutomorphismTables& group() const { return group_; }

private:
    infragauge::AutomorphismTables group_;
    uint16_t identity_{},r_{},s_{},z_{};
    std::array<uint16_t,8> element_{},bits_{};
    static uint16_t change(uint16_t code,uint16_t frame) {
        return code^((((frame&1)&((code>>1)&1))^(((frame>>1)&1)&(code&1)))<<2);
    }
};

struct CompleteLoopState {
    std::vector<std::vector<uint16_t>> based_loops;
    std::vector<std::vector<uint16_t>> normalized;
    std::vector<std::size_t> ranks;
};

// Forest compilation is topology-only. A probe computes each tree transport once,
// then one fundamental holonomy per chord: O(V+E), without materializing walks.
// This removes local frames on a fixed labeled graph, not vertex isomorphisms.
class D4LoopForest {
public:
    D4LoopForest(const infragauge::BaseGraph& base,const infragauge::AutomorphismTables& group)
        : coordinates_(group),vertices_(base.vertices().begin(),base.vertices().end()),
          edges_(base.edges().begin(),base.edges().end()) {
        std::map<infragauge::Vertex,std::size_t> ids;
        for(std::size_t i=0;i<vertices_.size();++i) ids.emplace(vertices_[i],i);
        struct Neighbor { std::size_t vertex,edge; bool reverse; };
        std::vector<std::vector<Neighbor>> adjacency(vertices_.size());
        for(std::size_t e=0;e<edges_.size();++e) {
            const auto [u,v]=edges_[e]; const auto a=ids.at(u),b=ids.at(v);
            endpoints_.push_back({a,b});
            adjacency[a].push_back({b,e,false}); adjacency[b].push_back({a,e,true});
        }
        for(auto& neighbors : adjacency) std::sort(neighbors.begin(),neighbors.end(),
            [](const auto& a,const auto& b){return a.vertex<b.vertex;});
        const auto absent=std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> component(vertices_.size(),absent);
        std::vector<bool> tree(edges_.size());
        for(std::size_t root=0;root<vertices_.size();++root) {
            if(component[root]!=absent) continue;
            const auto c=components_.size(); components_.push_back({root,{},{}});
            component[root]=c; std::vector<std::size_t> queue{root};
            for(std::size_t i=0;i<queue.size();++i) for(const auto next : adjacency[queue[i]])
                if(component[next.vertex]==absent) {
                    component[next.vertex]=c; tree[next.edge]=true; queue.push_back(next.vertex);
                    components_.back().tree.push_back({queue[i],next.vertex,next.edge,next.reverse});
                }
        }
        for(std::size_t e=0;e<edges_.size();++e)
            if(!tree[e]) components_[component[endpoints_[e].first]].chords.push_back(e);
    }

    CompleteLoopState probe(std::span<const uint16_t> links) const {
        if(links.size()!=edges_.size()) throw std::invalid_argument("wrong link vector length");
        for(const auto x : links) if(x>=8) throw std::out_of_range("invalid link value");
        const auto& group=coordinates_.group();
        std::vector<uint16_t> paths(vertices_.size(),coordinates_.identity());
        CompleteLoopState result;
        for(const auto& component : components_) {
            for(const auto [parent,child,edge,reverse] : component.tree)
                paths[child]=group.multiply(reverse ? group.inverse(links[edge]) : links[edge],paths[parent]);
            std::vector<uint16_t> loops; loops.reserve(component.chords.size());
            for(const auto edge : component.chords) {
                const auto [u,v]=endpoints_[edge];
                loops.push_back(group.multiply(group.inverse(paths[v]),group.multiply(links[edge],paths[u])));
            }
            const auto [normal,rank]=coordinates_.normalize(loops);
            result.based_loops.push_back(std::move(loops)); result.normalized.push_back(normal); result.ranks.push_back(rank);
        }
        return result;
    }

    std::vector<uint16_t> representative(const std::vector<std::vector<uint16_t>>& signature) const {
        if(signature.size()!=components_.size()) throw std::invalid_argument("wrong component count");
        std::vector<uint16_t> links(edges_.size(),coordinates_.identity());
        for(std::size_t c=0;c<components_.size();++c) {
            if(signature[c].size()!=components_[c].chords.size()) throw std::invalid_argument("wrong cycle rank");
            for(std::size_t i=0;i<signature[c].size();++i) {
                if(signature[c][i]>=8) throw std::out_of_range("invalid normalized holonomy");
                links[components_[c].chords[i]]=signature[c][i];
            }
        }
        return links;
    }

private:
    struct TreeEdge { std::size_t parent,child,edge; bool reverse; };
    struct Component { std::size_t root; std::vector<TreeEdge> tree; std::vector<std::size_t> chords; };
    D4CentralCoordinates coordinates_;
    std::vector<infragauge::Vertex> vertices_;
    std::vector<infragauge::Edge> edges_;
    std::vector<std::pair<std::size_t,std::size_t>> endpoints_;
    std::vector<Component> components_;
};

}  // namespace wgphysics::observables
