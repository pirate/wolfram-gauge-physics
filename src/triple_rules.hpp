#pragma once

#include <array>
#include "equivariant_pair_search.hpp"

namespace wgphysics::research {

struct TripleTable { std::vector<std::size_t> entries; };
using Triple = std::array<uint16_t,3>;

inline Triple decode_triple(std::size_t code,std::size_t n) {
    return {static_cast<uint16_t>(code/(n*n)),static_cast<uint16_t>((code/n)%n),static_cast<uint16_t>(code%n)};
}
inline std::size_t encode_triple(const Triple& t,std::size_t n) { return (t[0]*n+t[1])*n+t[2]; }
inline uint16_t triple_product(const AutomorphismTables& g,std::size_t code) {
    const auto t=decode_triple(code,g.order()); return g.multiply(g.multiply(t[0],t[1]),t[2]);
}
inline std::size_t conjugated_triple(const AutomorphismTables& g,std::size_t code,uint16_t frame) {
    auto t=decode_triple(code,g.order());
    for(auto& a : t) a=g.multiply(frame,g.multiply(a,g.inverse(frame)));
    return encode_triple(t,g.order());
}
inline std::size_t reversed_triple(const AutomorphismTables& g,std::size_t code) {
    const auto t=decode_triple(code,g.order());
    return encode_triple({g.inverse(t[2]),g.inverse(t[1]),g.inverse(t[0])},g.order());
}
inline void validate_triple_table(const AutomorphismTables& g,const TripleTable& rule) {
    const auto n=g.order(),id=g.index_of(Permutation::identity(g.fiber_size()));
    const auto& table=rule.entries; const std::size_t count=static_cast<std::size_t>(n)*n*n;
    if(table.size()!=count || std::any_of(table.begin(),table.end(),[&](auto x) { return x>=count; }))
        throw std::invalid_argument("triple table has invalid dimensions or entries");
    if(table[encode_triple({id,id,id},n)]!=encode_triple({id,id,id},n))
        throw std::invalid_argument("triple table moves the flat vacuum");
    for(std::size_t x=0;x<count;++x) {
        const auto y=table[x];
        if(table[y]!=x) throw std::invalid_argument("triple table is not an involution");
        if(triple_product(g,x)!=triple_product(g,y)) throw std::invalid_argument("triple table changes boundary product");
        if(table[reversed_triple(g,x)]!=reversed_triple(g,y))
            throw std::invalid_argument("triple table violates orientation inversion");
        for(uint16_t frame=0;frame<n;++frame)
            if(table[conjugated_triple(g,x,frame)]!=conjugated_triple(g,y,frame))
                throw std::invalid_argument("triple table violates gauge covariance");
    }
}

// Minimal transposition closures only, not every composite triple rule.
// Breadth-first closure is independent of the Python direct finite-action orbit.
inline std::vector<TripleTable> search_triple_involutions(const AutomorphismTables& g) {
    const auto n=g.order();
    if(n>8) throw std::invalid_argument("bounded triple census requires group order <=8");
    const std::size_t count=static_cast<std::size_t>(n)*n*n;
    const auto id=g.index_of(Permutation::identity(g.fiber_size()));
    const auto flat=encode_triple({id,id,id},n);
    std::vector<uint16_t> products(count);
    for(std::size_t x=0;x<count;++x) products[x]=triple_product(g,x);
    std::vector<std::vector<std::size_t>> action(n+1,std::vector<std::size_t>(count));
    for(std::size_t x=0;x<count;++x) {
        for(uint16_t frame=0;frame<n;++frame) action[frame][x]=conjugated_triple(g,x,frame);
        action[n][x]=reversed_triple(g,x);
    }
    using Exchanges=std::vector<std::pair<std::size_t,std::size_t>>;
    std::map<Exchanges,TripleTable> found;
    for(std::size_t x=0;x<count;++x) for(std::size_t y=x+1;y<count;++y) {
        if(x==flat || y==flat || products[x]!=products[y]) continue;
        TripleTable rule{std::vector<std::size_t>(count,count)};
        std::vector<std::pair<std::size_t,std::size_t>> queue{{x,y}};
        bool valid=true;
        for(std::size_t cursor=0;cursor<queue.size() && valid;++cursor) {
            const auto [a,b]=queue[cursor];
            if(rule.entries[a]!=count) { if(rule.entries[a]!=b) valid=false; continue; }
            rule.entries[a]=b; queue.emplace_back(b,a);
            for(const auto& op : action) queue.emplace_back(op[a],op[b]);
        }
        if(!valid) continue;
        Exchanges pairs;
        for(std::size_t i=0;i<count;++i) {
            if(rule.entries[i]==count) rule.entries[i]=i;
            if(i<rule.entries[i]) pairs.emplace_back(i,rule.entries[i]);
        }
        found.emplace(pairs,std::move(rule));
    }
    std::vector<TripleTable> result;
    for(auto& [pairs,rule] : found) { (void)pairs; validate_triple_table(g,rule); result.push_back(std::move(rule)); }
    return result;
}

} // namespace wgphysics::research
