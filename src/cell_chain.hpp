#pragma once

#include "equivariant_pair_search.hpp"

namespace wgphysics::research {

// Exact compiled holonomy dynamics on a fixed chain of triangular cells.
// Connector i transports from cell i+1 to cell i and never changes in this model.
// This fixture tests dynamics; its one-dimensional geometry is supplied.
class CellChain {
public:
    struct Event {
        std::size_t left{};
        std::vector<std::size_t> parents;
        std::size_t depth{};
        std::pair<uint16_t,uint16_t> before, after;
    };

    CellChain(const FiberGraph& fiber, std::vector<uint16_t> holonomies,
              std::vector<uint16_t> connectors, PairInteraction rule)
        : tables_(fiber.automorphisms()), values_(std::move(holonomies)),
          connectors_(std::move(connectors)), last_writer_(values_.size()) {
        if (values_.empty() || connectors_.size()+1 != values_.size())
            throw std::invalid_argument("chain requires n cells and n-1 connectors");
        for (const auto v : values_) validate(v);
        for (const auto t : connectors_) validate(t);
        for (const auto direction : {HurwitzDirection::Forward,HurwitzDirection::Inverse}) {
            auto& table = direction == HurwitzDirection::Forward ? forward_ : inverse_;
            for (const auto& a : tables_.elements()) for (const auto& b : tables_.elements()) {
                const auto result = pair_targets(a,b,rule,direction);
                table.emplace_back(tables_.index_of(result.first),tables_.index_of(result.second));
            }
        }
        const auto order = tables_.order();
        sectors_.resize(order);
        for (uint16_t a=0; a<order; ++a) {
            sectors_[a]=a;
            for (uint16_t g=0; g<order; ++g) sectors_[a]=std::min(sectors_[a],conjugate(g,a));
        }
    }

    const std::vector<uint16_t>& values() const { return values_; }

    CellChain(const FiberGraph& fiber,std::vector<uint16_t> holonomies,
              std::vector<uint16_t> connectors,const PairTable& rule)
        : CellChain(fiber,std::move(holonomies),std::move(connectors),PairInteraction::Hurwitz) {
        validate_pair_table(tables_,rule);
        for (std::size_t i=0;i<rule.size();++i) {
            forward_[i]={rule[i]/tables_.order(),rule[i]%tables_.order()};
            inverse_[rule[i]]={i/tables_.order(),i%tables_.order()};
        }
    }
    const std::vector<Event>& events() const { return events_; }
    const AutomorphismTables& tables() const { return tables_; }
    uint16_t sector(uint16_t value) const { validate(value); return sectors_[value]; }
    uint16_t identity() const {
        return tables_.index_of(Permutation::identity(tables_.fiber_size()));
    }

    void update(std::size_t left, HurwitzDirection direction = HurwitzDirection::Forward,
                bool record = true) {
        if (left+1 >= values_.size()) throw std::out_of_range("chain interaction is outside cells");
        const auto a=values_[left], b=values_[left+1], t=connectors_[left];
        const auto based_b=conjugate(t,b);
        const auto& table=direction == HurwitzDirection::Forward ? forward_ : inverse_;
        const auto [new_a,new_based_b]=table[a*tables_.order()+based_b];
        const auto new_b=conjugate(tables_.inverse(t),new_based_b);
        if (record) {
            if (unrecorded_) throw std::logic_error("cannot append provenance after unrecorded updates");
            Event event{left,{},0,{a,b},{new_a,new_b}};
            for (const auto cell : {left,left+1}) if (last_writer_[cell]) {
                const auto parent=*last_writer_[cell];
                if (std::find(event.parents.begin(),event.parents.end(),parent)==event.parents.end())
                    event.parents.push_back(parent);
                event.depth=std::max(event.depth,events_[parent].depth+1);
            }
            last_writer_[left]=last_writer_[left+1]=events_.size();
            events_.push_back(std::move(event));
        } else unrecorded_=true;
        values_[left]=new_a; values_[left+1]=new_b;
    }

    // The explicitly based ordered product transported to cell zero.
    uint16_t total() const {
        auto result=identity(), transport=identity();
        for (std::size_t i=0;i<values_.size();++i) {
            result=tables_.multiply(result,conjugate(transport,values_[i]));
            if (i<connectors_.size()) transport=tables_.multiply(transport,connectors_[i]);
        }
        return result;
    }

    // Word-rule evolution preserves this subgroup; general searched tables need
    // not. Changing the root frame conjugates the subgroup.
    std::set<uint16_t> generated_subgroup() const {
        std::set<uint16_t> generators;
        auto transport=identity();
        for (std::size_t i=0;i<values_.size();++i) {
            generators.insert(conjugate(transport,values_[i]));
            if (i<connectors_.size()) transport=tables_.multiply(transport,connectors_[i]);
        }
        std::set<uint16_t> result{identity()};
        std::vector<uint16_t> frontier{identity()};
        for (std::size_t i=0;i<frontier.size();++i) for (const auto generator : generators) {
            const auto next=tables_.multiply(frontier[i],generator);
            if (result.insert(next).second) frontier.push_back(next);
        }
        return result;
    }

private:
    void validate(uint16_t v) const {
        if (v>=tables_.order()) throw std::out_of_range("element is outside derived fiber group");
    }
    uint16_t conjugate(uint16_t g,uint16_t v) const {
        return tables_.multiply(g,tables_.multiply(v,tables_.inverse(g)));
    }
    AutomorphismTables tables_;
    std::vector<uint16_t> values_,connectors_,sectors_;
    std::vector<std::optional<std::size_t>> last_writer_;
    std::vector<std::pair<uint16_t,uint16_t>> forward_,inverse_;
    std::vector<Event> events_;
    bool unrecorded_{};
};

} // namespace wgphysics::research
