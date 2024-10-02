//
// Created by steffi on 05.01.22.
//

#ifndef STORM_HEURISTIC_H
#define STORM_HEURISTIC_H

#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/decision_tree/gini_gain.hpp>

namespace storm {
namespace modelchecker {
namespace dtstrat {
enum heuristicTypes {noneHeuristic, dtcontrolDT, moduleDT, mlpackDT, mlpackModule};

template<typename StateType>
class Heuristic {
   public:
    virtual bool useHeuristic()=0;
    virtual ~Heuristic() = default;
    virtual unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked)=0;
    virtual std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId)=0;
    virtual void printHeuristic(std::string const& filename) const=0;

    virtual void setExplorationEpsilon(double const& explEps)=0;
    virtual mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* getDT() {return nullptr;};
};
}
}
}

#endif  // STORM_HEURISTIC_H
