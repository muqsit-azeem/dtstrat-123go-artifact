//
// Created by steffi on 05.01.22.
//

#ifndef STORM_NONHEURISTIC_H
#define STORM_NONHEURISTIC_H

#include "storm/modelchecker/explorationDT/heuristics/Heuristic.h"
#include "storm/storage/sparse/StateValuations.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
class NonHeuristic : public Heuristic<StateType> {
   public:
    NonHeuristic() = default;
    bool useHeuristic() override;
    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(std::string const& filename) const override;
    void setExplorationEpsilon(double const& explEps);
    mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* getDT();
};

}
}
}

#endif  // STORM_NONHEURISTIC_H
