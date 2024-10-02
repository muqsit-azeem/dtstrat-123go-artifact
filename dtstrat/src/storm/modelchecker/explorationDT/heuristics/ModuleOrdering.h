//
// Created by steffi on 12.01.22.
//

#ifndef STORM_MODULEORDERING_H
#define STORM_MODULEORDERING_H

#include "storm/modelchecker/explorationDT/heuristics/ModuleHeuristic.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
class ModuleOrdering : public ModuleHeuristic<StateType> {
   public:
    ModuleOrdering();
    ModuleOrdering(uint const& numModules, bool const& random=false);
    ModuleOrdering(std::vector<uint> moduleOrdering);
    bool useHeuristic() override;

    uint getDecision() const override;

    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(std::string const& filename) const override;

   private:
    std::vector<uint> moduleOrdering;

};
}
}
}

#endif  // STORM_MODULEORDERING_H
