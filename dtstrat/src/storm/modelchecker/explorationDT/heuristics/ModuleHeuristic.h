//
// Created by steffi on 12.01.22.
//

#ifndef STORM_MODULEHEURISTIC_H
#define STORM_MODULEHEURISTIC_H

#include "storm/modelchecker/explorationDT/heuristics/Heuristic.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

enum moduleHeuristicTypes {ordering, decisiontree, orderingrandom, orderingtree};

template<typename StateType>
class ModuleHeuristic : public Heuristic<StateType> {
   public:
    virtual uint getDecision() const=0;
    void setExplorationEpsilon(double const& explEps) override {};
};
}
}
}
#endif  // STORM_MODULEHEURISTIC_H
