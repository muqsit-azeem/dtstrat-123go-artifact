//
// Created by steffi on 05.01.22.
//

#include "storm/modelchecker/explorationDT/learning/NonLearningConfig.h"
#include "storm/modelchecker/explorationDT/heuristics/NonHeuristic.h"

namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
NonLearningConfig<StateType>::NonLearningConfig(std::string const& modelName) : LearningConfig<StateType>(modelName){
    LearningConfig<StateType>::heuristic = new NonHeuristic<StateType>();
}

template<typename StateType>
NonLearningConfig<StateType>::NonLearningConfig(Heuristic<StateType>* heuristic, std::string const& modelName) : LearningConfig<StateType>(modelName){
    LearningConfig<StateType>::heuristic = heuristic;
}

template class NonLearningConfig<uint32_t>;
}
}
}
