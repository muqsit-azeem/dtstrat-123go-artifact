//
// Created by steffi on 05.01.22.
//

#ifndef STORM_NONLEARNINGCONFIG_H
#define STORM_NONLEARNINGCONFIG_H
#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"

namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
class NonLearningConfig : public LearningConfig<StateType> {
   public:
    NonLearningConfig(std::string const& modelName);
    NonLearningConfig(Heuristic<StateType>* heuristic, std::string const& modelName);
};
}
}
}
#endif  // STORM_NONLEARNINGCONFIG_H
