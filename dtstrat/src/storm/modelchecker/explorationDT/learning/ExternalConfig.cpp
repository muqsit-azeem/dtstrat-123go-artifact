//
// Created by steffi on 12.01.22.
//

#include "storm/modelchecker/explorationDT/learning/ExternalConfig.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleDecisionTreeHeuristic.h"


namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
ExternalConfig<StateType, ValueType>::ExternalConfig(const std::string &modelName, exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInfo, double const& explEpsilon) : LearningConfig<StateType>(modelName), brtdpInfo(&brtdpInfo) {
    LearningConfig<StateType>::heuristic = new ModuleDecisionTreeHeuristic<StateType, ValueType>(&brtdpInfo, explEpsilon, false);
}

template class ExternalConfig<uint32_t, double>;

}
}
}