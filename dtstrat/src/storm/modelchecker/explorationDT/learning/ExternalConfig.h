//
// Created by steffi on 12.01.22.
//

#ifndef STORM_EXTERNALCONFIG_H
#define STORM_EXTERNALCONFIG_H

#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"

namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class ExternalConfig : public LearningConfig<StateType> {
   public:
    ExternalConfig(std::string const& modelName, exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInfo, double const& eps);

   private:
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo;
};
}
}
}


#endif  // STORM_EXTERNALCONFIG_H
