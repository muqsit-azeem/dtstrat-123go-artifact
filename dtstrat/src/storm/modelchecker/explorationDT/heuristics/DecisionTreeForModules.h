//
// Created by steffi on 12.01.22.
//

#ifndef STORM_DECISIONTREEFORMODULES_H
#define STORM_DECISIONTREEFORMODULES_H

#include "storm/modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleHeuristic.h"
#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include "storm/storage/DecisionTree.h"


namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class DecisionTreeForModules : public ModuleHeuristic<StateType> {
   public:
    DecisionTreeForModules() = default;
    DecisionTreeForModules(exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInfo, std::string const& fileName);
    ~DecisionTreeForModules();

    DecisionTreeForModules& operator=(DecisionTreeForModules<StateType, ValueType> const& copy);

    virtual bool useHeuristic() override;
    /*!
     * should return the storm action id
     * @param inputValues
     * @return
     */
    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(const std::string &filename) const override;
    [[nodiscard]] uint getDecision() const override {return 0;};

   private:
    storm::storage::DecisionTree<uint>* decisionTree;
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo;
};

}
}
}

#endif  // STORM_DECISIONTREEFORMODULES_H
