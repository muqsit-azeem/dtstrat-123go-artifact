//
// Created by steffi on 12.01.22.
//

#ifndef STORM_MLDecisionTreeForModules_H
#define STORM_MLDecisionTreeForModules_H

#include "storm/modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleHeuristic.h"
#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include "storm/storage/DecisionTree.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/decision_tree/gini_gain.hpp>



namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class MLDecisionTreeForModules : public ModuleHeuristic<StateType> {
   public:
    MLDecisionTreeForModules() = default;
    MLDecisionTreeForModules(exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInfo, arma::Mat<int>& data, arma::Row<long unsigned int>& modules);
    ~MLDecisionTreeForModules();

    MLDecisionTreeForModules& operator=(MLDecisionTreeForModules<StateType, ValueType> const& copy);

    virtual bool useHeuristic() override;
    /*!
     * should return the storm action id
     * @param inputValues
     * @return
     */
    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(const std::string &filename) const override;
    uint getDecision() const override;

   private:
    mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* decisionTree;
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo;
};
template<typename StateType, typename ValueType>
uint MLDecisionTreeForModules<StateType, ValueType>::getDecision() const {
    return 0;
}
}
}
}

#endif  // STORM_MLDecisionTreeForModules_H
