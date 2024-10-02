//
// Created by steffi on 05.01.22.
//

#ifndef STORM_MLDecisionTreeHeuristic_H
#define STORM_MLDecisionTreeHeuristic_H
#include <string>

#include "storm/modelchecker/explorationDT/heuristics/Heuristic.h"
#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include "storm/storage/DecisionTree.h"
#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/decision_tree/gini_gain.hpp>
#include <armadillo>

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class MLDecisionTreeHeuristic : public Heuristic<StateType> {
   public:
    MLDecisionTreeHeuristic() = default;
    MLDecisionTreeHeuristic(exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation, double const& explorationEpsilon);
    MLDecisionTreeHeuristic(
        arma::mat& data, arma::Row<long unsigned int>& labels, uint numberOfActions,
        exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>& explorationInformation, double const& explorationEpsilon);
    ~MLDecisionTreeHeuristic();

    MLDecisionTreeHeuristic operator=(MLDecisionTreeHeuristic<StateType, ValueType> const& copy);

    virtual bool useHeuristic() override;
    /*!
     * should return the storm action id
     * @param inputValues
     * @return
     */
    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(const std::string &filename) const override;
    virtual void updateDT(arma::mat& data, arma::Row<long unsigned int>& labels, uint numberOfActions);

    void setExplorationEpsilon(double const& explEps) override;

    mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* getDT();

   protected:
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation;
    double explorationEpsilon;

   private:
    mlpack::tree::DecisionTree<mlpack::tree::GiniGain> *decisionTree;
};

}
}
}
#endif  // STORM_MLDecisionTreeHeuristic_H
