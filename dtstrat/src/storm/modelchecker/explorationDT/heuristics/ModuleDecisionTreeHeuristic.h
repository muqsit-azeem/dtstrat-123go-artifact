//
// Created by steffi on 08.01.22.
//

#ifndef STORM_MODULEDECISIONTREEHEURISTIC_H
#define STORM_MODULEDECISIONTREEHEURISTIC_H

#include "storm/modelchecker/explorationDT/heuristics/DecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleHeuristic.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class ModuleDecisionTreeHeuristic : public DecisionTreeHeuristic<StateType, ValueType> {
   public:
    ModuleDecisionTreeHeuristic() {};
    ModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo, double const& explorationEpsilon, bool const& useOnlyStateVariablesOfModuleItself);
    ModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo, double const& explorationEpsilon,
                                bool const& useOnlyStateVariablesOfModuleItself, moduleHeuristicTypes moduleHeuristic, std::vector<uint> moduleOrdering,
                                std::map<uint, uint> mapping, std::vector<const storm::storage::DecisionTree<uint>*>& treesOfThisRound);

    ~ModuleDecisionTreeHeuristic();

    ModuleDecisionTreeHeuristic& operator=(ModuleDecisionTreeHeuristic<StateType, ValueType> const& copy);

    static std::vector<const storm::storage::DecisionTree<uint>*> getTrees(uint const& numberOfModules, std::string const& modelName);

    ModuleHeuristic<StateType>* getModuleHeuristic(moduleHeuristicTypes moduleHeuristic, std::vector<uint> ordering={});
    void getModuleToDTMapping();

    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(const std::string &filename) const override;
    void updateDT(std::string const& someString, bool const& fromJson) override;
    bool useHeuristic() override;



   protected:
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo;
    std::vector<const storm::storage::DecisionTree<uint>*> decisionTrees;
    std::map<uint, uint> moduleDTMapping;
    ModuleHeuristic<StateType>* moduleHeuristic;
    bool useOnlyStateVariablesOfModuleItself = true;

   private:
    StateValuationDT reorder(StateValuationDT const& inputValues, uint const& module);
};

}
}
}

#endif  // STORM_MODULEDECISIONTREEHEURISTIC_H
