//
// Created by steffi on 05.01.22.
//

#ifndef STORM_DECISIONTREEHEURISTIC_H
#define STORM_DECISIONTREEHEURISTIC_H
#include <string>

#include "storm/modelchecker/explorationDT/heuristics/Heuristic.h"
#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include "storm/storage/DecisionTree.h"
#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class DecisionTreeHeuristic : public Heuristic<StateType> {
   public:
    DecisionTreeHeuristic() = default;
    DecisionTreeHeuristic(exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation, double const& explorationEpsilon);
    DecisionTreeHeuristic(std::string const& file, exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation, double const& explorationEpsilon, bool const& json=false);
    ~DecisionTreeHeuristic();

    DecisionTreeHeuristic operator=(DecisionTreeHeuristic<StateType, ValueType> const& copy);

    virtual bool useHeuristic() override;
    /*!
     * should return the storm action id
     * @param inputValues
     * @return
     */
    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(const std::string &filename) const override;
    virtual void updateDT(std::string const& someString, bool const& fromJson);

    void setExplorationEpsilon(double const& explEps) override;
    mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* getDT();


   protected:
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation;
    double explorationEpsilon;

   private:
    storm::storage::DecisionTree<uint>* decisionTree;
};

}
}
}
#endif  // STORM_DECISIONTREEHEURISTIC_H
