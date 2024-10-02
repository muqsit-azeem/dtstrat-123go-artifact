//
// Created by steffi on 05.09.22.
//

#ifndef STORM_DTINTERFACE_H
#define STORM_DTINTERFACE_H
#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/decision_tree/gini_gain.hpp>
#include <armadillo>
#include "storm/storage/DecisionTree.h"
#include "storm/modelchecker/explorationDT/StateValuationDT.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class DTInterface {
   public:
    DTInterface() = default;
    DTInterface(nlohmann2::json tree);
    DTInterface(std::string const &filename);
    DTInterface(arma::Mat<int> const& data, arma::Row<long unsigned int> const& labels, uint const& numberOfActions);
    DTInterface(mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* decisionTree);

    uint getDecision(StateValuationDT const& stateValuation, StateType const& currentStateId, bool& worked) const;

    void exportDT(std::string const& filename) const;

    uint size() const;

   private:
    mlpack::tree::DecisionTree<mlpack::tree::GiniGain> *decisionTreeML;
    storm::storage::DecisionTree<uint> *decisionTreeST;
};
}
}
}
#endif  // STORM_DTINTERFACE_H
