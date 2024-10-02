//
// Created by steffi on 05.09.22.
//
#include "storm/modelchecker/learning/DTInterface.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {
    template<typename StateType, typename ValueType>
    DTInterface<StateType, ValueType>::DTInterface(nlohmann2::json tree) {
        this->decisionTreeST = new storm::storage::DecisionTree<StateType>(tree);
    }

    template<typename StateType, typename ValueType>
    DTInterface<StateType, ValueType>::DTInterface(std::string const &filename) {
        this->decisionTreeST = new storm::storage::DecisionTree<StateType>(filename);
    }

    template<typename StateType, typename ValueType>
    DTInterface<StateType, ValueType>::DTInterface(arma::Mat<int> const& data, arma::Row<long unsigned int> const& labels, uint const& numberOfActions) {
        this->decisionTreeML =  new mlpack::tree::DecisionTree<mlpack::tree::GiniGain>(data, labels, numberOfActions, 1);
    }

    template<typename StateType, typename ValueType>
    DTInterface<StateType, ValueType>::DTInterface(mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* dt) {
        this->decisionTreeML = dt;
    }

    template<typename StateType, typename ValueType>
    uint DTInterface<StateType, ValueType>::getDecision(StateValuationDT const& stateValuation, StateType const& currentStateId, bool& worked) const {
        if (decisionTreeST) {
            auto allActions = decisionTreeST->getDecisions(stateValuation);
            if (allActions.empty()) return 0;
            else {
                worked = true;
                return allActions[0];
            }
        } else if (decisionTreeML) {
            auto matrix = stateValuation.transform();
            arma::Row<long unsigned int> predictions;
            arma::mat probabilities;
            decisionTreeML->Classify(matrix, predictions, probabilities);
            std::vector<double> probs(probabilities.n_rows);
            for(uint i=0; i<probs.size(); i++) {
                probs[i] = probabilities(i);
            }
            std::vector<uint> result(probs.size());
            std::iota(result.begin(), result.end(),0);
            sort(result.begin(), result.end(), [&](uint i, uint j){return probs[i]>probs[j];});
            if (result.empty()) return 0;
            else {
                worked = true;
                return result[0];
            }
        } else {
            return 0;
        }
    }

    template<typename StateType, typename ValueType>
    void DTInterface<StateType, ValueType>::exportDT(std::string const& filename) const {
        if (decisionTreeST) {
            std::ofstream ofs(filename, std::ofstream::trunc);
            ofs << decisionTreeST->getDotString();
            ofs.close();
        }
        else if (decisionTreeML) {
            std::ofstream ofs(filename);
            ofs << decisionTreeML->getDot();
            ofs.close();
        }
    }

    template<typename StateType, typename ValueType>
    uint DTInterface<StateType, ValueType>::size() const {
        if (decisionTreeST) {
            decisionTreeST->size();
        }
        else if (decisionTreeML) {
            decisionTreeML->size();
        }
        else return -1;
    }

    template class DTInterface<uint32_t, double>;
}
}
}