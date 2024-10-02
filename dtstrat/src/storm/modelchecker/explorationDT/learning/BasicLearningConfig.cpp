//
// Created by steffi on 05.01.22.
//

#include "storm/modelchecker/explorationDT/learning/BasicLearningConfig.h"
#include "storm/modelchecker/explorationDT/heuristics/NonHeuristic.h"
#include "storm/settings/modules/HeuristicLearningSettings.h"
#include "storm/settings/modules/IOSettings.h"
#include "storm/settings/SettingsManager.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/decision_tree/gini_gain.hpp>
#include <mlpack/methods/decision_tree/best_binary_numeric_split.hpp>
#include <mlpack/methods/decision_tree/all_categorical_split.hpp>
#include <mlpack/methods/decision_tree/all_dimension_select.hpp>
#include <armadillo>

#include "storm/exceptions/NotSupportedException.h"

namespace storm{
namespace modelchecker {
namespace dtstrat {

std::string exec(const char *cmd) {
    std::array<char, 128> buffer{};
    std::string result;
    std::unique_ptr < FILE, decltype(&pclose) > pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

template<typename StateType, typename ValueType>
BasicLearningConfig<StateType, ValueType>::BasicLearningConfig(double const& decisionEpsilon, int const& updateSteps, heuristicTypes const& heuristicType,
                                                               exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInformation)
    : decisionEpsilon(decisionEpsilon), updateSteps(updateSteps), set(false), heuristicType(heuristicType), brtdpInformation(&brtdpInformation), LearningConfig<StateType>(brtdpInformation.getModelName()) {
    if (heuristicType==heuristicTypes::dtcontrolDT) {
        LearningConfig<StateType>::heuristic = new DecisionTreeHeuristic<StateType, ValueType>();
    } else if (heuristicType==heuristicTypes::moduleDT) {
        LearningConfig<StateType>::heuristic = new ModuleDecisionTreeHeuristic<StateType, ValueType>();
    } else if (heuristicType==heuristicTypes::mlpackDT) {
        LearningConfig<StateType>::heuristic = new MLDecisionTreeHeuristic<StateType, ValueType>();
    } else if (heuristicType==heuristicTypes::mlpackModule) {
        LearningConfig<StateType>::heuristic = new MLModuleDecisionTreeHeuristic<StateType, ValueType>();
    } else {
        LearningConfig<StateType>::heuristic = new NonHeuristic<StateType>();
    }
}
template<typename StateType, typename ValueType>
Heuristic<StateType>* BasicLearningConfig<StateType, ValueType>::getHeuristic() {
    return this->heuristic;
}

template<typename StateType, typename ValueType>
bool BasicLearningConfig<StateType, ValueType>::updateNow() {
    if (updateSteps == -1) {
        return false;
    } else if (LearningConfig<StateType>::steps>0 && LearningConfig<StateType>::steps%updateSteps == 0) {
        //STORM_PRINT(LearningConfig<StateType>::steps << " simulations so far.\n");
        return true;
    } else {
        return false;
    }
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::update() {
    if (updateNow()) {
        if (heuristicType==heuristicTypes::dtcontrolDT) {
            updateDtcontrol();
        } else if (heuristicType==heuristicTypes::moduleDT) {
            LearningConfig<StateType>::heuristic = updateModuleDt();
        } else if (heuristicType==heuristicTypes::mlpackDT) {
            updateMLPackDT();
        } else if (heuristicType==heuristicTypes::mlpackModule) {
            updateMlPackModuleDTs();
        }
    }
    step();
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::resetSteps() {
    LearningConfig<StateType>::steps = 0;
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::step() {
    LearningConfig<StateType>::steps++;
}

void printVec(arma::colvec vec) {
    for (int i=0; i<vec.n_rows; i++) std::cout << vec(i) << ",";
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::updateMlPackModuleDTs() {
    auto configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    bool onlyReachable = configSettings.getOnlyReachable();
    arma::Row<long unsigned int> labels;
    arma::Mat<int> data;
    arma::Row<long unsigned int> modules;
    if (onlyReachable) {
        data = brtdpInformation->forwardPassStateValuations(labels, &modules);
    } else {
        data = brtdpInformation->getTransformedStateValuations(labels, &modules);
    }
    //arma::mat double_data = arma::conv_to<arma::mat>::from(data);
    LearningConfig<StateType>::heuristic = getNewMLModuleDTHeuristic(labels, data, modules);
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::updateMLPackDT() {
    auto configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    bool onlyReachable = configSettings.getOnlyReachable();
    arma::Row<long unsigned int> labels;
    arma::Mat<int> data;
    if (onlyReachable) {
        data = brtdpInformation->forwardPassStateValuations(labels);
    } else {
        data = brtdpInformation->getTransformedStateValuations(labels);
    }
    arma::mat double_data = arma::conv_to<arma::mat>::from(data);
    //arma::Row<long unsigned int> uniqueActions = arma::unique(labels);
    //mlpack::tree::DecisionTree<mlpack::tree::GiniGain,mlpack::tree::BestBinaryNumericSplit,mlpack::tree::AllCategoricalSplit,mlpack::tree::AllDimensionSelect, double,false> test_tree(double_data, labels, brtdpInformation->getNumberOfActions(), 1);
    LearningConfig<StateType>::heuristic = new MLDecisionTreeHeuristic<StateType, ValueType>(double_data, labels, brtdpInformation->getNumberOfActions(), brtdpInformation->getExplorationInformation(), decisionEpsilon);
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::updateDtcontrol() {
    auto configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    bool onlyReachable = configSettings.getOnlyReachable();
    // Output the file
    std::string fileName = brtdpInformation->getModelName() + ".csv";
    std::ofstream ofs(fileName, std::ofstream::trunc);
    brtdpInformation->printDTcontrolHeader(ofs, false);
    brtdpInformation->printData(ofs, nullptr, onlyReachable);
    ofs.close();

    // call DTcontrol
    std::string jsonString = callDTcontrol(fileName);
    ofs = std::ofstream("DTControl.json");
    ofs << jsonString;
    ofs.close();

    // get the new decision tree
    LearningConfig<StateType>::heuristic = getNewDT(jsonString, true);
}

template<typename StateType, typename ValueType>
std::string BasicLearningConfig<StateType, ValueType>::callDTcontrol(std::string const& fileName) {
    std::string cmd = ". venv/bin/activate; "
        "dtcontrol "
        "--input " + fileName +
        " --use-preset minnorm "
        "--output stdout:json "
        "--rerun "
        "2>\\dev\\null";
    const char *run_cmd = cmd.c_str();
    const std::string &result = exec(run_cmd);
    STORM_LOG_ASSERT(result != "", "dtcontrol not installed? : \nPlease run the following commands"
                     " to install dtcontrol in the \"Working directory\" of storm:\n"
                     "#   python -m venv venv\n"
                     "#   source venv/bin/activate\n"
                     "#   pip install dtcontrol\n"
    );
    auto start = result.find("START") + 5;
    auto end = result.find("END");
    const std::string jsonStringOfDT = result.substr(start, end - start);
    return jsonStringOfDT;
}

template<typename StateType, typename ValueType>
DecisionTreeHeuristic<StateType, ValueType>* BasicLearningConfig<StateType, ValueType>::getNewDT(std::string const& someString, bool const& fromJson) {
   return new DecisionTreeHeuristic<StateType, ValueType>(someString, &(brtdpInformation->getExplorationInformation()), decisionEpsilon ,fromJson);
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::prepareModuleDT() {
    auto configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    bool onlyReachable = configSettings.getOnlyReachable();
    // Output the file
    std::string fileName = brtdpInformation->getModelName() + ".csv";
    std::ofstream ofs(fileName, std::ofstream::trunc);
    brtdpInformation->printModuleDTHeader(ofs);
    brtdpInformation->printData(ofs, nullptr, onlyReachable, true, true);
    ofs.close();

    callSklearnMulti(fileName, useOnlyStateVariablesOfModuleItself);
}

template<typename StateType, typename ValueType>
std::string BasicLearningConfig<StateType, ValueType>::prepareModuleDTwoCall() {
    auto configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    bool onlyReachable = configSettings.getOnlyReachable();
    // Output the file
    std::string fileName = brtdpInformation->getModelName() + ".csv";
    std::ofstream ofs(fileName, std::ofstream::trunc);
    brtdpInformation->printModuleDTHeader(ofs);
    brtdpInformation->printData(ofs, nullptr, onlyReachable, true, true);
    ofs.close();
    return fileName;
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>* BasicLearningConfig<StateType, ValueType>::updateModuleDt() {
    prepareModuleDT();
    return getNewMultiDTHeuristic();
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::callSklearnMulti(std::string const& fileName, bool const& useOnlyStateVariablesOfModuleItself) const {
    // clear up -> delete all old module DTs
    for (uint i = 0; i < brtdpInformation->getNumberOfModules(); i++) {
        std::string filename = brtdpInformation->getModelName() + "_moduleDT_" + std::to_string(i) + ".dot";
        std::ifstream f(filename.c_str());
        if (f.good()) {
            remove(filename.c_str());
        }
    }
    std::string cmd = ". venv/bin/activate; "
        "python module_dt.py "
        + fileName;
    if (useOnlyStateVariablesOfModuleItself) cmd += " statesPerModule";
    const char *run_cmd = cmd.c_str();
    const std::string &result = exec(run_cmd);
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>* BasicLearningConfig<StateType, ValueType>::getNewMultiDTHeuristic() {
    return new ModuleDecisionTreeHeuristic<StateType, ValueType>(brtdpInformation, decisionEpsilon, useOnlyStateVariablesOfModuleItself);
}

template<typename StateType, typename ValueType>
MLModuleDecisionTreeHeuristic<StateType, ValueType>* BasicLearningConfig<StateType, ValueType>::getNewMLModuleDTHeuristic(arma::Row<long unsigned int>& labels,
                                                                                                                          arma::Mat<int>& data,
                                                                                                                          arma::Row<long unsigned int>& modules) {
    return new MLModuleDecisionTreeHeuristic<StateType, ValueType>(brtdpInformation, decisionEpsilon, useOnlyStateVariablesOfModuleItself, labels, data, modules);
}

template<typename StateType, typename ValueType>
std::string BasicLearningConfig<StateType, ValueType>::getModelName() const {
    return brtdpInformation->getModelName();
}

template<typename StateType, typename ValueType>
uint BasicLearningConfig<StateType, ValueType>::getNumberOfModules() const {
    return brtdpInformation->getNumberOfModules();
}

template<typename StateType, typename ValueType>
void BasicLearningConfig<StateType, ValueType>::setUpdateNow() {
    LearningConfig<StateType>::steps = updateSteps-1;
}

template class BasicLearningConfig<uint32_t, double>;
}
}
}
