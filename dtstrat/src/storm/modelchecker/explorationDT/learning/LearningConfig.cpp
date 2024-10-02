//
// Created by steffi on 05.01.22.
//
#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/HeuristicLearningSettings.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
bool LearningConfig<StateType>::useHeuristic() {
    if (heuristic) {
        return heuristic->useHeuristic();
    } else {
        return false;
    }
}

template<typename StateType>
uint LearningConfig<StateType>::getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) {
    if (heuristic) {
        return heuristic->getDecision(inputValues, currentStateId, worked);
    } else {
        return 0;
    }
}

template<typename StateType>
Heuristic<StateType>* LearningConfig<StateType>::getHeuristic() const {
    return heuristic;
}

template<typename StateType>
std::vector<uint> LearningConfig<StateType>::getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) {
    if (heuristic) {
        return heuristic->getAllDecisions(inputValues, currentStateId);
    } else {
        return {};
    }
}

template<typename StateType>
void LearningConfig<StateType>::printHeuristic(std::string const& fileName) {
    heuristic->printHeuristic(fileName);
}

template<typename StateType>
void LearningConfig<StateType>::update() {
    steps++;
}

template<typename StateType>
void LearningConfig<StateType>::firstLog() {
    bool addTimeoutLine = false;
    const std::string filename = "ExperimentsLog.csv";
    std::ifstream fs;
    fs.open(filename.c_str(), std::fstream::in);
    if(fs.is_open()) {
        //Got to the last character before EOF
        fs.seekg(-1, std::ios_base::end);
        if (fs.peek() == ',') {
            addTimeoutLine = true;
        }
    }
    fs.close();

    std::ofstream ofs("ExperimentsLog.csv",std::ios_base::app);
    if (addTimeoutLine) {
        ofs << "-1,-1";
    }
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
    std::string str(buffer);
    ofs << "\n" << str << "," << modelName << "," << learningType << "," << heuristicType << ",";
    ofs.close();
}

template<typename StateType>
void LearningConfig<StateType>::finalLog(double const& elapsed) {
    STORM_PRINT(steps << " simulations in total.\n");
    std::ofstream ofs("ExperimentsLog.csv",std::ios_base::app);
    //time_t rawtime;
    //struct tm * timeinfo;
    //char buffer[80];

    //time (&rawtime);
    //timeinfo = localtime(&rawtime);

    //strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
    //std::string str(buffer);
    //ofs << "\n" << str << "," << modelName << "," << learningType << "," << heuristicType << "," << steps << "," << elapsed ;
    ofs << steps << "," << elapsed;
    ofs.close();
}

template<typename StateType>
LearningConfig<StateType>::LearningConfig(std::string const& modelName) : modelName(modelName) {
    auto configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    configTypes configType = configSettings.getLearningConfigType();
    if (configType == configTypes::basic) {
        learningType = "basic";
    } else if (configType == configTypes::genetic) {
        learningType = "genetic";
    } else if (configType == configTypes::external) {
        learningType = "external";
    } else if (configType == configTypes::headstart) {
        learningType = "headstart";
    }
    else {
        learningType = "none";
    }
    auto heuristicType = configSettings.getHeuristicType();
    if (heuristicType==heuristicTypes::dtcontrolDT) {
        this->heuristicType = "dtcontrol";
    } else if (heuristicType==heuristicTypes::moduleDT) {
        this->heuristicType = "module";
    } else if (heuristicType==heuristicTypes::mlpackDT) {
        this->heuristicType = "mlpackDT";
    } else if (heuristicType==heuristicTypes::mlpackModule) {
        this->heuristicType = "mlpackModuleDt";
    } else {
        this->heuristicType = "none";
    }
}
template<typename StateType>
LearningConfig<StateType>::~LearningConfig() {
    delete heuristic;
}

template<typename StateType>
LearningConfig<StateType>& LearningConfig<StateType>::operator=(const LearningConfig<StateType>& copy) {
    if (this == &copy) return *this;
    else {
        delete heuristic;
        heuristic = copy.heuristic;
        steps = copy.steps;
        modelName = copy.modelName;
        learningType = copy.learningType;
        heuristicType = copy.heuristicType;
        return *this;
    }
}

template<typename StateType>
void LearningConfig<StateType>::prepareBoundLog() {
    std::ofstream ofs("bound.csv",std::ios_base::app);
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
    std::string str(buffer);
    ofs << "\n" << str << "," << modelName << "," << learningType << "," << heuristicType << ",";
    ofs.close();
}

template<typename StateType>
void LearningConfig<StateType>::writeBoundLog(double bound) {
    std::ofstream ofs("bound.csv",std::ios_base::app);
    ofs << bound << ",";
    ofs.close();
}

template<typename StateType>
void LearningConfig<StateType>::setUpdateNow() {
    // intentionally left empty
}

template class LearningConfig<uint32_t>;
}
}
}