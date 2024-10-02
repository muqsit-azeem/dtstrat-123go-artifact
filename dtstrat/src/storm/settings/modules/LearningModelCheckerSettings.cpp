//
// Created by steffi on 27.07.22
//
#include "storm/settings/modules/LearningModelCheckerSettings.h"

#include "storm/settings/Argument.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/SettingsManager.h"

#include <boost/algorithm/string/replace.hpp>

namespace storm {
namespace settings {
namespace modules {

const std::string LearningModelCheckerSettings::moduleName = "learningMC";
const std::string LearningModelCheckerSettings::decisiontree = "decisiontree";
const std::string LearningModelCheckerSettings::datafile = "datafile";
const std::string LearningModelCheckerSettings::learnparameter = "learnParameter";
const std::string LearningModelCheckerSettings::evluationmethod = "evaluationMethod";
const std::string LearningModelCheckerSettings::childrenfiles = "childrenFiles";
const std::string LearningModelCheckerSettings::treelearner = "dtLearning";
const std::string LearningModelCheckerSettings::learnparametername = "paramName";
const std::string LearningModelCheckerSettings::smcruns = "smcruns";
const std::string LearningModelCheckerSettings::maxsteps = "smcmaxsteps";
const std::string LearningModelCheckerSettings::onlyReachable = "reachable";
const std::string LearningModelCheckerSettings::numStepsDTStrat = "stepsDTStrat";
const std::string LearningModelCheckerSettings::subProgramsAndTasks = "childTasks";
const std::string LearningModelCheckerSettings::pacEpsilon = "pacEpsilon";
const std::string LearningModelCheckerSettings::pacDelta = "pacDelta";
const std::string LearningModelCheckerSettings::hyperparamType = "hypertype";
const std::string LearningModelCheckerSettings::hyperparamTime = "hypertime";
const std::string LearningModelCheckerSettings::semanticActions = "semanticAction";

LearningModelCheckerSettings::LearningModelCheckerSettings() : ModuleSettings(moduleName) {
    this->addOption(storm::settings::OptionBuilder(moduleName, decisiontree, false, "If set, gives a decision tree for evaluation (DT>data>parameter).")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "Path to the decision-tree file (.dot, .json).")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(
        storm::settings::OptionBuilder(moduleName, datafile, false, "If set, gives the data-file from which to learn a DT (DT>data>parameter).")
            .addArgument(
                storm::settings::ArgumentBuilder::createStringArgument("name", "Path to the data-file (.csv), states and action in the same line.")
                                        .setDefaultValueString("").build())
            .setIsAdvanced()
            .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, learnparameter, false, "If set, describes the parameters from which to learn from.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "Should be in form of x,y,z (seperated by commata).")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, evluationmethod, false, "Describes the way how to evaluate [statistical, building the MC].")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the heuristic, in [statistical,full,none]")
                                         .setDefaultValueString("statistical").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, childrenfiles, false, "If set, describes a dummy-file for the parameters to be placed in. File must have a @ instead of the parameter.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The path to the files.")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, treelearner, false, "Describes the method how to get the DT [sklearn,dtcontrol,mlpack].")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "")
                                         .setDefaultValueString("dtcontrol").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, learnparametername, false, "If set, describes the name of the parameter that should be changed.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "Should be the exact name of the parameter.")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, smcruns, false, "Number of runs for SMC.")
                        .addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("value","")
                                         .setDefaultValueUnsignedInteger(1000).build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, maxsteps, false, "Maximal number of steps per run for SMC.")
                        .addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("value","")
                                         .setDefaultValueUnsignedInteger(10000).build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, onlyReachable, false, "If set, use only reachable states under current strategy for export.")
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, numStepsDTStrat, false, "Number of simulations to learn form DTStrat.")
                        .addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("value","")
                                         .setDefaultValueUnsignedInteger(100).build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, subProgramsAndTasks, false, "Describes the children files and properties. Expected structure FILENAME#PROPERTY%FILENAME2#PROPERTY2.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "..")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, pacEpsilon, false, "Epsilon for PAC")
                        .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("value","")
                                         .setDefaultValueDouble(1e-1).build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, pacDelta, false, "Delta for PAC")
                        .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("value","")
                                         .setDefaultValueDouble(1e-1).build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, hyperparamType, false, "tbd")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "..")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, hyperparamTime, false, "tbd")
                        .addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("value","")
                                         .setDefaultValueUnsignedInteger(600).build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, semanticActions, false, "If set, try to use the semantic actions.")
                        .setIsAdvanced()
                        .build());
}

std::string LearningModelCheckerSettings::getDecisionTreeFile() const {
    return this->getOption(decisiontree).getArgumentByName("name").getValueAsString();
}


std::string LearningModelCheckerSettings::getDataFile() const {
    return this->getOption(datafile).getArgumentByName("name").getValueAsString();
}

std::string LearningModelCheckerSettings::getLearnHyperParamType() const {
    return this->getOption(hyperparamType).getArgumentByName("name").getValueAsString();
}

std::vector<int> LearningModelCheckerSettings::getLearnParameter() const {
    std::string params = this->getOption(learnparameter).getArgumentByName("name").getValueAsString();
    if (params.empty()) return {};
    else {
        std::stringstream ss(params);
        std::string number;
        std::vector<int> ret;
        while(std::getline(ss, number, ',')) {
            ret.push_back(std::stoi(number));
        }
        return ret;
    }
}

std::string LearningModelCheckerSettings::getEvaluationMethod() const {
    return this->getOption(evluationmethod).getArgumentByName("name").getValueAsString();
}

std::vector<std::string> LearningModelCheckerSettings::getChildrenFiles() const {
    std::string templateFile =  this->getOption(childrenfiles).getArgumentByName("name").getValueAsString();
    if (templateFile.empty()) return {};
    else {
        std::vector<std::string> ret;
        for(int param : this->getLearnParameter()) {
            ret.push_back(boost::replace_all_copy(templateFile, "@", std::to_string(param)));
        }
        return ret;
    }
}

std::string LearningModelCheckerSettings::getDTLearnAlgorithm() const {
    return this->getOption(treelearner).getArgumentByName("name").getValueAsString();
}

std::string LearningModelCheckerSettings::getParameterName() const {
    return this->getOption(learnparametername).getArgumentByName("name").getValueAsString();
}

uint64_t LearningModelCheckerSettings::getNumberOfSMCRuns() const {
    return this->getOption(smcruns).getArgumentByName("value").getValueAsUnsignedInteger();
}

uint64_t LearningModelCheckerSettings::getMaxStepsForSMC() const {
    return this->getOption(maxsteps).getArgumentByName("value").getValueAsUnsignedInteger();
}

bool LearningModelCheckerSettings::isOnlyReachableSet() const {
    return this->getOption(onlyReachable).getHasOptionBeenSet();
}

uint64_t LearningModelCheckerSettings::getNumDTStrat() const {
    return this->getOption(numStepsDTStrat).getArgumentByName("value").getValueAsUnsignedInteger();
}

std::vector<std::pair<std::string, std::string>> LearningModelCheckerSettings::getSubTasksAndProperties() const {
    std::string fileList =  this->getOption(subProgramsAndTasks).getArgumentByName("name").getValueAsString();
    if (fileList.empty()) return {};
    else {
        std::vector<std::pair<std::string, std::string>> fileAndTaskList;
        std::stringstream test(fileList);
        std::string segment;

        while(std::getline(test, segment, '%'))
        {
            std::stringstream subtest(segment);
            std::string first;
            std::string second;
            std::getline(subtest, first, '#');
            std::getline(subtest, second, '#');
            std::pair<std::string, std::string> fileAndProp(first, second);
            fileAndTaskList.push_back(fileAndProp);
        }
        return fileAndTaskList;
    }
}

double LearningModelCheckerSettings::getPACDelta() const {
    return this->getOption(pacDelta).getArgumentByName("value").getValueAsDouble();
}

double LearningModelCheckerSettings::getPACEpsilon() const {
    return this->getOption(pacEpsilon).getArgumentByName("value").getValueAsDouble();
}

uint64_t LearningModelCheckerSettings::getHyperparamLearnTime() const {
    return this->getOption(hyperparamTime).getArgumentByName("value").getValueAsUnsignedInteger();
}

bool LearningModelCheckerSettings::useSemanticActions() const {
    return this->getOption(semanticActions).getHasOptionBeenSet();
}

}  // namespace modules
}  // namespace settings
}  // namespace storm
