//
// Created by steffi on 05.01.22.
//
#include "storm/settings/modules/HeuristicLearningSettings.h"

#include "storm/settings/Argument.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/SettingsManager.h"

namespace storm {
namespace settings {
namespace modules {

const std::string HeuristicLearningSettings::moduleName = "heuristiclearning";
const std::string HeuristicLearningSettings::explorationEpsilon = "explEpsilon";
const std::string HeuristicLearningSettings::learningConfig = "learnConfig";
const std::string HeuristicLearningSettings::updateSteps = "updateSteps";
const std::string HeuristicLearningSettings::heuristicType = "heuristic";
const std::string HeuristicLearningSettings::onlyReachable = "onlyReachable";
const std::string HeuristicLearningSettings::moduleHeuristc = "moduleHeuristic";
const std::string HeuristicLearningSettings::externalFile = "filename";


HeuristicLearningSettings::HeuristicLearningSettings() : ModuleSettings(moduleName) {
    this->addOption(
        storm::settings::OptionBuilder(moduleName, explorationEpsilon, false, "If set, determines the probability of following the heuristic.")
            .addArgument(
                storm::settings::ArgumentBuilder::createDoubleArgument("value", "The value of the exploration epsilon.")
                    .setDefaultValueDouble(0).build())
            .setIsAdvanced()
            .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, learningConfig, false, "If set, describes the learning configuration of the heuristic.")
            .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the configuration, in [non, basic, genetic, external, headstart]")
                    .setDefaultValueString("non").build())
            .setIsAdvanced()
            .build());
    this->addOption(
        storm::settings::OptionBuilder(moduleName, updateSteps, false, "If set, determines the number of simulations before heuristic is updated.")
        .addArgument(
                storm::settings::ArgumentBuilder::createIntegerArgument("value", "Number of simulations before update.")
                .setDefaultValueInteger(-1).build())
            .setIsAdvanced()
            .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, heuristicType, false, "If set, describes the heuristic type.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the heuristic, in [non, dtcontrol, module, mlpack, mlpackmodule]")
                                         .setDefaultValueString("non").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, onlyReachable, false, "If set, use only reachable states under current strategy for export.")
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, moduleHeuristc, false, "If set, describes the heuristic type of the module decision.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the heuristic, in [orderingrandom, orderingtree, decisiontree]")
                                         .setDefaultValueString("orderingrandom").build())
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, externalFile, false, "If set, gives the external file from which to read the strategy from.")
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The path to the file.")
                                         .setDefaultValueString("").build())
                        .setIsAdvanced()
                        .build());
}

storm::modelchecker::dtstrat::configTypes HeuristicLearningSettings::getLearningConfigType() const {
    std::string configName = this->getOption(learningConfig).getArgumentByName("name").getValueAsString();
    if (configName=="basic") {
        return storm::modelchecker::dtstrat::configTypes::basic;
    } else if (configName=="genetic") {
        return storm::modelchecker::dtstrat::configTypes::genetic;
    } else if (configName=="external") {
        return storm::modelchecker::dtstrat::configTypes::external;
    } else if (configName=="headstart") {
        return storm::modelchecker::dtstrat::configTypes::headstart;
    }
    else {
        return storm::modelchecker::dtstrat::configTypes::noneConfig;
    }
}

storm::modelchecker::dtstrat::heuristicTypes HeuristicLearningSettings::getHeuristicType() const {
    std::string heuristicName = this->getOption(heuristicType).getArgumentByName("name").getValueAsString();
    if (heuristicName=="dtcontrol") {
        return storm::modelchecker::dtstrat::heuristicTypes::dtcontrolDT;
    } else if (heuristicName=="module") {
        return storm::modelchecker::dtstrat::heuristicTypes::moduleDT;
    } else if (heuristicName=="mlpack") {
        return storm::modelchecker::dtstrat::heuristicTypes::mlpackDT;
    } else if (heuristicName=="mlpackmodule") {
        return storm::modelchecker::dtstrat::heuristicTypes::mlpackModule;
    } else {
        return storm::modelchecker::dtstrat::heuristicTypes::noneHeuristic;
    }

}

int HeuristicLearningSettings::getNumberOfSimulationsBeforeUpdate() const {
    return this->getOption(updateSteps).getArgumentByName("value").getValueAsInteger();
}

double HeuristicLearningSettings::getHeuristicProbability() const {
    return this->getOption(explorationEpsilon).getArgumentByName("value").getValueAsDouble();
}

bool HeuristicLearningSettings::getOnlyReachable() const {
    return this->getOption(onlyReachable).getHasOptionBeenSet();
}

storm::modelchecker::dtstrat::moduleHeuristicTypes HeuristicLearningSettings::getModuleHeuristicType() const {
    std::string moduleHeuristicName = this->getOption(moduleHeuristc).getArgumentByName("name").getValueAsString();
    if (moduleHeuristicName == "decisiontree" || moduleHeuristicName == "dt") {
        return modelchecker::dtstrat::decisiontree;
    } else if (moduleHeuristicName == "orderingrandom") {
        return modelchecker::dtstrat::orderingrandom;
    } else if (moduleHeuristicName == "orderingtree") {
        return modelchecker::dtstrat::orderingtree;
    } else return modelchecker::dtstrat::ordering;
}

std::string HeuristicLearningSettings::getExternalFilename() const {
    return this->getOption(externalFile).getArgumentByName("name").getValueAsString();
}

}  // namespace modules
}  // namespace settings
}  // namespace storm
