//
// Created by steffi on 27.07.22.
//

#include "storm/modelchecker/learning/LearningStrategiesModelChecker.h"
#include "storm/modelchecker/learning/CandidateUntilDT.h"

#include "storm/modelchecker/explorationDT/BSCCCandidate.h"
#include "storm/modelchecker/explorationDT/BoundsDT.h"
#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"
#include "storm/modelchecker/explorationDT/StateGenerationDT.h"
#include "storm/modelchecker/explorationDT/StatisticsDT.h"

#include "storm/generator/CompressedState.h"

#include "storm/storage/prism/Program.h"

#include "storm/logic/FragmentSpecification.h"

#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"

#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"

#include "storm/exceptions/InvalidPropertyException.h"
#include "storm/modelchecker/explorationDT/learning/BasicLearningConfig.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/GeneralSettings.h"
#include "storm/settings/modules/HeuristicLearningSettings.h"
#include "storm/settings/modules/IOSettings.h"
#include "storm/settings/modules/LearningModelCheckerSettings.h"
#include "storm/storage/DecisionTree.h"
#include "storm/utility/graph.h"

#include <boost/algorithm/string/predicate.hpp>

#include "storm/api/builder.h"
#include "storm/api/properties.h"
#include "storm/environment/Environment.h"
#include "storm/modelchecker/explorationDT/heuristics/MLDecisionTreeHeuristic.h"
#include "storm/modelchecker/prctl/SparseMdpPrctlModelChecker.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/storage/SymbolicModelDescription.h"
#include "storm/utility/Engine.h"
#include "storm/utility/cli.h"

#include "storm/modelchecker/explorationDT/BRTDPInformation.h"

#include "storm-parsers/api/model_descriptions.h"
#include "storm-parsers/api/properties.h"

#include <chrono>

namespace storm {
namespace modelchecker {

using namespace exploration_detail_dt;

struct SymbolicInput {
    // The symbolic model description.
    boost::optional<storm::storage::SymbolicModelDescription> model;

    // The original properties to check.
    std::vector<storm::jani::Property> properties;

    // The preprocessed properties to check (in case they needed amendment).
    boost::optional<std::vector<storm::jani::Property>> preprocessedProperties;
};

struct ModelProcessingInformation {
    // The engine to use
    storm::utility::Engine engine;

    // If set, bisimulation will be applied.
    bool applyBisimulation;

    // If set, a transformation to Jani will be enforced
    bool transformToJani;

    // Which data type is to be used for numbers ...
    enum class ValueType { FinitePrecision, Exact, Parametric };
    ValueType buildValueType;         // ... during model building
    ValueType verificationValueType;  // ... during model verification

    // The Dd library to be used
    storm::dd::DdType ddType;

    // The environment used during model checking
    storm::Environment env;

    // A flag which is set to true, if the settings were detected to be compatible.
    // If this is false, it could be that the query can not be handled.
    bool isCompatible;
};

template<typename ModelType, typename StateType>
LearningStrategiesModelChecker<ModelType, StateType>::LearningStrategiesModelChecker(
    storm::prism::Program const& program, std::vector<storm::jani::Property>& properties,
    boost::optional<std::vector<storm::jani::Property>>& preprocessedProperties)
    : substprogram(program.substituteConstantsFormulas()),
      randomGenerator(std::chrono::system_clock::now().time_since_epoch().count()),
      comparator(storm::settings::getModule<storm::settings::modules::ExplorationSettings>().getPrecision()),
      properties(properties),
      preprocessedProperties(preprocessedProperties) {
    const auto& ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    std::stringstream test(ioSettings.getPrismInputFilename());
    std::string segment;
    std::vector<std::string> seglist;
    this->program = storm::api::parseProgram(ioSettings.getPrismInputFilename(), false, false).substituteConstantsFormulas();
    while (std::getline(test, segment, '/')) {
        seglist.push_back(segment);
    }
    this->modelName = seglist.back().substr(0, seglist.back().size() - 6);
    this->propertyName = properties[0].getName();
    const auto& learnSettings = storm::settings::getModule<storm::settings::modules::LearningModelCheckerSettings>();
    if (learnSettings.useSemanticActions())
        getAllActionsFromModel(program);
    STORM_PRINT("ACTIONS:\n");
    STORM_PRINT("indexToShortActionName\n");
    for(uint i=0; i<indexToShortActionName.size(); i++) STORM_PRINT(i<<" - "<<indexToShortActionName[i]<<'\n');
    STORM_PRINT("\ncompleteActionNameToIndexMapping\n");
    for(auto el : completeActionNameToIndexMapping) STORM_PRINT(el.first << " - " << el.second <<'\n');
    STORM_PRINT("\nshortActionNameToIndexMapping\n");
    for(auto el : shortActionNameToIndexMapping) STORM_PRINT(el.first << " - "<< el.second << '\n');
    STORM_PRINT("\nactionIndexToLabelMapping\n");
    for(uint i=0; i<actionIndexToLabelMapping.size(); i++) {
        STORM_PRINT(i<<" - [");
        for (auto el: actionIndexToLabelMapping[i]) STORM_PRINT(el<<" , ");
        STORM_PRINT("]\n");
    }
}

template<typename ModelType, typename StateType>
bool LearningStrategiesModelChecker<ModelType, StateType>::canHandle(CheckTask<storm::logic::Formula, ValueType> const& checkTask) const {
    storm::logic::Formula const& formula = checkTask.getFormula();
    storm::logic::FragmentSpecification fragment = storm::logic::reachability();
    return formula.isInFragment(fragment) && checkTask.isOnlyInitialStatesRelevantSet();
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getAllActionsFromModel(storm::prism::Program const& program) {
    // create a mapping of the actions in the program
    // each action (i.e. 'step') will have its unique index
    // actions of different modules, but the same action, shoudl have the same index
    // e.g. phil1_eat and phil2_eat should correspond to the same index

    std::vector<storm::prism::Module> allModules = program.getModules();
    for (storm::prism::Module module : allModules) {
        std::vector<storm::prism::Command> commands = module.getCommands();
        for (storm::prism::Command command : commands) {
            STORM_LOG_ASSERT(command.isLabeled(),"The actions must ALL have labels, otherwise it doesn't work.");
            std::string nameOfCommand = command.getActionName();
            // each of our actions is named as follows:
            // STEFFI_{modulename}_STEFFI_{linenumber}
            std::string trueNameOfCommand;
            if (nameOfCommand.rfind("STEFFI_", 0)==0) {
                // This action was named by us!
                size_t pos = nameOfCommand.rfind("STEFFI_");
                trueNameOfCommand = nameOfCommand.substr(pos+7);
                if (shortActionNameToIndexMapping.count(trueNameOfCommand)>0) {
                    uint actionIndex = shortActionNameToIndexMapping[trueNameOfCommand];
                    completeActionNameToIndexMapping[nameOfCommand] = actionIndex;
                    if (actionIndexToLabelMapping.size()<=actionIndex) {
                        std::vector<std::string> temp = {nameOfCommand};
                        actionIndexToLabelMapping.resize(actionIndex);
                        actionIndexToLabelMapping.push_back(temp);
                    } else {
                        actionIndexToLabelMapping[actionIndex].push_back(nameOfCommand);
                    }
                } else {
                    uint actionIndex = indexToShortActionName.size();
                    shortActionNameToIndexMapping[trueNameOfCommand] = actionIndex;
                    indexToShortActionName.push_back(trueNameOfCommand);
                    completeActionNameToIndexMapping[nameOfCommand] = actionIndex;
                    if (actionIndexToLabelMapping.size()<=actionIndex) {
                        std::vector<std::string> temp = {nameOfCommand};
                        actionIndexToLabelMapping.resize(actionIndex);
                        actionIndexToLabelMapping.push_back(temp);
                    } else {
                        actionIndexToLabelMapping[actionIndex].push_back(nameOfCommand);
                    }
                }
            }
            else {
                if (shortActionNameToIndexMapping.count(nameOfCommand)==0) {
                    uint actionIndex = indexToShortActionName.size();
                    shortActionNameToIndexMapping[nameOfCommand] = actionIndex;
                    indexToShortActionName.push_back(nameOfCommand);
                    completeActionNameToIndexMapping[nameOfCommand] = actionIndex;
                    if (actionIndexToLabelMapping.size()<=actionIndex) {
                        std::vector<std::string> temp = {nameOfCommand};
                        actionIndexToLabelMapping.resize(actionIndex);
                        actionIndexToLabelMapping.push_back(temp);
                    } else {
                        actionIndexToLabelMapping[actionIndex].push_back(nameOfCommand);
                    }
                }
            }
        }
    }
    uint actionIndex = indexToShortActionName.size();
    shortActionNameToIndexMapping["NEW-SELF-LOOP"] = actionIndex;
    indexToShortActionName.push_back("NEW-SELF-LOOP");
    completeActionNameToIndexMapping["NEW-SELF-LOOP"] = actionIndex;
    if (actionIndexToLabelMapping.size()<=actionIndex) {
        std::vector<std::string> temp = {"NEW-SELF-LOOP"};
        actionIndexToLabelMapping.resize(actionIndex);
        actionIndexToLabelMapping.push_back(temp);
    } else {
        actionIndexToLabelMapping[actionIndex].push_back("NEW-SELF-LOOP");
    }
}

template<typename ModelType, typename StateType>
std::unique_ptr<CheckResult> LearningStrategiesModelChecker<ModelType, StateType>::computeUntilProbabilities(
    Environment const& env, CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) {
    // GET THE DECISIONTREE
    auto start = std::chrono::high_resolution_clock::now();
    const auto& learnSettings = storm::settings::getModule<storm::settings::modules::LearningModelCheckerSettings>();
    if (!learnSettings.getDecisionTreeFile().empty()) {
        // The user has provided a decision-tree in a file (json or dot)
        // read the file and store it in a DT (methods exist already)
        std::string decisionTreeFilename = learnSettings.getDecisionTreeFile();
        getDTFromFile(decisionTreeFilename);
    } else {
        if (!learnSettings.getDataFile().empty()) {
            // The user has provided a file with the learning data
            // It should be a csv-file, with state-valuations and the chosen action as last entry per row
            std::string filename = learnSettings.getDataFile();
            auto method = learnSettings.getDTLearnAlgorithm();
            getDTFromDataFile(filename, method);
        } else {
            // We have to generate the data ourselves
            // check if we received parameters for the smaller models
            if (!learnSettings.getSubTasksAndProperties().empty() ||
                (!learnSettings.getLearnParameter().empty() && !learnSettings.getParameterName().empty()) ||
                (!learnSettings.getLearnHyperParamType().empty() && !learnSettings.getParameterName().empty())) {
                learnDTFromChildren(env, checkTask);
            } else {
                // Use DTStrat to get something
                useDTStrat(checkTask, learnSettings.getNumDTStrat(), learnSettings.useSemanticActions());
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    STORM_PRINT("Total time to get the DT: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000 << "s\n");
    if (decisionTree)
        STORM_PRINT("DT-size: " << decisionTree->size() << "\n");

    // EVALUATE THE STRATEGY
    ValueType returnValue;
    start = std::chrono::high_resolution_clock::now();
    if (learnSettings.getEvaluationMethod() == "full") {
        // do a similar thing as in prism (value iteration while iteratively building the MC)
        // change the returnValue here (=value of the strategy)
        program = program.substituteConstantsFormulas();
        substprogram = substprogram.substituteConstantsFormulas();
        evaluateVIMC(checkTask, learnSettings.useSemanticActions());
    } else if (learnSettings.getEvaluationMethod() == "statistical") {
        // do some statistical modelchecking
        // change the returnValue here (=value of the strategy)
        // changing the number of simulation to Hoeffding's inequality
        // numSMC = learnSettings.getNumberOfSMCRuns();
        storm::settings::modules::GeneralSettings gs;
        double epsilon = learnSettings.getPACEpsilon();
        double delta = learnSettings.getPACDelta();
        numSMC = (uint)1.0 / (2 * epsilon * epsilon) * log(2 / delta);
        STORM_PRINT("numSMC: " << numSMC << "\n");
        stopAfter = learnSettings.getMaxStepsForSMC();
        program = program.substituteConstantsFormulas();
        substprogram = substprogram.substituteConstantsFormulas();
        returnValue = evaluateSMC(checkTask, learnSettings.useSemanticActions());
    } else {
        auto const& ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
        if (ioSettings.isExportSchedulerSet()) {
            auto filename = ioSettings.getExportSchedulerFilename();
            this->decisionTree->exportDT(filename);
        }
        return std::make_unique<ExplicitQuantitativeCheckResult<ValueType>>(0, -1);
    }
    end = std::chrono::high_resolution_clock::now();
    STORM_PRINT("Time to evaluate the DT: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000 << "s\n");
    auto const& ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    if (ioSettings.isExportSchedulerSet()) {
        auto filename = ioSettings.getExportSchedulerFilename();
        this->decisionTree->exportDT(filename);
    }
    return std::make_unique<ExplicitQuantitativeCheckResult<ValueType>>(0, returnValue);
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::learnDTFromChildren(Environment const& env,
                                                                               CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) {
    const auto& learnSettings = storm::settings::getModule<storm::settings::modules::LearningModelCheckerSettings>();
    const auto& ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    bool storeToFile = (learnSettings.getDTLearnAlgorithm() != "mlpack");
    arma::Row<long unsigned int> labels;
    arma::Mat<int> data;
    bool first = true;
    std::string oldFileName = "";
    std::string learnParameterName = "";
    if (learnSettings.getSubTasksAndProperties().empty()) {
        learnParameterName = "%" + learnSettings.getParameterName();
    } else {
        auto childrenFilesAndProperties = learnSettings.getSubTasksAndProperties();
        for (auto fileAndProperty : childrenFilesAndProperties) {
            std::string file = fileAndProperty.first;
            std::stringstream test(file);
            std::string segment;
            std::vector<std::string> seglist;
            while (std::getline(test, segment, '.')) {
                seglist.push_back(segment);
            }
            learnParameterName += '%' + seglist[1];
        }
    }
    if (learnSettings.getSubTasksAndProperties().empty()) {
        auto datatimestart = std::chrono::high_resolution_clock::now();
        double totalTime = 0.0;
        // read the same file that is hidden somewhere in 'this->program' again with another constant given in val
        bool done = false;
        int val;
        if (learnSettings.getLearnHyperParamType() != "") {
            val = 1;
        } else {
            auto allVals = learnSettings.getLearnParameter();
            if (allVals.empty())
                done = true;
            else
                val = allVals[0];
        }
        uint count = 0;
        while (!done) {
            // for (auto val : learnSettings.getLearnParameter()) {
            auto innerstart = std::chrono::high_resolution_clock::now();
            storm::prism::Program changedProgram(this->program);

            // get all constants from the prism program
            std::vector<prism::Constant>& programConstants = changedProgram.getConstants();
            uint index = -1;

            // find the constant that we need to change (= the parameter)
            for (uint i = 0; i < programConstants.size(); i++) {
                auto& constant = programConstants[i];
                //STORM_PRINT(bla.getExpressionVariable().getExpression().toString() << " = " << bla.getExpression().toString() << "\n");
                // std::map<storm::expressions::Variable, storm::expressions::Expression> substitution;
                if (constant.getExpressionVariable().getExpression().toString() == learnSettings.getParameterName()) {
                    index = i;
                    break;
                }
            }
            auto& constant = programConstants[index];

            // generate a new expression that describes the value of this constant (far too complex, but storm was never meant to allow that..)
            // auto& expr = constant.getExpression();
            storm::expressions::IntegerType typ;
            auto valExpression = new storm::expressions::IntegerLiteralExpression(program.getManager(), val);
            std::shared_ptr<storm::expressions::BaseExpression> pt(valExpression);
            storm::expressions::Expression newExpression(pt);
            // constant = prism::Constant(constant.getExpressionVariable(), newExpression, constant.getFilename(), constant.getLineNumber());
            std::string constantDefinitionString = ioSettings.getConstantDefinitionString();
            std::map<storm::expressions::Variable, storm::expressions::Expression> constantDefinitions;
            constantDefinitions = storm::utility::cli::parseConstantDefinitionString(program.getManager(), constantDefinitionString);
            constantDefinitions[constant.getExpressionVariable()] = newExpression;
            changedProgram = changedProgram.defineUndefinedConstants(constantDefinitions).substituteNonStandardPredicates();
            std::string dataFileName;
            if (oldFileName.empty())
                dataFileName = this->modelName + "%" + this->propertyName + learnParameterName + "%" + std::to_string(val) + ".csv";
            else
                dataFileName = oldFileName.substr(0, oldFileName.size() - 4) + "_" + std::to_string(val) + ".csv";
            auto timeForChecking = checkSmallInstance(changedProgram, storeToFile, env, first, checkTask, learnSettings.isOnlyReachableSet(), labels, data,
                                                      this->properties, this->preprocessedProperties, dataFileName, oldFileName, learnSettings.useSemanticActions());
            oldFileName = dataFileName;
            totalTime += timeForChecking;
            auto innerend = std::chrono::high_resolution_clock::now();
            STORM_PRINT("Time to get the data for " << learnSettings.getParameterName() << "=" << val << " : "
                                                    << double(std::chrono::duration_cast<std::chrono::milliseconds>(innerend - innerstart).count()) / 1000
                                                    << "s\n");
            if (totalTime >= learnSettings.getHyperparamLearnTime()) {
                STORM_PRINT("Stop running because more than " << learnSettings.getHyperparamLearnTime() << "s.\n");
                break;
            } else if (totalTime + timeForChecking >= learnSettings.getHyperparamLearnTime()) {
                STORM_PRINT("Stop running because in the next step more than " << learnSettings.getHyperparamLearnTime() << "s.\n");
                break;
            };
            count++;
            if (learnSettings.getLearnHyperParamType() == "linear") {
                val += 1;
            } else if (learnSettings.getLearnHyperParamType() == "exponential") {
                val *= 2;
            } else {
                auto allVals = learnSettings.getLearnParameter();
                if (allVals.empty())
                    done = true;
                else if (count == allVals.size())
                    done = true;
                else
                    val = allVals[count];
            }
            if (count >= 30) {
                STORM_PRINT("More than 30 models verified -> stop here!\n");
                break;
            }
        }
        auto datatimeend = std::chrono::high_resolution_clock::now();
        STORM_PRINT("Total time to get the data: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(datatimeend - datatimestart).count()) / 1000
                                                   << "s\n");
    } else {
        double totalTime = 0.0;
        uint count = 0;
        auto datatimestart = std::chrono::high_resolution_clock::now();
        auto childrenFilesAndProperties = learnSettings.getSubTasksAndProperties();
        for (auto fileAndProperty : childrenFilesAndProperties) {
            std::string dataFileName;
            if (oldFileName.empty())
                dataFileName = this->modelName + "%" + this->propertyName + learnParameterName + "%" + std::to_string(count) + ".csv";
            else
                dataFileName = oldFileName.substr(0, oldFileName.size() - 4) + "_" + std::to_string(count) + ".csv";
            auto innerstart = std::chrono::high_resolution_clock::now();
            std::string file = fileAndProperty.first;
            std::string property = fileAndProperty.second;
            storm::prism::Program childProgram = storm::api::parseProgram(file, false, false);
            std::string constantDefinitionString = ioSettings.getConstantDefinitionString();
            childProgram =
                childProgram.defineUndefinedConstants(storm::utility::cli::parseConstantDefinitionString(childProgram.getManager(), constantDefinitionString))
                    .substituteNonStandardPredicates();
            std::vector<storm::jani::Property> newProperties;
            boost::optional<std::set<std::string>> propertyFilter;
            newProperties = storm::api::parsePropertiesForPrismProgram(property, childProgram, propertyFilter);
            newProperties = storm::api::substituteConstantsInProperties(newProperties, program.getConstantsFormulasSubstitution());
            auto formula = newProperties[0].getRawFormula();
            auto localCheckTask = storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>(*formula, true);
            // auto test = storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>()
            auto timeForChecking = checkSmallInstance(childProgram, storeToFile, env, first, checkTask, learnSettings.isOnlyReachableSet(), labels, data,
                                                      newProperties, this->preprocessedProperties, dataFileName, oldFileName, learnSettings.useSemanticActions());
            oldFileName = dataFileName;
            totalTime += timeForChecking;
            auto innerend = std::chrono::high_resolution_clock::now();
            STORM_PRINT("Time to get the data for " << file << " : "
                                                    << double(std::chrono::duration_cast<std::chrono::milliseconds>(innerend - innerstart).count()) / 1000
                                                    << "s\n");
            if (totalTime >= learnSettings.getHyperparamLearnTime()) {
                STORM_PRINT("Stop running because more than " << learnSettings.getHyperparamLearnTime() << "s.\n");
                break;
            } else if (totalTime + timeForChecking >= learnSettings.getHyperparamLearnTime()) {
                STORM_PRINT("Stop running because in the next step more than " << learnSettings.getHyperparamLearnTime() << "s.\n");
                break;
            };
            count++;
        }
        auto datatimeend = std::chrono::high_resolution_clock::now();
        STORM_PRINT("Total time to get the data: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(datatimeend - datatimestart).count()) / 1000
                                                   << "s\n");
    }
    auto method = learnSettings.getDTLearnAlgorithm();
    auto start2 = std::chrono::high_resolution_clock::now();
    // after collectin all the information, learn a decision-tree
    if (method == "none") {
        // do nothing
    } else if (method == "dtcontrol") {
        // Use the existing method to generate and load a dtcontrol-tree
        getDTControlTree(oldFileName);
    } else if (method == "sklearn") {
        // write a python-script to generate a sklearn-tree and store it as a dot (?)
        getSklearnTree(oldFileName);
    } else {
        // load the data in an armadillo-matrix to forward it to mlpack
        arma::Row<long unsigned int> possibleActions = arma::unique(labels);
        uint numberOfActions = possibleActions.size();
        getMLPackTree(data, labels, numberOfActions);
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    STORM_PRINT("Time to learn the DT: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count()) / 1000 << "s\n");
}

template<typename ModelType, typename StateType>
uint LearningStrategiesModelChecker<ModelType, StateType>::checkSmallInstance(storm::prism::Program& changedProgram, bool const& storeToFile,
                                                                              storm::Environment const& env, bool& first,
                                                                              CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask,
                                                                              bool const& onlyReachable, arma::Row<long unsigned int>& labels,
                                                                              arma::Mat<int>& data, std::vector<storm::jani::Property> const& properties,
                                                                              boost::optional<std::vector<storm::jani::Property>> const& preprocessedProperties,
                                                                              std::string const& dataFileName, std::string const& oldFileName,
                                                                              bool const& semanticActions) {
    auto start = std::chrono::high_resolution_clock::now();

    // get a model from the prism-program
    std::shared_ptr<storm::models::ModelBase> test = getModelFromProgram(changedProgram, properties, preprocessedProperties);
    std::shared_ptr<storm::models::sparse::Model<ValueType>> sparseModel = test->as<storm::models::sparse::Model<ValueType>>();
    auto mdp = sparseModel->template as<storm::models::sparse::Mdp<ValueType>>();
    auto cl = mdp->getChoiceLabeling();
    //STORM_PRINT("All actions have labels "<<cl.isComplete()<<"\n");
    /*uint count=0;
    for (auto row : mdp->getTransitionMatrix()) {
        std::ofstream dst("TEMP.txt", std::ios::binary);

        auto labels = cl.getLabelsOfChoice(count);
        std::string label;
        if (!labels.empty()) label = *(labels.begin());
        dst << count << " , " << label << "\n";
        count++;
    }*/
    // create a model-checking-instance
    auto const& properties_local = preprocessedProperties ? preprocessedProperties.get() : properties;
    auto property = properties_local[0];
    std::unique_ptr<storm::modelchecker::CheckResult> result;
    auto rawFormula = property.getRawFormula();
    bool filterForInitialStates = property.getFilter().getStatesFormula()->isInitialFormula();
    auto task = storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>(*rawFormula, filterForInitialStates);
    task.setProduceSchedulers(true);
    storm::modelchecker::SparseMdpPrctlModelChecker<storm::models::sparse::Mdp<ValueType>> modelchecker(*mdp);

    // check the property on the small model
    if (modelchecker.canHandle(task)) {
        result = modelchecker.check(env, task);
        if (first && storeToFile) {
            first = false;
            auto size = getNumberOfStateVariables(checkTask, changedProgram);
            printHeader(dataFileName, size);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    storm::modelchecker::ExplicitQualitativeCheckResult filter(sparseModel->getInitialStates());
    result->filter(filter);
    storm::storage::Scheduler<ValueType> const& scheduler = result->template asExplicitQuantitativeCheckResult<ValueType>().getScheduler();
    storm::storage::BitVector subsystem;

    // export the strategy in a file or a armadillo-matrix
    if (onlyReachable) {
        storm::storage::BitVector actionSelection =
            result->asExplicitQuantitativeCheckResult<ValueType>().getScheduler().computeActionSupport(mdp->getNondeterministicChoiceIndices());
        storm::storage::BitVector allStates(mdp->getNumberOfStates(), true);
        subsystem = storm::utility::graph::getReachableStates(mdp->getTransitionMatrix(), mdp->getInitialStates(), allStates,
                                                              storm::storage::BitVector(allStates.size(), false), false, 0, actionSelection);
    } else {
        subsystem = storm::storage::BitVector(mdp->getNumberOfStates(), true);
    }
    if (storeToFile) {
        if (!oldFileName.empty()) {
            std::ifstream src(oldFileName, std::ios::binary);
            std::ofstream dst(dataFileName, std::ios::binary);

            dst << src.rdbuf();
        }
        printToFile(dataFileName, scheduler, *mdp, subsystem, semanticActions);
    } else
        getValuationsAndLabels(scheduler, *mdp, data, labels, subsystem, semanticActions);
    return double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getDTFromFile(std::string const& decisionTreeFilename) {
    if (boost::algorithm::ends_with(decisionTreeFilename, ".dot")) {
        decisionTree = new dtstrat::DTInterface<StateType, ValueType>(decisionTreeFilename);
    } else if (boost::algorithm::ends_with(decisionTreeFilename, ".json")) {
        nlohmann2::json jsonTree = nlohmann2::json::parse(decisionTreeFilename);
        decisionTree = new dtstrat::DTInterface<StateType, ValueType>(jsonTree);
    }
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getDTFromDataFile(std::string const& datafilename, std::string const& method) {
    if (method == "dtcontrol") {
        // Use the existing method to generate and load a dtcontrol-tree
        getDTControlTree(datafilename);
    } else if (method == "sklearn") {
        // write a python-script to generate a sklearn-tree and store it as a dot (?)
        getSklearnTree(datafilename);
    } else {
        // load the data in an armadillo-matrix to forward it to mlpack
        getMLPackTree(datafilename);
    }
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::useDTStrat(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, uint const& numRuns, bool const& semanticActions) {
    STORM_LOG_WARN("THIS FUNCTION MAY BE DEPRECATED!");
    storm::logic::UntilFormula const& untilFormula = checkTask.getFormula();
    storm::logic::Formula const& conditionFormula = untilFormula.getLeftSubformula();
    storm::logic::Formula const& targetFormula = untilFormula.getRightSubformula();
    STORM_LOG_THROW(program.isDeterministicModel() || checkTask.isOptimizationDirectionSet(), storm::exceptions::InvalidPropertyException,
                    "For nondeterministic systems, an optimization direction (min/max) must be given in the property.");

    ExplorationInformationDT<StateType, ValueType> explorationInformation(
        checkTask.isOptimizationDirectionSet() ? checkTask.getOptimizationDirection() : storm::OptimizationDirection::Maximize, program);

    // The first row group starts at action 0.
    explorationInformation.newRowGroup(0);

    std::map<std::string, storm::expressions::Expression> labelToExpressionMapping = program.getLabelToExpressionMapping();
    generator::NextStateGeneratorOptions options;
    options.setBuildChoiceLabels();
    options.setBuildChoiceOrigins();
    options.setBuildStateValuations();
    StateGenerationDT<StateType, ValueType> stateGeneration(program, explorationInformation,
                                                            conditionFormula.toExpression(program.getManager(), labelToExpressionMapping),
                                                            targetFormula.toExpression(program.getManager(), labelToExpressionMapping), options);
    // Create a structure that holds the bounds for the states and actions.
    BoundsDT<StateType, ValueType> bounds;

    // Now perform the actual sampling.
    StatisticsDT<StateType, ValueType> stats;

    // Create the state valuations storage
    StateValuationsDT stateValuations;

    // Create the BRTDPInformation class that contains pointer to the important information of BRTDP
    BRTDPInformation<StateType, ValueType> brtdpInfo(explorationInformation, stateGeneration, bounds, stats, stateValuations,
                                                     storm::settings::getModule<storm::settings::modules::ExplorationSettings>().getPrecision());

    runBRTDPforNTimes(numRuns, brtdpInfo, semanticActions);
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::runBRTDPforNTimes(uint const& N, BRTDPInformation<StateType, ValueType>& brtdpInfo, bool const& semanticActions) {
    auto start = std::chrono::high_resolution_clock::now();

    StateGenerationDT<StateType, ValueType>& stateGeneration = brtdpInfo.getStateGeneration();
    ExplorationInformationDT<StateType, ValueType>& explorationInformation = brtdpInfo.getExplorationInformation();
    BoundsDT<StateType, ValueType>& bounds = brtdpInfo.getBounds();
    StatisticsDT<StateType, ValueType>& stats = brtdpInfo.getStats();

    StateValuationsDT& stateValuations = brtdpInfo.getStateValuations();

    // Generate the initial state so we know where to start the simulation.
    stateGeneration.computeInitialStates();
    STORM_LOG_THROW(stateGeneration.getNumberOfInitialStates() == 1, storm::exceptions::NotSupportedException,
                    "Currently only models with one initial state are supported by the exploration engine.");
    StateType initialStateIndex = stateGeneration.getFirstInitialState();

    // Create a stack that is used to track the path we sampled.
    std::vector<std::pair<StateType, ActionType>> stack;
    const auto& configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();

    double explorationEpsilon = configSettings.getHeuristicProbability();
    int updatesteps = configSettings.getNumberOfSimulationsBeforeUpdate();
    auto heuristicType = storm::modelchecker::dtstrat::heuristicTypes::mlpackDT;
    LearningConfig<StateType>* config = new BasicLearningConfig<StateType, ValueType>(explorationEpsilon, updatesteps, heuristicType, brtdpInfo);

    config->firstLog();
    config->prepareBoundLog();
    for (uint i = 0; i < N; i++) {
        bool result = brtdpInfo.samplePathFromInitialState(stateGeneration, explorationInformation, stack, bounds, stats, config, stateValuations, semanticActions);

        stats.sampledPath();
        stats.updateMaxPathLength(stack.size());

        // If a terminal state was found, we update the probabilities along the path contained in the stack.
        if (result) {
            // Update the bounds along the path to the terminal state.
            STORM_LOG_TRACE("Found terminal state, updating probabilities along path.");
            auto val = bounds.getLowerBoundForState(initialStateIndex, explorationInformation);
            brtdpInfo.updateProbabilityBoundsAlongSampledPath(stack, explorationInformation, bounds);
            STORM_LOG_TRACE(bounds.getLowerBoundForState(initialStateIndex, explorationInformation) << std::endl);
            val = bounds.getLowerBoundForState(initialStateIndex, explorationInformation) - val;
            // if (val==0) {config->setUpdateNow();}
            // config->writeBoundLog(bounds.getLowerBoundForState(initialStateIndex, explorationInformation));
        } else {
            // If not terminal state was found, the search aborted, possibly because of an EC-detection. In this
            // case, we cannot update the probabilities.
            STORM_LOG_TRACE("Did not find terminal state.");
        }
        config->update();

        STORM_LOG_DEBUG("Discovered states: " << explorationInformation.getNumberOfDiscoveredStates() << " (" << stats.numberOfExploredStates << " explored, "
                                              << explorationInformation.getNumberOfUnexploredStates() << " unexplored).");
        STORM_LOG_DEBUG("Value of initial state is in [" << bounds.getLowerBoundForState(initialStateIndex, explorationInformation) << ", "
                                                         << bounds.getUpperBoundForState(initialStateIndex, explorationInformation) << "].");
        ValueType difference = bounds.getDifferenceOfStateBounds(initialStateIndex, explorationInformation);
        STORM_LOG_DEBUG("Difference after iteration " << stats.pathsSampled << " is " << difference << ".");

        // If the number of sampled paths exceeds a certain threshold, do a precomputation.
        if (explorationInformation.performPrecomputationExcessiveSampledPaths(stats.pathsSampledSinceLastPrecomputation)) {
            brtdpInfo.performPrecomputation(stack, explorationInformation, bounds, stats);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    STORM_PRINT("Time to learn the DT from DTStrat: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000 << "s\n");

    auto mlHeuristic = config->getHeuristic();
    auto dt = mlHeuristic->getDT();
    decisionTree = new dtstrat::DTInterface<uint, ValueType>(dt);
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getDTControlTree(std::string const& dataFilename) {
    // Use a function that already exists to get a DTcontrol tree
    // Check out storm/modelchecker/explorationDT/heurstics/DecisionTreeHeuristic
    // this->decisionTree;
    const std::string fileName = dataFilename;
    std::fstream processedFile(fileName.c_str());
    std::string readString;
    std::getline(processedFile, readString);
    if (readString.empty() || readString[0] != '#') {
        uint size = std::count(readString.begin(), readString.end(), ',') + 1;
        processedFile.close();
        std::stringstream fileData;
        processedFile = std::fstream(fileName.c_str());
        fileData << "#NON-PERMISSIVE\n";
        fileData << "#BEGIN " << size << " 1\n";
        fileData << processedFile.rdbuf();
        processedFile.close();
        processedFile.open(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        processedFile << fileData.rdbuf();
    }
    std::string jsonString = dtstrat::BasicLearningConfig<StateType, ValueType>::callDTcontrol(dataFilename);
    std::ofstream ofs = std::ofstream("DTControl.json");
    ofs << jsonString;
    ofs.close();
    nlohmann2::json jsonTree = nlohmann2::json::parse(jsonString);
    this->decisionTree = new dtstrat::DTInterface<uint, ValueType>(jsonTree);
    STORM_PRINT("Generated DTCONTROL-tree\n");
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer{};
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getSklearnTree(std::string const& filename) {
    // Write a similar function as for dtcontrol to get a sklearn tree in dot-format
    // Sklearn.tree -> dot in module_dt.py
    // this->decisionTree = new dtstrat::DTInterface<uint, ValueType>(filename);
    std::string cmd =
        ". venv/bin/activate; "
        "python sklearn_dt.py " +
        filename;
    const char* run_cmd = cmd.c_str();
    const std::string& result = exec(run_cmd);
    STORM_LOG_ASSERT(result == "done\n", "ERROR! " + result);
    std::string dtfile = "sklearn.dot";
    this->decisionTree = new dtstrat::DTInterface<uint, ValueType>(dtfile);
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getMLPackTree(std::string const& filename) {
    auto res = readCSVFile(filename);
    auto data = res.first;
    auto labels = arma::conv_to<arma::Row<long unsigned int>>::from(res.second);
    arma::Row<long unsigned int> possibleActions = arma::unique(labels);
    uint numberOfActions = possibleActions.size();
    getMLPackTree(data, labels, numberOfActions);
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::storeData(std::string const& filename, arma::Mat<int> const& data,
                                                                     arma::Row<long unsigned int> const& labels) {
    arma::Mat<int> datat = data.t();
    datat.save(filename, arma::csv_ascii);
    std::fstream processesdFile(filename.c_str());
    std::stringstream fileData;
    fileData << "#NON-PERMISSIVE\n";
    fileData << "#BEGIN " << datat.n_cols << " 1\n";
    fileData << processesdFile.rdbuf();
    processesdFile.close();

    processesdFile.open(filename.c_str(), std::fstream::out | std::fstream::trunc);
    processesdFile << fileData.rdbuf();
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getMLPackTree(arma::Mat<int> const& data, arma::Row<long unsigned int> const& labels,
                                                                         uint const& numberOfActions) {
    this->decisionTree = new dtstrat::DTInterface<uint, ValueType>(data, labels, numberOfActions);
}

template<typename ModelType, typename StateType>
std::pair<arma::Mat<int>, arma::Row<int>> LearningStrategiesModelChecker<ModelType, StateType>::readCSVFile(std::string const& filename) {
    std::fstream processedFile(filename.c_str());
    std::string readString;
    std::getline(processedFile, readString);
    if (readString.empty() || readString[0] != '#') {
        arma::Mat<int> data;
        data.load(filename, arma::csv_ascii);
        data = data.t();
        arma::Row<int> labels;
        labels = data.row(data.n_rows - 1);
        data.shed_row(data.n_rows - 1);
        return {data, labels};
    } else {
        STORM_LOG_ERROR("NOT IMPLEMENTED: preprocess file with header. Please remove the header from the csv to continue.");
    }
}

template<typename ModelType, typename StateType>
float LearningStrategiesModelChecker<ModelType, StateType>::evaluateVIMC(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, bool const& semanticActions) {
    ExplorationInformationDT<StateType, ValueType> explorationInformation(
        checkTask.isOptimizationDirectionSet() ? checkTask.getOptimizationDirection() : storm::OptimizationDirection::Maximize, substprogram);

    // The first row group starts at action 0.
    explorationInformation.newRowGroup(0);
    std::map<std::string, storm::expressions::Expression> labelToExpressionMapping = substprogram.getLabelToExpressionMapping();
    generator::NextStateGeneratorOptions options;
    options.setBuildChoiceLabels();
    options.setBuildChoiceOrigins();
    options.setBuildStateValuations();
    storm::logic::UntilFormula const& untilFormula = checkTask.getFormula();
    storm::logic::Formula const& conditionFormula = untilFormula.getLeftSubformula();
    storm::logic::Formula const& targetFormula = untilFormula.getRightSubformula();
    StateGenerationDT<StateType, ValueType> stateGeneration(substprogram, explorationInformation,
                                                            conditionFormula.toExpression(substprogram.getManager(), labelToExpressionMapping),
                                                            targetFormula.toExpression(substprogram.getManager(), labelToExpressionMapping), options);

    StateValuationsDT stateValuations;
    stateGeneration.computeInitialStates();
    STORM_LOG_THROW(stateGeneration.getNumberOfInitialStates() == 1, storm::exceptions::NotSupportedException,
                    "Currently only models with one initial state are supported by the exploration engine.");
    createExplicitMCinDRN(stateGeneration, explorationInformation, stateValuations, semanticActions);
    return 0;
}

template<typename ModelType, typename StateType>
uint LearningStrategiesModelChecker<ModelType, StateType>::getNumberOfStateVariables(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask,
                                                                                     storm::prism::Program& program) {
    ExplorationInformationDT<StateType, ValueType> explorationInformation(
        checkTask.isOptimizationDirectionSet() ? checkTask.getOptimizationDirection() : storm::OptimizationDirection::Maximize, program);
    std::map<std::string, storm::expressions::Expression> labelToExpressionMapping = program.getLabelToExpressionMapping();
    generator::NextStateGeneratorOptions options;
    options.setBuildChoiceLabels();
    options.setBuildChoiceOrigins();
    options.setBuildStateValuations();
    storm::logic::UntilFormula const& untilFormula = checkTask.getFormula();
    storm::logic::Formula const& conditionFormula = untilFormula.getLeftSubformula();
    storm::logic::Formula const& targetFormula = untilFormula.getRightSubformula();
    StateGenerationDT<StateType, ValueType> stateGeneration(program, explorationInformation,
                                                            conditionFormula.toExpression(program.getManager(), labelToExpressionMapping),
                                                            targetFormula.toExpression(program.getManager(), labelToExpressionMapping), options);
    stateGeneration.computeInitialStates();
    StateType const currentStateId = stateGeneration.getFirstInitialState();
    auto unexploredIt = explorationInformation.findUnexploredState(currentStateId);
    storm::generator::CompressedState const& compressedState = unexploredIt->second;
    stateGeneration.load(compressedState, currentStateId);
    auto sv = stateGeneration.getStateValuation(currentStateId);
    return sv.size();
}

////Returns the candidate for the path.
// private Candidate<Pair<State, Integer>> updateCandidate(BRTDPInformation& brtdpInformation)
//{
//      brtdpInformation.
//      Candidate<Pair<State, Integer>> cnd = cndLst.size() > 0 ? cndLst.getLast() : null;
//     for (int i = lastIndx + 1; i <= path.size(); i++) {
//         State state = path.getState(i);
//         NatBitSet labels = computeLabels(state);
//         draState = dra.getEdgeDestByLabel(draState, labels);
//         Pair<State, Integer> st = new Pair<>(state, draState);
//
//         if (cnd != null && cnd.contains(st)) {
//             // state already in the current candidate
//             cnd.countState(st);
//             cnd.setNonTrivial();
//         } else {
//             Integer indx = cndMap.get(st);
//
//             if (indx == null) {
//                 // state seen for the first time
//                 cnd = new Candidate<>(i, st);
//                 cndLst.add(cnd);
//                 cndMap.put(st, cndLst.size() - 1);
//             } else {
//                 // state seen before - create a merge of all candidates from that point
//                 cnd = new Candidate<>(i, st);
//
//                 while (cndLst.size() > indx) {
//                     Candidate<Pair<State, Integer>> oldCnd = cndLst.removeLast();
//                     cnd.mergeWith(oldCnd);
//                 }
//
//                 cndLst.addLast(cnd);
//                 int v = cndLst.size() - 1;
//
//                 for (Pair<State, Integer> s : cnd.getStates()) {
//                     cndMap.put(s, v);
//                 }
//
//                 cnd.setNonTrivial();
//             }
//         }
//     }
//
//     cndTmr.stop();
//     return cnd;
// }

// template<typename ModelType, typename StateType>
// bool LearningStrategiesModelChecker<ModelType, StateType>::verifyCandidate(BSCCCandidate<StateType> cand)
//{
//     // check if the candidate is closed under transition relation
//     BSCCCandidate candidate;
//     set<StateType> states = cand.getStates();
//
//     for (StateType st : states) {
//         Collection<State> succs = getSucessors(st);
//
//         if (!states.containsAll(succs)) {
//             return false;
//         }
//     }
//     return true;
// }

template<typename ModelType, typename StateType>
double LearningStrategiesModelChecker<ModelType, StateType>::evaluateSMC(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, bool const& semanticActions) {
    ExplorationInformationDT<StateType, ValueType> explorationInformation(
        checkTask.isOptimizationDirectionSet() ? checkTask.getOptimizationDirection() : storm::OptimizationDirection::Maximize, substprogram);

    // The first row group starts at action 0.
    explorationInformation.newRowGroup(0);
    std::map<std::string, storm::expressions::Expression> labelToExpressionMapping = substprogram.getLabelToExpressionMapping();
    generator::NextStateGeneratorOptions options;
    options.setBuildChoiceLabels();
    options.setBuildChoiceOrigins();
    options.setBuildStateValuations();
    storm::logic::UntilFormula const& untilFormula = checkTask.getFormula();
    storm::logic::Formula const& conditionFormula = untilFormula.getLeftSubformula();
    storm::logic::Formula const& targetFormula = untilFormula.getRightSubformula();
    StateGenerationDT<StateType, ValueType> stateGeneration(substprogram, explorationInformation,
                                                            conditionFormula.toExpression(substprogram.getManager(), labelToExpressionMapping),
                                                            targetFormula.toExpression(substprogram.getManager(), labelToExpressionMapping), options);

    StateValuationsDT stateValuations;
    stateGeneration.computeInitialStates();
    STORM_LOG_THROW(stateGeneration.getNumberOfInitialStates() == 1, storm::exceptions::NotSupportedException,
                    "Currently only models with one initial state are supported by the exploration engine.");
    uint totalGood = 0;
    for (uint i = 0; i < numSMC; i++) {
        //STORM_PRINT("CURRENT SMCCOUNT: " << i << std::endl)
        bool result = samplePathFromInitialState(stateGeneration, explorationInformation, stateValuations, stopAfter, semanticActions);
        //STORM_PRINT("Good path:" << result << "\n");
        //STORM_PRINT("NUMBER OF SAMPLED PATHS SO FAR: " << i << "\n");
        totalGood += result;
        //STORM_PRINT("total Good paths so far:" << totalGood << "\n");
    }
    //STORM_PRINT("Evaluated Strategy: " << float(totalGood)/float(numSMC) << "\n");
    //STORM_PRINT("Total Good paths" << totalGood);
    return float(totalGood) / float(numSMC);
}

template<typename ModelType, typename StateType>
std::shared_ptr<storm::models::ModelBase> LearningStrategiesModelChecker<ModelType, StateType>::getModelFromProgram(
    prism::Program const& program, std::vector<storm::jani::Property> const& properties,
    boost::optional<std::vector<storm::jani::Property>> const& preprocessedProperties) {
    auto symbolicModel = storm::storage::SymbolicModelDescription(program);
    auto input = SymbolicInput();
    input.model = symbolicModel;
    input.properties = properties;
    input.preprocessedProperties = preprocessedProperties;
    auto mpi = ModelProcessingInformation();
    mpi.applyBisimulation = false;
    mpi.transformToJani = false;
    mpi.isCompatible = true;
    mpi.env = storm::Environment();
    mpi.ddType = dd::DdType::Sylvan;
    mpi.engine = storm::utility::Engine::Sparse;
    mpi.buildValueType = mpi.verificationValueType = storm::modelchecker::ModelProcessingInformation::ValueType::FinitePrecision;
    std::shared_ptr<storm::models::ModelBase> model;
    // copied from model-handling.h :: buildModelSparse()
    storm::builder::BuilderOptions options(createFormulasToRespect(input.properties), input.model.get());
    options.clearTerminalStates();
    options.setApplyMaximalProgressAssumption(false);
    options.setBuildAllRewardModels(true);
    options.setBuildChoiceLabels(true);
    options.setBuildStateValuations(true);
    options.setBuildAllLabels(true);
    options.setBuildChoiceOrigins(true);
    options.setReservedBitsForUnboundedVariables(32);
    options.setAddOutOfBoundsState(false);

    model = storm::api::buildSparseModel<ValueType>(input.model.get(), options, false, false);
    return model;
}

template<typename ModelType, typename StateType>
std::vector<std::shared_ptr<storm::logic::Formula const>> LearningStrategiesModelChecker<ModelType, StateType>::createFormulasToRespect(
    const std::vector<storm::jani::Property>& properties) {
    std::vector<std::shared_ptr<storm::logic::Formula const>> result = storm::api::extractFormulasFromProperties(properties);

    for (auto const& property : properties) {
        if (!property.getFilter().getStatesFormula()->isInitialFormula()) {
            result.push_back(property.getFilter().getStatesFormula());
        }
    }

    return result;
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::getValuationsAndLabels(storm::storage::Scheduler<ValueType> const& scheduler,
                                                                                  storm::models::sparse::Mdp<ValueType> const& model, arma::Mat<int>& data,
                                                                                  arma::Row<long unsigned int>& labels,
                                                                                  storm::storage::BitVector const& subsystem,
                                                                                  bool const& semanticActions) {
    storm::storage::sparse::StateValuations const& allSV = model.getStateValuations();
    int addValues = 0;
    bool adaptRows = (!data.empty());
    for (uint_fast64_t state = 0; state < model.getNumberOfStates(); state++) {
        if (subsystem[state]) {
            auto sv = allSV.getStateValuation(state);
            arma::Col<int> insertSV = sv.transform();
            if (adaptRows) {
                addValues = data.n_rows - insertSV.size();
                if (addValues < 0) {
                    data.insert_rows(data.n_rows, -addValues);
                }
                uint sizeBefore = insertSV.size();
                insertSV.resize(data.n_rows);
                for (uint i = 0; i < addValues; i++) {
                    insertSV(sizeBefore + i) = 0;
                }
            }
            data.insert_cols(data.n_cols, insertSV);
            uint_fast64_t chosenAction = scheduler.getChoice(state).getDeterministicChoice();
            if (semanticActions) {
                auto choiceLabels = model.getChoiceLabeling().getLabelsOfChoice(model.getTransitionMatrix().getRowGroupIndices()[state] + chosenAction);
                STORM_LOG_ASSERT(choiceLabels.size()<=1,"The action has too many labels!");
                std::string choiceLabel = *(choiceLabels.begin());
                chosenAction = completeActionNameToIndexMapping[choiceLabel];
            }
            labels.resize(data.n_cols);
            labels(data.n_cols - 1) = chosenAction;
        }
    }
}

template<typename ModelType, typename StateType>
bool LearningStrategiesModelChecker<ModelType, StateType>::samplePathFromInitialState(StateGenerationDT<StateType, ValueType>& stateGeneration,
                                                                                      ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                                      StateValuationsDT& stateValuations, int const& stopAfterNSteps, bool const& semanticActions) const {
    // todo change the actions
    // initialize a stack that will contain the sampled path
    std::vector<std::pair<StateType, ActionType>> stack;
    // initialize a set that contains all the state with unavailable DT suggested action
    std::set<StateType> statesWithUnavailDTAction;
    // Start the search from the initial state.
    stack.push_back(std::make_pair(stateGeneration.getFirstInitialState(), 0));
    // BSCCCandidate<StateType,ValueType> candidate;
    //  Initialize an updater for candidateDT
    CandidateUntilDT<StateType, ValueType> candidateUntilDt =
        storm::modelchecker::dtstrat::CandidateUntilDT<StateType, ValueType>(explorationInformation, stack, statesWithUnavailDTAction, decisionTree);

    // As long as we didn't find a terminal (accepting or rejecting) state in the search, sample a new successor.
    bool foundTerminalState = false;
    bool foundTargetState = false;
    //STORM_PRINT("New run\n")
    while (!foundTerminalState) {
        //        if (stopAfterNSteps > -1 && stack.size() >= stopAfterNSteps) {
        //            //STORM_PRINT("Stopped run\n");
        //            break;}
        //// the following condition is taken care in the update
        //        if(candidate.notTrivial() && verify(candidate)){
        //            //already hit a BSCC without reaching a terminal state
        //            return false;
        //        }
        StateType const& currentStateId = stack.back().first;
        //STORM_PRINT("Current state: "<<currentStateId << "\n");
        STORM_LOG_TRACE("State on top of stack is: " << currentStateId << ".");

        // If the state is not yet explored, we need to retrieve its behaviors.
        auto unexploredIt = explorationInformation.findUnexploredState(currentStateId);
        if (unexploredIt != explorationInformation.unexploredStatesEnd()) {
            STORM_LOG_TRACE("State was not yet explored.");

            // Explore the previously unexplored state.
            storm::generator::CompressedState const& compressedState = unexploredIt->second;
            std::pair<bool, bool> eval = exploreState(stateGeneration, currentStateId, compressedState, explorationInformation, stateValuations);
            foundTerminalState = eval.first;
            foundTargetState = eval.second;
            if (foundTerminalState) {
                //STORM_PRINT("TERMINAL STATE");
                auto stateValuation = stateValuations.getStateValuation(currentStateId);
                if (explorationInformation.isTarget(currentStateId)) {
                    //STORM_PRINT("TARGET STATE");
                    foundTargetState = true;
                }
            }
            explorationInformation.removeUnexploredState(unexploredIt);
        } else {
            // If the state was already explored, we check whether it is a terminal state or not.
            if (explorationInformation.isTerminal(currentStateId)) {
                auto stateValuation = stateValuations.getStateValuation(currentStateId);
                STORM_LOG_TRACE("Found already explored terminal state: " << currentStateId << ".");
                //STORM_PRINT(" is terminal: " << stateValuation << "\n");
                if (explorationInformation.isTarget(currentStateId)) {
                    //STORM_PRINT("TARGET STATE");
                    foundTargetState = true;
                }
                foundTerminalState = true;
            }
        }

        // If the state was not a terminal state, we continue the path search and sample the next state.
        if (!foundTerminalState) {
            // At this point, we can be sure that the state was expanded and that we can sample according to the
            // probabilities in the matrix.
            uint32_t chosenAction;
            auto stateValuation = stateValuations.getStateValuation(currentStateId);
            //STORM_PRINT(" is not terminal: " << stateValuation << "\n");
            bool worked = false;
            chosenAction = decisionTree->getDecision(stateValuation, currentStateId, worked);
            if (worked & semanticActions) {
                // get all action-names for this action
                std::vector<std::string> labels = actionIndexToLabelMapping[chosenAction];
                std::vector<uint> possibleActions;
                StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
                auto totalNumberOfActionsInState =
                    explorationInformation.getStartRowOfGroup(rowGroup + 1) - explorationInformation.getStartRowOfGroup(rowGroup);
                std::string finalLabel;
                for (uint32_t action = explorationInformation.getStartRowOfGroup(rowGroup); action < explorationInformation.getStartRowOfGroup(rowGroup + 1);
                     ++action) {
                    auto labelOfAction = explorationInformation.getLabel(action);
                    if (std::find(labels.begin(), labels.end(),labelOfAction)!=labels.end()) {
                        possibleActions.push_back(action-explorationInformation.getStartRowOfGroup(rowGroup));
                        if (finalLabel.empty()) finalLabel = labelOfAction;
                    }
                }
                if (possibleActions.empty()) {
                    worked = false;
                } else {
                    // currently just play the first action
                    // todo decide smartly
                    chosenAction = possibleActions[0];
                    //STORM_PRINT(" = " << chosenAction << " ("<< finalLabel << ")\n");
                }
            }
            if (!worked) {
                /// Sampling randomly for a lower bound
                chosenAction = sampleActionOfState(currentStateId, explorationInformation);
                statesWithUnavailDTAction.insert(currentStateId);
            }
            stack.back().second = chosenAction;
            STORM_LOG_TRACE("Sampled action " << chosenAction << " in state " << currentStateId << ".");
            //STORM_PRINT(" Chosen action: " << chosenAction << "\n");
            StateType successor = sampleSuccessorFromAction(currentStateId, chosenAction, explorationInformation, statesWithUnavailDTAction);
            STORM_LOG_TRACE("Sampled successor " << successor << " according to action " << chosenAction << " of state " << currentStateId << ".");
            //STORM_PRINT("Sampled successor " << successor << " according to action " << chosenAction << " of state " << currentStateId << ".\n");

            // initialize a BSCC candidate that is trivial initially (default)

            // CandidateUntilDT<StateType,ValueType>* candidateUntilDt = new storm::modelchecker::dtstrat::CandidateUntilDT<StateType,ValueType>();
            //  Update BSCCCandidate and check if we indeed reached a BSCC
            bool isACandidateBSCC = candidateUntilDt.update(explorationInformation, stack, statesWithUnavailDTAction, stateValuations, actionIndexToLabelMapping, semanticActions);

            // Put the successor state and a dummy action on top of the stack.
            stack.emplace_back(successor, 0);

            if (isACandidateBSCC) {
                //STORM_PRINT("***********************BSCC Reached************************** \n");
                return false;
            }
        }
    }
    //STORM_PRINT("foundTargetState" << foundTargetState << "\n");
    return foundTargetState;
}

// template<typename ModelType, typename StateType>
// void LearningStrategiesModelChecker<ModelType, StateType>::createExplicitMCinDRNOld(StateGenerationDT<StateType, ValueType>& stateGeneration,
//                                                                                       ExplorationInformationDT<StateType, ValueType>& explorationInformation,
//                                                                                       StateValuationsDT& stateValuations) {
//     // initialize a queue that will contain the states need to be considered
//     std::vector<StateType> queue;
//     // initialize set of marked states that has already been taken care
//     std::ofstream explicitMC_stream("markovchain-temp.drn", std::ofstream::trunc);
//     std::set<StateType> markedStates;
//     int numberOfStates= markedStates.size();
//     int seenStates = 0;
//     // corresponding state labeling
//     std::map<StateType, int> corrStateLabeling;
//     // Start the search from the initial state.
//     queue.push_back(stateGeneration.getFirstInitialState());
//     markedStates.insert(stateGeneration.getFirstInitialState());
//     // Before generating the behavior of the state, we need to determine whether it's a target state that
//     // does not need to be expanded.
//
//     //std::ostringstream explicitMC_stream;
//     bool foundTerminalState = false;
//     bool foundTargetState = false;
//     while (!queue.empty()) {
//         StateType currentStateId = queue.front();
//         numberOfStates+=1;
//         if (corrStateLabeling.find(currentStateId) == corrStateLabeling.end()){
//             corrStateLabeling[currentStateId]=seenStates; // storm2drnAdd()
//             seenStates+=1;
//         }
////        //STORM_PRINT("map values: ")
////        for(auto it = corrStateLabeling.cbegin(); it != corrStateLabeling.cend(); ++it)
////        {
////            std::cout << "(" << it->first << " " << it->second << "), \n";
////        }
////        //STORM_PRINT("picking the front of the queue: " << currentStateId << "\n");
//        queue.erase(queue.begin());
//        //STORM_PRINT("State on front of the queue: " << currentStateId << ".\n");
//
//        // If the state is not yet explored, we need to retrieve its behaviors.
//        auto unexploredIt = explorationInformation.findUnexploredState(currentStateId);
//        if (unexploredIt != explorationInformation.unexploredStatesEnd()) {
//            STORM_LOG_TRACE("State was not yet explored.");
//            // Explore the previously unexplored state.
//            storm::generator::CompressedState const& compressedState = unexploredIt->second;
//            std::pair<bool,bool> eval = exploreState(stateGeneration, currentStateId, compressedState, explorationInformation, stateValuations);
//            foundTerminalState = eval.first;
//            foundTargetState = eval.second;
//            if (foundTerminalState) {
//                if(foundTargetState){
//                    STORM_LOG_TRACE("A successor won't be pushed");
////                    //STORM_PRINT("A successor won't be pushed");
//                    if(stateGeneration.getFirstInitialState()== currentStateId){
////                        ofs << "state " << corrStateLabeling[currentStateId] << " done init\n"; ofs.close();
//                        explicitMC_stream << "state " << corrStateLabeling[currentStateId] << " done init\n";
//                    }else{
////                        ofs << result; ofs.close();
//                        explicitMC_stream << "state " << corrStateLabeling[currentStateId] << " done\n";
//                    }
//                    // ofs << result; ofs.close();
//                    explicitMC_stream << "\taction 0\n";
//                    explicitMC_stream << "\t\t"<< corrStateLabeling[currentStateId] <<" : 1\n";
//                }
//                else{
////                    ofs << result; ofs.close();
//                    explicitMC_stream << "state " << corrStateLabeling[currentStateId] << "\n";
//                    explicitMC_stream << "\taction 0\n";
//                    explicitMC_stream << "\t\t"<< corrStateLabeling[currentStateId] <<" : 1\n";
//                }
//            }
//            explorationInformation.removeUnexploredState(unexploredIt);
//        } else {
//            //STORM_PRINT("Not good!");
//            // If the state was already explored, we check whether it is a terminal state or not.
//            if (explorationInformation.isTerminal(currentStateId)) {
//                foundTerminalState = true;
//                auto stateValuation = stateValuations.getStateValuation(currentStateId);
//                STORM_LOG_TRACE("Found already explored terminal state: " << currentStateId << ".");
//                if (explorationInformation.isTarget(currentStateId)) {
//                    foundTargetState = true;
//                    STORM_LOG_TRACE("A successor won't be pushed");
////                    ofs << result; ofs.close();
//                    explicitMC_stream << "state " << corrStateLabeling[currentStateId];
//                    explicitMC_stream << "\taction 0";
//                    explicitMC_stream << "\t\t"<< corrStateLabeling[currentStateId] <<" : 1\n";
//                } else {
//                    STORM_LOG_TRACE("A successor won't be pushed");
//                    //STORM_PRINT("state " << currentStateId << "[0] not-done:terminal someNumber?");
//                    //STORM_PRINT("\taction 0 [0/1:what number is this?]\n");
//                    explicitMC_stream << "state " << corrStateLabeling[currentStateId];
//                    explicitMC_stream << "\taction 0";
//                    explicitMC_stream << "\t\t"<< corrStateLabeling[currentStateId] <<" : 1\n";
//                }
//            }
//        }
//
//        // If the state was not a terminal state, we continue the path search and sample the next state.
//        if (!foundTerminalState) {
//            // At this point, we can be sure that the state was expanded and that we can sample according to the
//            // probabilities in the matrix.
//            uint32_t chosenAction;
//            stateValuations.addStateValuation(currentStateId, stateGeneration.getStateValuation(currentStateId));
//            auto stateValuation = stateValuations.getStateValuation(currentStateId);
//            bool worked = false;
//            chosenAction = decisionTree->getDecision(stateValuation, currentStateId, worked);
//            if (!worked) { //and actionvalid?
//               /// this state will now be converted to state with single action with uniform probability of taking each action
//               StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
//               auto totalNumberOfActionsInState =  explorationInformation.getStartRowOfGroup(rowGroup + 1) -
//               explorationInformation.getStartRowOfGroup(rowGroup);
////               //STORM_PRINT("state " << currentStateId << " WITH NO DT ACTION");
////               //STORM_PRINT("\taction 0\n");
////               ofs << result; ofs.close();
//               explicitMC_stream << "state " << corrStateLabeling[currentStateId] << "\n";
//               explicitMC_stream << "\taction 0\n";
//
//               for (uint32_t action = explorationInformation.getStartRowOfGroup(rowGroup); action < explorationInformation.getStartRowOfGroup(rowGroup + 1);
//               ++action) {
//                   std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& row = explorationInformation.getRowOfMatrix(action);
////                   std::cout << "queue content before: " << "\n";
////                   for (auto x: queue)
////                       std::cout << x << ",";
////                   std::cout << "\n";
//                   for(auto entry: row){
////                     //STORM_PRINT("successor \t\t" << corrStateLabeling[entry.getColumn()] << " : " << entry.getValue()/totalNumberOfActionsInState <<
///"-probabilityOfReachingThisSucc / totalNumberOfActionsInState\n"); /                       ofs << result; ofs.close();
//                       explicitMC_stream << "\t\t" << corrStateLabeling[entry.getColumn()] << " : " << entry.getValue()/totalNumberOfActionsInState << "\n";
//                       if(markedStates.find(entry.getColumn())==markedStates.end()) {
//                           // std::cout << "the successor is not marked, so changes in the queue" << "\n";
//                           if (corrStateLabeling.find(entry.getColumn()) == corrStateLabeling.end()){
//                               corrStateLabeling[entry.getColumn()]=seenStates;
//                               seenStates+=1;
//                           }
//                           queue.push_back(entry.getColumn());
//                           markedStates.insert(entry.getColumn());
//                       } else {
////                           std::cout << "the successor is already marked, no change in the queue" << "\n";
//                       }
//                   }
////                   std::cout << "queue content after: " << "\n";
////                   for (auto x: queue)
////                       std::cout << x << ",";
////                   std::cout << "\n";
//               }
//            }
//            else{
////                //STORM_PRINT("state " << currentStateId << " WITH AVAILABLE DT ACTION");
//                if(stateGeneration.getFirstInitialState() == currentStateId){
////                    //STORM_PRINT("state " << currentStateId << "initial state");
////                    ofs << result; ofs.close();
//                    explicitMC_stream << "state " << corrStateLabeling[currentStateId] << " init\n";
//                }else{
////                    ofs << result; ofs.close();
////                    //STORM_PRINT("state " << currentStateId << "not a target");
//                    explicitMC_stream << "state " << corrStateLabeling[currentStateId] << "\n";
//                }
////                //STORM_PRINT("\taction "<< chosenAction << "\n");
//                // explicitMC_stream << "state " << currentStateId << "\n";
////                ofs << result; ofs.close();
//                explicitMC_stream << "\taction "<< chosenAction << "\n";
//
//                StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
//                auto initialAction = explorationInformation.getStartRowOfGroup(rowGroup);
//                auto lastPossibleAction = explorationInformation.getStartRowOfGroup(rowGroup+1);
//                std::vector<storm::storage::MatrixEntry<StateType, ValueType>> row;
//
//                if (initialAction+chosenAction>=lastPossibleAction) {
//                    ActionType act = sampleActionOfState(currentStateId, explorationInformation);
//                    row = explorationInformation.getRowOfMatrix(act);
//                } else {
//                    row = explorationInformation.getRowOfMatrix(initialAction + chosenAction);
//                }
//
//
////                std::cout << "queue content before: " << "\n";
////                for (auto x: queue)
////                    std::cout << x << ",";
////                std::cout << "\n";
//                for(auto entry: row){
////                    //STORM_PRINT("successor \t\t"<< entry.getColumn() <<" : "<< entry.getValue() << ": probabilityOfReachingThisSucc\n");
//                    if (corrStateLabeling.find(entry.getColumn()) == corrStateLabeling.end()){
//                        corrStateLabeling[entry.getColumn()]=seenStates; //stormTODRN()
//                        seenStates+=1;
//                    }
////                    ofs << result; ofs.close();
//                    explicitMC_stream << "\t\t"<< corrStateLabeling[entry.getColumn()] << " : "<< entry.getValue() <<"\n";
//                    if(markedStates.find(entry.getColumn())==markedStates.end()) {
////                        std::cout << "the successor is not marked, so changes in the queue" << "\n";
//                        queue.push_back(entry.getColumn());
//                        markedStates.insert(entry.getColumn());
//                    } else {
////                        std::cout << "the successor is already marked, no change in the queue" << "\n";
//                    }
//                }
////                std::cout << "queue content after: " << "\n";
////                for (auto x: queue)
////                    std::cout << x << ",";
////                std::cout << "\n";
//            }
//        }
//        stateValuations.addStateValuation(currentStateId, stateGeneration.getStateValuation(currentStateId));
////        //STORM_PRINT("current state valuations are: "<< std::endl);
////        for (int i = 0; i < stateValuations.getStateValuation(currentStateId).size(); ++i) {
////            //STORM_PRINT("i="  << i << ", SV: "<<stateValuations.getStateValuation(currentStateId)[i]<< std::endl);
////        }
////        //STORM_PRINT("numberOfStates:" << numberOfStates);
//    }
//    explicitMC_stream << "state " << numberOfStates << " done\n";
//    explicitMC_stream << "\taction 0\n";
//    explicitMC_stream << "\t\t"<< numberOfStates <<" : 1\n";
//    explicitMC_stream.close();
//    std::ifstream t("markovchain-temp.drn");
//    //"@type: DTMC\n@parameters\n\n@reward_models\n\n@nr_states\n" + std::to_string(numberOfStates+1) + "\n@nr_choices\n" + std::to_string(numberOfStates+1) +
//    "\n@model\n"; std::ofstream ofs("markovchain.drn", std::ofstream::trunc); ofs << "@type: DTMC\n@parameters\n\n@reward_models\n\n@nr_states\n" +
//    std::to_string(numberOfStates+1) + "\n@nr_choices\n" + std::to_string(numberOfStates+1) + "\n@model\n"; ofs << t.rdbuf(); ofs.close();
//}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::createExplicitMCinDRN(StateGenerationDT<StateType, ValueType>& stateGeneration,
                                                                                 ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                                 StateValuationsDT& stateValuations, bool const& semanticActions) {
    // todo change the actions
    //  initialize a queue that will contain the states need to be considered
    std::vector<StateType> queue;
    // initialize set of marked states that has already been taken care
    std::ofstream explicitMC_stream("markovchain-temp.drn", std::ofstream::trunc);

    // number of seen States
    int numberOfStatesInDRN = 0;
    int seenStates = 0;
    // get the initial state
    auto initState = stateGeneration.getFirstInitialState();
    // corresponding state labeling
    std::map<StateType, int> storm2drnStateLabeling;
    storm2drnStateLabeling[initState] = seenStates;

    // Start the search from the initial state.
    queue.push_back(initState);

    uint count_random_choices = 0;
    uint count_all_choices = 0;

    while (!queue.empty()) {
        StateType currentStateId = queue.front();
        // numberOfStatesExplored+=1;
        queue.erase(queue.begin());
        bool foundTerminalState = false;
        bool foundTargetState = false;
        STORM_LOG_TRACE("State should NOT have been yet explored.");
        // Explore the previously unexplored state.
        auto unexploredIt = explorationInformation.findUnexploredState(currentStateId);
        storm::generator::CompressedState const& compressedState = unexploredIt->second;
        std::pair<bool, bool> eval = exploreState(stateGeneration, currentStateId, compressedState, explorationInformation, stateValuations);
        stateValuations.addStateValuation(currentStateId, stateGeneration.getStateValuation(currentStateId));
        auto stateValuation = stateValuations.getStateValuation(currentStateId);
        //STORM_PRINT("Explore "<<queue.front()<< "=")
        // for (uint i=0; i<stateValuation.size(); i++) STORM_PRINT(stateValuation[i]<<",");
        //STORM_PRINT("\n");
        foundTerminalState = eval.first;
        foundTargetState = eval.second;
        if (foundTerminalState) {
            STORM_LOG_TRACE("A successor won't be pushed");
            if (foundTargetState) {
                if (stateGeneration.getFirstInitialState() == currentStateId) {
                    explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << " done init\n";
                } else {
                    explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << " done\n";
                }
                explicitMC_stream << "\taction 0\n";
                explicitMC_stream << "\t\t" << storm2drnStateLabeling[currentStateId] << " : 1\n";
            } else {
                if (stateGeneration.getFirstInitialState() == currentStateId) {
                    explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << " init\n";
                } else {
                    explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << " \n";
                }
                explicitMC_stream << "\taction 0\n";
                explicitMC_stream << "\t\t" << storm2drnStateLabeling[currentStateId] << " : 1\n";
            }
        }

        explorationInformation.removeUnexploredState(unexploredIt);

        if (!foundTerminalState) {
            uint32_t chosenAction;
            bool worked = false;
            count_all_choices++;
            chosenAction = decisionTree->getDecision(stateValuation, currentStateId, worked);
            if (!worked) count_random_choices++;
            //STORM_PRINT(" action from DT: " <<chosenAction << " =" << indexToShortActionName[chosenAction] << "\n ");
            if (worked & semanticActions) {
                // get all action-names for this action
                std::vector<std::string> labels = actionIndexToLabelMapping[chosenAction];
                std::vector<uint> possibleActions;
                StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
                auto totalNumberOfActionsInState =
                    explorationInformation.getStartRowOfGroup(rowGroup + 1) - explorationInformation.getStartRowOfGroup(rowGroup);
                std::string finalLabel;
                for (uint32_t action = explorationInformation.getStartRowOfGroup(rowGroup); action < explorationInformation.getStartRowOfGroup(rowGroup + 1);
                     ++action) {
                    auto labelOfAction = explorationInformation.getLabel(action);
                    if (std::find(labels.begin(), labels.end(),labelOfAction)!=labels.end()) {
                        possibleActions.push_back(action-explorationInformation.getStartRowOfGroup(rowGroup));
                        if (finalLabel.empty()) finalLabel = labelOfAction;
                    }
                }
                if (possibleActions.empty()) {
                    worked = false;
                } else {
                    // currently just play the first action
                    // todo decide smartly
                    chosenAction = possibleActions[0];
                    //STORM_PRINT(" = " << chosenAction << " ("<< finalLabel << ")\n");
                }
            }
            if (stateGeneration.getFirstInitialState() == currentStateId) {
                explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << " init\n";
            } else {
                explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << "\n";
            }
            explicitMC_stream << "\taction " << chosenAction << "\n";

            /// To add the queue to the successors we need to know what successors we have
            /// preparation to check if the action is valid or not
            StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
            auto initialAction = explorationInformation.getStartRowOfGroup(rowGroup);
            auto lastPossibleAction = explorationInformation.getStartRowOfGroup(rowGroup + 1);
            std::vector<storm::storage::MatrixEntry<StateType, ValueType>> row;

            bool isActionValid = initialAction + chosenAction < lastPossibleAction;
            //STORM_PRINT(" Is action valid: " <<isActionValid);
            //if (!isActionValid) STORM_PRINT("not valid action\n");
            if (worked && isActionValid) {
                row = explorationInformation.getRowOfMatrix(initialAction + chosenAction);
                //STORM_PRINT(" -> valid action, possible successors ");
                double prob = 0;
                uint count =0;
                std::vector<double> values;
                for (auto entry : row) {
                    //STORM_PRINT(entry.getColumn() << " ");
                    if (storm2drnStateLabeling.find(entry.getColumn()) == storm2drnStateLabeling.end()) {
                        seenStates += 1;
                        storm2drnStateLabeling[entry.getColumn()] = seenStates;  // stormTODRN()
                        queue.push_back(entry.getColumn());
                    } else {
                        /// the successor is already marked, no change in the queue
                    }
                    if (count == row.size()-1) {
                        explicitMC_stream << "\t\t" << storm2drnStateLabeling[entry.getColumn()] << " : " <<  (double)((int)((1.0 - prob)*pow(10,5)))/pow(10,5) << "\n";
                        //STORM_PRINT("1-prob " << (1.0 - prob) << "\n");
                    }
                    else {
                        explicitMC_stream << "\t\t" << storm2drnStateLabeling[entry.getColumn()] << " : " << (double)((int)((entry.getValue()) *pow(10,5)))/pow(10,5)<< "\n";
                        //STORM_PRINT("value " << entry.getValue() << "\n");
                    }
                    prob += (double)((int)((entry.getValue()) *pow(10,5)))/pow(10,5);
                    values.push_back(entry.getValue());
                    count++;
                    double v = 1.0/3.0;
                }
            }
            else {
                /// we don't want to sample, instead for a lower bound, we want to add
                ///  all the successors uniformly at random
                ///  this state will now be converted to state with single action with uniform probability of taking each action
                StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
                auto totalNumberOfActionsInState =
                    explorationInformation.getStartRowOfGroup(rowGroup + 1) - explorationInformation.getStartRowOfGroup(rowGroup);

                // explicitMC_stream << "state " << storm2drnStateLabeling[currentStateId] << "\n";
                // explicitMC_stream << "\taction 0\n";

                for (uint32_t action = explorationInformation.getStartRowOfGroup(rowGroup); action < explorationInformation.getStartRowOfGroup(rowGroup + 1);
                     ++action) {
                    std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& row = explorationInformation.getRowOfMatrix(action);
                    double total=0;
                    uint count=0;
                    std::vector<double> values;
                    for (auto entry : row) {
                        if (storm2drnStateLabeling.find(entry.getColumn()) == storm2drnStateLabeling.end()) {
                            seenStates += 1;
                            storm2drnStateLabeling[entry.getColumn()] = seenStates;
                            //STORM_PRINT("Column: " << entry.getColumn() <<"seenStates: " << seenStates << std::endl);
                            queue.push_back(entry.getColumn());
                        } else {
                            /// the successor is already marked, no change in the queue
                        }

                        if (count==row.size()-1)
                            explicitMC_stream << "\t\t" << storm2drnStateLabeling[entry.getColumn()] << " : " << (double)((int)(((1.0-total) / totalNumberOfActionsInState) *pow(10,5)))/pow(10,5)
                                              << "\n";
                        else
                        explicitMC_stream << "\t\t" << storm2drnStateLabeling[entry.getColumn()] << " : " << (double)((int)((entry.getValue()/ totalNumberOfActionsInState) *pow(10,5)))/pow(10,5)
                                          << "\n";
                        total +=(double)((int)((entry.getValue()) *pow(10,6)))/pow(10,6);
                        values.push_back(entry.getValue());
                        count++;
                    }
                }
            }
            //STORM_PRINT("\n");
        }
    }
    numberOfStatesInDRN = storm2drnStateLabeling.size();
    explicitMC_stream << "state " << numberOfStatesInDRN << " done\n";
    explicitMC_stream << "\taction 0\n";
    explicitMC_stream << "\t\t" << numberOfStatesInDRN << " : 1\n";
    explicitMC_stream.close();
    numberOfStatesInDRN += 1;

    std::ifstream t("markovchain-temp.drn");
    std::ofstream ofs("markovchain.drn", std::ofstream::trunc);
    ofs << "@type: DTMC\n@parameters\n\n@reward_models\n\n@nr_states\n" + std::to_string(numberOfStatesInDRN) + "\n@nr_choices\n" +
               std::to_string(numberOfStatesInDRN) + "\n@model\n";
    ofs << t.rdbuf();
    ofs.close();
    STORM_PRINT("Random-choices: "<< count_random_choices << " of total " << count_all_choices << " choices\n");
}

template<typename ModelType, typename StateType>
typename LearningStrategiesModelChecker<ModelType, StateType>::ActionType LearningStrategiesModelChecker<ModelType, StateType>::sampleActionOfState(
    StateType const& currentStateId, ExplorationInformationDT<StateType, ValueType> const& explorationInformation) const {
    // Get all the actions
    StateType rowGroup = explorationInformation.getRowGroup(currentStateId);

    // Check for cases in which we do not need to perform more work.
    if (explorationInformation.onlyOneActionAvailable(rowGroup)) {
        return explorationInformation.getStartRowOfGroup(rowGroup);
    }

    // If there are more choices to consider, choose randomly.
    STORM_LOG_TRACE("Sampling from actions leaving the state.");
    std::vector<ActionType> actions;

    for (uint32_t row = explorationInformation.getStartRowOfGroup(rowGroup); row < explorationInformation.getStartRowOfGroup(rowGroup + 1); ++row) {
        actions.push_back(row);
    }
    int randomIndex = rand() % actions.size();

    return actions[randomIndex];
}

template<typename ModelType, typename StateType>
std::pair<bool, bool> LearningStrategiesModelChecker<ModelType, StateType>::exploreState(StateGenerationDT<StateType, ValueType>& stateGeneration,
                                                                                         const StateType& currentStateId,
                                                                                         const generator::CompressedState& currentState,
                                                                                         ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                                         StateValuationsDT& stateValuations) const {
    bool isTerminalState = false;
    bool isTargetState = false;

    // Finally, map the unexplored state to the row group.
    explorationInformation.assignStateToNextRowGroup(currentStateId);
    STORM_LOG_TRACE("Assigning row group " << explorationInformation.getRowGroup(currentStateId) << " to state " << currentStateId << ".");

    // Before generating the behavior of the state, we need to determine whether it's a target state that
    // does not need to be expanded.
    stateGeneration.load(currentState, currentStateId);
    auto st = stateGeneration.getStateValuation(currentStateId);
    stateValuations.addStateValuation(currentStateId, stateGeneration.getStateValuation(currentStateId));
    if (stateGeneration.isTargetState()) {
        isTargetState = true;
        isTerminalState = true;
    } else if (stateGeneration.isConditionState()) {
        STORM_LOG_TRACE("Exploring state.");

        // If it needs to be expanded, we use the generator to retrieve the behavior of the new state.
        storm::generator::StateBehavior<ValueType, StateType> behavior = stateGeneration.expand();
        STORM_LOG_TRACE("State has " << behavior.getNumberOfChoices() << " choices.");

        // Clumsily check whether we have found a state that forms a trivial BMEC.
        bool otherSuccessor = false;
        for (auto const& choice : behavior) {
            for (auto const& entry : choice) {
                if (entry.first != currentStateId) {
                    otherSuccessor = true;
                    break;
                }
            }
        }
        isTerminalState = !otherSuccessor;

        // If the state was neither a trivial (non-accepting) terminal state nor a target state, we
        // need to store its behavior.
        if (!isTerminalState) {
            // Next, we insert the behavior into our matrix structure.
            StateType startAction = explorationInformation.getActionCount();
            explorationInformation.addActionsToMatrix(behavior.getNumberOfChoices());

            ActionType localAction = 0;

            for (auto const& choice : behavior) {
                explorationInformation.setDTStratIndex(startAction + localAction, choice);
                explorationInformation.setLabel(startAction + localAction, choice);
                for (auto const& entry : choice) {
                    explorationInformation.getRowOfMatrix(startAction + localAction).emplace_back(entry.first, entry.second);
                    STORM_LOG_TRACE("Found transition " << currentStateId << "-[" << (startAction + localAction) << ", " << entry.second << "]-> "
                                                        << entry.first << ".");
                }
                ++localAction;
            }

            // Terminate the row group.
            explorationInformation.terminateCurrentRowGroup();
        }
    } else {
        // In this case, the state is neither a target state nor a condition state and therefore a rejecting
        // terminal state.
        isTerminalState = true;
    }

    if (isTerminalState) {
        STORM_LOG_TRACE("State does not need to be explored, because it is " << (isTargetState ? "a target state" : "a rejecting terminal state") << ".");
        explorationInformation.addTerminalState(currentStateId);
        // Increase the size of the matrix, but leave the row empty.
        explorationInformation.addActionsToMatrix(1);
        explorationInformation.setDTStratIndex(explorationInformation.getActionCount() - 1);
        explorationInformation.setInternalLabel(explorationInformation.getActionCount() -1, "TERMINAL");
        // Terminate the row group.
        explorationInformation.newRowGroup();
        if (isTargetState)
            explorationInformation.addTargetState(currentStateId);
    }
    return {isTerminalState, isTargetState};
}

template<typename ModelType, typename StateType>
StateType LearningStrategiesModelChecker<ModelType, StateType>::sampleSuccessorFromAction(
    const StateType& currentStateId, const ActionType& chosenAction, const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
    std::set<StateType>& statesWithUnavailDTAction) const {
    StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
    auto initialAction = explorationInformation.getStartRowOfGroup(rowGroup);
    auto lastPossibleAction = explorationInformation.getStartRowOfGroup(rowGroup + 1);
    std::vector<storm::storage::MatrixEntry<StateType, ValueType>> row;
    if (initialAction + chosenAction >= lastPossibleAction) {
        statesWithUnavailDTAction.insert(currentStateId);
        ActionType act = sampleActionOfState(currentStateId, explorationInformation);
        row = explorationInformation.getRowOfMatrix(act);
    } else {
        row = explorationInformation.getRowOfMatrix(initialAction + chosenAction);
    }
    if (row.size() == 1) {
        return row.front().getColumn();
    }
    std::vector<ValueType> probabilities(row.size());
    std::transform(row.begin(), row.end(), probabilities.begin(),
                   [](storm::storage::MatrixEntry<StateType, ValueType> const& entry) { return entry.getValue(); });

    // Now sample according to the probabilities.
    std::discrete_distribution<StateType> distribution(probabilities.begin(), probabilities.end());
    //STORM_PRINT("probabilities: ");
    // for (auto val : probabilities) STORM_PRINT( val << " ");
    //STORM_PRINT("distr: " << distribution);
    auto chosen = distribution(randomGenerator);
    //STORM_PRINT(" -> chosen "<< chosen << "\n");
    return row[chosen].getColumn();
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::printHeader(const std::string& filename, uint const& size) {
    //std::ofstream ofs(filename);
    //ofs << "#NON-PERMISSIVE\n";
    //ofs << "#BEGIN " << size << " 1\n";
    //ofs.close();
}

template<typename ModelType, typename StateType>
void LearningStrategiesModelChecker<ModelType, StateType>::printToFile(std::string const& filename, storm::storage::Scheduler<ValueType> const& scheduler,
                                                                       storm::models::sparse::Mdp<ValueType> const& model,
                                                                       storm::storage::BitVector const& subsystem, bool const& semanticActions) {
    // todo change the actions here
    storm::storage::sparse::StateValuations const& allSV = model.getStateValuations();
    uint currentsize = 0;
    for (uint i = 0; i < allSV.getNumberOfStates(); i++) {
        // find a non-empty state-valuation
        auto firstSV = allSV.getStateValuation(i);
        if (firstSV.size() == 0)
            continue;
        else
            currentsize = firstSV.size();
    }
    uint addZeros = 0;
    std::fstream processedFile(filename.c_str());
    if (processedFile.good()) {
        std::string readString;
        std::getline(processedFile, readString);
        // STORM_LOG_THROW((readString.empty() || readString[0] = '#'),storm::exceptions::InvalidOperationException,"Something weird happened with file " +
        // filename);
        std::getline(processedFile, readString);
        std::stringstream test(readString);
        std::string sizestr;
        std::getline(test, sizestr, ' ');
        std::getline(test, sizestr, ' ');
        uint size = std::stoi(sizestr);
        if (size < currentsize) {
            std::stringstream fileData;
            fileData << "#NON-PERMISSIVE\n";
            fileData << "#BEGIN " << currentsize << " 1\n";
            // reread the whole file and add a column
            while (std::getline(processedFile, readString)) {
                // change this line
                // get last comma
                std::string firstpart;
                std::string secondpart;
                for (int i = readString.size() - 1; i >= 0; i--) {
                    if (readString[i] == ',') {
                        firstpart = readString.substr(0, i + 1);
                        secondpart = readString.substr(i + 1, readString.size() - (i + 1));
                        break;
                    }
                }
                fileData << firstpart;
                for (uint i = 0; i < (currentsize - size); i++) fileData << "0,";
                fileData << secondpart << "\n";
            }
            processedFile.close();
            processedFile.open(filename.c_str(), std::fstream::out | std::fstream::trunc);
            processedFile << fileData.rdbuf();
        } else if (size > currentsize) {
            // add values to the state-valuations before printing
            addZeros = size - currentsize;
        }
        processedFile.close();
    }
    std::ofstream ofs(filename, std::ios::app);
    for (uint_fast64_t state = 0; state < model.getNumberOfStates(); state++) {
        if (subsystem[state]) {
            auto sv = allSV.getStateValuation(state);
            for (uint i = 0; i < sv.size(); i++) ofs << sv[i] << ",";
            for (uint i = 0; i < addZeros; i++) ofs << "0,";
            uint_fast64_t chosenAction = scheduler.getChoice(state).getDeterministicChoice();
            if (semanticActions) {
                auto choiceLabels = model.getChoiceLabeling().getLabelsOfChoice(model.getTransitionMatrix().getRowGroupIndices()[state] + chosenAction);
                if (choiceLabels.empty()) {
                    STORM_LOG_WARN("There are no labels");
                }
                STORM_LOG_ASSERT(choiceLabels.size()<=1,"The action has too many labels!");
                std::string choiceLabel = *(choiceLabels.begin());
                chosenAction = completeActionNameToIndexMapping[choiceLabel];
            }
            ofs << chosenAction << "\n";
        }
    }
    ofs.close();
}

template class LearningStrategiesModelChecker<storm::models::sparse::Dtmc<double>, uint32_t>;
template class LearningStrategiesModelChecker<storm::models::sparse::Mdp<double>, uint32_t>;
}  // namespace modelchecker
}  // namespace storm
