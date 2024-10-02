#include "storm/modelchecker/explorationDT/SparseExplorationModelCheckerDT.h"

#include "storm/modelchecker/explorationDT/BoundsDT.h"
#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"
#include "storm/modelchecker/explorationDT/StateGenerationDT.h"
#include "storm/modelchecker/explorationDT/StatisticsDT.h"

#include "storm/generator/CompressedState.h"

#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/storage/SparseMatrix.h"

#include "storm/storage/prism/Program.h"

#include "storm/logic/FragmentSpecification.h"

#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"

#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/CoreSettings.h"
#include "storm/settings/modules/ExplorationSettings.h"

#include "storm/utility/constants.h"
#include "storm/utility/graph.h"
#include "storm/utility/macros.h"
#include "storm/utility/prism.h"

#include "storm/exceptions/InvalidOperationException.h"
#include "storm/exceptions/InvalidPropertyException.h"
#include "storm/exceptions/NotSupportedException.h"

#include "storm/modelchecker/explorationDT/learning/NonLearningConfig.h"
#include "storm/modelchecker/explorationDT/learning/BasicLearningConfig.h"
#include "storm/modelchecker/explorationDT/learning/GeneticLearningConfig.h"
#include "storm/modelchecker/explorationDT/learning/ExternalConfig.h"
#include "storm/modelchecker/explorationDT/learning/BasicLearningConfigHeadstart.h"

#include "storm/settings/modules/HeuristicLearningSettings.h"
#include "storm/settings/modules/IOSettings.h"

#include <chrono>

namespace storm {
namespace modelchecker {

template<typename ModelType, typename StateType>
SparseExplorationDTStratModelChecker<ModelType, StateType>::SparseExplorationDTStratModelChecker(storm::prism::Program const& program)
    : program(program.substituteConstantsFormulas()),
      randomGenerator(std::chrono::system_clock::now().time_since_epoch().count()),
      comparator(storm::settings::getModule<storm::settings::modules::ExplorationSettings>().getPrecision()) {
    // Intentionally left empty.
}

template<typename ModelType, typename StateType>
bool SparseExplorationDTStratModelChecker<ModelType, StateType>::canHandle(CheckTask<storm::logic::Formula, ValueType> const& checkTask) const {
    storm::logic::Formula const& formula = checkTask.getFormula();
    storm::logic::FragmentSpecification fragment = storm::logic::reachability();
    return formula.isInFragment(fragment) && checkTask.isOnlyInitialStatesRelevantSet();
}

template<typename ModelType, typename StateType>
std::unique_ptr<CheckResult> SparseExplorationDTStratModelChecker<ModelType, StateType>::computeUntilProbabilities(
    Environment const& env, CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) {
    storm::logic::UntilFormula const& untilFormula = checkTask.getFormula();
    storm::logic::Formula const& conditionFormula = untilFormula.getLeftSubformula();
    storm::logic::Formula const& targetFormula = untilFormula.getRightSubformula();
    STORM_LOG_THROW(program.isDeterministicModel() || checkTask.isOptimizationDirectionSet(), storm::exceptions::InvalidPropertyException,
                    "For nondeterministic systems, an optimization direction (min/max) must be given in the property.");

    ExplorationInformationDT<StateType, ValueType> explorationInformation(checkTask.isOptimizationDirectionSet() ? checkTask.getOptimizationDirection()
                                                                                                               : storm::OptimizationDirection::Maximize,
                                                                        program);

    // The first row group starts at action 0.
    explorationInformation.newRowGroup(0);

    std::map<std::string, storm::expressions::Expression> labelToExpressionMapping = program.getLabelToExpressionMapping();
    generator::NextStateGeneratorOptions options;
    options.setBuildChoiceLabels();
    options.setBuildChoiceOrigins();
    options.setBuildStateValuations();
    StateGenerationDT<StateType, ValueType> stateGeneration(program, explorationInformation,
                                                          conditionFormula.toExpression(program.getManager(), labelToExpressionMapping),
                                                          targetFormula.toExpression(program.getManager(), labelToExpressionMapping),
                                                          options);
    // Create a structure that holds the bounds for the states and actions.
    BoundsDT<StateType, ValueType> bounds;

    // Now perform the actual sampling.
    StatisticsDT<StateType, ValueType> stats;

    // Create the state valuations storage
    StateValuationsDT stateValuations;

    // Create the BRTDPInformation class that contains pointer to the important information of BRTDP
    BRTDPInformation<StateType, ValueType> brtdpInfo(explorationInformation, stateGeneration, bounds, stats, stateValuations, storm::settings::getModule<storm::settings::modules::ExplorationSettings>().getPrecision());

    // Compute and return re`sult.
    std::tuple<StateType, ValueType, ValueType> boundsForInitialState = performExploration(brtdpInfo);
    return std::make_unique<ExplicitQuantitativeCheckResult<ValueType>>(std::get<0>(boundsForInitialState), std::get<1>(boundsForInitialState));
}

template<typename ModelType, typename StateType>
std::tuple<StateType, typename ModelType::ValueType, typename ModelType::ValueType>
SparseExplorationDTStratModelChecker<ModelType, StateType>::performExploration(BRTDPInformation<StateType, ValueType>& brtdpInfo) const {
    // start timer
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
    StateActionStack stack;

    LearningConfig<StateType>* config = getLearningConfig(brtdpInfo);
    config->firstLog();
    config->prepareBoundLog();

    bool convergenceCriterionMet = false;
    while (!convergenceCriterionMet) {
        bool result = brtdpInfo.samplePathFromInitialState(stateGeneration, explorationInformation, stack, bounds, stats, config, stateValuations);

        stats.sampledPath();
        stats.updateMaxPathLength(stack.size());

        // If a terminal state was found, we update the probabilities along the path contained in the stack.
        if (result) {
            // Update the bounds along the path to the terminal state.
            STORM_LOG_TRACE("Found terminal state, updating probabilities along path.");
            auto val = bounds.getLowerBoundForState(initialStateIndex, explorationInformation);
            brtdpInfo.updateProbabilityBoundsAlongSampledPath(stack, explorationInformation, bounds);
            STORM_LOG_TRACE(bounds.getLowerBoundForState(initialStateIndex, explorationInformation) <<std::endl);
            val = bounds.getLowerBoundForState(initialStateIndex, explorationInformation) - val;
            //if (val==0) {config->setUpdateNow();}
            //config->writeBoundLog(bounds.getLowerBoundForState(initialStateIndex, explorationInformation));
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
        convergenceCriterionMet = comparator.isZero(difference);

        // If the number of sampled paths exceeds a certain threshold, do a precomputation.
        if (!convergenceCriterionMet && explorationInformation.performPrecomputationExcessiveSampledPaths(stats.pathsSampledSinceLastPrecomputation)) {
            brtdpInfo.performPrecomputation(stack, explorationInformation, bounds, stats);
        }
    }

    // Show statistics if required.
    if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isShowStatisticsSet()) {
        stats.printToStream(std::cout, explorationInformation);
    }

    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    if (ioSettings.isExportSchedulerSet()) {
        config->printHeuristic(ioSettings.getExportSchedulerFilename());
    }
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = (end - start) / std::chrono::milliseconds(1);
    config->finalLog(elapsed);
    return std::make_tuple(initialStateIndex, bounds.getLowerBoundForState(initialStateIndex, explorationInformation),
                           bounds.getUpperBoundForState(initialStateIndex, explorationInformation));
}

template<typename ModelType, typename StateType>
LearningConfig<StateType>* SparseExplorationDTStratModelChecker<ModelType, StateType>::getLearningConfig(BRTDPInformation<StateType, ValueType>& brtdpInformation) const {
    const auto& configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    configTypes configType = configSettings.getLearningConfigType();
    if (configType == configTypes::basic) {
        double explorationEpsilon = configSettings.getHeuristicProbability();
        int updatesteps = configSettings.getNumberOfSimulationsBeforeUpdate();
        auto heuristicType = configSettings.getHeuristicType();
        return new BasicLearningConfig<StateType, ValueType>(explorationEpsilon, updatesteps, heuristicType, brtdpInformation);
    } else if (configType == configTypes::genetic) {
        double explorationEpsilon = configSettings.getHeuristicProbability();
        int updatesteps = configSettings.getNumberOfSimulationsBeforeUpdate();
        auto heuristicType = configSettings.getHeuristicType();
        return new GeneticLearningConfig<StateType, ValueType>(explorationEpsilon, updatesteps, heuristicType, brtdpInformation);
    } else if (configType == configTypes::external) {
        double explorationEpsilon = configSettings.getHeuristicProbability();
        return new ExternalConfig<StateType, ValueType>(brtdpInformation.getModelName(), brtdpInformation, explorationEpsilon);
    }
    else if (configType == configTypes::headstart) {
        double explorationEpsilon = configSettings.getHeuristicProbability();
        int updatesteps = configSettings.getNumberOfSimulationsBeforeUpdate();
        auto heuristicType = configSettings.getHeuristicType();
        auto filename = configSettings.getExternalFilename();
        return new BasicLearningConfigHeadstart<StateType, ValueType>(explorationEpsilon, updatesteps, heuristicType, brtdpInformation, filename);
    }
    else if (configType == configTypes::noneConfig) {
        return new NonLearningConfig<StateType>(brtdpInformation.getModelName());
    }
    else {
        return new NonLearningConfig<StateType>(brtdpInformation.getModelName());
    }
}

template class SparseExplorationDTStratModelChecker<storm::models::sparse::Dtmc<double>, uint32_t>;
template class SparseExplorationDTStratModelChecker<storm::models::sparse::Mdp<double>, uint32_t>;
}  // namespace modelchecker
}  // namespace storm
