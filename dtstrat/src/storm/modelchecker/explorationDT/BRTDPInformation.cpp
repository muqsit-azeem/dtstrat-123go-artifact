//
// Created by steffi on 07.01.22.
//
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/utility/graph.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"

namespace storm{
namespace modelchecker {
namespace exploration_detail_dt {
template<typename StateType, typename ValueType>
BRTDPInformation<StateType, ValueType>::BRTDPInformation(ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                         StateGenerationDT<StateType, ValueType>& stateGeneration, BoundsDT<StateType, ValueType>& bounds,
                                                         StatisticsDT<StateType, ValueType>& stats,
                                                         StateValuationsDT& stateValuations,
                                                         double const& precision)
    : explorationInformation(&explorationInformation),
      stateGeneration(&stateGeneration),
      bounds(&bounds),
      stats(&stats),
      stateValuations(&stateValuations),
      comparator(precision),
      randomGenerator(std::chrono::system_clock::now().time_since_epoch().count()){
          // Intentionally left empty
      };

template<typename StateType, typename ValueType>
StateGenerationDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getStateGeneration() {
    return *stateGeneration;
}

template<typename StateType, typename ValueType>
ExplorationInformationDT<unsigned int, double>& BRTDPInformation<StateType, ValueType>::getExplorationInformation() {
    return *explorationInformation;
}

template<typename StateType, typename ValueType>
StateGenerationDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getStateGeneration() const {
    return *stateGeneration;
}

template<typename StateType, typename ValueType>
ExplorationInformationDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getExplorationInformation() const {
    return *explorationInformation;
}

template<typename StateType, typename ValueType>
BoundsDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getBounds() {
    return *bounds;
}

template<typename StateType, typename ValueType>
BoundsDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getBounds() const {
    return *bounds;
}

template<typename StateType, typename ValueType>
StatisticsDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getStats() {
    return *stats;
}

template<typename StateType, typename ValueType>
StatisticsDT<StateType, ValueType>& BRTDPInformation<StateType, ValueType>::getStats() const {
    return *stats;
}

template<typename StateType, typename ValueType>
StateValuationsDT& BRTDPInformation<StateType, ValueType>::getStateValuations()  {
    return *stateValuations;
}

template<typename StateType, typename ValueType>
StateValuationsDT& BRTDPInformation<StateType, ValueType>::getStateValuations() const {
    return *stateValuations;
}

template<typename StateType, typename ValueType>
uint BRTDPInformation<StateType, ValueType>::getBRTDPActionOfState(StateType const& currentStateId) {
    std::vector<std::pair<StateType, std::pair<ValueType, ValueType>>> actionValues;
    StateType rowGroup = explorationInformation->getRowGroup(currentStateId);

    // Check for cases in which we do not need to perform more work.
    if (explorationInformation->onlyOneActionAvailable(rowGroup)) {
        return explorationInformation->getStartRowOfGroup(rowGroup);
    }
    for (uint32_t row = explorationInformation->getStartRowOfGroup(rowGroup);
         row < explorationInformation->getStartRowOfGroup(rowGroup + 1); ++row) {
        std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const &Actualrow = explorationInformation->getRowOfMatrix(
            row);
        if (Actualrow.size() == 1) {
            if (Actualrow.front().getColumn() == currentStateId) {
                // Self-loop action
                continue;
            }
        }
        //auto upperBoundState = bounds.getUpperBoundForState(currentStateId, explorationInformation);
        auto diff1 = 0;//1.0 - upperBoundState;
        auto diff2 = 0;//bounds.getDifferenceOfStateBounds(currentStateId, explorationInformation);
        if (comparator.isZero(diff1) && comparator.isZero(diff2))
            actionValues.push_back(std::make_pair(row, std::make_pair(bounds->getLowerBoundForAction(row),
                                                                      bounds->getUpperBoundForAction(row))));
        else
            actionValues.push_back(std::make_pair(row, std::make_pair(
                                                           bounds->getBoundForAction(explorationInformation->getOptimizationDirection(), row),
                                                           bounds->getLowerBoundForAction(row))));
    }
    if (explorationInformation->maximize()) {
        std::sort(actionValues.begin(), actionValues.end(),
                  [](std::pair<StateType, std::pair<ValueType, ValueType>> const &a,
                     std::pair<StateType, std::pair<ValueType, ValueType>> const &b) {
                      return a.second.first > b.second.first;
                  });
    } else {
        std::sort(actionValues.begin(), actionValues.end(),
                  [](std::pair<StateType, std::pair<ValueType, ValueType>> const &a,
                     std::pair<StateType, std::pair<ValueType, ValueType>> const &b) {
                      return a.second.first < b.second.first;
                  });
    }

    // Determine the first elements of the sorted range that agree on their value.
    auto end = ++actionValues.begin();
    while (end != actionValues.end() &&
           comparator.isEqual(actionValues.begin()->second.first, end->second.first)) {
        ++end;
    }
    if (explorationInformation->maximize()) {
        std::sort(actionValues.begin(), end, [](std::pair<StateType, std::pair<ValueType, ValueType>> const &a,
                                                std::pair<StateType, std::pair<ValueType, ValueType>> const &b) {
            return a.second.second > b.second.second;
        });
    } else {
        std::sort(actionValues.begin(), end, [](std::pair<StateType, std::pair<ValueType, ValueType>> const &a,
                                                std::pair<StateType, std::pair<ValueType, ValueType>> const &b) {
            return a.second.second < b.second.second;
        });
    }

    return actionValues.begin()->first;
}

template<typename StateType, typename ValueType>
std::string BRTDPInformation<StateType, ValueType>::printValuation(StateValuationDT const& valuation, bool const& printStateVariableAssignment, char const& seperator) {
    std::string result;
    auto assignment = stateGeneration->getStateVariableModuleAssociation();
    STORM_LOG_THROW(assignment.size()==valuation.size(),exceptions::UnexpectedException,"Number of state variables does not match the size of the state valuation.");
    for (uint i=0; i<valuation.size(); i++) {
        if (printStateVariableAssignment) {
            result += std::to_string(assignment[i]) + seperator;
        }
        result += std::to_string(valuation[i]) + seperator;
    }
    result.pop_back();
    return result;
}

template<typename StateType, typename ValueType>
std::vector<std::pair<uint, uint>> BRTDPInformation<StateType, ValueType>::getModuleAndActions(uint const& chosenAction, std::vector<storm::prism::Module> const& allModules, bool& warning) {
    dtstrat::Action dtstratAction = explorationInformation->getAction(chosenAction);
    auto modules = dtstratAction.getModules();
    auto cmdindices = dtstratAction.getCommandIndices();
    if (!warning & (modules.size()>1)) {
        STORM_LOG_WARN("Only module 0 is printed, but there are more.");
        warning = true;
    }
    std::vector<std::pair<uint,uint>> res;
    for (uint i = 0; i < 1; i++) {
        int modIndx = -1;
        auto modi = modules[i];
        auto it = find_if(allModules.begin(), allModules.end(),
                          [&modi](const storm::prism::Module &obj) {
                              return obj.getName() == modi;
                          });
        if (it != allModules.end()) {
            modIndx = it - allModules.begin();
        }
        res.push_back(std::make_pair(modIndx,cmdindices[i]));
    }
    return res;
}


template<typename StateType, typename ValueType>
std::string BRTDPInformation<StateType, ValueType>::getModuleAndAction(uint const& chosenAction, std::vector<storm::prism::Module> const& allModules, bool& warning, char const& seperator) {
    auto modulesAndActions = getModuleAndActions(chosenAction, allModules, warning);
    std::string res;
    for (auto pair : modulesAndActions) {
        res += std::to_string(pair.first) + seperator + std::to_string(pair.second) + seperator;
    }
    res.pop_back();
    return res;
}

template<typename StateType, typename ValueType>
void BRTDPInformation<StateType, ValueType>::printData(std::ofstream& ofs, dtstrat::LearningConfig<StateType>* config, bool const& onlyReachableStates,
                                                       bool const& printModules, bool const& printStateVariableAssignment, char const& seperator) {
    bool warning = false;
    storm::storage::BitVector reachableStates(stateValuations->getNumberOfStates(), true);
    if (onlyReachableStates) {
        reachableStates = forwardPass(config);
    }
    std::vector<storm::prism::Module> allModules = explorationInformation->getAllModules();
    for (StateType state =0; state<stateValuations->getNumberOfStates(); state++) {
        if (explorationInformation->isTerminal(state)) continue;
        if (!stateValuations->isEmpty(state) && reachableStates[state]) {
            // print state and action
            uint chosenAction = explorationInformation->getDTStratIndex(getBestActionOfState(state),state);
            ofs << printValuation(stateValuations->getStateValuation(state), printStateVariableAssignment, seperator);
            ofs << seperator;
            if (printModules) {
                ofs << getModuleAndAction(chosenAction, allModules, warning, seperator);
            } else {
                ofs << chosenAction;
            }
            ofs << "\n";
        }
    }
}

template<typename StateType, typename ValueType>
void BRTDPInformation<StateType, ValueType>::printDTcontrolHeader(std::ofstream &ofs, bool const& printModules) {
    STORM_LOG_ASSERT(!printModules,"Can't use dtcontrol and print modules!");
    stateGeneration->printDTcontrolHeader(ofs);
}

template<typename StateType, typename ValueType>
void BRTDPInformation<StateType, ValueType>::printModuleDTHeader(std::ofstream& ofs) {
    stateGeneration->printModuleDTHeader(ofs);
}

template<typename StateType>
void addState(StateType const& state, storm::storage::BitVector& states) {
    if (states.size()<=state) {
        states.resize(state+1);
    }
    states.set(state, true);
}

template<typename StateType, typename ValueType>
    void BRTDPInformation<StateType, ValueType>::addStateValuation(StateType const& state, arma::Mat<int>& data, arma::Row<long unsigned int>& chosenActions, std::pair<uint, uint>& chosenAction, arma::Row<long unsigned int>* modules) {
    auto stateValuation = stateValuations->getStateValuation(state);
    arma::Col<int> sV = stateValuation.transform();
    if (data.size()==0) data.resize(sV.size(),0);
    data.insert_cols(data.n_cols, sV);
    chosenActions.resize(chosenActions.size()+1);
    chosenActions(chosenActions.size()-1) = chosenAction.second;
    if (modules) {
        modules->resize(modules->size()+1);
        (*modules)(modules->size()-1) = chosenAction.first;
    }
}
template<typename StateType, typename ValueType>
uint BRTDPInformation<StateType, ValueType>::getDTActionOfState(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config) {
    auto stateValuation = stateValuations->getStateValuation(currentStateId);
    bool worked = false;
    return config->getDecision(stateValuation, currentStateId, worked);

}


template<typename StateType, typename ValueType>
uint BRTDPInformation<StateType, ValueType>::getBestActionOfState(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config) {
    if (config && config->useHeuristic()) {
        auto stateValuation = stateValuations->getStateValuation(currentStateId);
        bool worked = false;
        return config->getDecision(stateValuation, currentStateId, worked);
    } else {
        return getBRTDPActionOfState(currentStateId);
    }
}

template<typename StateType, typename ValueType>
std::vector<StateType> BRTDPInformation<StateType, ValueType>::sampleActionsVectorOfState(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config) {
    if (config && config->useHeuristic()) {
        auto stateValuation = stateValuations->getStateValuation(currentStateId);
        return config->getAllDecisions(stateValuation, currentStateId);
    } else {
        auto rowGroup = explorationInformation->getRowGroup(currentStateId);
        std::vector<StateType> actions;
        for(auto action=explorationInformation->getStartRowOfGroup(rowGroup); action<explorationInformation->getStartRowOfGroup(rowGroup+1);action++) {
            actions.push_back(action);
        }
        return actions;
    }
}

template<typename StateType, typename ValueType>
std::vector<std::vector<std::pair<StateType, ValueType>>> BRTDPInformation<StateType, ValueType>::getSuccessors(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config) {
    bool goalState = false;
    std::vector<std::vector<std::pair<StateType, ValueType>>>  successors;
    auto unexploredIt = explorationInformation->findUnexploredState(currentStateId);
    if (unexploredIt != explorationInformation->unexploredStatesEnd()) {
        return successors;
    } else {
        if (explorationInformation->isTerminal(currentStateId)) {
            return successors;
        } else {
            std::vector<std::pair<StateType, ValueType>> successorsForOneAction;
            uint chosenAction = getBestActionOfState(currentStateId, config);
            auto tempSucc = explorationInformation->getRowOfMatrix(chosenAction);
            bool selfLoopAction = true;
            for (auto v: tempSucc) {
                successorsForOneAction.push_back(std::make_pair(v.getColumn(),v.getValue()));
                if (v.getColumn() != currentStateId) {
                    selfLoopAction = false;
                }
            }
            successors.push_back(successorsForOneAction);
            if (successors.empty()) {
                //STORM_PRINT("Only self-loop actions for this state" << std::endl);
                // We only have one action that is self-looping, or we don't have an action at all
                // then add all possible successors for all actions
                std::vector<StateType> actionVector = sampleActionsVectorOfState(currentStateId, config);
                for (auto chosenAction: actionVector) {
                    auto tempSucc = explorationInformation->getRowOfMatrix(chosenAction);
                    std::vector<std::pair<StateType, ValueType>> successorsForOneAction;
                    for (auto v: tempSucc)
                        successorsForOneAction.push_back(std::make_pair(v.getColumn(),v.getValue()));
                    successors.push_back(successorsForOneAction);
                }
            }
        }
        return successors;
        // we should not have seen this state, should we?
    }
}


template<typename StateType, typename ValueType>
std::vector<StateType> BRTDPInformation<StateType, ValueType>::getSuccessorsOfStateDTActionPair(StateType const& currentStateId, ActionType chosenAction) {
    bool goalState = false;
    std::vector<StateType> successors;
    auto unexploredIt = explorationInformation->findUnexploredState(currentStateId);
    if (unexploredIt != explorationInformation->unexploredStatesEnd()) {
        return successors;
    } else {
        if (explorationInformation->isTerminal(currentStateId)) {
            return successors;
        } else {
            std::vector<StateType> successorsForChosenAction;
            auto succRow = explorationInformation->getRowOfMatrix(chosenAction);
            for (auto v: succRow) {
                successorsForChosenAction.push_back(v.getColumn());
            }
        }
        return successors;
    }
}

template<typename StateType, typename ValueType>
storm::storage::BitVector BRTDPInformation<StateType, ValueType>::forwardPass(dtstrat::LearningConfig<StateType>* config) {
    stateGeneration->computeInitialStates();
    STORM_LOG_THROW(stateGeneration->getNumberOfInitialStates() == 1, storm::exceptions::NotSupportedException,
                    "Currently only models with one initial state are supported by the exploration engine.");
    StateType initialStateIndex = stateGeneration->getFirstInitialState();
    storm::storage::BitVector seenStates(1,false);
    std::list<StateType> stack;
    stack.push_back(initialStateIndex);
    while (stack.size() > 0) {
        StateType currentState = stack.front();
        addState(currentState, seenStates);
        stack.pop_front();
        bool terminalState = false;
        std::vector<std::vector<std::pair<StateType, ValueType>>> successors = getSuccessors(currentState, config);
        for (const std::vector<std::pair<StateType, ValueType>>& successorsPerAction : successors) {
            for (std::pair<StateType, ValueType> stateAndProbability: successorsPerAction) {
                // we have not seen this state before: either because its index is greater than what we have stored or because it is marked as "false" in seen states
                if (stateAndProbability.first>=seenStates.size() || !(seenStates[stateAndProbability.first])) {
                    // don't put state several times in the stack
                    if (std::find(stack.begin(), stack.end(), stateAndProbability.first) == stack.end()) {
                        stack.push_back(stateAndProbability.first);
                    }
                }
            }
        }
    }
    return seenStates;
}



template<typename StateType, typename ValueType>
arma::Mat<int> BRTDPInformation<StateType, ValueType>::forwardPassStateValuations(arma::Row<long unsigned int>& chosenActions, arma::Row<long unsigned int>* modules, dtstrat::LearningConfig<StateType>* config) {
    std::vector<storm::prism::Module> allModules = explorationInformation->getAllModules();
    stateGeneration->computeInitialStates();
    STORM_LOG_THROW(stateGeneration->getNumberOfInitialStates() == 1, storm::exceptions::NotSupportedException,
                    "Currently only models with one initial state are supported by the exploration engine.");
    arma::Mat<int> data;
    StateType initialStateIndex = stateGeneration->getFirstInitialState();
    storm::storage::BitVector seenStates(1,false);
    std::list<StateType> stack;
    stack.push_back(initialStateIndex);
    while (stack.size() > 0) {
        StateType currentState = stack.front();
        addState(currentState, seenStates);
        uint chosenAction = explorationInformation->getDTStratIndex(getBestActionOfState(currentState), currentState);
        bool warning;
        auto moduleAndAction = getModuleAndActions(chosenAction, allModules, warning);
        addStateValuation(currentState, data, chosenActions, moduleAndAction[0], modules);
        stack.pop_front();
        bool terminalState = false;
        std::vector<std::vector<std::pair<StateType, ValueType>>> successors = getSuccessors(currentState, config);
        for (const std::vector<std::pair<StateType, ValueType>>& successorsPerAction : successors) {
            for (std::pair<StateType, ValueType> stateAndProbability: successorsPerAction) {
                // we have not seen this state before: either because its index is greater than what we have stored or because it is marked as "false" in seen states
                if (stateAndProbability.first>=seenStates.size() || !(seenStates[stateAndProbability.first])) {
                    if (explorationInformation->isUnexplored(stateAndProbability.first)) {
                        continue;
                    }
                    if (explorationInformation->isTerminal(stateAndProbability.first)) {
                        continue;
                    }
                    // don't put state several times in the stack
                    if (std::find(stack.begin(), stack.end(), stateAndProbability.first) == stack.end()) {
                        stack.push_back(stateAndProbability.first);
                    }
                }
            }
        }
    }
    return data;
}

template<typename StateType, typename ValueType>
uint BRTDPInformation<StateType, ValueType>::getNumberOfModules() const {
    return explorationInformation->getNumberOfModules();
}

template<typename StateType, typename ValueType>
std::string BRTDPInformation<StateType, ValueType>::getModelName() const {
    return explorationInformation->getModelName();
}

template<typename StateType, typename ValueType>
bool BRTDPInformation<StateType, ValueType>::samplePathFromInitialState(StateGenerationDT<StateType, ValueType>& stateGeneration,
                                                                        ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                        BRTDPInformation::StateActionStack& stack, BoundsDT<StateType, ValueType>& bounds,
                                                                        StatisticsDT<StateType, ValueType>& stats, LearningConfig<StateType>* config,
                                                                        StateValuationsDT& stateValuations,
                                                                        int const& stopAfterNSteps) const {
    uint countBTRDP=0;
    uint countDTStrat=0;
    // Start the search from the initial state.
    stack.push_back(std::make_pair(stateGeneration.getFirstInitialState(), 0));

    // As long as we didn't find a terminal (accepting or rejecting) state in the search, sample a new successor.
    bool foundTerminalState = false;
    while (!foundTerminalState) {
        if (stopAfterNSteps>-1 && stack.size()>=stopAfterNSteps) break;
        StateType const& currentStateId = stack.back().first;
        STORM_LOG_TRACE("State on top of stack is: " << currentStateId << ".");

        // If the state is not yet explored, we need to retrieve its behaviors.
        auto unexploredIt = explorationInformation.findUnexploredState(currentStateId);
        if (unexploredIt != explorationInformation.unexploredStatesEnd()) {
            STORM_LOG_TRACE("State was not yet explored.");

            // Explore the previously unexplored state.
            storm::generator::CompressedState const& compressedState = unexploredIt->second;
            foundTerminalState = exploreState(stateGeneration, currentStateId, compressedState, explorationInformation, bounds, stats, stateValuations);
            if (foundTerminalState) {
                auto stateValuation = stateValuations.getStateValuation(currentStateId);
                STORM_LOG_TRACE("Aborting sampling of path, because a terminal state was reached.");
            }
            explorationInformation.removeUnexploredState(unexploredIt);
        } else {
            // If the state was already explored, we check whether it is a terminal state or not.
            if (explorationInformation.isTerminal(currentStateId)) {
                auto stateValuation = stateValuations.getStateValuation(currentStateId);
                STORM_LOG_TRACE("Found already explored terminal state: " << currentStateId << ".");
                foundTerminalState = true;
            }
        }

        // Notify the stats about the performed exploration step.
        stats.explorationStep();

        // If the state was not a terminal state, we continue the path search and sample the next state.
        if (!foundTerminalState) {
            // At this point, we can be sure that the state was expanded and that we can sample according to the
            // probabilities in the matrix.
            uint32_t chosenAction;
            if (config->useHeuristic()) {
                auto stateValuation = stateValuations.getStateValuation(currentStateId);
                bool worked = false;
                chosenAction = config->getDecision(stateValuation, currentStateId, worked);
                if (!worked) {
                    if (stopAfterNSteps>-1) return false;
                    chosenAction = sampleActionOfState(currentStateId, explorationInformation, bounds);
                    countBTRDP++;
                } else countDTStrat++;
            } else {
                chosenAction = sampleActionOfState(currentStateId, explorationInformation, bounds);
                countBTRDP++;
            }
            stack.back().second = chosenAction;
            STORM_LOG_TRACE("Sampled action " << chosenAction << " in state " << currentStateId << ".");

            StateType successor = sampleSuccessorFromAction(chosenAction, explorationInformation, bounds);
            STORM_LOG_TRACE("Sampled successor " << successor << " according to action " << chosenAction << " of state " << currentStateId << ".");

            // Put the successor state and a dummy action on top of the stack.
            stack.emplace_back(successor, 0);

            // If the number of exploration steps exceeds a certain threshold, do a precomputation.
            if (explorationInformation.performPrecomputationExcessiveExplorationSteps(stats.explorationStepsSinceLastPrecomputation)) {
                performPrecomputation(stack, explorationInformation, bounds, stats);

                STORM_LOG_TRACE("Aborting the search after precomputation.");
                stack.clear();
                break;
            }
        }
    }
    STORM_LOG_TRACE("BRTDP: "<<countBTRDP<< " DTStrat: " << countDTStrat << "\n");
    return foundTerminalState;
}

template<typename StateType, typename ValueType>
bool BRTDPInformation<StateType, ValueType>::exploreState(StateGenerationDT<StateType, ValueType>& stateGeneration, const StateType& currentStateId,
                                                          const generator::CompressedState& currentState,
                                                          ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                          BoundsDT<StateType, ValueType>& bounds, StatisticsDT<StateType, ValueType>& stats,
                                                          StateValuationsDT& stateValuations) const {
        bool isTerminalState = false;
        bool isTargetState = false;

        ++stats.numberOfExploredStates;

        // Finally, map the unexplored state to the row group.
        explorationInformation.assignStateToNextRowGroup(currentStateId);
        STORM_LOG_TRACE("Assigning row group " << explorationInformation.getRowGroup(currentStateId) << " to state " << currentStateId << ".");

        // Initialize the bounds, because some of the following computations depend on the values to be available for
        // all states that have been assigned to a row-group.
        bounds.initializeBoundsForNextState();

        // Before generating the behavior of the state, we need to determine whether it's a target state that
        // does not need to be expanded.
        stateGeneration.load(currentState, currentStateId);
        stateValuations.addStateValuation(currentStateId, stateGeneration.getStateValuation(currentStateId));
        if (stateGeneration.isTargetState()) {
            ++stats.numberOfTargetStates;
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

                // Retrieve the lowest state bounds (wrt. to the current optimization direction).
                std::pair<ValueType, ValueType> stateBounds = getLowestBounds(explorationInformation.getOptimizationDirection());

                for (auto const& choice : behavior) {
                    explorationInformation.setDTStratIndex(startAction + localAction, choice);
                    for (auto const& entry : choice) {
                        explorationInformation.getRowOfMatrix(startAction + localAction).emplace_back(entry.first, entry.second);
                        STORM_LOG_TRACE("Found transition " << currentStateId << "-[" << (startAction + localAction) << ", " << entry.second << "]-> "
                                                            << entry.first << ".");
                    }

                    std::pair<ValueType, ValueType> actionBounds = computeBoundsOfAction(startAction + localAction, explorationInformation, bounds);
                    bounds.initializeBoundsForNextAction(actionBounds);
                    stateBounds = combineBounds(explorationInformation.getOptimizationDirection(), stateBounds, actionBounds);

                    STORM_LOG_TRACE("Initializing bounds of action " << (startAction + localAction) << " to "
                                                                     << bounds.getLowerBoundForAction(startAction + localAction) << " and "
                                                                     << bounds.getUpperBoundForAction(startAction + localAction) << ".");

                    ++localAction;
                }

                // Terminate the row group.
                explorationInformation.terminateCurrentRowGroup();

                bounds.setBoundsForState(currentStateId, explorationInformation, stateBounds);
                STORM_LOG_TRACE("Initializing bounds of state " << currentStateId << " to " << bounds.getLowerBoundForState(currentStateId, explorationInformation)
                                                                << " and " << bounds.getUpperBoundForState(currentStateId, explorationInformation) << ".");
            }
        } else {
            // In this case, the state is neither a target state nor a condition state and therefore a rejecting
            // terminal state.
            isTerminalState = true;
        }

        if (isTerminalState) {
            STORM_LOG_TRACE("State does not need to be explored, because it is " << (isTargetState ? "a target state" : "a rejecting terminal state") << ".");
            explorationInformation.addTerminalState(currentStateId);

            if (isTargetState) {
                bounds.setBoundsForState(currentStateId, explorationInformation,
                                         std::make_pair(storm::utility::one<ValueType>(), storm::utility::one<ValueType>()));
                bounds.initializeBoundsForNextAction(std::make_pair(storm::utility::one<ValueType>(), storm::utility::one<ValueType>()));
            } else {
                bounds.setBoundsForState(currentStateId, explorationInformation,
                                         std::make_pair(storm::utility::zero<ValueType>(), storm::utility::zero<ValueType>()));
                bounds.initializeBoundsForNextAction(std::make_pair(storm::utility::zero<ValueType>(), storm::utility::zero<ValueType>()));
            }

            // Increase the size of the matrix, but leave the row empty.
            explorationInformation.addActionsToMatrix(1);
            explorationInformation.setDTStratIndex(explorationInformation.getActionCount()-1);
            // Terminate the row group.
            explorationInformation.newRowGroup();
        }

        return isTerminalState;
}
template<typename StateType, typename ValueType>
typename BRTDPInformation<StateType, ValueType>::ActionType BRTDPInformation<StateType, ValueType>::sampleActionOfState(const StateType& currentStateId,
                                                                       const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                       BoundsDT<StateType, ValueType>& bounds) const {
        // Determine the values of all available actions.
        std::vector<std::pair<ActionType, ValueType>> actionValues;
        StateType rowGroup = explorationInformation.getRowGroup(currentStateId);

        // Check for cases in which we do not need to perform more work.
        if (explorationInformation.onlyOneActionAvailable(rowGroup)) {
            return explorationInformation.getStartRowOfGroup(rowGroup);
        }

        // If there are more choices to consider, start by gathering the values of relevant actions.
        STORM_LOG_TRACE("Sampling from actions leaving the state.");

        for (uint32_t row = explorationInformation.getStartRowOfGroup(rowGroup); row < explorationInformation.getStartRowOfGroup(rowGroup + 1); ++row) {
            actionValues.push_back(std::make_pair(row, bounds.getBoundForAction(explorationInformation.getOptimizationDirection(), row)));
        }

        STORM_LOG_ASSERT(!actionValues.empty(), "Values for actions must not be empty.");

        // Sort the actions wrt. to the optimization direction.
        if (explorationInformation.maximize()) {
            std::sort(actionValues.begin(), actionValues.end(),
                      [](std::pair<ActionType, ValueType> const& a, std::pair<ActionType, ValueType> const& b) { return a.second > b.second; });
        } else {
            std::sort(actionValues.begin(), actionValues.end(),
                      [](std::pair<ActionType, ValueType> const& a, std::pair<ActionType, ValueType> const& b) { return a.second < b.second; });
        }

        // Determine the first elements of the sorted range that agree on their value.
        auto end = ++actionValues.begin();
        while (end != actionValues.end() && comparator.isEqual(actionValues.begin()->second, end->second)) {
            ++end;
        }

        // Now sample from all maximizing actions.
        std::uniform_int_distribution<ActionType> distribution(0, std::distance(actionValues.begin(), end) - 1);
        return actionValues[distribution(randomGenerator)].first;
}

template<typename StateType, typename ValueType>
StateType BRTDPInformation<StateType, ValueType>::sampleSuccessorFromAction(const ActionType& chosenAction,
                                                                            const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                            const BoundsDT<StateType, ValueType>& bounds) const {
        std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& row = explorationInformation.getRowOfMatrix(chosenAction);
        if (row.size() == 1) {
            return row.front().getColumn();
        }

        // Depending on the selected next-state heuristic, we give the states other likelihoods of getting chosen.
        if (explorationInformation.useDifferenceProbabilitySumHeuristic() || explorationInformation.useProbabilityHeuristic()) {
            std::vector<ValueType> probabilities(row.size());
            if (explorationInformation.useDifferenceProbabilitySumHeuristic()) {
                std::transform(row.begin(), row.end(), probabilities.begin(),
                               [&bounds, &explorationInformation](storm::storage::MatrixEntry<StateType, ValueType> const& entry) {
                                   return entry.getValue() + bounds.getDifferenceOfStateBounds(entry.getColumn(), explorationInformation);
                               });
            } else if (explorationInformation.useProbabilityHeuristic()) {
                std::transform(row.begin(), row.end(), probabilities.begin(),
                               [](storm::storage::MatrixEntry<StateType, ValueType> const& entry) { return entry.getValue(); });
            }

            // Now sample according to the probabilities.
            std::discrete_distribution<StateType> distribution(probabilities.begin(), probabilities.end());
            return row[distribution(randomGenerator)].getColumn();
        } else {
            STORM_LOG_ASSERT(explorationInformation.useUniformHeuristic(), "Illegal next-state heuristic.");
            std::uniform_int_distribution<ActionType> distribution(0, row.size() - 1);
            return row[distribution(randomGenerator)].getColumn();
        }
}

template<typename StateType, typename ValueType>
bool BRTDPInformation<StateType, ValueType>::performPrecomputation(const BRTDPInformation::StateActionStack& stack,
                                                                   ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                   BoundsDT<StateType, ValueType>& bounds, StatisticsDT<StateType, ValueType>& stats) const {
    ++stats.numberOfPrecomputations;

    // Outline:
    // 1. construct a sparse transition matrix of the relevant part of the state space.
    // 2. use this matrix to compute states with probability 0/1 and an MEC decomposition (in the max case).
    // 3. use MEC decomposition to collapse MECs.
    STORM_LOG_TRACE("Starting " << (explorationInformation.useLocalPrecomputation() ? "local" : "global") << " precomputation.");

    // Construct the matrix that represents the fragment of the system contained in the currently sampled path.
    storm::storage::SparseMatrixBuilder<ValueType> builder(0, 0, 0, false, true, 0);

    // Determine the set of states that was expanded.
    std::vector<StateType> relevantStates;
    if (explorationInformation.useLocalPrecomputation()) {
        for (auto const& stateActionPair : stack) {
            if (explorationInformation.maximize() || !storm::utility::isOne(bounds.getLowerBoundForState(stateActionPair.first, explorationInformation))) {
                relevantStates.push_back(stateActionPair.first);
            }
        }
        std::sort(relevantStates.begin(), relevantStates.end());
        auto newEnd = std::unique(relevantStates.begin(), relevantStates.end());
        relevantStates.resize(std::distance(relevantStates.begin(), newEnd));
    } else {
        for (StateType state = 0; state < explorationInformation.getNumberOfDiscoveredStates(); ++state) {
            // Add the state to the relevant states if they are not unexplored.
            if (!explorationInformation.isUnexplored(state)) {
                relevantStates.push_back(state);
            }
        }
    }
    StateType sink = relevantStates.size();

    // Create a mapping for faster look-up during the translation of flexible matrix to the real sparse matrix.
    // While doing so, record all target states.
    std::unordered_map<StateType, StateType> relevantStateToNewRowGroupMapping;
    storm::storage::BitVector targetStates(sink + 1);
    for (StateType index = 0; index < relevantStates.size(); ++index) {
        relevantStateToNewRowGroupMapping.emplace(relevantStates[index], index);
        if (storm::utility::isOne(bounds.getLowerBoundForState(relevantStates[index], explorationInformation))) {
            targetStates.set(index);
        }
    }

    // Do the actual translation.
    StateType currentRow = 0;
    for (auto const& state : relevantStates) {
        builder.newRowGroup(currentRow);
        StateType rowGroup = explorationInformation.getRowGroup(state);
        for (auto row = explorationInformation.getStartRowOfGroup(rowGroup); row < explorationInformation.getStartRowOfGroup(rowGroup + 1); ++row) {
            ValueType unexpandedProbability = storm::utility::zero<ValueType>();
            for (auto const& entry : explorationInformation.getRowOfMatrix(row)) {
                auto it = relevantStateToNewRowGroupMapping.find(entry.getColumn());
                if (it != relevantStateToNewRowGroupMapping.end()) {
                    // If the entry is a relevant state, we copy it over (and compensate for the offset change).
                    builder.addNextValue(currentRow, it->second, entry.getValue());
                } else {
                    // If the entry is an unexpanded state, we gather the probability to later redirect it to an unexpanded sink.
                    unexpandedProbability += entry.getValue();
                }
            }
            if (unexpandedProbability != storm::utility::zero<ValueType>()) {
                builder.addNextValue(currentRow, sink, unexpandedProbability);
            }
            ++currentRow;
        }
    }
    // Then, make the unexpanded state absorbing.
    builder.newRowGroup(currentRow);
    builder.addNextValue(currentRow, sink, storm::utility::one<ValueType>());
    storm::storage::SparseMatrix<ValueType> relevantStatesMatrix = builder.build();
    storm::storage::SparseMatrix<ValueType> transposedMatrix = relevantStatesMatrix.transpose(true);
    STORM_LOG_TRACE("Successfully built matrix for precomputation.");

    storm::storage::BitVector allStates(sink + 1, true);
    storm::storage::BitVector statesWithProbability0;
    storm::storage::BitVector statesWithProbability1;
    if (explorationInformation.maximize()) {
        // If we are computing maximal probabilities, we first perform a detection of states that have
        // probability 01 and then additionally perform an MEC decomposition. The reason for this somewhat
        // duplicate work is the following. Optimally, we would only do the MEC decomposition, because we need
        // it anyway. However, when only detecting (accepting) MECs, we do not infer which of the other states
        // (not contained in MECs) also have probability 0/1.
        targetStates.set(sink, true);
        statesWithProbability0 = storm::utility::graph::performProb0A(transposedMatrix, allStates, targetStates);
        targetStates.set(sink, false);
        statesWithProbability1 =
            storm::utility::graph::performProb1E(relevantStatesMatrix, relevantStatesMatrix.getRowGroupIndices(), transposedMatrix, allStates, targetStates);

        storm::storage::MaximalEndComponentDecomposition<ValueType> mecDecomposition(relevantStatesMatrix, relevantStatesMatrix.transpose(true));
        ++stats.ecDetections;
        STORM_LOG_TRACE("Successfully computed MEC decomposition. Found " << (mecDecomposition.size() > 1 ? (mecDecomposition.size() - 1) : 0) << " MEC(s).");

        // If the decomposition contains only the MEC consisting of the sink state, we count it as 'failed'.
        STORM_LOG_ASSERT(mecDecomposition.size() > 0, "Expected at least one MEC (the trivial sink MEC).");
        if (mecDecomposition.size() == 1) {
            ++stats.failedEcDetections;
        } else {
            stats.totalNumberOfEcDetected += mecDecomposition.size() - 1;

            // 3. Analyze the MEC decomposition.
            for (auto const& mec : mecDecomposition) {
                // Ignore the (expected) MEC of the sink state.
                if (mec.containsState(sink)) {
                    continue;
                }

                collapseMec(mec, relevantStates, relevantStatesMatrix, explorationInformation, bounds);
            }
        }
    } else {
        // If we are computing minimal probabilities, we do not need to perform an EC-detection. We rather
        // compute all states (of the considered fragment) that have probability 0/1. For states with
        // probability 0, we have to mark the sink as being a target. For states with probability 1, however,
        // we must treat the sink as being rejecting.
        targetStates.set(sink, true);
        statesWithProbability0 =
            storm::utility::graph::performProb0E(relevantStatesMatrix, relevantStatesMatrix.getRowGroupIndices(), transposedMatrix, allStates, targetStates);
        targetStates.set(sink, false);
        statesWithProbability1 =
            storm::utility::graph::performProb1A(relevantStatesMatrix, relevantStatesMatrix.getRowGroupIndices(), transposedMatrix, allStates, targetStates);
    }

    // Set the bounds of the identified states.
    STORM_LOG_ASSERT((statesWithProbability0 & statesWithProbability1).empty(), "States with probability 0 and 1 overlap.");
    for (auto state : statesWithProbability0) {
        // Skip the sink state as it is not contained in the original system.
        if (state == sink) {
            continue;
        }

        StateType originalState = relevantStates[state];
        bounds.setUpperBoundForState(originalState, explorationInformation, storm::utility::zero<ValueType>());
        explorationInformation.addTerminalState(originalState);
    }
    for (auto state : statesWithProbability1) {
        // Skip the sink state as it is not contained in the original system.
        if (state == sink) {
            continue;
        }

        StateType originalState = relevantStates[state];
        bounds.setLowerBoundForState(originalState, explorationInformation, storm::utility::one<ValueType>());
        explorationInformation.addTerminalState(originalState);
    }
    return true;
}

template<typename StateType, typename ValueType>
void BRTDPInformation<StateType, ValueType>::collapseMec(const storage::MaximalEndComponent& mec, const std::vector<StateType>& relevantStates,
                                                         const storage::SparseMatrix<ValueType>& relevantStatesMatrix,
                                                         ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                         BoundsDT<StateType, ValueType>& bounds) const {
    bool containsTargetState = false;

    // Now we record all actions leaving the EC.
    std::vector<ActionType> leavingActions;
    std::vector<StateType> correspStates;
    std::vector<StateType> states;
    for (auto const& stateAndChoices : mec) {
        // Compute the state of the original model that corresponds to the current state.
        StateType originalState = relevantStates[stateAndChoices.first];
        StateType originalRowGroup = explorationInformation.getRowGroup(originalState);
        if (std::find(states.begin(), states.end(), originalState) == states.end()) {
            states.push_back(originalState);
        }

        // Check whether a target state is contained in the MEC.
        if (!containsTargetState && storm::utility::isOne(bounds.getLowerBoundForRowGroup(originalRowGroup))) {
            containsTargetState = true;
        }

        // For each state, compute the actions that leave the MEC.
        auto includedChoicesIt = stateAndChoices.second.begin();
        auto includedChoicesIte = stateAndChoices.second.end();
        for (auto action = explorationInformation.getStartRowOfGroup(originalRowGroup);
             action < explorationInformation.getStartRowOfGroup(originalRowGroup + 1); ++action) {
            if (includedChoicesIt != includedChoicesIte) {
                STORM_LOG_TRACE("Next (local) choice contained in MEC is "
                                << (*includedChoicesIt - relevantStatesMatrix.getRowGroupIndices()[stateAndChoices.first]));
                STORM_LOG_TRACE("Current (local) choice iterated is " << (action - explorationInformation.getStartRowOfGroup(originalRowGroup)));
                if (action - explorationInformation.getStartRowOfGroup(originalRowGroup) !=
                    *includedChoicesIt - relevantStatesMatrix.getRowGroupIndices()[stateAndChoices.first]) {
                    STORM_LOG_TRACE("Choice leaves the EC.");
                    leavingActions.push_back(action);
                    correspStates.push_back(originalState);
                } else {
                    STORM_LOG_TRACE("Choice stays in the EC.");
                    ++includedChoicesIt;
                }
            } else {
                STORM_LOG_TRACE("Choice leaves the EC, because there is no more choice staying in the EC.");
                leavingActions.push_back(action);
                correspStates.push_back(originalState);
            }
        }
    }

    // Now, we need to collapse the EC only if it does not contain a target state and the leaving actions are
    // non-empty, because only then have the states a (potentially) non-zero, non-one probability.
    if (!containsTargetState && !leavingActions.empty()) {
        // In this case, no target state is contained in the MEC, but there are actions leaving the MEC. To
        // prevent the simulation getting stuck in this MEC again, we replace all states in the MEC by a new
        // state whose outgoing actions are the ones leaving the MEC. We do this, by assigning all states in the
        // MEC to a new row group, which will then hold all the outgoing choices.

        // Remap all contained states to the new row group.
        StateType nextRowGroup = explorationInformation.getNextRowGroup();
        std::vector<std::vector<std::pair<StateType, int>>> predecessors;
        predecessors.resize(states.size());
        for (auto const& stateAndChoices : mec) {
            StateType originalState = relevantStates[stateAndChoices.first];
            StateType originalRowGroup = explorationInformation.getRowGroup(originalState);
            explorationInformation.assignStateToRowGroup(originalState, nextRowGroup);
            for (uint action = explorationInformation.getStartRowOfGroup(originalRowGroup);
                 action < explorationInformation.getStartRowOfGroup(originalRowGroup + 1); action++) {
                auto row = explorationInformation.getRowOfMatrix(action);
                for (auto entry: row) {
                    StateType succState = entry.getColumn();
                    int statesIndex;
                    auto it = find(states.begin(), states.end(), succState);
                    if (it != states.end()) {
                        statesIndex = it - states.begin();
                        std::pair<StateType, int> seachPair = {originalState, action};
                        auto it2 = std::find_if(predecessors[statesIndex].begin(),
                                                predecessors[statesIndex].end(),
                                                [&](const auto &pair) { return pair.first == originalState; });
                        if (it2 != predecessors[statesIndex].end()) {
                            // the state is already predecessor
                            int existantAction = predecessors[statesIndex][it2 -
                                                                           predecessors[statesIndex].begin()].second;
                            float beforeProb = 0;
                            auto row = explorationInformation.getRowOfMatrix(existantAction);
                            for (auto entry: row) {
                                if (entry.getColumn() == succState) {
                                    beforeProb = entry.getValue();
                                }
                            }
                            auto row2 = explorationInformation.getRowOfMatrix(action);
                            float nowProb = 0;
                            for (auto entry: row2) {
                                if (entry.getColumn() == succState) {
                                    nowProb = entry.getValue();
                                }
                            }
                            if (nowProb > beforeProb) {
                                predecessors[statesIndex][it2 -
                                                          predecessors[statesIndex].begin()].second = action;
                            }
                        } else {
                            predecessors[statesIndex].push_back(seachPair);
                        }
                    } else {
                        // state is not in states list (probably a state outside the MEC)
                    }
                }
            }
        }

        bounds.initializeBoundsForNextState();

        // Add to the new row group all leaving actions of contained states and set the appropriate bounds for
        // the actions and the new state.
        std::pair<ValueType, ValueType> stateBounds = getLowestBounds(explorationInformation.getOptimizationDirection());
        uint i = 0;
        for (auto const& action : leavingActions) {
            explorationInformation.moveActionToBackOfMatrix(action);
            std::pair<ValueType, ValueType> const& actionBounds = bounds.getBoundsForAction(action);
            bounds.initializeBoundsForNextAction(actionBounds);
            stateBounds = combineBounds(explorationInformation.getOptimizationDirection(), stateBounds, actionBounds);
            uint stormActionIndex = explorationInformation.getActionCount() - 1;
            // Store the additional information
            std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> liste;
            for (auto state: states) {
                std::vector<std::pair<StateType, int>> actions = explorationInformation.getPath(state, action, correspStates[i], predecessors, states);
                liste.push_back(std::make_pair(state, actions));
            }
            explorationInformation.setDTStratIndex(stormActionIndex, liste);

            i++;
        }
        bounds.setBoundsForRowGroup(nextRowGroup, stateBounds);

        // Terminate the row group of the newly introduced state.
        explorationInformation.terminateCurrentRowGroup();
    }
}

template<typename StateType, typename ValueType>
ValueType BRTDPInformation<StateType, ValueType>::getLowestBound(const OptimizationDirection& direction) const {
    if (direction == storm::OptimizationDirection::Maximize) {
        return storm::utility::zero<ValueType>();
    } else {
        return storm::utility::one<ValueType>();
    }
}

template<typename StateType, typename ValueType>
std::pair<ValueType, ValueType> BRTDPInformation<StateType, ValueType>::getLowestBounds(const OptimizationDirection& direction) const {
    ValueType val = getLowestBound(direction);
    return std::make_pair(val, val);
}

template<typename StateType, typename ValueType>
std::pair<ValueType, ValueType> BRTDPInformation<StateType, ValueType>::computeBoundsOfAction(
    const ActionType& action, const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
    const BoundsDT<StateType, ValueType>& bounds) const {
    // TODO: take into account self-loops?
    std::pair<ValueType, ValueType> result = std::make_pair(storm::utility::zero<ValueType>(), storm::utility::zero<ValueType>());
    for (auto const& element : explorationInformation.getRowOfMatrix(action)) {
        result.first += element.getValue() * bounds.getLowerBoundForState(element.getColumn(), explorationInformation);
        result.second += element.getValue() * bounds.getUpperBoundForState(element.getColumn(), explorationInformation);
    }
    return result;
}

template<typename StateType, typename ValueType>
std::pair<ValueType, ValueType> BRTDPInformation<StateType, ValueType>::combineBounds(const OptimizationDirection& direction,
                                                                                      const std::pair<ValueType, ValueType>& bounds1,
                                                                                      const std::pair<ValueType, ValueType>& bounds2) const {
    if (direction == storm::OptimizationDirection::Maximize) {
        return std::make_pair(std::max(bounds1.first, bounds2.first), std::max(bounds1.second, bounds2.second));
    } else {
        return std::make_pair(std::min(bounds1.first, bounds2.first), std::min(bounds1.second, bounds2.second));
    }
}
template<typename StateType, typename ValueType>
ValueType BRTDPInformation<StateType, ValueType>::updateProbabilityBoundsAlongSampledPath(
    BRTDPInformation::StateActionStack& stack, const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
    BoundsDT<StateType, ValueType>& bounds) const {
    ValueType result = 0.0;
    stack.pop_back();
    while (!stack.empty()) {
        result = updateProbabilityOfAction(stack.back().first, stack.back().second, explorationInformation, bounds);
        stack.pop_back();
    }
    return result;
}

template<typename StateType, typename ValueType>
ValueType BRTDPInformation<StateType, ValueType>::updateProbabilityOfAction(const StateType& state, const ActionType& action,
                                                                       const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                       BoundsDT<StateType, ValueType>& bounds) const {
    double result = 0.0;
    ValueType oldValue = bounds.getLowerBoundForState(state, explorationInformation);
    // Compute the new lower/upper values of the action.
    std::pair<ValueType, ValueType> newBoundsForAction = computeBoundsOfAction(action, explorationInformation, bounds);

    // And set them as the current value.
    bounds.setBoundsForAction(action, newBoundsForAction);

    // Check if we need to update the values for the states.
    if (explorationInformation.maximize()) {
        bounds.setLowerBoundOfStateIfGreaterThanOld(state, explorationInformation, newBoundsForAction.first);

        StateType rowGroup = explorationInformation.getRowGroup(state);
        if (newBoundsForAction.second < bounds.getUpperBoundForRowGroup(rowGroup)) {
            if (explorationInformation.getRowGroupSize(rowGroup) > 1) {
                newBoundsForAction.second = std::max(newBoundsForAction.second, computeBoundOverAllOtherActions(storm::OptimizationDirection::Maximize, state,
                                                                                                                action, explorationInformation, bounds));
            }

            bounds.setUpperBoundForRowGroup(rowGroup, newBoundsForAction.second);
        }
        ValueType newValue = bounds.getLowerBoundForState(state, explorationInformation);
        result = newValue - oldValue;
    } else {
        bounds.setUpperBoundOfStateIfLessThanOld(state, explorationInformation, newBoundsForAction.second);

        StateType rowGroup = explorationInformation.getRowGroup(state);
        if (bounds.getLowerBoundForRowGroup(rowGroup) < newBoundsForAction.first) {
            if (explorationInformation.getRowGroupSize(rowGroup) > 1) {
                ValueType min = computeBoundOverAllOtherActions(storm::OptimizationDirection::Minimize, state, action, explorationInformation, bounds);
                newBoundsForAction.first = std::min(newBoundsForAction.first, min);
            }

            bounds.setLowerBoundForRowGroup(rowGroup, newBoundsForAction.first);
        }
        ValueType newValue = bounds.getLowerBoundForState(state, explorationInformation);
        result = newValue - oldValue;
    }
    return result;
}
template<typename StateType, typename ValueType>
ValueType BRTDPInformation<StateType, ValueType>::computeBoundOverAllOtherActions(const OptimizationDirection& direction, const StateType& state,
                                                                                  const ActionType& action,
                                                                                  const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                                  const BoundsDT<StateType, ValueType>& bounds) const {
        ValueType bound = getLowestBound(explorationInformation.getOptimizationDirection());

        ActionType group = explorationInformation.getRowGroup(state);
        for (auto currentAction = explorationInformation.getStartRowOfGroup(group); currentAction < explorationInformation.getStartRowOfGroup(group + 1);
             ++currentAction) {
            if (currentAction == action) {
                continue;
            }

            if (direction == storm::OptimizationDirection::Maximize) {
                bound = std::max(bound, computeUpperBoundOfAction(currentAction, explorationInformation, bounds));
            } else {
                bound = std::min(bound, computeLowerBoundOfAction(currentAction, explorationInformation, bounds));
            }
        }
        return bound;
}
template<typename StateType, typename ValueType>
ValueType BRTDPInformation<StateType, ValueType>::computeLowerBoundOfAction(const ActionType& action,
                                                                            const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                            const BoundsDT<StateType, ValueType>& bounds) const {
    ValueType result = storm::utility::zero<ValueType>();
    for (auto const& element : explorationInformation.getRowOfMatrix(action)) {
        result += element.getValue() * bounds.getLowerBoundForState(element.getColumn(), explorationInformation);
    }
    return result;
}

template<typename StateType, typename ValueType>
ValueType BRTDPInformation<StateType, ValueType>::computeUpperBoundOfAction(const ActionType& action,
                                                                            const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                                            const BoundsDT<StateType, ValueType>& bounds) const {
    ValueType result = storm::utility::zero<ValueType>();
    for (auto const& element : explorationInformation.getRowOfMatrix(action)) {
        result += element.getValue() * bounds.getUpperBoundForState(element.getColumn(), explorationInformation);
    }
    return result;
}
template<typename StateType, typename ValueType>
std::pair<ValueType, ValueType> BRTDPInformation<StateType, ValueType>::computeBoundsOfState(
    const StateType& currentStateId, const ExplorationInformationDT<StateType, ValueType>& explorationInformation,
    const BoundsDT<StateType, ValueType>& bounds) const {
        StateType group = explorationInformation.getRowGroup(currentStateId);
        std::pair<ValueType, ValueType> result = getLowestBounds(explorationInformation.getOptimizationDirection());
        for (ActionType action = explorationInformation.getStartRowOfGroup(group); action < explorationInformation.getStartRowOfGroup(group + 1); ++action) {
            std::pair<ValueType, ValueType> actionValues = computeBoundsOfAction(action, explorationInformation, bounds);
            result = combineBounds(explorationInformation.getOptimizationDirection(), result, actionValues);
        }
        return result;
}
template<typename StateType, typename ValueType>
uint BRTDPInformation<StateType, ValueType>::getMaxNumberOfSteps() const {
    return stats->maxPathLength;
}

template<typename StateType, typename ValueType>
int BRTDPInformation<StateType, ValueType>::getRandomNum(int const& start, int const& end) {
        std::uniform_int_distribution<uint> distribution(start, end);
        return distribution(randomGenerator);
}
template<typename StateType, typename ValueType>
BRTDPInformation<StateType, ValueType>::~BRTDPInformation() {
}

template<typename StateType, typename ValueType>
arma::Mat<int> BRTDPInformation<StateType, ValueType>::getTransformedStateValuations(arma::Row<long unsigned int>& labels, arma::Row<long unsigned int>* modules) {
    std::vector<storm::prism::Module> allModules;
    if (modules)  allModules = explorationInformation->getAllModules();
    arma::Mat<int> result;
    for (uint i=0; i<stateValuations->getNumberOfStates(); i++){
        if (!stateValuations->isEmpty(i)) {
            if (!(explorationInformation->isTerminal(i))) {
                auto sv = stateValuations->getStateValuation(i);
                uint chosenAction = explorationInformation->getDTStratIndex(getBestActionOfState(i), i);
                result.insert_cols(result.n_cols, sv.transform());
                labels.resize(result.n_cols);
                if (modules) {
                    bool warning;
                    auto moduleAndAction = getModuleAndActions(chosenAction, allModules, warning)[0];
                    modules->resize(modules->size()+1);
                    (*modules)(modules->size()-1) = moduleAndAction.first;
                    labels(result.n_cols-1) = moduleAndAction.second;
                } else {
                    labels(result.n_cols-1) = chosenAction;
                }
            }
        }
    }
    return result;
}

template<typename StateType, typename ValueType>
uint BRTDPInformation<StateType, ValueType>::getNumberOfActions() const {
    return explorationInformation->getNumberOfActions();
}

template<typename StateType, typename ValueType>
arma::Row<int> BRTDPInformation<StateType, ValueType>::getAssignment() const {
    auto assignment = stateGeneration->getStateVariableModuleAssociation();
    arma::Row<int> res(assignment);
    return res;
}


template class BRTDPInformation<uint32_t, double>;
}
}
}
