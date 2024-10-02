#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/ExplorationSettings.h"
#include "storm/settings/modules/IOSettings.h"

#include "storm/utility/macros.h"

#include "storm/exceptions/NotImplementedException.h"

namespace storm {
namespace modelchecker {
namespace exploration_detail_dt {

template<typename StateType, typename ValueType>
ExplorationInformationDT<StateType, ValueType>::ExplorationInformationDT(storm::OptimizationDirection const& direction, storm::prism::Program& model,
                                                                     ActionType const& unexploredMarker)
    : unexploredMarker(unexploredMarker),
      optimizationDirection(direction),
      localPrecomputation(false),
      numberOfExplorationStepsUntilPrecomputation(100000),
      numberOfSampledPathsUntilPrecomputation(),
      nextStateHeuristic(storm::settings::modules::ExplorationSettings::NextStateHeuristic::DifferenceProbabilitySum),
      model(&model) {
    storm::settings::modules::ExplorationSettings const& settings = storm::settings::getModule<storm::settings::modules::ExplorationSettings>();
    localPrecomputation = settings.isLocalPrecomputationSet();
    numberOfExplorationStepsUntilPrecomputation = settings.getNumberOfExplorationStepsUntilPrecomputation();
    if (settings.isNumberOfSampledPathsUntilPrecomputationSet()) {
        numberOfSampledPathsUntilPrecomputation = settings.getNumberOfSampledPathsUntilPrecomputation();
    }

    nextStateHeuristic = settings.getNextStateHeuristic();
    getAllActions();
    moduleActionToDTStrat = getModuleActionToDTStratMapping();
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::getAllActions() {
    dtstrat::getAllActionsFromModel(allActions, *model);
}

template<typename StateType, typename ValueType>
typename ExplorationInformationDT<StateType, ValueType>::const_iterator ExplorationInformationDT<StateType, ValueType>::findUnexploredState(
    StateType const& state) const {
    return unexploredStates.find(state);
}

template<typename StateType, typename ValueType>
typename ExplorationInformationDT<StateType, ValueType>::const_iterator ExplorationInformationDT<StateType, ValueType>::unexploredStatesEnd() const {
    return unexploredStates.end();
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::removeUnexploredState(const_iterator it) {
    unexploredStates.erase(it);
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::addUnexploredState(StateType const& stateId, storm::generator::CompressedState const& compressedState) {
    stateToRowGroupMapping.push_back(unexploredMarker);
    unexploredStates[stateId] = compressedState;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::assignStateToRowGroup(StateType const& state, ActionType const& rowGroup) {
    stateToRowGroupMapping[state] = rowGroup;
}

template<typename StateType, typename ValueType>
StateType ExplorationInformationDT<StateType, ValueType>::assignStateToNextRowGroup(StateType const& state) {
    stateToRowGroupMapping[state] = rowGroupIndices.size() - 1;
    return stateToRowGroupMapping[state];
}

template<typename StateType, typename ValueType>
StateType ExplorationInformationDT<StateType, ValueType>::getNextRowGroup() const {
    return rowGroupIndices.size() - 1;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::newRowGroup(ActionType const& action) {
    rowGroupIndices.push_back(action);
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::newRowGroup() {
    newRowGroup(matrix.size());
}

template<typename StateType, typename ValueType>
int ExplorationInformationDT<StateType, ValueType>::findActionInList(dtstrat::Action const &searchAction) {
    int foundIndex = -1;
    auto it = std::find(allActions.begin(), allActions.end(), searchAction);
    if (it != allActions.end()) {
        foundIndex = it - allActions.begin();
    }
    if (foundIndex == -1) {
        STORM_PRINT("Action not in actionList");
    };
    return foundIndex;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::setDTStratIndex(const uint& action) {
    if (rowToDTStratActionMapping.size()<action) {
        rowToDTStratActionMapping.resize(action);
    }
    if (rowToDTStratActionMapping.size() == action) {
        rowToDTStratActionMapping.push_back(dtstrat::ActionRepresentation<StateType>(-1));
    } else {
        rowToDTStratActionMapping[action] = dtstrat::ActionRepresentation<StateType>(-1);
    }
}

template<typename StateType, typename ValueType>
std::vector<dtstrat::ActionRepresentation<StateType>> ExplorationInformationDT<StateType, ValueType>::getRowToDTStratActionMapping() const {
    return rowToDTStratActionMapping;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::setDTStratIndex(const uint& action, storm::generator::Choice<ValueType, StateType> const& choice) {
    if (rowToDTStratActionMapping.size()<action) {
        rowToDTStratActionMapping.resize(action);
    }
    if (rowToDTStratActionMapping.size() == action) {
        rowToDTStratActionMapping.push_back(dtstrat::ActionRepresentation<StateType>(getDTStratIndex(choice)));
    } else {
        rowToDTStratActionMapping[action] = dtstrat::ActionRepresentation<StateType>(getDTStratIndex(choice));
    }
}


template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::setInternalLabel(const uint& action, std::string const& label) {
    if (actionLabelling.size()<action) {
        actionLabelling.resize(action);
    }
    if (actionLabelling.size() == action) {
        actionLabelling.push_back(label);
    } else {
        actionLabelling[action] = label;
    }
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::setLabel(const uint& action, storm::generator::Choice<ValueType, StateType> const& choice) {
    storm::storage::FlatSet<uint_fast64_t> currentCommandSet = choice.getOriginData().empty()
                                                                   ? storm::storage::FlatSet<uint_fast64_t>()
                                                                   : boost::any_cast<storm::storage::FlatSet<uint_fast64_t>>(std::move(choice.getOriginData()));
    std::vector<std::string> modules;
    std::vector<int> cmdIndices;
    std::string label;
    STORM_LOG_ASSERT(!currentCommandSet.empty(), "The choice has no commands!");
    std::set<std::string> labels;
    for (auto const& entry : currentCommandSet) {
        auto moduleAndCommand = model->getModuleCommandIndexByGlobalCommandIndex(entry);
        auto command = model->getModule(moduleAndCommand.first).getCommand(moduleAndCommand.second);
        label = command.getActionName();
        labels.insert(label);
    }
    STORM_LOG_ASSERT(labels.size()==1,"There are too many or too few labels for this action!");
    setInternalLabel(action, *(labels.begin()));
}

template<typename StateType, typename ValueType>
std::string ExplorationInformationDT<StateType, ValueType>::getLabel(uint const& actionIndex) const {
    STORM_LOG_ASSERT(actionIndex<actionLabelling.size(),"The action is out of scope for the action labelling!");
    return actionLabelling[actionIndex];
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::setDTStratIndex(const uint& action, std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> pathsPerStates) {
    if (rowToDTStratActionMapping.size()<action) {
        rowToDTStratActionMapping.resize(action);
    }
    if (rowToDTStratActionMapping.size() == action) {
        rowToDTStratActionMapping.push_back(dtstrat::ActionRepresentation<StateType>(pathsPerStates));
    } else {
        rowToDTStratActionMapping[action] = dtstrat::ActionRepresentation<StateType>(pathsPerStates);
    }
}

template<typename StateType, typename ValueType>
int ExplorationInformationDT<StateType, ValueType>::getDTStratIndex(storm::generator::Choice<ValueType, StateType> const& choice) {
        storm::storage::FlatSet<uint_fast64_t> currentCommandSet = choice.getOriginData().empty()
                                                                       ? storm::storage::FlatSet<uint_fast64_t>()
                                                                       : boost::any_cast<storm::storage::FlatSet<uint_fast64_t>>(
                                                                             std::move(choice.getOriginData()));
        std::vector<std::string> modules;
        std::vector<int> cmdIndices;
        std::string label;
        if (!currentCommandSet.empty()) {
            for (auto const &entry: currentCommandSet) {
                auto moduleAndCommand = model->getModuleCommandIndexByGlobalCommandIndex(entry);
                auto module = model->getModule(moduleAndCommand.first);
                std::string moduleName = module.getName();
                modules.push_back(moduleName);
                cmdIndices.push_back(moduleAndCommand.second);

                auto command = model->getModule(moduleAndCommand.first).getCommand(moduleAndCommand.second);
                label = command.getActionName();
            }
            dtstrat::Action newAction(modules, cmdIndices, label);
            return findActionInList(newAction);
        }
        STORM_PRINT("Action not found");
        return -1;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::terminateCurrentRowGroup() {
    rowGroupIndices.push_back(matrix.size());
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::moveActionToBackOfMatrix(ActionType const& action) {
    matrix.emplace_back(std::move(matrix[action]));
}

template<typename StateType, typename ValueType>
StateType ExplorationInformationDT<StateType, ValueType>::getActionCount() const {
    return matrix.size();
}

template<typename StateType, typename ValueType>
std::size_t ExplorationInformationDT<StateType, ValueType>::getNumberOfUnexploredStates() const {
    return unexploredStates.size();
}

template<typename StateType, typename ValueType>
std::size_t ExplorationInformationDT<StateType, ValueType>::getNumberOfDiscoveredStates() const {
    return stateToRowGroupMapping.size();
}

template<typename StateType, typename ValueType>
StateType const& ExplorationInformationDT<StateType, ValueType>::getRowGroup(StateType const& state) const {
    return stateToRowGroupMapping[state];
}

template<typename StateType, typename ValueType>
StateType const& ExplorationInformationDT<StateType, ValueType>::getUnexploredMarker() const {
    return unexploredMarker;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::isUnexplored(StateType const& state) const {
    return stateToRowGroupMapping[state] == unexploredMarker;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::isTerminal(StateType const& state) const {
    return terminalStates.find(state) != terminalStates.end();
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::isTarget(StateType const& state) const {
    return targetStates.find(state) != targetStates.end();
}

template<typename StateType, typename ValueType>
typename ExplorationInformationDT<StateType, ValueType>::ActionType const& ExplorationInformationDT<StateType, ValueType>::getStartRowOfGroup(
    StateType const& group) const {
    return rowGroupIndices[group];
}

template<typename StateType, typename ValueType>
std::size_t ExplorationInformationDT<StateType, ValueType>::getRowGroupSize(StateType const& group) const {
    return rowGroupIndices[group + 1] - rowGroupIndices[group];
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::onlyOneActionAvailable(StateType const& group) const {
    return getRowGroupSize(group) == 1;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::addTerminalState(StateType const& state) {
    terminalStates.insert(state);
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::addTargetState(StateType const& state) {
    targetStates.insert(state);
}

template<typename StateType, typename ValueType>
std::vector<storm::storage::MatrixEntry<StateType, ValueType>>& ExplorationInformationDT<StateType, ValueType>::getRowOfMatrix(ActionType const& row) {
    return matrix[row];
}

template<typename StateType, typename ValueType>
std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& ExplorationInformationDT<StateType, ValueType>::getRowOfMatrix(
    ActionType const& row) const {
    return matrix[row];
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::addActionsToMatrix(std::size_t const& count) {
    matrix.resize(matrix.size() + count);
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::maximize() const {
    return optimizationDirection == storm::OptimizationDirection::Maximize;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::minimize() const {
    return !maximize();
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::performPrecomputationExcessiveExplorationSteps(
    std::size_t& numberExplorationStepsSinceLastPrecomputation) const {
    bool result = numberExplorationStepsSinceLastPrecomputation > numberOfExplorationStepsUntilPrecomputation;
    if (result) {
        numberExplorationStepsSinceLastPrecomputation = 0;
    }
    return result;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::performPrecomputationExcessiveSampledPaths(std::size_t& numberOfSampledPathsSinceLastPrecomputation) const {
    if (!numberOfSampledPathsUntilPrecomputation) {
        return false;
    } else {
        bool result = numberOfSampledPathsSinceLastPrecomputation > numberOfSampledPathsUntilPrecomputation.get();
        if (result) {
            numberOfSampledPathsSinceLastPrecomputation = 0;
        }
        return result;
    }
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::useLocalPrecomputation() const {
    return localPrecomputation;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::useGlobalPrecomputation() const {
    return !useLocalPrecomputation();
}

template<typename StateType, typename ValueType>
storm::settings::modules::ExplorationSettings::NextStateHeuristic const& ExplorationInformationDT<StateType, ValueType>::getNextStateHeuristic() const {
    return nextStateHeuristic;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::useDifferenceProbabilitySumHeuristic() const {
    return nextStateHeuristic == storm::settings::modules::ExplorationSettings::NextStateHeuristic::DifferenceProbabilitySum;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::useProbabilityHeuristic() const {
    return nextStateHeuristic == storm::settings::modules::ExplorationSettings::NextStateHeuristic::Probability;
}

template<typename StateType, typename ValueType>
bool ExplorationInformationDT<StateType, ValueType>::useUniformHeuristic() const {
    return nextStateHeuristic == storm::settings::modules::ExplorationSettings::NextStateHeuristic::Uniform;
}

template<typename StateType, typename ValueType>
storm::OptimizationDirection const& ExplorationInformationDT<StateType, ValueType>::getOptimizationDirection() const {
    return optimizationDirection;
}

template<typename StateType, typename ValueType>
void ExplorationInformationDT<StateType, ValueType>::setOptimizationDirection(storm::OptimizationDirection const& direction) {
    optimizationDirection = direction;
}

template<typename StateType, typename ValueType>
std::vector<std::pair<StateType, int>> ExplorationInformationDT<StateType, ValueType>::getPath(StateType& state, ActionType leavingActionST,
                                                                                             StateType& correspState,
                                                                                             std::vector<std::vector<std::pair<StateType, int>>> predecessors,
                                                                                             std::vector<StateType>& states) const {
    // Given a state, we want to find a path to the correspState of the leaving action
    std::vector<std::pair<StateType, int>> actions;
    if (correspState == state) {
        // if the state is already the one where we can directly take the leaving action, we're done
    } else
    {
        // we have to find a path to the leaving state
        std::vector<std::pair<StateType, int>> seen;
        std::list<std::pair<StateType, int>> stack;
        stack.push_back(std::make_pair(correspState, leavingActionST));
        bool found = false;
        while (!found && stack.size() > 0) {
            // get a state from the stack
            std::pair<StateType, int> currentStateAndActionSt = stack.front();
            StateType currentState = currentStateAndActionSt.first;
            stack.pop_front();
            if (currentState == state) {
                // we have backtracked from the corresponding state back to the state from which we're searching for a path
                found = true;
                // we have to find the path now (state and action pairs)
                std::list<std::pair<StateType, int>> newStack;
                newStack.push_back(currentStateAndActionSt);
                while (!newStack.empty()) {
                    std::pair<StateType, int> traceBack = newStack.front();
                    if (traceBack.first == correspState) break;
                    newStack.pop_front();
                    // for each step in the path, add a pair (state, DT-action) to the list
                    dtstrat::ActionRepresentation<StateType> actionRep = rowToDTStratActionMapping[traceBack.second];
                    if (!actionRep.isOne()) {
                        // we access a collapsed MEC (apparently the MEC was not a MAXIMAL EC before...)
                        // get the whole path in the original MEC
                        // vector of pairs
                        // each pair is a state and another vector. The state tells you from which state the action was called
                        // the vector contains again pairs. Each pair consists of a state and an corresponding action.
                        // The vector contains the path from the state to the leaving action of the mec.
                        std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> liste = actionRep.getPathForStates();
                        auto it2 = std::find_if(liste.begin(), liste.end(),
                                                [&](const auto &pair) {
                                                    return pair.first == traceBack.first;
                                                });
                        if (it2 != liste.end()) {
                            int listIndex = it2 - liste.begin();
                            std::vector<std::pair<StateType, int>> path = liste[listIndex].second;
                            for (auto paar: path) {
                                actions.push_back(paar);
                            }
                        }

                    } else {
                        actions.push_back(std::make_pair(traceBack.first, actionRep.getAction()));
                    }
                    auto row = getRowOfMatrix(traceBack.second);
                    for (auto entry: row) {
                        StateType nextState = entry.getColumn();
                        for (uint i = 0; i < seen.size(); i++) {
                            if (seen[i].first == nextState) {
                                newStack.push_back(seen[i]);
                            }
                        }
                    }
                }
                break;
            }
            seen.push_back(currentStateAndActionSt);
            int indxState;
            auto it = find(states.begin(), states.end(), currentState);
            if (it != states.end()) indxState = it - states.begin();
            else {
                continue;
            }
            // get all the predecessors of the current state
            // predecessors contain pairs of states and (storm-) actions
            std::vector<std::pair<StateType, int>> preds = predecessors[indxState];
            for (auto pandA: preds) {
                StateType X = pandA.first;
                if (!std::any_of(seen.begin(), seen.end(),
                                 [&X](const std::pair<StateType, int> &p) { return p.first == X; })
                    && !std::any_of(stack.begin(), stack.end(),
                                    [&X](const std::pair<StateType, int> &p) { return p.first == X; })) {
                    // only visit the state if we have not seen it before
                    stack.push_back(pandA);
                }
            }
        }
    }
    ActionRep dtstratAction = rowToDTStratActionMapping[leavingActionST];
    if (dtstratAction.isOne()) {
        actions.push_back(std::make_pair(correspState, dtstratAction.getAction()));
    } else {
        actions.push_back(std::make_pair(correspState, dtstratAction.getActionOfState(correspState)));
    }
    return actions;
}

template<typename StateType, typename ValueType>
dtstrat::Action ExplorationInformationDT<StateType, ValueType>::getAction(const uint& action) {
    STORM_LOG_ASSERT(action<allActions.size(), "Action is not in action list!");
    return allActions[action];
}

template<typename StateType, typename ValueType>
int ExplorationInformationDT<StateType, ValueType>::getDTStratIndex(const uint& stormAction, StateType const& currentStateId) {
    STORM_LOG_ASSERT(stormAction<rowToDTStratActionMapping.size(), "Action is not in storm-actions!");
    ActionRep action = rowToDTStratActionMapping[stormAction];
    if (action.isOne()) {
        return action.getAction();
    }
    else {
        return action.getActionOfState(currentStateId);
    }
}

template<typename StateType, typename ValueType>
std::vector<storm::prism::Module> ExplorationInformationDT<StateType, ValueType>::getAllModules() const {
    return model->getModules();
}

template<typename StateType, typename ValueType>
std::vector<std::pair<StateType, ValueType>> ExplorationInformationDT<StateType, ValueType>::getSuccessorsOfAction(const uint& chosenAction) {
    std::vector<std::pair<StateType, ValueType>> result;
    std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& row = getRowOfMatrix(chosenAction);
    for (auto element: row) {
        result.push_back(std::make_pair(element.getColumn(), element.getValue()));
    }
    return result;
}

template<typename StateType, typename ValueType>
uint ExplorationInformationDT<StateType, ValueType>::getNumberOfModules() const {
    return getAllModules().size();
}

template<typename StateType, typename ValueType>
std::string ExplorationInformationDT<StateType, ValueType>::getModelName() const{
    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    auto prismName = ioSettings.getPrismInputFilename();
    std::size_t found = prismName.find_last_of('.');
    return prismName.substr(0,found);
}

template<typename StateType, typename ValueType>
std::vector<uint> ExplorationInformationDT<StateType, ValueType>::getStormActionFromDTStrat(std::vector<uint> const& possibleActionIds, StateType const& currentStateId) const {
    std::vector<uint> result;
    for (auto action : possibleActionIds) {
        auto startRowOfGroup = getStartRowOfGroup(getRowGroup(currentStateId));
        auto nextRowGroupStart = getStartRowOfGroup(getRowGroup(currentStateId)+1);
        for (uint stormAction=startRowOfGroup; stormAction<nextRowGroupStart; stormAction++) {
            ActionRep stormActionRep = rowToDTStratActionMapping[stormAction];
            if ((stormActionRep.isOne() && stormActionRep.getAction() == action) || (!stormActionRep.isOne() && stormActionRep.getActionOfState(currentStateId) == action)) {
                result.push_back(stormAction);
            }
        }
    }
    return result;
}

template<typename StateType, typename ValueType>
std::vector<uint> ExplorationInformationDT<StateType, ValueType>::getDTStratActionFromModuleAction(std::string const &module, std::vector<uint> const &actions) {
    std::vector<uint> result;
    for (uint const &action: actions) {
        if (moduleActionToDTStrat.find(module) != moduleActionToDTStrat.end()) {
            std::map<uint, std::vector<uint>>& cmdIndicesActionMapping = moduleActionToDTStrat.at(module);
            if (cmdIndicesActionMapping.find(action) != cmdIndicesActionMapping.end()) {
                auto possibleCmdIndices = cmdIndicesActionMapping.at(action);
                result.insert(result.end(), possibleCmdIndices.begin(), possibleCmdIndices.end());
            } else {
                //STORM_LOG_THROW(false, exceptions::UnexpectedException,"This module '" + module + "' has no action " + std::to_string(action) + "!");
            }
        } else {
            STORM_LOG_THROW(false,exceptions::UnexpectedException,"This module '" + module + "' does not exist!");
        }
    }
    return result;
}

template<typename StateType, typename ValueType>
std::map<std::string, std::map<uint, std::vector<uint>>>
ExplorationInformationDT<StateType, ValueType>::getModuleActionToDTStratMapping() {
    std::map<std::string, std::map<uint, std::vector<uint>>> moduleActionToDTStratMapping;
    // iterate over all actions
    for (uint j = 0; j < allActions.size(); j++) {
        auto const &action = allActions[j];
        auto modules = action.getModules();
        auto cmdIndices = action.getCommandIndices();

        // iterate over modules in which the action is performed
        for (uint i = 0; i < modules.size(); i++) {
            auto module = modules[i];
            auto cmdIndex = cmdIndices[i];
            // enter the cmdIndex for this module in the mapping
            if (moduleActionToDTStratMapping.find(module) != moduleActionToDTStratMapping.end()) {
                std::map<uint, std::vector<uint>> &cmdIndicesActionMapping = moduleActionToDTStratMapping[module];
                if (cmdIndicesActionMapping.find(cmdIndex) != cmdIndicesActionMapping.end()) {
                    cmdIndicesActionMapping[cmdIndex].push_back(j);
                } else {
                    std::vector<uint> actionsForCmdIndex = {j};
                    cmdIndicesActionMapping[cmdIndex] = actionsForCmdIndex;
                }
            } else {
                std::map<uint, std::vector<uint>> cmdIndicesActionMapping;
                std::vector<uint> actionsForCmdIndex = {j};
                cmdIndicesActionMapping[cmdIndex] = actionsForCmdIndex;
                moduleActionToDTStratMapping[module] = cmdIndicesActionMapping;
            }
        }
    }
    return moduleActionToDTStratMapping;
}

template<typename StateType, typename ValueType>
uint ExplorationInformationDT<StateType, ValueType>::getNumberOfActions() const {
    return allActions.size();
}

template class ExplorationInformationDT<uint32_t, double>;
}  // namespace exploration_detail
}  // namespace modelchecker
}  // namespace storm
