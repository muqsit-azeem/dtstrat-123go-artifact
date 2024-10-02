#ifndef STORM_MODELCHECKER_EXPLORATIONDT_EXPLORATION_DETAIL_EXPLORATIONINFORMATION_H_
#define STORM_MODELCHECKER_EXPLORATIONDT_EXPLORATION_DETAIL_EXPLORATIONINFORMATION_H_

#include <limits>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "storm/solver/OptimizationDirection.h"

#include "storm/generator/CompressedState.h"

#include "storm/storage/BoostTypes.h"
#include "storm/storage/SparseMatrix.h"

#include "storm/settings/modules/ExplorationSettings.h"
#include "storm/modelchecker/explorationDT/ActionDT.h"
#include "storm/generator/Choice.h"
#include "storm/modelchecker/explorationDT/StateGenerationDT.h"

namespace storm {
namespace modelchecker {
namespace exploration_detail_dt {

template<typename StateType, typename ValueType>
class ExplorationInformationDT {
    typedef dtstrat::ActionRepresentation<StateType> ActionRep;
   public:
    typedef StateType ActionType;
    typedef storm::storage::FlatSet<StateType> StateSet;
    typedef std::unordered_map<StateType, storm::generator::CompressedState> IdToStateMap;
    typedef typename IdToStateMap::const_iterator const_iterator;
    typedef std::vector<std::vector<storm::storage::MatrixEntry<StateType, ValueType>>> MatrixType;

    ExplorationInformationDT(storm::OptimizationDirection const& direction, storm::prism::Program& model,
                           ActionType const& unexploredMarker = std::numeric_limits<ActionType>::max());

    void getAllActions();

    const_iterator findUnexploredState(StateType const& state) const;

    const_iterator unexploredStatesEnd() const;

    void removeUnexploredState(const_iterator it);

    void addUnexploredState(StateType const& stateId, storm::generator::CompressedState const& compressedState);

    void assignStateToRowGroup(StateType const& state, ActionType const& rowGroup);

    StateType assignStateToNextRowGroup(StateType const& state);

    StateType getNextRowGroup() const;

    void newRowGroup(ActionType const& action);

    void newRowGroup();

    void terminateCurrentRowGroup();

    void moveActionToBackOfMatrix(ActionType const& action);

    StateType getActionCount() const;

    std::size_t getNumberOfUnexploredStates() const;

    std::size_t getNumberOfDiscoveredStates() const;

    StateType const& getRowGroup(StateType const& state) const;

    StateType const& getUnexploredMarker() const;

    bool isUnexplored(StateType const& state) const;

    bool isTerminal(StateType const& state) const;

    bool isTarget(StateType const& state) const;

    ActionType const& getStartRowOfGroup(StateType const& group) const;

    std::size_t getRowGroupSize(StateType const& group) const;

    bool onlyOneActionAvailable(StateType const& group) const;

    void addTerminalState(StateType const& state);

    void addTargetState(StateType const& state);

    std::vector<storm::storage::MatrixEntry<StateType, ValueType>>& getRowOfMatrix(ActionType const& row);

    std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& getRowOfMatrix(ActionType const& row) const;

    void addActionsToMatrix(std::size_t const& count);

    bool maximize() const;

    bool minimize() const;

    bool performPrecomputationExcessiveExplorationSteps(std::size_t& numberExplorationStepsSinceLastPrecomputation) const;

    bool performPrecomputationExcessiveSampledPaths(std::size_t& numberOfSampledPathsSinceLastPrecomputation) const;

    bool useLocalPrecomputation() const;

    bool useGlobalPrecomputation() const;

    storm::settings::modules::ExplorationSettings::NextStateHeuristic const& getNextStateHeuristic() const;

    bool useDifferenceProbabilitySumHeuristic() const;

    bool useProbabilityHeuristic() const;

    bool useUniformHeuristic() const;

    storm::OptimizationDirection const& getOptimizationDirection() const;

    void setOptimizationDirection(storm::OptimizationDirection const& direction);

    int getDTStratIndex(storm::generator::Choice<ValueType, StateType> const& choice);
    int getDTStratIndex(uint const& stormAction, StateType const& currentStateId);

    int findActionInList(dtstrat::Action const &searchAction);

    void setDTStratIndex(uint const& action);
    void setDTStratIndex(uint const& action, storm::generator::Choice<ValueType, StateType> const& choice);
    void setLabel(uint const& action, storm::generator::Choice<ValueType, StateType> const& choice);

    void setInternalLabel(uint const& action, std::string const& label);

    void setDTStratIndex(const uint& action, std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> pathsPerStates);

    std::string getLabel(uint const& actionIndex) const;

    std::vector<ActionRep> getRowToDTStratActionMapping() const;

    std::vector<std::pair<StateType, int>> getPath(StateType& state, ActionType leavingAction, StateType& correspState,
                                                   std::vector<std::vector<std::pair<StateType, int>>> predecessors, std::vector<StateType>& states) const;

    dtstrat::Action getAction(uint const& action);

    uint getNumberOfActions() const;

    std::vector<storm::prism::Module> getAllModules() const;
    uint getNumberOfModules() const;

    std::string getModelName() const;

    std::vector<std::pair<StateType, ValueType>> getSuccessorsOfAction(uint const& action);

    std::vector<uint> getStormActionFromDTStrat(std::vector<uint> const& possibleActionIds, StateType const& currentStateId) const;
    std::vector<uint> getDTStratActionFromModuleAction(std::string const &module, std::vector<uint> const &actions);

   private:
    std::map<std::string, std::map<uint, std::vector<uint>>> getModuleActionToDTStratMapping();
    std::vector<std::string> actionLabelling;

    std::vector<dtstrat::Action> allActions;
    MatrixType matrix;
    std::vector<ActionRep> rowToDTStratActionMapping;
    std::map<std::string, std::map<uint, std::vector<uint>>> moduleActionToDTStrat;
    std::vector<StateType> rowGroupIndices;

    std::vector<StateType> stateToRowGroupMapping;
    StateType unexploredMarker;
    IdToStateMap unexploredStates;

    storm::OptimizationDirection optimizationDirection;
    StateSet terminalStates;
    StateSet targetStates;

    bool localPrecomputation;
    std::size_t numberOfExplorationStepsUntilPrecomputation;
    boost::optional<std::size_t> numberOfSampledPathsUntilPrecomputation;

    storm::settings::modules::ExplorationSettings::NextStateHeuristic nextStateHeuristic;

    storm::prism::Program* model;
};
}  // namespace exploration_detail
}  // namespace modelchecker
}  // namespace storm

#endif /* STORM_MODELCHECKER_EXPLORATIONDT_EXPLORATION_DETAIL_EXPLORATIONINFORMATION_H_ */
