#include "storm/modelchecker/explorationDT/BoundsDT.h"

#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"

namespace storm {
namespace modelchecker {
namespace exploration_detail_dt {

template<typename StateType, typename ValueType>
std::pair<ValueType, ValueType> BoundsDT<StateType, ValueType>::getBoundsForState(
    StateType const& state, ExplorationInformationDT<StateType, ValueType> const& explorationInformation) const {
    ActionType index = explorationInformation.getRowGroup(state);
    if (index == explorationInformation.getUnexploredMarker()) {
        return std::make_pair(storm::utility::zero<ValueType>(), storm::utility::one<ValueType>());
    } else {
        return boundsPerState[index];
    }
}

template<typename StateType, typename ValueType>
ValueType BoundsDT<StateType, ValueType>::getLowerBoundForState(StateType const& state,
                                                              ExplorationInformationDT<StateType, ValueType> const& explorationInformation) const {
    ActionType index = explorationInformation.getRowGroup(state);
    if (index == explorationInformation.getUnexploredMarker()) {
        return storm::utility::zero<ValueType>();
    } else {
        return getLowerBoundForRowGroup(index);
    }
}

template<typename StateType, typename ValueType>
ValueType const& BoundsDT<StateType, ValueType>::getLowerBoundForRowGroup(StateType const& rowGroup) const {
    return boundsPerState[rowGroup].first;
}

template<typename StateType, typename ValueType>
ValueType BoundsDT<StateType, ValueType>::getUpperBoundForState(StateType const& state,
                                                              ExplorationInformationDT<StateType, ValueType> const& explorationInformation) const {
    ActionType index = explorationInformation.getRowGroup(state);
    if (index == explorationInformation.getUnexploredMarker()) {
        return storm::utility::one<ValueType>();
    } else {
        return getUpperBoundForRowGroup(index);
    }
}

template<typename StateType, typename ValueType>
ValueType const& BoundsDT<StateType, ValueType>::getUpperBoundForRowGroup(StateType const& rowGroup) const {
    return boundsPerState[rowGroup].second;
}

template<typename StateType, typename ValueType>
std::pair<ValueType, ValueType> const& BoundsDT<StateType, ValueType>::getBoundsForAction(ActionType const& action) const {
    return boundsPerAction[action];
}

template<typename StateType, typename ValueType>
ValueType const& BoundsDT<StateType, ValueType>::getLowerBoundForAction(ActionType const& action) const {
    return boundsPerAction[action].first;
}

template<typename StateType, typename ValueType>
ValueType const& BoundsDT<StateType, ValueType>::getUpperBoundForAction(ActionType const& action) const {
    return boundsPerAction[action].second;
}

template<typename StateType, typename ValueType>
ValueType const& BoundsDT<StateType, ValueType>::getBoundForAction(storm::OptimizationDirection const& direction, ActionType const& action) const {
    if (direction == storm::OptimizationDirection::Maximize) {
        return getUpperBoundForAction(action);
    } else {
        return getLowerBoundForAction(action);
    }
}

template<typename StateType, typename ValueType>
ValueType BoundsDT<StateType, ValueType>::getDifferenceOfStateBounds(StateType const& state,
                                                                   ExplorationInformationDT<StateType, ValueType> const& explorationInformation) const {
    std::pair<ValueType, ValueType> bounds = getBoundsForState(state, explorationInformation);
    return bounds.second - bounds.first;
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::initializeBoundsForNextState(std::pair<ValueType, ValueType> const& vals) {
    boundsPerState.push_back(vals);
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::initializeBoundsForNextAction(std::pair<ValueType, ValueType> const& vals) {
    boundsPerAction.push_back(vals);
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setLowerBoundForState(StateType const& state, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                         ValueType const& value) {
    setLowerBoundForRowGroup(explorationInformation.getRowGroup(state), value);
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setLowerBoundForRowGroup(StateType const& group, ValueType const& value) {
    boundsPerState[group].first = value;
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setUpperBoundForState(StateType const& state, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                         ValueType const& value) {
    setUpperBoundForRowGroup(explorationInformation.getRowGroup(state), value);
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setUpperBoundForRowGroup(StateType const& group, ValueType const& value) {
    boundsPerState[group].second = value;
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setBoundsForAction(ActionType const& action, std::pair<ValueType, ValueType> const& values) {
    boundsPerAction[action] = values;
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setBoundsForState(StateType const& state, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                     std::pair<ValueType, ValueType> const& values) {
    StateType const& rowGroup = explorationInformation.getRowGroup(state);
    setBoundsForRowGroup(rowGroup, values);
}

template<typename StateType, typename ValueType>
void BoundsDT<StateType, ValueType>::setBoundsForRowGroup(StateType const& rowGroup, std::pair<ValueType, ValueType> const& values) {
    boundsPerState[rowGroup] = values;
}

template<typename StateType, typename ValueType>
bool BoundsDT<StateType, ValueType>::setLowerBoundOfStateIfGreaterThanOld(StateType const& state,
                                                                        ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                                        ValueType const& newLowerValue) {
    StateType const& rowGroup = explorationInformation.getRowGroup(state);
    if (boundsPerState[rowGroup].first < newLowerValue) {
        boundsPerState[rowGroup].first = newLowerValue;
        return true;
    }
    return false;
}

template<typename StateType, typename ValueType>
bool BoundsDT<StateType, ValueType>::setUpperBoundOfStateIfLessThanOld(StateType const& state,
                                                                     ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                                     ValueType const& newUpperValue) {
    StateType const& rowGroup = explorationInformation.getRowGroup(state);
    if (newUpperValue < boundsPerState[rowGroup].second) {
        boundsPerState[rowGroup].second = newUpperValue;
        return true;
    }
    return false;
}

template class BoundsDT<uint32_t, double>;

}  // namespace exploration_detail
}  // namespace modelchecker
}  // namespace storm
