//
// Created by steffi on 04.01.22.
//
#ifndef STORM_ACTIONDT_H
#define STORM_ACTIONDT_H

#include <vector>
#include <string>

#include "storm/storage/prism/Program.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {
class Action {
   public:
    Action(std::vector<std::string> mods, std::vector<int> cmdIndices, std::string lab);

    bool operator==(Action const& one) const;

    uint size() const;

    std::string getLabel() const;

    std::vector<std::string> getModules() const;

    std::vector<int> getCommandIndices() const;

    void appendModule(std::string newModule);

    void appendCommand(int newCommand);

   private:
    // An action is characterized by all modules that it appears in:
    std::vector<std::string> modules;
    // and the corresponding command indices from storm (basically "lines" in the prism-model)
    std::vector<int> commandIndices;
    // and the label (e.g. 'reset')
    std::string label;
};

/*!
 * Class that represents either one action (simple case) or a set of actions for different states
 * This is needed because SCC are collapsed and we have to choose an action that is not actually possible in the original model
 * This is reflected in the "set-of-actions"-case: each state in the SCC has an action that finally leads to the "one" action that STORM has
 * @tparam StateType
 */
template<typename StateType>
class ActionRepresentation {
   private:
    // if we look at an SCC, there are actions for all states:
    std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> pathsForStates;
    // if there's just one action:
    uint action;
   public:
    ActionRepresentation() = default;
    ActionRepresentation(uint action) : action(action) {};
    ActionRepresentation(std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> vec) : pathsForStates(vec) {};

    ActionRepresentation& operator=(const ActionRepresentation<StateType>& other);
    /*!
     * checks if the representation represents one action (true) or a set for an SCC
     * @return
     */
    bool isOne() const;

    /*!
     * returns the action (if set)
     * @return dtstrat-id of the action
     */
    uint getAction() const;

    /*!
     * returns the action for a set of states in the SCC (if set)
     * @return vector of pairs (State + Vector (of pairs (state + action))
     */
    std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> getPathForStates() const;

    /*!
     * returns the action of the state if it is contained in the SCC of this action
     * @param state id of the state of interest
     * @return dtstrat-id of the action
     */
    uint getActionOfState(StateType const& state) const;

    /*!
     * returns the set of actions that represent this action starting from a state
     * @param state id of the state of interest
     * @return vector of pairs (state + action)
     */
    std::vector<std::pair<StateType, int>> getPathForActionOfSTate(StateType const& state) const;
};

// This function reads out the model and returns all possible actions
void getAllActionsFromModel(std::vector<Action> &allActions, storm::prism::Program &model);
}
}
}


#endif  // STORM_ACTIONDT_H
