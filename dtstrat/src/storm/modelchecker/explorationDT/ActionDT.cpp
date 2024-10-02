//
// Created by steffi on 04.01.22.
//
#include "storm/modelchecker/explorationDT/ActionDT.h"


namespace storm {
namespace modelchecker {
namespace dtstrat {

Action::Action(std::vector<std::string> mods, std::vector<int> cmdIndices, std::string lab) {
    modules = mods;
    commandIndices = cmdIndices;
    label = lab;
}

bool Action::operator==(Action const& one) const {
    bool ret = true;
    if (one.size() != this->size()) {return false;}
    if (one.getLabel() != this->getLabel()) {return false;}
    std::vector<std::string> oneModules = one.getModules();
    std::vector<int> oneCmdIndices = one.getCommandIndices();
    for (uint i = 0; i<modules.size(); i++) {
        if (modules[i] != oneModules[i]) {
            ret = false;
            break;
        }
        if (commandIndices[i] != oneCmdIndices[i]) {
            ret = false;
            break;
        }
    }
    return ret;
}

uint Action::size() const {
    return modules.size();
}

std::string Action::getLabel() const {
    return label;
}

std::vector<std::string> Action::getModules() const {
    return modules;
}

std::vector<int> Action::getCommandIndices() const {
    return commandIndices;
}

void Action::appendModule(std::string newModule) {
    modules.push_back(newModule);
    return;
}

void Action::appendCommand(int newCommand) {
    commandIndices.push_back(newCommand);
    return;
}

template<typename StateType>
bool ActionRepresentation<StateType>::isOne() const {
    if (!pathsForStates.empty()) return false;
    else return true;
};

template<typename StateType>
uint ActionRepresentation<StateType>::getAction() const {
    STORM_LOG_ASSERT(pathsForStates.empty(),"This action is not a one-action (but an SCC-action).");
    return action;
};

template<typename StateType>
std::vector<std::pair<StateType, std::vector<std::pair<StateType, int>>>> ActionRepresentation<StateType>::getPathForStates() const {
    STORM_LOG_ASSERT(!pathsForStates.empty(), "This action is a one-action, not an SCC-action.");
    return pathsForStates;
}

template<typename StateType>
ActionRepresentation<StateType>& ActionRepresentation<StateType>::operator=(const ActionRepresentation<StateType>& other) {
    if (other.isOne()) {
        this->action = other.getAction();
    } else {
        this->pathsForStates.clear();
        auto otherPFS = other.getPathForStates();
        std::copy(otherPFS.begin(), otherPFS.end(),this->pathsForStates.begin());
    }
    return *this;
}

// This function reads out the model and returns all possible actions
void getAllActionsFromModel(std::vector<Action> &allActions, storm::prism::Program &model) {
    std::vector<storm::prism::Module> allModules = model.getModules();
    //  We iterate over all modules and all commands (i.e. actions) in the model
    for (storm::prism::Module module : allModules) {
        std::map<std::string, std::vector<int>> idicti;
        std::vector<storm::prism::Command> commands = module.getCommands();
        int countCommands = 0;
        for (storm::prism::Command command : commands) {
            if (command.isLabeled()) {
                // Labelled commands may be used in different modules.
                // They are thus just put into a dictionary.
                // idicti consists of lists for each action-name, where the list contains all commands that
                // correspond to the action.
                std::string name_of_command = command.getActionName();
                if (idicti.find(name_of_command) == idicti.end()) {
                    std::vector<int> newV;
                    newV.push_back(countCommands);
                    idicti[name_of_command] = newV;
                } else {
                    idicti[name_of_command].push_back(countCommands);
                }
            } else {
                // Unlabeled commands are local commands, so they can be added to the action-list without any
                // further complication
                std::vector<std::string> modules;
                modules.push_back(module.getName());
                std::vector<int> cmdIndices;
                cmdIndices.push_back(countCommands);
                Action newAction(modules, cmdIndices, "");
                allActions.push_back(newAction);
            }
            countCommands++;
        }

        // After iterating over all commands, we still need to evluate the labelled actions.
        // Each combination of labelled actions in different modules should get its own index.
        std::map<std::string, std::vector<int>>::iterator it;
        for (it = idicti.begin(); it != idicti.end(); it++) {
            bool found = false;
            std::vector<Action> newActions;
            for (Action &action : allActions) {
                if (action.getLabel() == it->first) {
                    // If the action-name already exists, we need a copy of the already existing action
                    // for each new command that we add with the same action-name
                    // e.g. we have action "send" that consitst of: [('module one',2),('module two',3)]
                    //      we want two other possible commands for module three namely command 4, and 6
                    //      so, we need two actions in total [('module one',2),('module two',3),('module three',4)]
                    //                                       [('module one',2),('module two',3),('module three',6)]
                    found = true;
                    if ((it->second).size() > 1) {
                        for (uint i = 1; i < (it->second).size(); i++) {
                            std::vector<std::string> mods(action.getModules());
                            std::vector<int> cmds(action.getCommandIndices());
                            mods.push_back(module.getName());
                            cmds.push_back((it->second)[i]);
                            Action newAction(mods, cmds, action.getLabel());
                            newActions.push_back(newAction);
                        }
                    }
                    action.appendModule(module.getName());
                    action.appendCommand((it->second)[0]);
                }
            }
            if (!found) {
                for (int cmdIndx : it->second) {
                    std::vector<std::string> mods;
                    std::vector<int> cmdIdc;
                    mods.push_back(module.getName());
                    cmdIdc.push_back(cmdIndx);
                    Action newAction(mods, cmdIdc, it->first);
                    allActions.push_back(newAction);
                }
            } else if (newActions.size() > 0) {
                for (auto newAction : newActions) {
                    allActions.push_back(newAction);
                }
            }
        }
    }
    return;
}

template<typename StateType>
uint ActionRepresentation<StateType>::getActionOfState(const StateType &state) const {
    STORM_LOG_ASSERT(!pathsForStates.empty(), "This is not an SCC-action.");
    auto pathForThisState = getPathForActionOfSTate(state);
    if (!pathForThisState.empty()) {
        return pathForThisState[0].second;
    }
    else return -1;
}

template<typename StateType>
std::vector<std::pair<StateType, int>> ActionRepresentation<StateType>::getPathForActionOfSTate(StateType const& state) const {
    STORM_LOG_ASSERT(!pathsForStates.empty(), "This is not an SCC-action.");
    auto it = std::find_if(pathsForStates.begin(), pathsForStates.end(), [state] (const std::pair<StateType, std::vector<std::pair<StateType, int>>>& a) { return a.first == state;});
    if (!(it==pathsForStates.end())) {
        return it->second;
    } else {
        return {};
    }
}

template class ActionRepresentation<uint32_t>;
}
}
}