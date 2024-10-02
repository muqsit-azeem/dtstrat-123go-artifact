//
// Created by steffi on 07.01.22.
//

#ifndef STORM_BRTDPINFORMATION_H
#define STORM_BRTDPINFORMATION_H

#include <random>

#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"
#include "storm/modelchecker/explorationDT/StateGenerationDT.h"
#include "storm/modelchecker/explorationDT/BoundsDT.h"
#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"
#include "storm/modelchecker/explorationDT/StatisticsDT.h"
#include "storm/storage/MaximalEndComponent.h"
#include <mlpack/core.hpp>
#include <armadillo>

namespace storm{
namespace modelchecker {
namespace exploration_detail_dt {
using namespace dtstrat;

template<typename StateType, typename ValueType>
/*!
 * The main class that gathers all the information and can run simulations
 */
class BRTDPInformation {
   public:
    typedef StateType ActionType;
    typedef std::vector<std::pair<StateType, ActionType>> StateActionStack;

    BRTDPInformation(ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                     StateGenerationDT<StateType, ValueType>& stateGeneration,
                     BoundsDT<StateType, ValueType>& bounds, StatisticsDT<StateType, ValueType>& stats,
                     StateValuationsDT& stateValuations,
                     double const& precision);

    ~BRTDPInformation();

    /*!
     * returns the state generation object
     * @return state generation object
     */
    StateGenerationDT<StateType, ValueType>& getStateGeneration();

    /*!
     * returns the exploration information object
     * @return exploration information object
     */
    ExplorationInformationDT<unsigned int, double>& getExplorationInformation();

    /*!
     * returns the bounds
     * @return bounds
     */
    BoundsDT<StateType, ValueType>& getBounds();

    /*!
     * returns the statistics
     * @return statistics
     */
    StatisticsDT<StateType, ValueType>& getStats();

    /*!
     * returns the state valuations
     * @return state valuations
     */
    StateValuationsDT& getStateValuations();

    /*!
     * returns the constant state generation object
     * @return constant state generation object
     */
    StateGenerationDT<StateType, ValueType>& getStateGeneration() const;

    /*!
     * returns a constant exploration information object
     * @return constant exploration information object
     */
    ExplorationInformationDT<StateType, ValueType>& getExplorationInformation() const;

    /*!
     * returns constant bounds
     * @return constant bounds
     */
    BoundsDT<StateType, ValueType>& getBounds() const;

    /*!
     * returns constant statistics
     * @return constant statistics
     */
    StatisticsDT<StateType, ValueType>& getStats() const;

    /*!
     * returns the constant state valuations
     * @return constant state valuations
     */
     StateValuationsDT& getStateValuations() const;


    /*!
     * Outputs the data for learning a decision tree into a file
     * It reads all states that are in the state valuations and outputs them including the best action as given from the configuration (if set) otherwise
     * from BRTDP.
     * @param ofs the output stream
     * @param config (optional) the configuration, if set, it will be used to output the actions, otherwise BRTDP bounds will be used
     * @param onlyReachableStates determines whether only the reachable states under the current strategy should be printed
     * @param printModules determines whether to output the module of each action
     * @param printStateVariableAssignment determines whether to output the association of each state variable to which module (-1 is global)
     * @param seperator character for the seperation
     */
    void printData(std::ofstream& ofs, dtstrat::LearningConfig<StateType>* config= nullptr, bool const& onlyReachableStates = false, bool const& printModules = false,
                   bool const& printStateVariableAssignment = false, char const& seperator = ',');

    /*!
     * @param currentStateId current state
     * @param config  learning configuration with the heuristic DT
     * @return dtstrat-action id
     */
    uint getDTActionOfState(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config);

    /*!
     * Uses BRTP to determine which action is best in the current state.
     * @param currentStateId current state
     * @param config (optional) learning configuration with a heuristic
     * @return dtstrat-action id
     */
    uint getBestActionOfState(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config=nullptr);

    /*!
     * add the state valuation of "state" to the matrix in "data"
     * @param state storm state index
     * @param data armadillo double matrix
     */
    void addStateValuation(StateType const& state, arma::Mat<int>& data, arma::Row<long unsigned int>& chosenActions, std::pair<uint, uint>& chosenAction, arma::Row<long unsigned int>* modules=nullptr);

    /*!
     * Print the header in the file for later calling dtcontrol.
     * @param ofs the output stream
     * @param printModules whether or not to print the modules
     */
    void printDTcontrolHeader(std::ofstream& ofs, bool const& printModules);

    /*!
     * Print the header in the file for later calling sklearn via a python script.
     * @param ofs the output stream
     */
    void printModuleDTHeader(std::ofstream& ofs);

    /*!
     * Get the module and corresponding action to print into the file.
     * @param chosenAction DTStrat action index
     * @param allModules list of all modules (prism modules)
     * @param seperator character for seperation (optional, default is ,)
     * @return string containing module and action
     */
    std::string getModuleAndAction(uint const& chosenAction, std::vector<storm::prism::Module> const& allModules, bool& warning, char const& seperator=',');

    /*!
     * Get the module and corresponding action
     * @param chosenAction  DTStrat action index
     * @param allModules list of all modules (prism modules)
     * @param warning
     * @return list of pairs (module and action)
     */
    std::vector<std::pair<uint, uint>> getModuleAndActions(uint const& chosenAction, std::vector<storm::prism::Module> const& allModules, bool& warning);

    /*!
     * Generates a string out of a state valuation, to be printed in a file.
     * @param valuation state-valuation
     * @param printStateVariableAssignment whether or not to print the assignment of the state variables to a module
     * @param seperator character for seperation (optional, default is ,)
     * @return
     */
    std::string printValuation(StateValuationDT const& valuation, bool const& printStateVariableAssignment=false, char const& seperator=',');

    /*!
     * Simple forward pass to get the reachable states under a current strategy (given in the config resp. none)
     * @param config if set, the strategy of the config is used to get the reachable states, unless the BRTDP choice is used
     * @return a set of reachable states
     */
    storm::storage::BitVector forwardPass(dtstrat::LearningConfig<StateType>* config= nullptr);

    /*!
     * Forward pass that automatically generates a matrix for the decision tree
     * @param config
     * @return
     */
    arma::Mat<int> forwardPassStateValuations(arma::Row<long unsigned int>& chosenActions, arma::Row<long unsigned int>* modules=nullptr, dtstrat::LearningConfig<StateType>* config=nullptr);

    /*!
     * Retrieves the successors of a state under a given strategy
     * @param currentStateId state id
     * @param config if set, the strategy is used to get the successor, unless the BRTDP choice is used
     * @return
     */
    std::vector<std::vector<std::pair<StateType, ValueType>>> getSuccessors(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config);


    /*!
     * Retrieves the successors of a state under a given strategy
     * @param currentStateId state id
     * @param chosenAction action id
     * @return
     */
    std::vector<StateType> getSuccessorsOfStateDTActionPair(StateType const& currentStateId, ActionType chosenAction);

    /*!
     * Retrieves the number of modules in the model from the exploration information (that contains the model)
     * @return number of modules
     */
    uint getNumberOfModules() const;

    /*!
     * Retrieves the model name (as a string)
     * @return model name
     */
    std::string getModelName() const;

    /*!
     * Retrieves the maximum number of steps that were needed to reach a terminal state, so far
     * @return maximum number of steps to reach a terminal state
     */
    uint getMaxNumberOfSteps() const;

    /*!
     * Given a state, returns the action that BRTDP would suggest (the one with the highest upper bound)
     * @param currentStateId current state index
     * @return storm action-index
     */
    uint getBRTDPActionOfState(StateType const& currentStateId);

    /*!
     * Retrieves all actions of a current state. If the config is set, it returns the actions ordered according to the heuristic of the config, else just all
     * actions as given in BRTDP
     * @param currentStateId current state index
     * @param config optional, learning configuration
     * @return vector of storm-action indices
     */
    std::vector<StateType> sampleActionsVectorOfState(StateType const& currentStateId, dtstrat::LearningConfig<StateType>* config=nullptr);

    /*!
     * Runs one simulation
     * @param stateGeneration State generation to get new states in exploration
     * @param explorationInformation storage of information during exploration
     * @param stack stack of states and actions that are visited in a simulation
     * @param bounds bounds of states and actions
     * @param stats statistics
     * @param config learning heuristic
     * @param stopAfterNSteps optional, default=-1, if set to something bigger than -1, the simulation stops after this number of steps
     * @return whether or not a terminal state is found
     */
    bool samplePathFromInitialState(StateGenerationDT<StateType, ValueType>& stateGeneration,
                                    ExplorationInformationDT<StateType, ValueType>& explorationInformation, StateActionStack& stack,
                                    BoundsDT<StateType, ValueType>& bounds, StatisticsDT<StateType, ValueType>& stats,
                                    LearningConfig<StateType>* config, StateValuationsDT& stateValuations,
                                    int const& stopAfterNSteps = -1) const;

    /*!
     * Explores a state, i.e. loads it in the state generation to get its possible actions and successors, store them in the exploration information, and store
     * the state-valuation
     * @param stateGeneration State generation to get new states in exploration
     * @param currentStateId state index
     * @param currentState compressed state (stored in the unexplored states of the exploration information)
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @param stats statistics
     * @return whether or not the state is a terminal state
     */
    bool exploreState(StateGenerationDT<StateType, ValueType>& stateGeneration, StateType const& currentStateId,
                      storm::generator::CompressedState const& currentState, ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                      BoundsDT<StateType, ValueType>& bounds, StatisticsDT<StateType, ValueType>& stats, StateValuationsDT& stateValuations) const;

    /*!
     * Based on the upper bounds of the actions, sample one of the actions with the highest upper bound
     * @param currentStateId state index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return storm action index
     */
    ActionType sampleActionOfState(StateType const& currentStateId, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                   BoundsDT<StateType, ValueType>& bounds) const;

    /*!
     * Given (a state, implicitly, and) an action, randomly sample a successor, based on their probability
     * @param chosenAction storm action index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return state index of successor
     */
    StateType sampleSuccessorFromAction(ActionType const& chosenAction, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                        BoundsDT<StateType, ValueType> const& bounds) const;

    /*!
     *  Outline:
     * 1. construct a sparse transition matrix of the relevant part of the state space.
     * 2. use this matrix to compute states with probability 0/1 and an MEC decomposition (in the max case).
     * 3. use MEC decomposition to collapse MECs.
     * @param stack state and action stack
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @param stats statistics
     * @return
     */
    bool performPrecomputation(StateActionStack const& stack, ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                               BoundsDT<StateType, ValueType>& bounds, StatisticsDT<StateType, ValueType>& stats) const;

    /*!
     * Collapse a found MEC
     * @param mec MEC
     * @param relevantStates vector of states appearing in the MEC
     * @param relevantStatesMatrix sub-transitionmatrix only for the relevant states
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     */
    void collapseMec(storm::storage::MaximalEndComponent const& mec, std::vector<StateType> const& relevantStates,
                     storm::storage::SparseMatrix<ValueType> const& relevantStatesMatrix, ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                     BoundsDT<StateType, ValueType>& bounds) const;

    /*!
     * Based on the bounds of the successor states of this action, calculate the bounds of the action itself
     * @param action storm action index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return lower and upper bound of the action
     */
    std::pair<ValueType, ValueType> computeBoundsOfAction(ActionType const& action, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                          BoundsDT<StateType, ValueType> const& bounds) const;

    /*!
     * Based on the optimization direction, return 0 or 1
     * @param direction optimization direction (maximize or minimize)
     * @return 0 or 1
     */
    ValueType getLowestBound(storm::OptimizationDirection const& direction) const;

    /*!
     * Get two bounds, set to the lowest possible bound, given the optimization direction
     * either (0,0) or (1,1)
     * @param direction optimization direction (maximize or minimize)
     * @return pair of bounds
     */
    std::pair<ValueType, ValueType> getLowestBounds(storm::OptimizationDirection const& direction) const;

    /*!
     * Given the optimization direction, combine two bounds (max of both, resp. min of both)
     * @param direction optimization direction (maximize or minimize)
     * @param bounds1
     * @param bounds2
     * @return bounds
     */
    std::pair<ValueType, ValueType> combineBounds(storm::OptimizationDirection const& direction, std::pair<ValueType, ValueType> const& bounds1,
                                                  std::pair<ValueType, ValueType> const& bounds2) const;

    /*!
     * Backpropagation of BRTDP. After a simulation, call that function to update the lower and upper bounds
     * @param stack Stack of states and actions
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions (old, will be udpated in this function)
     * @return change of the lower bound of the initial state (aka "convergence progress")
     */
    ValueType updateProbabilityBoundsAlongSampledPath(StateActionStack& stack, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                 BoundsDT<StateType, ValueType>& bounds) const;

    /*!
     * Update the lower and upper bound of an action, given the bounds of its successors
     * @param state state index
     * @param action storm action-index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions (old, will be udpated in this function)
     * @return change in the lower bound
     */
    ValueType updateProbabilityOfAction(StateType const& state, ActionType const& action, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                   BoundsDT<StateType, ValueType>& bounds) const;

    /*!
     * Given an action, compute the lower resp. upper bounds of all other actions (given the optimization direction)
     * @param direction optimization direction (either maximize or minimize)
     * @param state state index
     * @param action storm action index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return returns the biggest upper (max) resp. smallest lower (min) bound of all other actions
     */
    ValueType computeBoundOverAllOtherActions(storm::OptimizationDirection const& direction, StateType const& state, ActionType const& action,
                                              ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                              BoundsDT<StateType, ValueType> const& bounds) const;

    /*!
     * Given an action, calculate its lower bound based on the successors
     * @param action storm action-index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return lower bound of the action
     */
    ValueType computeLowerBoundOfAction(ActionType const& action, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                        BoundsDT<StateType, ValueType> const& bounds) const;

    /*!
     * Given an action, calculate its upper bound based on the successors
     * @param action storm action-index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return upper bound of the action
     */
    ValueType computeUpperBoundOfAction(ActionType const& action, ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                        BoundsDT<StateType, ValueType> const& bounds) const;

    /*!
     * Based on the possible actions in this state, calculate its bounds
     * @param currentStateId state action index
     * @param explorationInformation storage of information during exploration
     * @param bounds bounds of states and actions
     * @return lower and upper bound of the state
     */
    std::pair<ValueType, ValueType> computeBoundsOfState(StateType const& currentStateId,
                                                         ExplorationInformationDT<StateType, ValueType> const& explorationInformation,
                                                         BoundsDT<StateType, ValueType> const& bounds) const;

    /*!
     * returns a random number in a range
     * @param start lower bound of the range
     * @param end upper bound of the range
     * @return random number
     */
    int getRandomNum(int const& start, int const& end);

    /*!
     *
     */
    arma::Mat<int> getTransformedStateValuations(arma::Row<long unsigned int>& labels, arma::Row<long unsigned int>* modules=nullptr);

    /*!
     *
     */
    uint getNumberOfActions() const;

    arma::Row<int> getAssignment() const;

   private:
    // stores the information of BRTDP

    // Exploration Information: information that we gather during exploration (number of states, actions, ...)
    ExplorationInformationDT<StateType, ValueType>* explorationInformation;

    // State Generation: uses the underlying prism program to generate new states and their behavior
    StateGenerationDT<StateType, ValueType>* stateGeneration;

    // Bounds: stores the bounds of states and actions, and calculates new ones
    BoundsDT<StateType, ValueType>* bounds;

    // Comparator: compares values up to the set precision
    storm::utility::ConstantsComparator<ValueType> comparator;

    // Statistics: struct that stores important statistics (number of simulations, steps, ...)
    StatisticsDT<StateType, ValueType>* stats;

    // The random number generator.
    mutable std::default_random_engine randomGenerator;

    // The state valuation
    StateValuationsDT* stateValuations;
};
}
}
}

#endif  // STORM_BRTDPINFORMATION_H
