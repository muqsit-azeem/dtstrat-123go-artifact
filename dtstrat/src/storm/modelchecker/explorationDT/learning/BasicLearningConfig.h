//
// Created by steffi on 05.01.22.
//

#ifndef STORM_BASICLEARNING_H
#define STORM_BASICLEARNING_H

#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"
#include "storm/modelchecker/explorationDT/heuristics/DecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleDecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/heuristics/MLDecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/heuristics/MLModuleDecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"

namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class BasicLearningConfig : public LearningConfig<StateType> {
   public:
    BasicLearningConfig(double const& decisionEpsilon, int const& updateSteps, heuristicTypes const& heuristicType,
                        exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInformation);

    /*!
     * Update function that overrides the LearningConfig-update
     * Updates the underlying heuristic
     */
    virtual void update() override;

    /*!
     * Calls a python-script using sklearn to generate several trees for each module
     * @param fileName the csv-file that contains the state-valuations and the actions
     * @param useOnlyStateVariablesOfModuleItself learn a module tree only by using the state-variables of the module itself (not of the others)
     */
    void callSklearnMulti(std::string const& fileName, bool const& useOnlyStateVariablesOfModuleItself) const;

    /*!
     * Retrieve the model name from the BRTDP-information
     */
    std::string getModelName() const;

    /*!
     * Retrieve the number of modules in the model by looking at the BRTDP-information
     */
    uint getNumberOfModules() const;

    /*!
     * Resets the variable that counts when to do an update. It will now be performed in the next iteration.
     */
    void setUpdateNow();

    /*!
     * Calls dtcontrol from command line and retrieves the json-output
     * @param fileName file containing the state-valuations and actions
     * @return decision-tree in json-format
     */
    static std::string callDTcontrol(std::string const& fileName);

   protected:
    /*!
     * Returns if the update should be performed now
     */
    bool updateNow();

    /*!
     * Performs a step (counting the number of steps)
     */
    void step();

    /*!
     * Updates a DTControl-heuristic
     * First: gets all state-valuations and actions and writes them in a file
     * Second: call DTControl in Python
     * Third: retrieve a decision-tree from the DTControl-output (=json-format)
     */
    void updateDtcontrol();

    /*!
     * Updates a MLPackDecisionTree-heuristic
     * First: gets all state-valuations and actions
     * Second: calls MLPack to generate a new decision tree
     */
    void updateMLPackDT();

    /*!
     * Updates the MLPack-Module-dt_heuristic
     */
     void updateMlPackModuleDTs();

    /*!
     * Generates a new ModuleDT-Heuristic
     * First: calls prepareModuleDT (to write the file and call python)
     * Second: calls getNewMultiDTHeuristic (generates a new moduleDT-heuristic)
     */
    ModuleDecisionTreeHeuristic<StateType, ValueType>* updateModuleDt();

    /*!
     * Generates a new DecisionTree-Heuristic based on a string (dot-file or json)
     * @param someString either dot-file or json
     * @param fromJson determines whether dot-file (False) or json-format (True)
     */
    DecisionTreeHeuristic<StateType, ValueType>* getNewDT(std::string const& someString, bool const& fromJson);

    /*!
     * Gets all state-valuations and actions and writes them in a file
     */
    void prepareModuleDT();

    /*!
     * Generates a new MultiDecisionTree-Heuristic by just calling the initializer
     */
    ModuleDecisionTreeHeuristic<StateType, ValueType>* getNewMultiDTHeuristic();

    MLModuleDecisionTreeHeuristic<StateType, ValueType>* getNewMLModuleDTHeuristic(arma::Row<long unsigned int>& labels,
                                                                                   arma::Mat<int>& data,
                                                                                   arma::Row<long unsigned int>& modules);

    /*!
     * Prints the state-valuations and actions into a file and returns the filename
     */
    std::string prepareModuleDTwoCall();

    /*!
     * probability-threshold when to use exploration and when explotation
     */
    double decisionEpsilon;

    /*!
     * Pointer to the BRTDP-information
     */
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInformation;

    /*!
     * defines the heuristic-type that should be learned
     */
    heuristicTypes heuristicType;

    /*
     * Was probably meant to check if the Config already contains a valid heuristic, but I think it is never used
     */
    // TODO see if the variable is actually necessary
    bool set;

    /*!
     * Counts the steps since the last update
     */
    int updateSteps;

    Heuristic<StateType>* getHeuristic();


private:
    /*!
     * resets the update-step-counter to 0
     */
    void resetSteps();

    /*!
     * determines if the module-DT should only be learned by using the state-variables that belong to the module itself
     */
    bool useOnlyStateVariablesOfModuleItself = true;
};
}
}
}

#endif  // STORM_BASICLEARNING_H
