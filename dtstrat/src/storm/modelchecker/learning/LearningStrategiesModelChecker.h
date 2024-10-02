//
// Created by steffi on 27.07.22.
//

#ifndef STORM_LEARNINGSTRATEGIESMODELCHECKER_H
#define STORM_LEARNINGSTRATEGIESMODELCHECKER_H

#include <random>

#include "storm/modelchecker/AbstractModelChecker.h"

#include "storm/storage/prism/Program.h"

#include "storm/generator/CompressedState.h"
#include "storm/generator/VariableInformation.h"

#include "storm/utility/ConstantsComparator.h"
#include "storm/storage/DecisionTree.h"
#include "storm/modelchecker/learning/DTInterface.h"
#include "storm/storage/jani/Property.h"
#include "storm/models/ModelBase.h"
#include "storm/storage/Scheduler.h"
#include "storm/models/sparse/Mdp.h"

#include "modelchecker/explorationDT/BSCCCandidate.h"
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"

#include <armadillo>

namespace storm {

class Environment;

namespace storage {
class MaximalEndComponent;
}
namespace prism {
class Program;
}

namespace modelchecker {

using namespace exploration_detail_dt;

template<typename ModelType, typename StateType = uint32_t>
class LearningStrategiesModelChecker : public AbstractModelChecker<ModelType> {
   public:
    typedef typename ModelType::ValueType ValueType;
    typedef StateType ActionType;
    typedef std::vector<std::pair<StateType, ActionType>> StateActionStack;

    LearningStrategiesModelChecker(storm::prism::Program const& program,
                                   std::vector<storm::jani::Property>& properties,
                                   boost::optional<std::vector<storm::jani::Property>>& preprocessedProperties);

    virtual bool canHandle(CheckTask<storm::logic::Formula, ValueType> const& checkTask) const override;

    void getAllActionsFromModel(storm::prism::Program const& program);

    /*!
     * main function, starts all the computation
     * @param env environment
     * @param checkTask task that should be checked
     * @return check result (usually two values, and (optionally) the scheduler
     */
    virtual std::unique_ptr<CheckResult> computeUntilProbabilities(Environment const& env,
                                                                   CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) override;

    /*!
     * Evaluate a small instance of a possibly bigger program
     * @param changedProgram Prims-program with a changed constant
     * @param storeToFile Whether the states and actions should be stored in a file
     * @param env storm-environment (has no apparent use, but has to be handed over)
     * @param first necessary for DTControl (then we need to print a header). If this is the first model that is evaluated
     * @param checkTask The property that should be checked
     * @param onlyReachable Whether or not to print only reachable states under the strategy.
     * @param labels actions (if they should not be stored in a file, they end up in that variable)
     * @param data state-valuations for the actions (if they should not be stored in a file, they end up in that variable)
     * @param properties also the formula that is checked
     * @param preprocessedProperties some more formula stuff that was just copy-pasted and no one actually knows why we need it
     */
    uint checkSmallInstance(storm::prism::Program& changedProgram, bool const& storeToFile, storm::Environment const& env, bool& first,
                            CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, bool const& onlyReachable, arma::Row<long unsigned int>& labels,
                            arma::Mat<int>& data, std::vector<storm::jani::Property> const& properties,
                            boost::optional<std::vector<storm::jani::Property>> const& preprocessedProperties, std::string const& dataFileName="data.csv",
                            std::string const& oldFileName="", bool const& semanticActions=false);

    /*!
     * Given a file that contains a decision-tree, load it and store it in the class-attribute "decision-tree"
     * Can either be a json (from dtcontrol) or a dot (in a specific expected format)
     * @param filename the filename
     */
    void getDTFromFile(std::string const& filename);

    /*!
     * Learn a decision tree on data that is stored in a csv-file
     * It should be a plain csv-file, that works for all
     * If it contains a dtcontrol-header, it can be used for dtcontrol but for nothing else
     * @param datafilename filename
     * @param method method (dtcontrol, sklearn, or mlpack)
     */
    void getDTFromDataFile(std::string const& datafilename, std::string const& method);

    void learnDTFromChildren(Environment const& env, CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask);

    void useDTStrat(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, uint const& numRuns, bool const& semanticActions);

    void getDTControlTree(std::string const& filename);

    void getSklearnTree(std::string const& filename);

    void getMLPackTree(std::string const& filename);

    void getMLPackTree(arma::Mat<int> const& data,
                       arma::Row<long unsigned int> const& labels,
                        uint const& numberOfActions);

    std::pair<arma::Mat<int>, arma::Row<int>> readCSVFile(std::string const& filename);

    void runBRTDPforNTimes(uint const& N, BRTDPInformation<StateType, ValueType>& brtdpInfo, bool const& semanticActions);

    float evaluateVIMC(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, bool const& semanticActions);

//    void createExplicitMCinDRNOld(StateGenerationDT<StateType, ValueType>& stateGeneration,
//                                      ExplorationInformationDT<StateType, ValueType>& explorationInformation,
//                                      StateValuationsDT& stateValuations);
    void createExplicitMCinDRN(StateGenerationDT<StateType, ValueType>& stateGeneration,
                               ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                               StateValuationsDT& stateValuations, bool const& semanticActions);

    double evaluateSMC(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, bool const& semanticActions);

    std::shared_ptr<storm::models::ModelBase> getModelFromProgram(prism::Program const& program, std::vector<storm::jani::Property> const& properties,
                                                                  boost::optional<std::vector<storm::jani::Property>> const& preprocessedProperties);

    std::vector<std::shared_ptr<storm::logic::Formula const>> createFormulasToRespect(std::vector<storm::jani::Property> const& properties);

    void storeData(std::string const& filename, arma::Mat<int> const& data, arma::Row<long unsigned int> const& labels);

    void printHeader(std::string const& filename, uint const& size);

    void printToFile(std::string const& filename, storm::storage::Scheduler<ValueType> const& scheduler, storm::models::sparse::Mdp<ValueType> const& model,
                     storm::storage::BitVector const& subsystem, bool const& semanticActions);

    uint getNumberOfStateVariables(CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask, storm::prism::Program& program);

    void getValuationsAndLabels(storm::storage::Scheduler<ValueType> const& scheduler, storm::models::sparse::Mdp<ValueType> const& model, arma::Mat<int>& data,
                                arma::Row<long unsigned int>& labels, storm::storage::BitVector const& subsystem, bool const& semanticAction);

    bool samplePathFromInitialState(StateGenerationDT<StateType, ValueType>& stateGeneration,
                                    ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                    StateValuationsDT& stateValuations,
                                    int const& stopAfterNSteps, bool const& semanticActions) const;

    typename LearningStrategiesModelChecker<ModelType, StateType>::ActionType sampleActionOfState(
        StateType const& currentStateId, ExplorationInformationDT<StateType, ValueType> const& explorationInformation) const;

    std::pair<bool, bool>  exploreState(StateGenerationDT<StateType, ValueType>& stateGeneration, const StateType& currentStateId,
                 const generator::CompressedState& currentState,
                 ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                 StateValuationsDT& stateValuations) const;

    StateType sampleSuccessorFromAction(const StateType& currentStateId, const ActionType& chosenAction,
                                        const ExplorationInformationDT<StateType, ValueType>& explorationInformation, std::set<StateType>& statesWithUnavailDTAction) const;

   private:
    // The program that defines the model to check.
    storm::prism::Program program;
    storm::prism::Program substprogram;
    std::string modelName;
    std::string propertyName;

    std::map<std::string, uint> completeActionNameToIndexMapping;
    std::vector<std::vector<std::string>> actionIndexToLabelMapping;
    std::map<std::string, uint> shortActionNameToIndexMapping;
    std::vector<std::string> indexToShortActionName;

    // The random number generator.
    mutable std::default_random_engine randomGenerator;

    // A comparator used to determine whether values are equal.
    storm::utility::ConstantsComparator<ValueType> comparator;

    // DecisionTree
    //storm::storage::DecisionTree<uint>* decisionTree;
    dtstrat::DTInterface<uint, ValueType>* decisionTree;

    std::vector<storm::jani::Property> properties;
    boost::optional<std::vector<storm::jani::Property>> preprocessedProperties;

    uint stopAfter;
    uint numSMC;

};
}  // namespace modelchecker
}  // namespace storm

#endif  // STORM_LEARNINGSTRATEGIESMODELCHECKER_H
