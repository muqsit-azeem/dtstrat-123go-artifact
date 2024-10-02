#ifndef STORM_MODELCHECKER_EXPLORATIONDT_SPARSEEXPLORATIONMODELCHECKERDT_H_
#define STORM_MODELCHECKER_EXPLORATIONDT_SPARSEEXPLORATIONMODELCHECKERDT_H_

#include <random>

#include "storm/modelchecker/AbstractModelChecker.h"

#include "storm/storage/prism/Program.h"

#include "storm/generator/CompressedState.h"
#include "storm/generator/VariableInformation.h"

#include "storm/utility/ConstantsComparator.h"
#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"

namespace storm {

class Environment;

namespace storage {
class MaximalEndComponent;
}
namespace prism {
class Program;
}

namespace modelchecker {
namespace exploration_detail_dt {
template<typename StateType, typename ValueType>
class StateGenerationDT;
template<typename StateType, typename ValueType>
class ExplorationInformationDT;
template<typename StateType, typename ValueType>
class BoundsDT;
template<typename StateType, typename ValueType>
struct StatisticsDT;
}  // namespace exploration_detail

using namespace exploration_detail_dt;
using namespace dtstrat;

template<typename ModelType, typename StateType = uint32_t>
class SparseExplorationDTStratModelChecker : public AbstractModelChecker<ModelType> {
   public:
    typedef typename ModelType::ValueType ValueType;
    typedef StateType ActionType;
    typedef std::vector<std::pair<StateType, ActionType>> StateActionStack;

    SparseExplorationDTStratModelChecker(storm::prism::Program const& program);

    virtual bool canHandle(CheckTask<storm::logic::Formula, ValueType> const& checkTask) const override;

    /*!
     * main function, starts all the computation
     * @param env environment
     * @param checkTask task that should be checked
     * @return check result (usually two values, and (optionally) the scheduler
     */
    virtual std::unique_ptr<CheckResult> computeUntilProbabilities(Environment const& env,
                                                                   CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) override;

   private:
    /*!
     * Depending on the settings, return the corresponding learning configuration.
     * @param brtdpInformation contains pointer to all information of BRTDP
     * @return a pointer to the learning config (pointer because then we can have any config)
     */
    LearningConfig<StateType>* getLearningConfig(BRTDPInformation<StateType, ValueType>& brtdpInformation) const;

    /*!
     * Starts off the exploration = the algorithm
     * @param stateGeneration the class that generates states and information about them
     * @param explorationInformation the information we gather during runtime
     * @return inital state, lower bound, upper bound of the initial state
     */
    std::tuple<StateType, ValueType, ValueType> performExploration(BRTDPInformation<StateType, ValueType>& brtdpInfo) const;

    // The program that defines the model to check.
    storm::prism::Program program;

    // The random number generator.
    mutable std::default_random_engine randomGenerator;

    // A comparator used to determine whether values are equal.
    storm::utility::ConstantsComparator<ValueType> comparator;
};
}  // namespace modelchecker
}  // namespace storm

#endif /* STORM_MODELCHECKER_EXPLORATIONDT_SPARSEEXPLORATIONMODELCHECKERDT_H_ */
