#ifndef STORM_MODELCHECKER_EXPLORATION_EXPLORATION_DETAIL_STATEGENERATION_H_
#define STORM_MODELCHECKER_EXPLORATION_EXPLORATION_DETAIL_STATEGENERATION_H_

#include "storm/generator/CompressedState.h"
#include "storm/generator/PrismNextStateGenerator.h"

#include "storm/storage/sparse/StateStorage.h"

namespace storm {
namespace generator {
template<typename ValueType, typename StateType>
class PrismNextStateGenerator;
}

namespace modelchecker {
namespace exploration_detail_dt {

template<typename StateType, typename ValueType>
class ExplorationInformationDT;

template<typename StateType, typename ValueType>
class StateGenerationDT {
   public:
    StateGenerationDT(storm::prism::Program const& program, ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                    storm::expressions::Expression const& conditionStateExpression, storm::expressions::Expression const& targetStateExpression,
                    generator::NextStateGeneratorOptions const& options);

    StateGenerationDT(StateGenerationDT<StateType, ValueType>& old, storm::prism::Program& program, generator::NextStateGeneratorOptions& options, ExplorationInformationDT<StateType, ValueType>& explorationInformation);

    void load(storm::generator::CompressedState const& state, uint const& stateId);

    std::vector<StateType> getInitialStates();

    storm::generator::StateBehavior<ValueType, StateType> expand();

    void computeInitialStates();

    StateType getFirstInitialState() const;

    std::size_t getNumberOfInitialStates() const;

    bool isConditionState() const;

    bool isTargetState() const;

    void addStateValuation(storm::storage::sparse::state_type const& currentStateIndex, storm::storage::sparse::StateValuationsBuilder& valuationsBuilder) const;

    std::vector<int> getStateVariableModuleAssociation() const;

    storm::storage::sparse::StateValuations::StateValuation getStateValuation(StateType const& stateId);

    storm::storage::sparse::StateValuations getAllStateValuations();

    void printDTcontrolHeader(std::ofstream& ofs);
    void printModuleDTHeader(std::ofstream& ofs);

    void getProgramOptions(storm::prism::Program& program, generator::NextStateGeneratorOptions& options);

   private:
    storm::generator::PrismNextStateGenerator<ValueType, StateType> generator;
    std::function<StateType(storm::generator::CompressedState const&)> stateToIdCallback;
    storm::storage::sparse::StateValuationsBuilder stateValuationsBuilder;

    storm::storage::sparse::StateStorage<StateType> stateStorage;

    storm::expressions::Expression conditionStateExpression;
    storm::expressions::Expression targetStateExpression;
    generator::NextStateGeneratorOptions options;

    storm::storage::sparse::StateStorage<StateType> getAllData(storm::prism::Program& program, generator::NextStateGeneratorOptions& options,
                                                               storm::storage::sparse::StateValuationsBuilder& stateValuationsBuilder,
                                                               storm::expressions::Expression& conditionStateExpression, storm::expressions::Expression& targetStateExpression, ExplorationInformationDT<StateType, ValueType>& explorationInformation);
};

}  // namespace exploration_detail
}  // namespace modelchecker
}  // namespace storm

#endif /* STORM_MODELCHECKER_EXPLORATION_EXPLORATION_DETAIL_STATEGENERATION_H_ */
