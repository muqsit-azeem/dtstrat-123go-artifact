#include "storm/modelchecker/explorationDT/StateGenerationDT.h"

#include "storm/modelchecker/explorationDT/ExplorationInformationDT.h"

namespace storm {
namespace modelchecker {
namespace exploration_detail_dt {

template<typename StateType, typename ValueType>
StateGenerationDT<StateType, ValueType>::StateGenerationDT(storm::prism::Program const& program,
                                                       ExplorationInformationDT<StateType, ValueType>& explorationInformation,
                                                       storm::expressions::Expression const& conditionStateExpression,
                                                       storm::expressions::Expression const& targetStateExpression,
                                                       generator::NextStateGeneratorOptions const& options)
    : generator(program, options),
      options(options),
      stateStorage(generator.getStateSize()),
      conditionStateExpression(conditionStateExpression),
      targetStateExpression(targetStateExpression),
      stateValuationsBuilder(generator.initializeStateValuationsBuilder())
{
    stateToIdCallback = [&explorationInformation, this](storm::generator::CompressedState const& state) -> StateType {
        StateType newIndex = stateStorage.getNumberOfStates();

        // Check, if the state was already registered.
        std::pair<StateType, std::size_t> actualIndexBucketPair = stateStorage.stateToId.findOrAddAndGetBucket(state, newIndex);

        if (actualIndexBucketPair.first == newIndex) {
            explorationInformation.addUnexploredState(newIndex, state);
        }
        return actualIndexBucketPair.first;
    };
}

template<typename StateType, typename ValueType>
void StateGenerationDT<StateType, ValueType>::load(storm::generator::CompressedState const& state, uint const& stateId) {
    generator.load(state);
    generator.addStateValuation(stateId, stateValuationsBuilder);
}

template<typename StateType, typename ValueType>
std::vector<StateType> StateGenerationDT<StateType, ValueType>::getInitialStates() {
    return stateStorage.initialStateIndices;
}

template<typename StateType, typename ValueType>
storm::generator::StateBehavior<ValueType, StateType> StateGenerationDT<StateType, ValueType>::expand() {
    return generator.expand(stateToIdCallback);
}

template<typename StateType, typename ValueType>
bool StateGenerationDT<StateType, ValueType>::isConditionState() const {
    return generator.satisfies(conditionStateExpression);
}

template<typename StateType, typename ValueType>
bool StateGenerationDT<StateType, ValueType>::isTargetState() const {
    return generator.satisfies(targetStateExpression);
}

template<typename StateType, typename ValueType>
void StateGenerationDT<StateType, ValueType>::computeInitialStates() {
    stateStorage.initialStateIndices = generator.getInitialStates(stateToIdCallback);
}

template<typename StateType, typename ValueType>
StateType StateGenerationDT<StateType, ValueType>::getFirstInitialState() const {
    return stateStorage.initialStateIndices.front();
}

template<typename StateType, typename ValueType>
std::size_t StateGenerationDT<StateType, ValueType>::getNumberOfInitialStates() const {
    return stateStorage.initialStateIndices.size();
}

template <typename StateType, typename ValueType>
void StateGenerationDT<StateType, ValueType>::addStateValuation(storm::storage::sparse::state_type const& currentStateIndex, storm::storage::sparse::StateValuationsBuilder& valuationsBuilder) const {
    generator.addStateValuation(currentStateIndex, valuationsBuilder);
}

template <typename StateType, typename ValueType>
std::vector<int> StateGenerationDT<StateType, ValueType>::getStateVariableModuleAssociation() const {
    return generator.getVariableInformation().getStateVariableModuleAssociation();
}

template <typename StateType, typename ValueType>
storm::storage::sparse::StateValuations::StateValuation StateGenerationDT<StateType, ValueType>::getStateValuation(StateType const& stateId) {
    //storm::storage::sparse::StateValuations currentStateValuations = stateValuationsBuilder.build(0); // 0 is not correct, but could work
    return stateValuationsBuilder.getStateValuation(stateId);
}

template<typename StateType, typename ValueType>
void StateGenerationDT<StateType, ValueType>::printDTcontrolHeader(std::ofstream& ofs) {
    ofs << "#NON-PERMISSIVE\n";
    auto stateValuation = getStateValuation(0);
    ofs << "#BEGIN " << stateValuation.size() << " 1\n";
}

template<typename StateType, typename ValueType>
void StateGenerationDT<StateType, ValueType>::printModuleDTHeader(std::ofstream& ofs) {
    auto stateValuation = getStateValuation(0);
    ofs << "#BEGIN " << stateValuation.size()*2 << " 1";
    for (uint l = 0; l < stateValuation.size() * 2; l++) ofs << ",";
    ofs << ",\n";
}

template<typename StateType, typename ValueType>
storm::storage::sparse::StateValuations StateGenerationDT<StateType, ValueType>::getAllStateValuations() {
    storm::storage::sparse::StateValuations currentStateValuations = stateValuationsBuilder.build(0); // 0 is not correct, but could work
    return currentStateValuations;
}

template<typename StateType, typename ValueType>
StateGenerationDT<StateType, ValueType>::StateGenerationDT(StateGenerationDT<StateType, ValueType>& old, storm::prism::Program& program, generator::NextStateGeneratorOptions& options, ExplorationInformationDT<StateType, ValueType>& explorationInformation) : generator(program, options), stateStorage(generator.getStateSize()) {
    stateStorage = old.getAllData(program, options, stateValuationsBuilder, conditionStateExpression, targetStateExpression, explorationInformation);
    stateToIdCallback = [&explorationInformation, this](storm::generator::CompressedState const& state) -> StateType {
        StateType newIndex = stateStorage.getNumberOfStates();

        // Check, if the state was already registered.
        std::pair<StateType, std::size_t> actualIndexBucketPair = stateStorage.stateToId.findOrAddAndGetBucket(state, newIndex);

        if (actualIndexBucketPair.first == newIndex) {
            explorationInformation.addUnexploredState(newIndex, state);
        }
        return actualIndexBucketPair.first;
    };
}

template<typename StateType, typename ValueType>
storm::storage::sparse::StateStorage<StateType> StateGenerationDT<StateType, ValueType>::getAllData(prism::Program& program, generator::NextStateGeneratorOptions& optionsCopy,
                                                         storm::storage::sparse::StateValuationsBuilder& stateValuationsBuilderC,
                                                         storm::expressions::Expression& conditionStateExpressionC, storm::expressions::Expression& targetStateExpressionC, ExplorationInformationDT<StateType, ValueType>& explorationInformation) {
    program = generator.getProgram();
    optionsCopy = options;
    stateValuationsBuilderC = stateValuationsBuilder;
    conditionStateExpressionC = conditionStateExpression;
    targetStateExpressionC = targetStateExpression;
    return stateStorage;
}

template<typename StateType, typename ValueType>
void StateGenerationDT<StateType, ValueType>::getProgramOptions(prism::Program& program, generator::NextStateGeneratorOptions& optionsC) {
    program = generator.getProgram();
    optionsC = options;
}

template class StateGenerationDT<uint32_t, double>;
}  // namespace exploration_detail
}  // namespace modelchecker
}  // namespace storm
