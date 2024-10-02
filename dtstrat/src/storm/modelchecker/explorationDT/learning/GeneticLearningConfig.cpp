//
// Created by steffi on 11.01.22.
//

#include "storm/modelchecker/explorationDT/learning/GeneticLearningConfig.h"
#include "storm/modelchecker/explorationDT/learning/NonLearningConfig.h"
#include "storm/modelchecker/explorationDT/heuristics/NonHeuristic.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/ExplorationSettings.h"
#include <omp.h>

namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
GeneticLearningConfig<StateType, ValueType>::GeneticLearningConfig(double const& decisionEpsilon, int const& updateSteps, heuristicTypes const& heuristicType,
                                                                   exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInformation)
    : BasicLearningConfig<StateType, ValueType>(decisionEpsilon, updateSteps, heuristicType, brtdpInformation) {
    LearningConfig<StateType>::heuristic = new NonHeuristic<StateType>();
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::update() {
    if (BasicLearningConfig<StateType, ValueType>::updateNow()) {
        if (population.empty()) population = Population(*this, 10, decisionTreeSet);
        else {
            population.nextGeneration(decisionTreeSet);
        }
        prepareDecisionTrees();
        population.evovle(3);
        LearningConfig<StateType>::heuristic = population.getBestHeuristic(decisionTreeSet);
    }
    BasicLearningConfig<StateType, ValueType>::step();
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::prepareDecisionTrees() {
    auto existingGenerations = population.getGenerations();

    // write down the state valuations and choices in a file
    std::string fileName = BasicLearningConfig<StateType, ValueType>::prepareModuleDTwoCall();

    decisionTreeSet.newGeneration();
    // fill the decision trees based on which individuals we have
    for (auto const& individual : population.individuals) {
        if (decisionTreeSet.isSet(individual.isUseOnlySVOfModule())) continue;
        else {
            decisionTreeSet.getTrees(individual.isUseOnlySVOfModule(), fileName, this);
        }
    }
}

template<typename StateType, typename ValueType>
double GeneticLearningConfig<StateType, ValueType>::evaluateHeuristic(Heuristic<StateType>* heuristic) {
    ValueType result = 0;
    std::vector<std::pair<StateType, StateType>> stack;
    auto* config = new NonLearningConfig<StateType>(heuristic, BasicLearningConfig<StateType, ValueType>::brtdpInformation->getModelName());
    bool foundTerminalState = BasicLearningConfig<StateType, ValueType>::brtdpInformation->samplePathFromInitialState(
        BasicLearningConfig<StateType, ValueType>::brtdpInformation->getStateGeneration(),
        BasicLearningConfig<StateType, ValueType>::brtdpInformation->getExplorationInformation(), stack,
        BasicLearningConfig<StateType, ValueType>::brtdpInformation->getBounds(), BasicLearningConfig<StateType, ValueType>::brtdpInformation->getStats(),
        config, BasicLearningConfig<StateType, ValueType>::brtdpInformation->getStateValuations(),
        BasicLearningConfig<StateType, ValueType>::brtdpInformation->getMaxNumberOfSteps()*5);
    if (foundTerminalState) {
        result = 1 + BasicLearningConfig<StateType, ValueType>::brtdpInformation->updateProbabilityBoundsAlongSampledPath(
            stack, BasicLearningConfig<StateType, ValueType>::brtdpInformation->getExplorationInformation(),
            BasicLearningConfig<StateType, ValueType>::brtdpInformation->getBounds());
    }
    return result;
}

template<typename StateType, typename ValueType>
int GeneticLearningConfig<StateType, ValueType>::randomNum(int const& start, int const& end)
{
    return BasicLearningConfig<StateType,ValueType>::brtdpInformation->getRandomNum(start, end);
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::Population::evovle(uint const& numGenerations) {
    for (uint i=0; i<numGenerations; i++) {
        //#pragma omp parallel for
            for (auto& individual : individuals) {
                individual.getFitness();
            }
        // sort the population in decreasing order of fitness score (highest value - greatest fitness)
        std::sort(individuals.begin(), individuals.end(), [](Individual const& a, Individual const& b) {return a.getFitness()>b.getFitness();});
        // generate new offsprings for new generation
        std::vector<Individual> new_generation;

        // Perform Elitism, that mean 10% of fittest population
        // goes to the next generation
        int s = int(std::ceil(float(populationSize)/10.0));
        new_generation.reserve(s);
        for(int k = 0;k<s;k++)
            new_generation.push_back(individuals[k]);

        // From 50% of fittest population, Individuals
        // will mate to produce offspring
        s = populationSize-s;
        for(int k = 0;k<s;k++)
        {
            int len = individuals.size();
            int r = config->randomNum(0, int(len/2));
            Individual parent1 = individuals[r];
            r = config->randomNum(0, int(len/2));
            Individual parent2 = individuals[r];
            Individual offspring = parent1.mate(parent2, numIndividual);
            new_generation.push_back(offspring);
        }
        individuals = new_generation;
    }
}

template<typename StateType, typename ValueType>
GeneticLearningConfig<StateType, ValueType>::Population::Population(GeneticLearningConfig<StateType, ValueType>& config, uint const& size,
                                                                    DecisionTreeSet& decisionTreeSet)
    : populationSize(size), config(&config), currentGeneration(0) {
    std::vector<uint> bounds;
    uint numModules = config.brtdpInformation->getNumberOfModules();
    for (uint i=0; i<numModules; i++) {
        bounds.push_back(numModules-1);
    }
    bounds.push_back(1);
    bounds.push_back(1);
    for (uint i=0; i<numModules; i++) {
        bounds.push_back(numModules-1);
    }
    this->bounds = bounds;
    for (uint i=0; i<size; i++) individuals.push_back(Individual(config, bounds, decisionTreeSet, currentGeneration, numIndividual));
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>* GeneticLearningConfig<StateType, ValueType>::Population::getBestHeuristic(DecisionTreeSet& decisionTreeSet) {
    for (auto individual : individuals) individual.getFitness();
    std::sort(individuals.begin(), individuals.end(), [](Individual& a, Individual& b) {return a.getFitness()>b.getFitness();});
    STORM_LOG_TRACE("Best heuristic has fitness: " << individuals[0].getFitness() << std::endl);
    return individuals[0].transformToHeuristic(*(config->BasicLearningConfig<StateType, ValueType>::brtdpInformation));
}

template<typename StateType, typename ValueType>
GeneticLearningConfig<StateType, ValueType>::Population::Population() : config(nullptr), populationSize(0) {}

template<typename StateType, typename ValueType>
uint GeneticLearningConfig<StateType, ValueType>::Population::size() {
    return populationSize;
}

template<typename StateType, typename ValueType>
bool GeneticLearningConfig<StateType, ValueType>::Population::empty() {
    return (populationSize<1);
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::Population::nextGeneration(DecisionTreeSet& decisionTreeSet) {
    currentGeneration++;
    individuals.erase(individuals.begin()+int(populationSize/2), individuals.end());
    while(individuals.size()<populationSize) individuals.push_back(Individual(*config, bounds, decisionTreeSet, currentGeneration, numIndividual));
    for (auto& indi : individuals) indi.eraseFitness();
}
template<typename StateType, typename ValueType>
storm::storage::BitVector GeneticLearningConfig<StateType, ValueType>::Population::getGenerations() {
    storm::storage::BitVector result;
    for (auto const& individual : individuals) {
        uint generation = individual.getGeneration();
        if (result.size() <= generation) {
            auto test = generation+1;
            result.resize(test);
        }
        result.set(generation);
    }
    return result;
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>* GeneticLearningConfig<StateType, ValueType>::Individual::transformToHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInfo, double const& givenExplEps) {
    std::vector<const storm::storage::DecisionTree<uint>*> fittingDecisionTrees = decisionTreeSet->getDecisionTrees(isUseOnlySVOfModule());
    double explEps;
    if (givenExplEps>-1) explEps = givenExplEps;
    else {
        if (fitness>0) explEps = 0;
        else explEps = 0.5;
    }
    return new ModuleDecisionTreeHeuristic<StateType, ValueType>(&brtdpInfo, explEps, isUseOnlySVOfModule(), getModuleHeuristicType(), getOrdering(), getModuleToDTMapping(), fittingDecisionTrees);
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::Individual::eraseFitness() {
    fitness = -1;
}

template<typename StateType, typename ValueType>
double GeneticLearningConfig<StateType, ValueType>::Individual::getFitness() {
    STORM_LOG_TRACE("individual " << number);
    if (fitness>-1) {
        STORM_LOG_TRACE(" fitness: " << fitness << "\n");
        return fitness;
    }
    else {
        //double result = evaluateModuleHeuristic(config->brtdpInformation->getStateGeneration(), config->brtdpInformation->getExplorationInformation(),
        //                                        config->brtdpInformation->getBounds(), config->brtdpInformation->getStats());
        double result = evaluateModuleHeuristicWOCopy();
        fitness = result;
        STORM_LOG_TRACE(" fitness: " << fitness << "\n");
        return fitness;
    }
}

template<typename StateType, typename ValueType>
GeneticLearningConfig<StateType, ValueType>::Individual::Individual(GeneticLearningConfig<StateType, ValueType>& config, std::vector<uint>& bounds, DecisionTreeSet& decisionTreeSet, uint const& generationNum, long& numIndividual)
    : bounds(bounds), geneLength(bounds.size()), fitness(-1), config(&config), decisionTreeSet(&decisionTreeSet), generationNum(generationNum), number(numIndividual++) {
    for(uint i=0; i<geneLength/2+1; i++) {
        genes.push_back(config.randomNum(0, bounds[i]));
    }
    std::vector<uint> ordering(geneLength/2-1);
    std::iota(ordering.begin(), ordering.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(ordering.begin(), ordering.end(), g);
    genes.insert(genes.end(), ordering.begin(), ordering.end());
}

template<typename StateType, typename ValueType>
typename GeneticLearningConfig<StateType, ValueType>::Individual GeneticLearningConfig<StateType, ValueType>::Individual::mate(
    const GeneticLearningConfig::Individual& parent2, long& numIndividual) {
    uint numModules = config->brtdpInformation->getNumberOfModules();
    std::vector<uint> childGenes;
    childGenes.reserve(2*numModules+2);
    // only inherit this part randomly, don't mess up the ordering
    for (uint i=0; i<numModules+2; i++) {
        int prob = config->randomNum(0,100);
        if (prob<45) {
            childGenes.push_back(genes[i]);
        } else if (prob<90){
            childGenes.push_back(parent2.genes[i]);
        } else {
            childGenes.push_back(config->randomNum(0,bounds[i]));
        }
    }
    STORM_LOG_ASSERT(childGenes.size()==numModules+2,"Messed up genes 1...");
    int prob = config->randomNum(0,100);
    if (prob<45) {
        childGenes.insert(childGenes.end(), genes.begin()+numModules+2, genes.end());
        STORM_LOG_ASSERT(childGenes.size()==geneLength,"Messed up genes 2.1...");
    } else if (prob<90){
        childGenes.insert(childGenes.end(), parent2.genes.begin()+numModules+2, parent2.genes.end());
        STORM_LOG_ASSERT(childGenes.size()==geneLength,"Messed up genes 2.2...");
    } else {
        bool prob = bool(config->randomNum(0,1));
        std::vector<uint> ordering;
        if (prob) ordering=mutateOrdering(*this);
        else ordering=mutateOrdering(parent2);
        childGenes.insert(childGenes.end(), ordering.begin(), ordering.end());
        STORM_LOG_ASSERT(childGenes.size()==geneLength,"Messed up genes 2.3...");
    }
    STORM_LOG_ASSERT(childGenes.size()==geneLength,"Messed up genes 2...");
    return {*config, bounds, genes, *decisionTreeSet, generationNum, numIndividual};
}

template<typename StateType, typename ValueType>
GeneticLearningConfig<StateType, ValueType>::Individual::Individual(GeneticLearningConfig<StateType, ValueType>& config, std::vector<uint>& bounds,
                                                                    std::vector<uint>& genes,
                                                                    DecisionTreeSet& decisionTreeSet, uint const& generationNum, long& numIndividual)
    : bounds(bounds), geneLength(bounds.size()), fitness(-1), config(&config), genes(genes), decisionTreeSet(&decisionTreeSet), generationNum(generationNum), number(numIndividual++) {}

template<typename StateType, typename ValueType>
std::vector<uint> GeneticLearningConfig<StateType, ValueType>::Individual::mutateOrdering(const GeneticLearningConfig::Individual& parent) {
    uint numModules = config->brtdpInformation->getNumberOfModules();
    int index1 = config->randomNum(numModules+2,2*numModules+1);
    int index2 = config->randomNum(numModules+2,2*numModules+1);
    std::vector<uint> result = parent.genes;
    result[index1] = parent.genes[index2];
    result[index2] = parent.genes[index1];;
    return {result.begin()+numModules+2, result.end()};
}

template<typename StateType, typename ValueType>
bool GeneticLearningConfig<StateType, ValueType>::Individual::isUseOnlySVOfModule() const {
    uint numModules = config->brtdpInformation->getNumberOfModules();
    return genes[numModules];
}

template<typename StateType, typename ValueType>
std::map<uint, uint> GeneticLearningConfig<StateType, ValueType>::Individual::getModuleToDTMapping() const {
    uint numModules = config->brtdpInformation->getNumberOfModules();
    std::map<uint,uint> result;
    for (uint i=0; i<numModules; i++) result[i]=genes[i];
    return result;
}

template<typename StateType, typename ValueType>
moduleHeuristicTypes GeneticLearningConfig<StateType, ValueType>::Individual::getModuleHeuristicType() const {
    uint numModules = config->brtdpInformation->getNumberOfModules();
    uint val = genes[numModules+1];
    return static_cast<moduleHeuristicTypes>(val);
}

template<typename StateType, typename ValueType>
std::vector<uint> GeneticLearningConfig<StateType, ValueType>::Individual::getOrdering() const {
    uint numModules = config->brtdpInformation->getNumberOfModules();
    return {genes.begin()+numModules+2, genes.end()};
}

template<typename StateType, typename ValueType>
double GeneticLearningConfig<StateType, ValueType>::Individual::evaluateModuleHeuristicWOCopy() {
    ValueType result = 0;
    std::vector<std::pair<StateType, StateType>> stack;
    auto heuristic = transformToHeuristic(*(config->brtdpInformation), 0);
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo = config->brtdpInformation;
    auto* configNew = new NonLearningConfig<StateType>(heuristic, brtdpInfo->getModelName());
    bool foundTerminalState =
            config->brtdpInformation->samplePathFromInitialState(brtdpInfo->getStateGeneration(), brtdpInfo->getExplorationInformation(), stack, brtdpInfo->getBounds(), brtdpInfo->getStats(), config, brtdpInfo->getStateValuations(), brtdpInfo->getMaxNumberOfSteps() * 5);
    if (foundTerminalState) {
        // Target or Sink state?
        auto boundTerminalState = brtdpInfo->getBounds().getLowerBoundForState(stack.back().first, brtdpInfo->getExplorationInformation());
        STORM_LOG_TRACE("Bound of Terminal State: " << boundTerminalState << "\n");
        StateType initialStateIndex = brtdpInfo->getStateGeneration().getFirstInitialState();
        auto val = brtdpInfo->getBounds().getLowerBoundForState(initialStateIndex, brtdpInfo->getExplorationInformation());
        STORM_LOG_TRACE(" path: ");
        for (auto bla : stack) STORM_LOG_TRACE(bla.first << ":" << bla.second << "   ");
        result = 1 + brtdpInfo->updateProbabilityBoundsAlongSampledPath(stack, brtdpInfo->getExplorationInformation(), brtdpInfo->getBounds());
        val = brtdpInfo->getBounds().getLowerBoundForState(initialStateIndex, brtdpInfo->getExplorationInformation()) - val;
        double test = 0;
    }
    return result;
}

template<typename StateType, typename ValueType>
double GeneticLearningConfig<StateType, ValueType>::Individual::evaluateModuleHeuristic(exploration_detail_dt::StateGenerationDT<StateType, ValueType>& stateGenerationOld,
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> explorationInformationOld,
    exploration_detail_dt::BoundsDT<StateType, ValueType> boundsOld, exploration_detail_dt::StatisticsDT<StateType, ValueType> statsOld, StateValuationsDT stateValuationsOld) {
    ValueType result = 0;
    std::vector<std::pair<StateType, StateType>> stack;
    prism::Program program;
    generator::NextStateGeneratorOptions optionsC;
    stateGenerationOld.getProgramOptions(program, optionsC);
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> explorationInformation(explorationInformationOld);
    exploration_detail_dt::StateGenerationDT<StateType, ValueType> stateGeneration(stateGenerationOld, program, optionsC, explorationInformation);
    exploration_detail_dt::BoundsDT<StateType, ValueType> boundsBRTDP(boundsOld);
    exploration_detail_dt::StatisticsDT<StateType, ValueType> stats(statsOld);
    StateValuationsDT stateValuations(stateValuationsOld);
    exploration_detail_dt::BRTDPInformation<StateType, ValueType> brtdpInfo(
        explorationInformation, stateGeneration, boundsBRTDP, stats, stateValuations, storm::settings::getModule<storm::settings::modules::ExplorationSettings>().getPrecision());
    auto heuristic = transformToHeuristic(brtdpInfo, 0);
    auto* config = new NonLearningConfig<StateType>(heuristic, brtdpInfo.getModelName());
    bool foundTerminalState =
        brtdpInfo.samplePathFromInitialState(stateGeneration, explorationInformation, stack, boundsBRTDP, stats, config, stateValuations,brtdpInfo.getMaxNumberOfSteps() * 5);
    if (foundTerminalState) {
        // Target or Sink state?
        auto boundTerminalState = boundsBRTDP.getLowerBoundForState(stack.back().first, explorationInformation);
        STORM_LOG_TRACE("Bound of Terminal State: " << boundTerminalState << "\n");
        StateType initialStateIndex = stateGeneration.getFirstInitialState();
        auto val = boundsBRTDP.getLowerBoundForState(initialStateIndex, explorationInformation);
        STORM_LOG_TRACE(" path: ");
        for (auto bla : stack) STORM_LOG_TRACE(bla.first << ":" << bla.second << "   ");
        result = 1 + brtdpInfo.updateProbabilityBoundsAlongSampledPath(stack, explorationInformation, boundsBRTDP);
        val = boundsBRTDP.getLowerBoundForState(initialStateIndex, explorationInformation) - val;
        double test = 0;
    }
    return result;
}

template<typename StateType, typename ValueType>
double GeneticLearningConfig<StateType, ValueType>::Individual::getFitness() const {
    return fitness;
}
template<typename StateType, typename ValueType>
uint GeneticLearningConfig<StateType, ValueType>::Individual::getGeneration() const {
    return generationNum;
}

template<typename StateType, typename ValueType>
uint GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::size() const {
    return decisionTrees.size();
}

template<typename StateType, typename ValueType>
bool GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::isSet(const bool& isUseOnlySVOfModule) const {
    if (onlySVForModule.size() <= generationBegin.back()) return false;
    else if ((isUseOnlySVOfModule && onlySVForModule.getNextSetIndex(generationBegin.back())!=onlySVForModule.size())
             || (!isUseOnlySVOfModule && onlySVForModule.getNextUnsetIndex(generationBegin.back())!=onlySVForModule.size())) return true;
    else return false;
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::getTrees(const bool& isUseOnlySVOfModule, std::string const& filename, BasicLearningConfig<StateType, ValueType>* config) {
    config->callSklearnMulti(filename, isUseOnlySVOfModule);
    std::vector<const storm::storage::DecisionTree<uint>*> trees = ModuleDecisionTreeHeuristic<StateType, ValueType>::getTrees(config->getNumberOfModules(), config->getModelName());
    for (auto dtptr : trees) {
        auto tree = new storm::storage::DecisionTree<uint>(*dtptr);
        decisionTrees.push_back(tree);
    }
    onlySVForModule.resize(decisionTrees.size(), isUseOnlySVOfModule);
}

template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::newGeneration() {
    generationBegin.push_back(decisionTrees.size());
}

template<typename StateType, typename ValueType>
std::vector<const storm::storage::DecisionTree<uint>*> GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::getDecisionTrees(const bool& useOnlySV) const {
    std::vector<const storm::storage::DecisionTree<uint>*> result;
    for (uint i=generationBegin.back(); i<decisionTrees.size(); i++) {
        if (useOnlySV && onlySVForModule.get(i) || !useOnlySV && !onlySVForModule.get(i)) {
            const storm::storage::DecisionTree<uint>* tree = decisionTrees[i];
            result.push_back(tree);
        }
    }
    return result;
}
template<typename StateType, typename ValueType>
void GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::removeUnneccessaryTrees(const storage::BitVector& appearingGenerations) {
    std::vector<uint> remove;
    for (uint i=0; i<decisionTrees.size(); i++) {
        if (appearingGenerations.get(getGeneration(i))) {
            decisionTrees[i]->clear();
            delete decisionTrees[i];
            remove.push_back(i);
        }
    }
}
template<typename StateType, typename ValueType>
uint GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::getGeneration(const uint& i) {
    for (uint j=generationBegin.size()-1; j>=0; j--) {
        if (i>=generationBegin[j]) return j;
    }
    return -1;
}

template<typename StateType, typename ValueType>
GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::~DecisionTreeSet() {
    for (auto tree : decisionTrees) delete tree;
}

template<typename StateType, typename ValueType>
typename GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet& GeneticLearningConfig<StateType, ValueType>::DecisionTreeSet::operator=(
    const GeneticLearningConfig::DecisionTreeSet& copy) {
    if (this == &copy) return *this;
    else {
        for (auto tree : decisionTrees) delete tree;
        decisionTrees = copy.decisionTrees;
        onlySVForModule = copy.onlySVForModule;
        generationBegin = copy.generationBegin;
        return *this;
    }
}

    template<typename StateType, typename ValueType>
    void GeneticLearningConfig<StateType, ValueType>::setUpdateNow() {
        STORM_LOG_TRACE("Set ExplEpsilon to 0.5!\n");
        LearningConfig<StateType>::heuristic->setExplorationEpsilon(0.5);
    }

template class GeneticLearningConfig<uint32_t, double>;

}
}
}
