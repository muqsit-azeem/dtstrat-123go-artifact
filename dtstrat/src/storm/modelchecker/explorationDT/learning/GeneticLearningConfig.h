//
// Created by steffi on 11.01.22.
//

#ifndef STORM_GENETICLEARNINGCONFIG_H
#define STORM_GENETICLEARNINGCONFIG_H

#include "storm/modelchecker/explorationDT/learning/BasicLearningConfig.h"

namespace storm{
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class GeneticLearningConfig : public BasicLearningConfig<StateType,ValueType> {
    class DecisionTreeSet{
       public:
        DecisionTreeSet()=default;
        uint size() const;
        bool isSet(bool const& isUseOnlySVOfModule) const;
        void getTrees(bool const& isUseOnlySVOfModule, std::string const& filename, BasicLearningConfig<StateType, ValueType>* config);
        void newGeneration();
        std::vector<const storm::storage::DecisionTree<uint>*> getDecisionTrees(bool const& useOnlySV) const;
        void removeUnneccessaryTrees(storm::storage::BitVector const& appearingGenerations);
        ~DecisionTreeSet();
        DecisionTreeSet& operator=(DecisionTreeSet const& copy);

       private:
        std::vector<storm::storage::DecisionTree<uint>*> decisionTrees;
        storm::storage::BitVector onlySVForModule;
        std::vector<uint> generationBegin;

        uint getGeneration(uint const& i);
    };

    class Individual {
       public:
        Individual(GeneticLearningConfig<StateType, ValueType>& config, std::vector<uint>& bounds, DecisionTreeSet& decisionTreeSet, uint const& generationNum, long& numIndividual);
        Individual(GeneticLearningConfig<StateType, ValueType>& config, std::vector<uint>& bounds, std::vector<uint>& genes, DecisionTreeSet& decisionTreeSet, uint const& generationNum, long& numIndividual);

        Individual mate(Individual const& parent2, long& numIndividual);
        std::vector<uint> mutateOrdering(Individual const& parent);
        double getFitness();
        double getFitness() const;
        ModuleDecisionTreeHeuristic<StateType, ValueType>* transformToHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInfo, double const& givenExplEps=-1.0);
        double evaluateModuleHeuristic(exploration_detail_dt::StateGenerationDT<StateType, ValueType>& stateGenerationOld,
                                       exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> explorationInformationOld,
                                       exploration_detail_dt::BoundsDT<StateType, ValueType> boundsOld,
                                       exploration_detail_dt::StatisticsDT<StateType, ValueType> statsOld,
                                       StateValuationsDT stateValuationsOld);
        double evaluateModuleHeuristicWOCopy();
        bool isUseOnlySVOfModule() const;
        uint getGeneration() const;

        void eraseFitness();

       private:
        std::vector<uint> bounds;
        std::vector<uint> genes;
        uint geneLength;
        GeneticLearningConfig<StateType, ValueType>* config;
        std::vector<storm::storage::DecisionTree<uint>*> decisionTrees;
        double fitness;
        DecisionTreeSet* decisionTreeSet;
        uint generationNum;

        uint number;


        std::map<uint, uint> getModuleToDTMapping() const;
        moduleHeuristicTypes getModuleHeuristicType() const;
        std::vector<uint> getOrdering() const;
    };

    class Population {
       public:
        std::vector<Individual> individuals;
        Population();
        Population(GeneticLearningConfig<StateType, ValueType>& config, uint const& size, DecisionTreeSet& decisionTreeSet);
        void evovle(uint const& numGenerations);
        ModuleDecisionTreeHeuristic<StateType, ValueType>* getBestHeuristic(DecisionTreeSet& decisionTreeSet);
        uint size();
        bool empty();
        void nextGeneration(DecisionTreeSet& decisionTreeSet);
        storm::storage::BitVector getGenerations();

       private:
        uint populationSize;
        GeneticLearningConfig<StateType, ValueType>* config;
        std::vector<uint> bounds;
        uint currentGeneration;

        long numIndividual=0;
    };

   public:
    GeneticLearningConfig(double const& decisionEpsilon, int const& updateSteps, heuristicTypes const& heuristicType,
                          exploration_detail_dt::BRTDPInformation<StateType, ValueType>& brtdpInformation);

    void update() override;

    void setUpdateNow() override;

   private:
    double evaluateHeuristic(Heuristic<StateType>* heuristic); // fitness function
    int randomNum(int const& start, int const& end);
    void prepareDecisionTrees();

    DecisionTreeSet decisionTreeSet;
    Population population;


};



}
}
}
#endif  // STORM_GENETICLEARNINGCONFIG_H
