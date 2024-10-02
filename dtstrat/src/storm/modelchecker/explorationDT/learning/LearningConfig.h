//
// Created by steffi on 05.01.22.
//

#ifndef STORM_LEARNINGCONFIG_H
#define STORM_LEARNINGCONFIG_H

#include "storm/modelchecker/explorationDT/heuristics/Heuristic.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

enum configTypes {noneConfig, basic, genetic, external, headstart};

template<typename StateType>
class LearningConfig {
   public:
    explicit LearningConfig(std::string const& modelName);

    virtual ~LearningConfig();

    LearningConfig& operator=(const LearningConfig<StateType>& copy);

    /*!
     * Update function of the configuration. Will be called after each simulation.
     */
    virtual void update();

    /*!
     * Get-function for the heuristic stored in the configuration.
     */
    Heuristic<StateType>* getHeuristic() const;

    /*!
     * Asks the heuristic whether it should be used or not.
     * @return
     */
    virtual bool useHeuristic();

    /*!
     * Asks the heuristic what to do next.
     * @param inputValues state observation
     * @param currentStateId current state index
     * @param worked output if we actually found a possible action
     * @return storm-action index
     */
    virtual uint getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked);

    /*!
     * Asks the heuristic for a set of possible actions.
     * @param inputValues state observation
     * @param currentStateId current state index
     * @return vector of storm-action-indices
     */
    virtual std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId);

    void printHeuristic(std::string const& fileName);

    virtual void firstLog();

    virtual void finalLog(double const& elapsed);

    void prepareBoundLog();
    void writeBoundLog(double bound);

    virtual void setUpdateNow();

   protected:
    Heuristic<StateType>* heuristic;
    uint steps=0;

   private:
    std::string modelName;
    std::string learningType;
    std::string heuristicType;
};
}
}
}

#endif  // STORM_LEARNINGCONFIG_H
