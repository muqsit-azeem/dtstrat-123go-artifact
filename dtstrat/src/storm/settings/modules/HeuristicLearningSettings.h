//
// Created by steffi on 05.01.22.
//

#ifndef STORM_HEURISTICLEARNINGSETTINGS_H
#define STORM_HEURISTICLEARNINGSETTINGS_H
#include "storm-config.h"
#include "storm/settings/modules/ModuleSettings.h"

#include "storm/builder/ExplorationOrder.h"

#include "storm/modelchecker/explorationDT/learning/LearningConfig.h"
#include "storm/modelchecker/explorationDT/heuristics/Heuristic.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleHeuristic.h"

namespace storm {
namespace settings {
namespace modules {

/*!
 * This class represents the general settings.
 */
class HeuristicLearningSettings : public ModuleSettings {
   public:
    /*!
     * Creates a new set of general settings.
     */
    HeuristicLearningSettings();

    /*!
     * Retrieves the name of the configuration to be used.
     *
     * @return The name of the configuration
     */
    storm::modelchecker::dtstrat::configTypes getLearningConfigType() const;

    /*!
     * Retrieves the number of simulations after which an update of the config should be performed.
     *
     * @return Number of simulations
     */
     int getNumberOfSimulationsBeforeUpdate() const;

     /*!
      * Retrieves the probability of following the heuristic instead of BRTDP.
      *
      * @return probability of following the heuristic
      */
      double getHeuristicProbability() const;

      /*!
       * Retrieves the type of the heuristic
       */
       storm::modelchecker::dtstrat::heuristicTypes getHeuristicType() const;

       /*!
        * Retrieves whether to use only reachable states under the current strategy for export
        */
       bool getOnlyReachable() const;

       /*!
        * Retrieves the type of heuristic for the module
        */
        storm::modelchecker::dtstrat::moduleHeuristicTypes getModuleHeuristicType() const;

    /*!
     * Retrieve the filename for an external config
     */
     std::string getExternalFilename() const;

    // The name of the module.
    static const std::string moduleName;

   private:
    // Define the string names of the options as constants.
    static const std::string learningConfig;
    static const std::string explorationEpsilon;
    static const std::string updateSteps;
    static const std::string heuristicType;
    static const std::string onlyReachable;
    static const std::string moduleHeuristc;
    static const std::string externalFile;
};

}  // namespace modules
}  // namespace settings
}  // namespace storm
#endif  // STORM_HEURISTICLEARNINGSETTINGS_H
