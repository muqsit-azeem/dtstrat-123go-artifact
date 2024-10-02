//
// Created by steffi on 27.07.22.
//

#ifndef STORM_LEARNINGMODELCHECKERSETTINGS_H
#define STORM_LEARNINGMODELCHECKERSETTINGS_H
#include "storm-config.h"
#include "storm/settings/modules/ModuleSettings.h"

namespace storm {
namespace settings {
namespace modules {

/*!
 * This class represents the general settings.
 */
class LearningModelCheckerSettings : public ModuleSettings {
   public:
    /*!
     * Creates a new set of general settings.
     */
    LearningModelCheckerSettings();

    std::string getDecisionTreeFile() const;

    std::string getDataFile() const;

    std::vector<int> getLearnParameter() const;

    std::string getEvaluationMethod() const;

    std::vector<std::string> getChildrenFiles() const;

    std::string getDTLearnAlgorithm() const;

    std::string getParameterName() const;

    uint64_t getNumberOfSMCRuns() const;

    uint64_t getMaxStepsForSMC() const;

    bool isOnlyReachableSet() const;

    uint64_t getNumDTStrat() const;

    std::vector<std::pair<std::string, std::string>> getSubTasksAndProperties() const;

    double getPACEpsilon() const;

    double getPACDelta() const;

    std::string getLearnHyperParamType() const;

    uint64_t getHyperparamLearnTime() const;

    bool useSemanticActions() const;

    // The name of the module.
    static const std::string moduleName;

   private:
    // Define the string names of the options as constants.
    static const std::string decisiontree;
    static const std::string datafile;
    static const std::string learnparameter;
    static const std::string evluationmethod;
    static const std::string childrenfiles;
    static const std::string treelearner;
    static const std::string learnparametername;
    static const std::string smcruns;
    static const std::string maxsteps;
    static const std::string onlyReachable;
    static const std::string numStepsDTStrat;
    static const std::string subProgramsAndTasks;
    static const std::string pacEpsilon;
    static const std::string pacDelta;
    static const std::string hyperparamType;
    static const std::string hyperparamTime;
    static const std::string semanticActions;
};

}  // namespace modules
}  // namespace settings
}  // namespace storm

#endif  // STORM_LEARNINGMODELCHECKERSETTINGS_H
