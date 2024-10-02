//
// Created by muqsit on 29.09.22.
//
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include "DTInterface.h"
#include "modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/modelchecker/explorationDT/BSCCCandidate.h"

#ifndef STORM_CANDIDATEUNTILDT_H
#define STORM_CANDIDATEUNTILDT_H



namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class CandidateUntilDT {
    typedef StateType ActionType;

   private:
    // CHECKING CANDIDATES
    // last inspected index	(-1 if none)
    int lastIndx=-1;
    //private Candidate candidate;
    BSCCCandidate<StateType, ValueType> lastNonTrivialCnd;
    // sequence of non-subsumed candidates on the path
    std::vector<BSCCCandidate<StateType,ValueType>> cndLst;
    // state -> index in cndLst
    std::map<StateType,int> cndStateIndexMap;


    //Statistics
    int candLen;
    int transLen;
    int cndFound;
    int largestCnd;

  public:
   const exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformationDt;
   std::vector<std::pair<StateType, ActionType>> stack;
   std::set<StateType> statesWithUnavailDtActions;
   // DT for a strategy
   dtstrat::DTInterface<uint, ValueType>* decisionTree;


   CandidateUntilDT(exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>const& explorationInformationDt, std::vector<std::pair<StateType, ActionType>> const& stack, std::set<StateType> const& statesWithUnavailDtActions, dtstrat::DTInterface<uint, ValueType>* decisionTree);
//    void resetStats();
   bool verifyCandidate(BSCCCandidate<StateType, ValueType> const& cand,
                        std::set<StateType> const& statesWithUnavailDtActions,
                        exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> const& explorationInformationDt,
                        StateValuationsDT& stateValuations, dtstrat::DTInterface<uint, ValueType>* decisionTree,
                        std::vector<std::vector<std::string>> const& actionIndexToLabelMapping,
                        bool const& semanticActions);
   BSCCCandidate<StateType, ValueType> updateCandidate(std::vector<std::pair<StateType, ActionType>> const& stack);
   bool update(exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> const& explorationInformationDt,
               std::vector<std::pair<StateType, ActionType>> const& stack,
               std::set<StateType> const& statesWithUnavailDtActions,
               StateValuationsDT& stateValuations, std::vector<std::vector<std::string>> const& actionIndexToLabelMapping,
               bool const& semanticActions);
};

}  // namespace dtstrat
}  // namespace modelchecker
}  // namespace storm

#endif  // STORM_CANDIDATEUNTILDT_H
