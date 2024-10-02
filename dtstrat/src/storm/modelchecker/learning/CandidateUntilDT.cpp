//
// Created by muqsit on 29.09.22.
//

#include "CandidateUntilDT.h"
#include "modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/modelchecker/explorationDT/BSCCCandidate.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {
    template <typename StateType, typename ValueType>
    CandidateUntilDT<StateType,ValueType>::CandidateUntilDT(exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>const& explorationInformationDt, std::vector<std::pair<StateType, ActionType>> const& stack, std::set<StateType> const& statesWithUnavailDtActions, dtstrat::DTInterface<uint, ValueType>* decisionTree):
      explorationInformationDt(&explorationInformationDt) {
        this->stack = stack;
        this->statesWithUnavailDtActions = statesWithUnavailDtActions;
        this->decisionTree = decisionTree;
    }

    template<typename StateType, typename ValueType>
    bool CandidateUntilDT<StateType,ValueType>::verifyCandidate(BSCCCandidate<StateType, ValueType> const& cand,
                                                                 std::set<StateType> const& statesWithUnavailDtActions,
                                                                 exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> const& explorationInformationDt,
                                                                 StateValuationsDT& stateValuations, dtstrat::DTInterface<uint, ValueType>* decisionTree,
                                                                 std::vector<std::vector<std::string>> const& actionIndexToLabelMapping,
                                                                 bool const& semanticActions) {

        // std::set<ActionType> AllSuccOfAllActionsofAllStatesInCandidate;
        // need to look at all the successors under the strategy not for all action
        // Determine all available actions.
//        STORM_PRINT(std::endl << "verifying BSCC Candidate" << std::endl)
        std::set<StateType> candidateStates;
        candidateStates = cand.states;
        for(StateType currentStateId: candidateStates) {
            //StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
            // For every action to consider, find all of its successors.
            const bool isDTActionUnavailable = statesWithUnavailDtActions.find(currentStateId)!=statesWithUnavailDtActions.end();
            if(!isDTActionUnavailable){
            //                STORM_PRINT("Successors from DT-Action leaving the state."<< std::endl);
            //if a successor is outside the candidate, it's not a BSCC
            auto stateValuation = stateValuations.getStateValuation(currentStateId);
            bool worked = false;
            ActionType chosenAction = decisionTree->getDecision(stateValuation, currentStateId, worked);

            if (worked & semanticActions) {
                // get all action-names for this action
                std::vector<std::string> labels = actionIndexToLabelMapping[chosenAction];
                std::vector<uint> possibleActions;
                StateType rowGroup = explorationInformationDt.getRowGroup(currentStateId);
                auto totalNumberOfActionsInState =
                    explorationInformationDt.getStartRowOfGroup(rowGroup + 1) - explorationInformationDt.getStartRowOfGroup(rowGroup);
                std::string finalLabel;
                for (uint32_t action = explorationInformationDt.getStartRowOfGroup(rowGroup); action < explorationInformationDt.getStartRowOfGroup(rowGroup + 1);
                     ++action) {
                    auto labelOfAction = explorationInformationDt.getLabel(action);
                    if (std::find(labels.begin(), labels.end(),labelOfAction)!=labels.end()) {
                        possibleActions.push_back(action-explorationInformationDt.getStartRowOfGroup(rowGroup));
                        if (finalLabel.empty()) finalLabel = labelOfAction;
                    }
                }
                if (possibleActions.empty()) {
                    worked = false;
                } else {
                    // currently just play the first action
                    // todo decide smartly
                    chosenAction = possibleActions[0];
                    //STORM_PRINT(" = " << chosenAction << " ("<< finalLabel << ")\n");
                }
            }

            StateType rowGroup = explorationInformationDt.getRowGroup(currentStateId);
            auto initialAction = explorationInformationDt.getStartRowOfGroup(rowGroup);
            auto lastPossibleAction = explorationInformationDt.getStartRowOfGroup(rowGroup+1);
            std::vector<storm::storage::MatrixEntry<StateType, ValueType>> row;
            //                STORM_PRINT(std::endl<< "Successors from DT-Action leaving the state.");
//            STORM_PRINT("CHOSEN ACTION" <<chosenAction <<std::endl);
            //STORM_PRINT("Count of Available Actions" << lastPossibleAction - initialAction <<std::endl);
            STORM_LOG_ASSERT(lastPossibleAction > initialAction+chosenAction, "ChosenAction is Unavailable in the State");
            row = explorationInformationDt.getRowOfMatrix(initialAction + chosenAction);

            //                auto startRow = explorationInformationDt.getRowOfMatrix(initialAction);
            //                auto nxtRow = explorationInformationDt.getRowOfMatrix(initialAction+1);
//            row = explorationInformationDt.getRowOfMatrix(initialAction + chosenAction);

            //                STORM_PRINT(std::endl<< "Successors from DT-Action leaving the state.");
            std::vector<StateType> successorsForChosenAction;
            //                auto succRow = explorationInformationDt.getRowOfMatrix(chosenAction);
            for (auto v: row) {
                successorsForChosenAction.push_back(v.getColumn());
            }
            for(StateType succ: successorsForChosenAction){
                const bool isin= candidateStates.find(succ)!=candidateStates.end();
                if(!isin){
                    // if new state ret false
                    //                    if(std::count(candidateStates.begin(), candidateStates.end(), succ)){
                    //                        STORM_PRINT("Candidate States ret false :")
                    //                        for(StateType currentStateId: candidateStates){
                    //                            STORM_PRINT(currentStateId<<",");
                    //                        }
                    //                        STORM_PRINT(std::endl<<"a DT successor (getColumn) from the row ret false: " << succ << std::endl);
                    return false;
                }
                //                    else{
                //                        STORM_PRINT("Candidate States:")
                //                        for(StateType currentStateId: candidateStates){
                //                            STORM_PRINT(currentStateId<<",");
                //                        }
                //                        STORM_PRINT(std::endl<<"a DT successor (getColumn) from the row: " << succ << std::endl);
                //                    }
            }
            } else{
                //                STORM_PRINT("Successors from a random leaving the state.");
                StateType rowGroup = explorationInformationDt.getRowGroup(currentStateId);
                for (uint32_t action = explorationInformationDt.getStartRowOfGroup(rowGroup); action < explorationInformationDt.getStartRowOfGroup(rowGroup + 1); ++action) {
                    std::vector<storm::storage::MatrixEntry<StateType, ValueType>> const& row = explorationInformationDt.getRowOfMatrix(action);
                    for(auto entry: row){
                        //                            STORM_PRINT("a successor (getColumn) from the row " << entry.getColumn());
                        const bool isin = candidateStates.find(entry.getColumn())!=candidateStates.end();
                        if(!isin){
                            return false;
                        }
                    }
                }

            }
        }
        return true;
    }


    template<typename StateType, typename ValueType>
    BSCCCandidate<StateType, ValueType> CandidateUntilDT<StateType,ValueType>::updateCandidate(std::vector<std::pair<StateType, ActionType>> const& path)
    {
        //for each path there is at most one candidate
        BSCCCandidate<StateType, ValueType> candidate;

        if(cndLst.size() > 0)
            candidate=cndLst[cndLst.size()-1];

        for (int i = lastIndx + 1; i < path.size(); i++) {
            StateType st = path.at(i).first;

            if (candidate.states.size()!=0 && candidate.contains(st)) {
                // state already in the current candidate
                candidate.countState(st);
                candidate.setNonTrivial();
            } else {
                // find the index of the candidate,
                // if this state is already in a candidate
                if (cndStateIndexMap.find(st)== cndStateIndexMap.end()) {
                    // state seen for the first time, and,
                    // i is the path index where candidate was born
                    candidate = BSCCCandidate<StateType, ValueType>(st, i);
                    cndLst.push_back(candidate);

                    // remembering the index in the candidate list
                    cndStateIndexMap[st] = cndLst.size()-1;
                } else {
                    // state seen before - create a merge of all candidates from that point
                    int indx = cndStateIndexMap[st];
                    candidate = BSCCCandidate<StateType, ValueType>(st, i);
                    while (cndLst.size() > indx) {
                        BSCCCandidate<StateType, ValueType> oldCnd = cndLst.back();
                        cndLst.pop_back();
                        candidate.mergeWith(oldCnd);
                    }

                    cndLst.push_back(candidate);
                    int v = cndLst.size() - 1;

                    for (auto s : candidate.getStates()) {
                        // update the indices according to candidate list
                        cndStateIndexMap[s]= v;
                    }
                    candidate.setNonTrivial();
                }
            }
        }
        return candidate;
    }

    template <typename StateType, typename ValueType>
    bool CandidateUntilDT<StateType,ValueType>::update(exploration_detail_dt::ExplorationInformationDT<StateType, ValueType> const& explorationInformationDt,
                                                        std::vector<std::pair<StateType, ActionType>> const& stack,
                                                        std::set<StateType> const& statesWithUnavailDtActions,
                                                        StateValuationsDT& stateValuations, std::vector<std::vector<std::string>> const& actionIndexToLabelMapping,
                                                        bool const& semanticActions)
    {
        transLen = stack.size();
        BSCCCandidate<StateType, ValueType> candidate = updateCandidate(stack);

//        if (candidate.!=lastNonTrivialCnd && !candidate.isTrivial()) {
//            // a new non-trivial candidate is found
//            lastNonTrivialCnd = candidate;
//        }

        if (!candidate.isTrivial() && !candidate.isIncomplete()) {
            if (verifyCandidate(candidate, statesWithUnavailDtActions, explorationInformationDt, stateValuations, decisionTree, actionIndexToLabelMapping, semanticActions)) {
                // we hit a BSCC, without reaching the target state
                return true;
            } else {
                //
                candidate.setIncomplete();
            }

            transLen = candidate.firstentry;
            candLen = stack.size() - transLen;
        }
        lastIndx = stack.size()-1;
        return false;
    }

    template class CandidateUntilDT<uint32_t, double>;
}  // namespace dtstrat
}  // namespace modelchecker
}  // namespace storm