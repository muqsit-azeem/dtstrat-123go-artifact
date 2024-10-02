//
// Created by muqsit on 27.09.22.
//
#include <set>
#include <map>
#ifndef STORM_BSCCCANDIDATE_H
#define STORM_BSCCCANDIDATE_H

namespace storm{
namespace modelchecker{
namespace dtstrat{

    template<typename StateType, typename ValueType>
    class BSCCCandidate{
       private:
        // path index where the candidate was born
        int birth;
        // the candidate is known not be the bottom
        bool incomplete = false;
        //the candidate is not an SCC
        bool trivial = true;
        // how many times a state occurs after the birth, null means that all entries are zero
        std::map<StateType,int> counter;


       public:
        //states in the candidate
        std::set<StateType> states;
        // the point where the state appeared first
        int firstentry;
        BSCCCandidate()=default;
        explicit BSCCCandidate(StateType st, int idx);

        // Count occurrence of candidate's state.
       void countState(StateType st);

       bool isTrivial();

       void setNonTrivial();

       int getCount(StateType st);

       bool contains(StateType st) const;

       std::set<StateType> getStates();

       int getBirthIdx();

        //Merge with another candidate.
       void mergeWith(BSCCCandidate old);


       // Returns true iff the candidate is k-strong.
       bool isStrong(int k);

       void setIncomplete();

       bool isIncomplete();

    };
}
}
}


#endif  // STORM_BSCCCANDIDATE_H
