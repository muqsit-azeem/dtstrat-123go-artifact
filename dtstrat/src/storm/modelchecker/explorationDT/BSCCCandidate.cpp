//
// Created by muqsit on 27.09.22.
//

#include "storm/modelchecker/explorationDT/BSCCCandidate.h"
#include <algorithm>
#include <vector>

namespace storm{
namespace modelchecker{
namespace dtstrat{
    // template<typename StateType, typename ValueType>
    //    BSCCCandidate<StateType, ValueType>::BSCCCandidate() {
    //        //intentionally left blank
    //    }
    template<typename StateType, typename ValueType>
    BSCCCandidate<StateType, ValueType>::BSCCCandidate(StateType st, int idx) {
        birth = idx;
        states.insert(st);
        firstentry = idx;
    }

    template<typename StateType, typename ValueType>
    void BSCCCandidate<StateType,ValueType>::countState(StateType st)
    {
        auto it = counter.find(st);
        if (it == counter.end()) {
            counter.insert({st, 0});
        } else{
            it->second = it->second+1;
        }
    }

    template<typename StateType, typename ValueType>
    bool BSCCCandidate<StateType, ValueType>::isTrivial()
    {
        return trivial;
    }

    template<typename StateType, typename ValueType>
    void BSCCCandidate<StateType, ValueType>::setNonTrivial()
    {
        trivial = false;
    }

    template<typename StateType, typename ValueType>
    int BSCCCandidate<StateType, ValueType>::getCount(StateType st)
    {
        auto it = counter.find(st);
        if (it != counter.end()) {
            return it->second;
        }
        return -1;
    }


    template<typename StateType, typename ValueType>
    bool BSCCCandidate<StateType, ValueType>::contains(StateType st) const
    {
        return states.find(st) != states.end();
    }

    template<typename StateType, typename ValueType>
    std::set<StateType> BSCCCandidate<StateType, ValueType>::getStates()
    {
        return states;
    }

    template<typename StateType, typename ValueType>
    int BSCCCandidate<StateType, ValueType>::getBirthIdx()
    {
        return birth;
    }


    //Merge with another candidate.
    template<typename StateType, typename ValueType>
    void BSCCCandidate<StateType, ValueType>::mergeWith(BSCCCandidate old)
    {

        std::set_union(old.states.begin(), old.states.end(), states.begin(), states.end(),
                       std::inserter(states, states.begin()));
    }

    // Returns true iff the candidate is k-strong.
    template<typename StateType, typename ValueType>
    bool BSCCCandidate<StateType, ValueType>::isStrong(int k)
    {
        if (counter.empty()) {
            return false;
        }

        for (int st : states) {
            if (counter.find(st)->second < k) {
                return false;
            }
        }

        return true;
    }

    template<typename StateType, typename ValueType>
    void BSCCCandidate<StateType, ValueType>::setIncomplete()
    {
        incomplete = true;
    }

    template<typename StateType, typename ValueType>
    bool BSCCCandidate<StateType, ValueType>::isIncomplete()
    {
        return incomplete;
    }

    template class BSCCCandidate<u_int32_t, double>;
}
}
}

