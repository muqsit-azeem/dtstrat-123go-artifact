//
// Created by steffi on 21.03.22.
//

#ifndef STORM_STATEVALUATIONDT_H
#define STORM_STATEVALUATIONDT_H
#include<vector>
#include "storm/storage/sparse/StateValuations.h"
#include <mlpack/core.hpp>
#include <armadillo>

namespace storm {
namespace modelchecker {
namespace dtstrat {

    class StateValuationDT {
    public:
        StateValuationDT()=default;
        explicit StateValuationDT(std::vector<int> const& vals);
        explicit StateValuationDT(uint const& size);
        explicit StateValuationDT(const storage::sparse::StateValuations::StateValuation & inputValues);

        int& operator[] (uint const& index);

        const int& operator[] (uint const& index) const;

        uint size() const;

        arma::Col<int> transform();

        arma::Col<int> transform() const;

        friend std::ostream& operator<<(std::ostream& os, const StateValuationDT& dt) {
            for (uint i=0; i<dt.size(); i++) os << dt[i] << ",";
            return os;
        }

    private:
        std::vector<int> values;
        uint sizeOfV;
    };

    class StateValuationsDT {

    public:
        StateValuationsDT() = default;
        void addStateValuation(uint const& currentStateId, storage::sparse::StateValuations::StateValuation const& inputValues);
        StateValuationDT getStateValuation(uint const& currentStateId);
        uint getNumberOfStates();
        bool isEmpty(uint const& stateId);

        arma::Mat<int> transform();

    private:
        std::vector <StateValuationDT> valuations;

    };
}
}
}


#endif //STORM_STATEVALUATIONDT_H
