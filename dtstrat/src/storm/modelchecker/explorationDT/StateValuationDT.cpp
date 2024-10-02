//
// Created by steffi on 21.03.22.
//

#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include "storm/utility/macros.h"

namespace storm {
    namespace modelchecker {
        namespace dtstrat {
            StateValuationDT::StateValuationDT(std::vector<int> const& vals) : values(vals), sizeOfV(vals.size()) {
                // intentionally left empty
            }

            StateValuationDT::StateValuationDT(uint const& size) : values(std::vector<int>(size)), sizeOfV(size) {
                // intentionally left empty
            }

            int &StateValuationDT::operator[](uint const& index) {
                STORM_LOG_ASSERT(index<=sizeOfV,"index is out of bounds");
                return values[index];
            }

            const int &StateValuationDT::operator[](uint const& index) const {
                STORM_LOG_ASSERT(index<=sizeOfV,"index is out of bounds");
                return values[index];
            }

            StateValuationDT::StateValuationDT(const storage::sparse::StateValuations::StateValuation &inputValues) {
                sizeOfV = inputValues.size();
                values.resize(sizeOfV);
                for (uint i=0; i<sizeOfV; i++) values[i] = inputValues[i];
            }

            uint StateValuationDT::size() const {
                return values.size();
            }

            arma::Col<int> StateValuationDT::transform() {
                return values;
            }

            arma::Col<int> StateValuationDT::transform() const {
                return values;
            }

            void StateValuationsDT::addStateValuation(const uint &currentStateId,
                                                      const storage::sparse::StateValuations::StateValuation &inputValues) {
                if (currentStateId>=valuations.size()){
                    valuations.resize(currentStateId);
                    valuations.emplace_back(inputValues);
                }
                else {
                    valuations[currentStateId] = StateValuationDT(inputValues);
                }
            }

            StateValuationDT StateValuationsDT::getStateValuation(uint const& currentStateId) {
                if (currentStateId>=valuations.size())
                    return StateValuationDT(0);
                else {
                    return valuations[currentStateId];
                }
            }

            uint StateValuationsDT::getNumberOfStates() {
                return valuations.size();
            }

            bool StateValuationsDT::isEmpty(uint const& stateId) {
                return valuations[stateId].size()==0;
            }

            arma::Mat<int> StateValuationsDT::transform() {
                arma::Mat<int> result;
                for (uint i=0; i<valuations.size(); i++){
                    result.insert_cols(i, valuations[i].transform());
                }
                return result;
            }
        }
    }
}
