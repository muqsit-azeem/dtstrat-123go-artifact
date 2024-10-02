#ifndef UPDATE_DT_H
#define UPDATE_DT_H

#include "storm/dtstrat/decisionTree.h"
#include <vector>
#include <tuple>
#include <math.h> 
#include <iostream>
#include <string>
#include <random>
//
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include "storm/dtstrat/json.hpp"
//
#include "storm/modelchecker/exploration/ExplorationInformation.h"
#include "storm/modelchecker/exploration/Bounds.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/settings/modules/ExplorationSettings.h"
#include "storm/settings/SettingsManager.h"
#include "storm/storage/sparse/StateValuations.h"
#include "storm/dtstrat/writeToFile.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "storm/modelchecker/exploration/SparseExplorationModelChecker.h"
#include "storm/storage/prism/Program.h"



namespace dtstrat {
    namespace update {
        enum moduleChoice { random, treeSize, smallest };
        // A class that remmebers all important hyperparameter for the updating
        template<typename ModelType, typename ValueType, typename StateType = uint32_t>
        class UpdaterDT {
            private:
                //if storm is to  be run in parallel
                bool parallel;
                // learning rate
                double alpha;
                // discount factor, always 1 in reachability setting
                double gamma;
                // P&H: how often does a node to be visited before we want to update it
                uint historyMin;
                // probability to determine how likely we want to do "random exploration" (=BRTDP)
                double epsilon;
                // stores if the DTStrat-algorithm has to be used at all
                bool dtIsUsed;
                // stores if the DTStrat2-algorithm has to be used at all
                bool dt2IsUsed;
                // determines how often the update should be performed
                uint updateSteps;
                // determines which calculation of QValues to use
                bool useTrueQValue;

                // lookup-Table
                boost::numeric::ublas::mapped_matrix<int> lookUpTable;

                // target states
                std::vector<dtstrat::dt::Observation<ValueType>> targetStates;
                // batchwise update
                bool updateBatchwise;

                void addTargetState(dtstrat::dt::Observation<ValueType> targetState) {
                    bool found = false;
                    for (dtstrat::dt::Observation<ValueType> state : this->targetStates) {
                        found = (targetState == state);
                    }
                    if (!found) {
                        this->targetStates.push_back(targetState);
                    }
                }

                moduleChoice moduleMethodChoice;

            public:
            typedef StateType ActionType;
                // creates an instance
                UpdaterDT(bool dtIsUse = false, bool dt2IsUse = false, double alph=0.9, double gamm=0.9, uint histmin=500, double epsi=0.5, uint updateStep=0, bool uTrueQValue=false, std::string splitDec="tstatistic", bool updateBatchwis=false, bool deleteTree=false, bool parallel=false, bool printModules=false) {
                    this->alpha = alph;
                    this->gamma = gamm;
                    this->historyMin = histmin;
                    this->epsilon = epsi;
                    this->dtIsUsed = dtIsUse;
                    this->dt2IsUsed = dt2IsUse;
                    this->updateSteps = updateStep;
                    this->useTrueQValue = uTrueQValue;
                    this->splitDecision = splitDec;
                    this->updateBatchwise = updateBatchwis;
                    this->updateCounter = 0;
                    this->deleteTree = deleteTree;
                    this->useDT2Counter = 0;
                    this->parallel = parallel;
                    this->printModules = printModules;
                }

                // counts how often no update was performed
                int updateCounter;
                // counts after how many simulations start using DT2
                uint useDT2Counter;
                // stores the property as a string
                std::string property = "";
                // stores the model name as a string
                std::string modelName = "";
                // action id mapping from dtstrat to storm
                std::vector<int> stormDtstratActionMapping;
                // mapping from stateIds to chosen action (fixed)- gets updated after k steps
                std::map<StateType, ActionType> stateIdChosenActionMapping;
                // State valuations
                dt::StateObservations<ValueType> allStateObservations;
                // updateBatch
                std::vector<std::tuple<dtstrat::dt::Observation<ValueType>,int,int,bool>> updateBatch;
                // for brtdp update
                std::vector<uint> updateBatchStormStates;
                // which statistic to use
                std::string splitDecision;
                // dt should be deleted and then created new from scratch
                bool deleteTree = false;
                //std::string jsonStringOfDT = "";
                dtstrat::dt::DecisionTree<ValueType, int> jsonDT;
                dtstrat::dt::DecisionTree<ValueType, int> jsonDT2;
                std::vector<dtstrat::dt::DecisionTree<ValueType,int>> multiDts;
                // debugging option
                bool debug = false;
                bool printStates = false;

                std::vector<std::vector<std::pair<StateType, std::vector<std::pair<StateType,int>>>>> mecActions;
                std::vector<int> stormActionMecListMapping;
                bool printModules;

                bool isLookUpSet() {
                    return (this->lookUpTable.size1()>0);
                }

            // get functions for the values
                bool getUseDT() const {
                    return this->dtIsUsed;
                }

                bool getUseDT2() const {
                    return this->dt2IsUsed;
                }

                double getEpsilon() const {
                    return this->epsilon;
                }

                double getAlpha() const {
                    return this->alpha;
                }

                double getGamma() const {
                    return this->gamma;
                }

                uint getHistoryMin() const {
                    return this->historyMin;
                }

                uint getUpdateSteps() const {
                    return this->updateSteps;
                }

                bool getUseTrueQ() const {
                    return this->useTrueQValue;
                }

                bool getUpdateBatchwise() const {
                    return this->updateBatchwise;
                }

                bool isParallel() const {
                    return this->parallel;
                }

                bool isPrintModules() const {
                    return this->printModules;
                }

                // calculate the mean of a list of float-values
                // used to calculate the mean of changes in the qvalues
                float mean(std::vector<float> history) const {
                    float m = 0.0;
                    for (float val : history) {
                        m += val;
                    }
                    return m/history.size();
                }

                // calculate the standard deviation of a list of float-values
                // used to calculate the standard deviation of changes in the qvalues
                float stdev(std::vector<float> history) const {
                    float me = mean(history);
                    float s = 0;
                    for (float val : history) {
                        s += (val-me)*(val-me);
                    }
                    s /= history.size();
                    return sqrt(s);
                }

                // calculate the mean of a list of Observations
                // used to calculate the mean of the dimensions of the observations
                float* mean_observation(std::vector<dtstrat::dt::Observation<ValueType>> hist) const {
                    uint size = hist[0].size;
                    float* m = new float[size];
                    for (uint i=0; i<size; i++) {m[i] = 0.0;}
                    for (auto val : hist) {
                        for (uint i = 0; i<size; i++) {
                            m[i] += (val)[i];
                        }
                    }
                    for (uint i = 0; i<size; i++) {
                        m[i] /= hist.size();
                    }
                    return m;
                }

                // calculate the standard deviation of a list of Observations
                // used to calculate the standard deviation of the dimensions of the observations
                float* stdev_observation(std::vector<dtstrat::dt::Observation<ValueType>> his, float* meanVals=nullptr, float* minVals=nullptr, float* maxVals=nullptr) const {
                    uint size = his[0].size;
                    float* s = new float[size];
                    for (uint i=0; i<size; i++) {
                        s[i] = 0.0;
                        if (minVals) {
                            minVals[i] = std::numeric_limits<float>::max();
                        }
                        if (maxVals) {
                            maxVals[i] = std::numeric_limits<float>::min();
                        }
                    }
                    float* m = mean_observation(his);
                    for (auto val : his) {
                        for (uint i = 0; i<size; i++) {
                            s[i] += ((val)[i]-m[i])*((val)[i]-m[i]);
                            if (minVals) {
                                if (minVals[i]>val[i]) minVals[i] = val[i];
                            }
                            if (maxVals) {
                                if (maxVals[i]<val[i]) maxVals[i] = val[i];
                            }
                        }
                    }
                    for (uint i = 0; i<size; i++) {
                        s[i] /= his.size();
                        s[i] = sqrt(s[i]);
                    }
                    if (meanVals) {
                         std::memcpy(meanVals, m, size*sizeof(float));
                    }
                    return s;
                }

                // decide if a split should be performed in a node
                bool decide_split(dtstrat::dt::Node<ValueType, float>& node) const {
                    // check if the node was visited at least history_min times
                    if (node.history.size()<this->historyMin) {
                        return false;
                    }
                    else {
                        uint size = node.observation_history[0].size;
                        // calculate mu and stdev of the changes of the qvalues
                        float mu = mean(node.history);
                        float sigma = stdev(node.history);
                        if (abs(mu)<2*sigma) {
                            // the splitting criterion by Pyeatt&Howe is fullfilled
                            bool all_zero = true;
                            float* std_obs = stdev_observation(node.observation_history);
                            for (uint i = 0; i<size; i++) {
                                if (std_obs[i] != 0.0) {
                                    all_zero = false;
                                    break;
                                }
                            }
                            // if the observations are all zero, do not perform a split
                            // NOTE: Learning from experience, the algorithm tended to split a node, even if there was no reason for it
                            if (all_zero) { 
                                //std::cout << " All zero -> no split "<< std::endl;
                                return false;
                            }
                            else {
                                return true;
                            }
                        }
                        else{
                            //std::cout << " abs(mu) " << abs(mu) << " >= 2*sigma " << 2*sigma << " -> no split" << std::endl;
                            return false;
                        }
                    }
                }

                // checks whether there are no positive values in the list
                bool no_positive(std::vector<float> inp) const {
                    bool alzero = false;
                    for (auto val : inp) {
                        if (val < 0) {
                            alzero = true;
                            break;
                        }
                    }
                    return alzero;
                }

                // calculates the t-statistic
                // NOTE: I am not sure if this is correct. The T-statstic is a bit hard to understand in 2min
                float* t_statistic(std::vector<dtstrat::dt::Observation<ValueType>>& inp) const {
                    float* mean = mean_observation(inp);
                    float* stdv = stdev_observation(inp);
                    int n = inp.size();
                    uint size = inp[0].size;

                    float* res = new float[size];
                    for (uint i=0; i<size; i++) {
                        res[i] = mean[i] / (stdv[i] / sqrt(n));
                    }
                    return res;
                }

                /*
                std::vector<dtstrat::dt::HistoryElement> goodChoicesElements (const dtstrat::dt::Node<ValueType>& currentNode) {
                    std::vector<dtstrat::dt::HistoryElement> result;
                    for (dtstrat::dt::HistoryElement element : currentNode.complete_history) {
                        if (element.deltaQ>=0 && element.qValue>0) {
                            result.push_back(element);
                        }
                    }
                    return resuconstlt;
                }*/

                std::vector<uint> getGoodActions(const dtstrat::dt::Node<ValueType, float>& currentNode, std::vector<std::vector<dtstrat::dt::Observation<ValueType>>>& observationPerAction) const {
                    // only get good choices
                    std::vector<uint> goodActions;
                    for (dtstrat::dt::HistoryElement<ValueType> element : currentNode.complete_history) {
                        if (element.deltaQ>=0 && element.qValue>0) {
                            auto it = std::find(goodActions.begin(), goodActions.end(), element.action); 
                            uint index = 0;
                            if (it == goodActions.end()) {
                                goodActions.push_back(element.action);
                                index = goodActions.size()-1;
                            } else {
                                index =  std::distance(goodActions.begin(), it);
                            }
                            if (observationPerAction.size()<goodActions.size()) {
                                observationPerAction.resize(goodActions.size());
                            }
                            observationPerAction[index].push_back(element.state);
                        }
                    }
                    return goodActions;
                }

                float giniVal(float num_total, float num_left, float num_right, float num_good_left, float num_good_right) const {
                    float giniValue = (num_left/num_total) * (1 - ((num_good_left/num_left)*(num_good_left/num_left)) - ((1-num_good_left/num_left) * (1-num_good_left/num_left)))
                                        + (num_right/num_total) * (1 - ((num_good_right/num_right)*(num_good_right/num_right)) - ( (1-num_good_right/num_right)*(1-num_good_right/num_right) ));
                    return giniValue;
                }

                float gini(uint action, dtstrat::dt::Node<ValueType, float>& currentNode, uint dimension, float mean) const {
                    uint num_left = 0;
                    uint num_right = 0;
                    uint num_good_left = 0;
                    uint num_good_right = 0;
                    for (dtstrat::dt::HistoryElement<ValueType> element : currentNode.complete_history) {
                        // decide if element is left or right
                        if (element.state[dimension]>=mean) {
                            //right
                            num_right++;
                            if (element.action == action) {
                                num_good_right++;
                            }
                        }
                        else {
                            //left
                            num_left++;
                            if (element.action == action) {
                                num_good_left++;
                            }
                        }
                    }
                    uint num_total = currentNode.complete_history.size();
                    if (num_total == 0 || num_left == 0 || num_right == 0) {
                        return 1;
                    }
                    else {
                        float giniValue = giniVal((float) num_total, (float) num_left, (float) num_right, (float) num_good_left, (float) num_good_right);
                        return giniValue;     
                    }              
                }

                std::tuple<uint, float> getSplitBasedOnGoodActionAndGini(dtstrat::dt::Node<ValueType, float>& currentNode) const {
                    // get all good actions
                    std::vector<std::vector<dtstrat::dt::Observation<ValueType>>> observationPerAction;
                    std::vector<uint> goodActions = getGoodActions(currentNode, observationPerAction);

                    float bestGini = 1;
                    //int bestAction = -1;
                    uint bestDimension = 0;
                    float bestSplitValue = 0.0;
                    
                    // for each action in the good choices
                    if (goodActions.size()>0) {
                        uint sizeStateObservation = observationPerAction[0][0].size;
                        float* means = new float[sizeStateObservation];
                        float* deviation = stdev_observation(currentNode.observation_history, means);
                        for (uint i = 0; i<goodActions.size(); i++) {
                            //STORM_PRINT("/SPLIT ON "<<goodActions[i] <<std::endl);
                            // for each dimension to split on
                            // check the variance (if == 0, then you can't split), if >0, split at mean and calculate Gini or IG1
                            for (uint dimension = 0; dimension<sizeStateObservation; dimension++){  
                                if (deviation[dimension]>0) {
                                    float giniIndex = gini(goodActions[i], currentNode, dimension, means[dimension]);
                                    //STORM_PRINT("  Splitting on " << dimension << " at " << means[dimension] << " has gini " << giniIndex << std::endl);
                                    if (giniIndex < bestGini) {
                                        bestGini = giniIndex;
                                        //bestAction = goodActions[i];
                                        bestDimension = dimension;
                                        bestSplitValue = means[dimension];
                                    }
                                }
                            }
                        }
                    }
                    return std::make_tuple(bestDimension, bestSplitValue);
                }

                std::tuple<uint, float> getSplitBasedOnGoodActionAndGiniExhaustive(dtstrat::dt::Node<ValueType, float>& currentNode) const {
                    // get all good actions
                    std::vector<std::vector<dtstrat::dt::Observation<ValueType>>> observationPerAction;
                    std::vector<uint> goodActions = getGoodActions(currentNode, observationPerAction);

                    float bestGini = 1;
                    //int bestAction = -1;
                    uint bestDimension = 0;
                    float bestSplitValue = 0.0;
                    
                    // for each action in the good choices
                    if (goodActions.size()>0) {
                        uint sizeStateObservation = observationPerAction[0][0].size;
                        float* minVals = new float[sizeStateObservation];
                        float* maxVals = new float[sizeStateObservation];
                        float* deviation = stdev_observation(currentNode.observation_history, nullptr, minVals, maxVals);
                        for (uint i = 0; i<goodActions.size(); i++) {
                            //STORM_PRINT("/SPLIT ON "<<goodActions[i] <<std::endl);
                            // for each dimension to split on
                            // check the variance (if == 0, then you can't split), if >0, try to split on every possible splitvalue
                            for (uint dimension = 0; dimension<sizeStateObservation; dimension++){  
                                if (deviation[dimension]>0) {
                                    std::vector<int> possibleSplitValues((int)(maxVals[dimension]-minVals[dimension])+1);
                                    // generates a range in possibleSplitValues, starting from (int)minVals[dimension] and ending when the vector ends
                                    std::iota(std::begin(possibleSplitValues), std::end(possibleSplitValues), (int)minVals[dimension]); 
                                    for (int val : possibleSplitValues) {
                                        float giniIndex = gini(goodActions[i], currentNode, dimension, (float)val);
                                        //STORM_PRINT("  Splitting on " << dimension << " at " << val << " has gini " << giniIndex << std::endl);
                                        if (giniIndex < bestGini) {
                                            bestGini = giniIndex;
                                            //bestAction = goodActions[i];
                                            bestDimension = dimension;
                                            bestSplitValue = val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    return std::make_tuple(bestDimension, bestSplitValue);
                }

                std::vector<std::vector<dtstrat::dt::HistoryElement<ValueType>>> partitionCompleteHistory(dtstrat::dt::Node<ValueType, float>& currentNode) const {
                    std::vector<dtstrat::dt::HistoryElement<ValueType>> left;
                    std::vector<dtstrat::dt::HistoryElement<ValueType>> right;
                    for (uint i = 0; i<currentNode.complete_history.size(); i++) {
                        dtstrat::dt::HistoryElement<ValueType> historyElement = currentNode.complete_history[i];
                        if (true) {
                            left.push_back(historyElement);
                        }
                        else {
                            right.push_back(historyElement);
                        }
                    }
                    std::vector<std::vector<dtstrat::dt::HistoryElement<ValueType>>> ret{left, right};
                    return ret;
                }

                uint sampleActionOfState(StateType const& currentStateId, storm::modelchecker::exploration_detail::ExplorationInformation<StateType, ValueType> const& explorationInformation, storm::modelchecker::exploration_detail::Bounds<StateType, ValueType>& bounds) const {
                    // Determine the values of all available actions.
                    storm::utility::ConstantsComparator<ValueType> comparator(storm::settings::getModule<storm::settings::modules::ExplorationSettings>().getPrecision());
                    std::vector<std::pair<uint, ValueType>> actionValues;
                    StateType rowGroup = explorationInformation.getRowGroup(currentStateId);
                    
                    // Check for cases in which we do not need to perform more work.
                    if (explorationInformation.onlyOneActionAvailable(rowGroup)) {
                        return explorationInformation.getStartRowOfGroup(rowGroup);
                    }
                                    
                    for (uint32_t row = explorationInformation.getStartRowOfGroup(rowGroup); row < explorationInformation.getStartRowOfGroup(rowGroup + 1); ++row) {
                        actionValues.push_back(std::make_pair(row, bounds.getBoundForAction(explorationInformation.getOptimizationDirection(), row)));
                    }
                    
                    // Sort the actions wrt. to the optimization direction.
                    if (explorationInformation.maximize()) {
                        std::sort(actionValues.begin(), actionValues.end(), [] (std::pair<uint, ValueType> const& a, std::pair<uint, ValueType> const& b) { return a.second > b.second; } );
                    } else {
                        std::sort(actionValues.begin(), actionValues.end(), [] (std::pair<uint, ValueType> const& a, std::pair<uint, ValueType> const& b) { return a.second < b.second; } );
                    }
                    
                    // Determine the first elements of the sorted range that agree on their value.
                    auto end = ++actionValues.begin();
                    while (end != actionValues.end() && comparator.isEqual(actionValues.begin()->second, end->second)) {
                        ++end;
                    }
                    
                    // Now sample from all maximizing actions.
                    std::uniform_int_distribution<uint> distribution(0, std::distance(actionValues.begin(), end) - 1);
                    std::default_random_engine randomGenerator( std::chrono::system_clock::now().time_since_epoch().count() );
                    return actionValues[distribution(randomGenerator)].first;
                }

                std::vector<uint> getBestActionsFromBRTDP(dtstrat::dt::Node<ValueType, float>& node, storm::modelchecker::exploration_detail::ExplorationInformation<StateType, ValueType> const& explorationInformation, storm::modelchecker::exploration_detail::Bounds<StateType, ValueType>& bounds, std::vector<std::vector<dtstrat::dt::Observation<ValueType>>>& observationPerAction, std::vector<uint>& uniqueDTActions) const {
                    std::vector<uint> stormActions;
                    for (uint i =0; i<node.observation_history.size(); i++) {
                        dtstrat::dt::Observation<ValueType> stateObservation = node.observation_history[i];
                        // this takes longer, but less memory
                        // uint stormStateId = getCurrentStateIdInStorm(state);
                        uint stormStateId = node.stormStateList[i];
                        uint stormAction = sampleActionOfState(stormStateId, explorationInformation, bounds);
                        uint dtAction = stormDtstratActionMapping[stormAction];
                        stormActions.push_back(dtAction);
                        auto it = std::find(uniqueDTActions.begin(), uniqueDTActions.end(), dtAction); 
                        uint index =0;
                        if (it == uniqueDTActions.end()) {
                            uniqueDTActions.push_back(dtAction);
                            index = uniqueDTActions.size()-1;
                        } else {
                            index = std::distance(uniqueDTActions.begin(), it);
                        }
                        if (observationPerAction.size()<uniqueDTActions.size()) {
                            observationPerAction.resize(uniqueDTActions.size());
                        }
                        observationPerAction[index].push_back(stateObservation);
                    }
                    return stormActions;
                }

                float entropy(float total, float good, float bad) const {
                    float val = 0.0;
                    if (good>0) {
                        val += -good/total*log2(good/total);
                    }
                    if (bad>0) {
                        val += -bad/total*log2(bad/total);
                    }
                    return val;
                }

                float informationGainBRTDP(uint action, dtstrat::dt::Node<ValueType, float>& node, uint dimension, float mean, std::vector<std::vector<dtstrat::dt::Observation<ValueType>>>& observationsPerAction, uint index, std::vector<uint>& goodActions) const {
                    uint numStates = node.observation_history.size();
                    uint numWithAction = observationsPerAction[index].size();
                    uint numLeft = 0;
                    uint numRight = 0;
                    uint numLeftWithAction = 0;
                    uint numRightWithAction = 0;
                    for (uint i=0; i<node.observation_history.size(); i++) {
                        dtstrat::dt::Observation<ValueType> stateObservation = node.observation_history[i];
                        if (stateObservation[dimension]>=mean) {
                            numRight++;
                            if (goodActions[i]==action) {
                                numRightWithAction++;
                            }
                        }
                        else {
                            numLeft++;
                            if (goodActions[i]==action) {
                                numLeftWithAction++;
                            }
                        }
                    }
                    if (numStates == 0 || numLeft == 0 || numRight == 0) {
                        return 0;
                    }
                    // entropy before the split
                    float ig = entropy((float) numStates, (float) numWithAction, (float) (numStates-numWithAction));
                    ig -= ( (float)numLeft/(float)numStates*entropy((float) numLeft, (float) numLeftWithAction, (float) (numLeft-numLeftWithAction))
                            +(float)numRight/(float)numStates*entropy((float) numRight, (float) numRightWithAction, (float) (numRight-numRightWithAction)) );
                    return ig;
                };

                std::tuple<uint,float> getSplitBasedOnBRTDP(dtstrat::dt::Node<ValueType, float>& node, storm::modelchecker::exploration_detail::ExplorationInformation<StateType, ValueType> const& explorationInformation, storm::modelchecker::exploration_detail::Bounds<StateType, ValueType>& bounds) const {
                    std::vector<std::vector<dtstrat::dt::Observation<ValueType>>> observationsPerAction;
                    std::vector<uint> uniqueDTActions;
                    std::vector<uint> goodActions = getBestActionsFromBRTDP(node, explorationInformation, bounds, observationsPerAction, uniqueDTActions);
                    // sort by amount maybe? or do that later
                    // similar to the above variant
                    float bestIG = 0;
                    //int bestAction = -1;
                    uint bestDimension = 0;
                    float bestSplitValue = 0.0;
                    
                    // for each action in the good choices
                    if (uniqueDTActions.size()>0) {
                        uint sizeStateObservation = observationsPerAction[0][0].size;
                        float* means = new float[sizeStateObservation];
                        float* deviation = stdev_observation(node.observation_history, means);
                        for (uint i = 0; i<uniqueDTActions.size(); i++) {
                            //STORM_PRINT("/SPLIT ON "<<uniqueDTActions[i] <<std::endl);
                            // for each dimension to split on
                            // check the variance (if == 0, then you can't split), if >0, split at mean and calculate Gini or IG1
                            for (uint dimension = 0; dimension<sizeStateObservation; dimension++){  
                                if (deviation[dimension]>0) {
                                    float igB = informationGainBRTDP(uniqueDTActions[i], node, dimension, means[dimension], observationsPerAction, i, goodActions);
                                    //STORM_PRINT("  Splitting on " << dimension << " at " << means[dimension] << " has IG " << igB << std::endl);
                                    if (igB > bestIG) {
                                        bestIG = igB;
                                        //bestAction = uniqueDTActions[i];
                                        bestDimension = dimension;
                                        bestSplitValue = means[dimension];
                                    }
                                }
                            }
                        }
                    }
                    return std::make_tuple(bestDimension, bestSplitValue);
                }

                std::tuple<uint,float> getSplitBasedOnBRTDPExhaustive(dtstrat::dt::Node<ValueType, float>& node, storm::modelchecker::exploration_detail::ExplorationInformation<StateType, ValueType> const& explorationInformation, storm::modelchecker::exploration_detail::Bounds<StateType, ValueType>& bounds) const {
                    std::vector<std::vector<dtstrat::dt::Observation<ValueType>>> observationsPerAction;
                    std::vector<uint> uniqueDTActions;
                    std::vector<uint> goodActions = getBestActionsFromBRTDP(node, explorationInformation, bounds, observationsPerAction, uniqueDTActions);
                    // sort by amount maybe? or do that later
                    // similar to the above variant
                    float bestIG = 0;
                    //int bestAction = -1;
                    uint bestDimension = 0;
                    float bestSplitValue = 0.0;
                    
                    // for each action in the good choices
                    if (uniqueDTActions.size()>0) {
                        uint sizeStateObservation = observationsPerAction[0][0].size;
                        float* maxVals = new float[sizeStateObservation];
                        float* minVals = new float[sizeStateObservation];
                        float* deviation = stdev_observation(node.observation_history, nullptr, minVals, maxVals);
                        for (uint i = 0; i<uniqueDTActions.size(); i++) {
                            //STORM_PRINT("/SPLIT ON "<<uniqueDTActions[i] <<std::endl);
                            // for each dimension to split on
                            // check the variance (if == 0, then you can't split), if >0, split at mean and calculate Gini or IG1
                            for (uint dimension = 0; dimension<sizeStateObservation; dimension++){  
                                if (deviation[dimension]>0) {
                                    std::vector<int> possibleSplitValues((int)(maxVals[dimension]-minVals[dimension])+1);
                                    // generates a range in possibleSplitValues, starting from (int)minVals[dimension] and ending when the vector ends
                                    std::iota(std::begin(possibleSplitValues), std::end(possibleSplitValues), (int)minVals[dimension]); 
                                    for (int val : possibleSplitValues) {
                                        float igB = informationGainBRTDP(uniqueDTActions[i], node, dimension, val, observationsPerAction, i, goodActions);
                                        //STORM_PRINT("  Splitting on " << dimension << " at " << val << " has IG " << igB << std::endl);
                                        if (igB > bestIG) {
                                            bestIG = igB;
                                            //bestAction = uniqueDTActions[i];
                                            bestDimension = dimension;
                                            bestSplitValue = val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    return std::make_tuple(bestDimension, bestSplitValue);
                }

                std::vector<std::tuple<uint,float>> getActionsWithQvalue(const dtstrat::dt::Node<ValueType, float>& currentNode, std::vector<std::vector<dtstrat::dt::Observation<ValueType>>>& observationPerAction) const {
                    // only get good choices
                    std::vector<std::tuple<uint,float>> goodActions;
                    for (dtstrat::dt::HistoryElement<ValueType> element : currentNode.complete_history) {
                            std::tuple<uint,float> tup = std::make_tuple(element.action, element.qValue);
                            auto it = std::find(goodActions.begin(), goodActions.end(), tup); 
                            uint index = 0;
                            if (it == goodActions.end()) {
                                goodActions.push_back(tup);
                                index = goodActions.size()-1;
                            } else {
                                index =  std::distance(goodActions.begin(), it);
                            }
                            if (observationPerAction.size()<goodActions.size()) {
                                observationPerAction.resize(goodActions.size());
                            }
                            observationPerAction[index].push_back(element.state);
                    }
                    return goodActions;
                }

                 float gini2(std::tuple<uint,float> actionAndQVal, dtstrat::dt::Node<ValueType, float>& currentNode, uint dimension, float mean) const {
                    uint num_left = 0;
                    uint num_right = 0;
                    uint num_good_left = 0;
                    uint num_good_right = 0;
                    for (dtstrat::dt::HistoryElement<ValueType> element : currentNode.complete_history) {
                        // decide if element is left or right
                        if (element.state[dimension]>=mean) {
                            //right
                            num_right++;
                            if (element.action == std::get<0>(actionAndQVal) && std::abs(element.qValue-std::get<1>(actionAndQVal))<1e-6 ) {
                                num_good_right++;
                            }
                        }
                        else {
                            //left
                            num_left++;
                            if (element.action == std::get<0>(actionAndQVal) && std::abs(element.qValue-std::get<1>(actionAndQVal))<1e-6 ) {
                                num_good_left++;
                            }
                        }
                    }
                    uint num_total = currentNode.complete_history.size();
                    if (num_total == 0 || num_left == 0 || num_right == 0) {
                        return 1;
                    }
                    else {
                        float giniValue = giniVal((float) num_total, (float) num_left, (float) num_right, (float) num_good_left, (float) num_good_right);
                        return giniValue;     
                    }              
                }

                std::tuple<uint,float> getSplitBasedOnGoodActionByQValue(dtstrat::dt::Node<ValueType, float>& currentNode) const {
                    // get all good actions
                    std::vector<std::vector<dtstrat::dt::Observation<ValueType>>> observationPerAction;
                    std::vector<std::tuple<uint,float>> goodActions = getActionsWithQvalue(currentNode, observationPerAction);

                    float bestGini = 1;
                    uint bestDimension = 0;
                    float bestSplitValue = 0.0;
                    
                    // for each action in the good choices
                    if (goodActions.size()>0) {
                        uint sizeStateObservation = observationPerAction[0][0].size;
                        float* means = new float[sizeStateObservation];
                        float* deviation = stdev_observation(currentNode.observation_history, means);
                        for (uint i = 0; i<goodActions.size(); i++) {
                            //STORM_PRINT("/SPLIT ON "<<goodActions[i] <<std::endl);
                            // for each dimension to split on
                            // check the variance (if == 0, then you can't split), if >0, split at mean and calculate Gini or IG1
                            for (uint dimension = 0; dimension<sizeStateObservation; dimension++){  
                                if (deviation[dimension]>0) {
                                    float giniIndex = gini2(goodActions[i], currentNode, dimension, means[dimension]);
                                    //STORM_PRINT("  Splitting on " << dimension << " at " << means[dimension] << " has gini " << giniIndex << std::endl);
                                    if (giniIndex < bestGini) {
                                        bestGini = giniIndex;
                                        //bestAction = goodActions[i];
                                        bestDimension = dimension;
                                        bestSplitValue = means[dimension];
                                    }
                                }
                            }
                        }
                    }
                    return std::make_tuple(bestDimension, bestSplitValue);
                }

                std::tuple <int, float> get_split(dtstrat::dt::Node<ValueType, float>& node, storm::modelchecker::exploration_detail::ExplorationInformation<StateType, ValueType> const& explorationInformation, storm::modelchecker::exploration_detail::Bounds<StateType, ValueType>& bounds) const {
                    //std::ofstream logger;
                    //logger.open ("DEBUG.OUT", std::ios_base::app);
                    //for (auto it = node.observation_history.begin(); it!=node.observation_history.end(); it++) {
                    //    logger << (*it).printToCSV() <<"\n";
                    //}
                    //logger.close();

                    int dim;
                    int splitvalue;
                    if (this->splitDecision=="tstatistic") {
                        float* tstats = t_statistic(node.observation_history);
                        uint size = node.observation_history[0].size;
                        if (no_positive(node.history)) {
                            float* deviations = stdev_observation(node.observation_history);
                            dim = -1;
                            float maxdev = -numeric_limits<float>::max();
                            for (uint i = 0; i<size; i++) {
                                if (deviations[i]>maxdev) {
                                    dim = i;
                                    maxdev = deviations[i];
                                }
                            }
                        }
                        // calculate the t-statistic of each parameter
                        else {
                            dim = -1;
                            float maxdev = -numeric_limits<float>::max();;
                            for (uint i = 0; i<size; i++) {
                                if (tstats[i]>maxdev) {
                                    dim = i;
                                    maxdev = tstats[i];
                                }
                            }
                        }
                        // NOTE: Is this correct? There is nothing in the paper about actually getting the splitvalue...
                        //float splitvalue = mean_observation(node.observation_history)[dim];
                        splitvalue = int(mean_observation(node.observation_history)[dim]);
                    } 
                    else if (this->splitDecision=="groupingActionGini") {
                        // New
                        std::tuple<uint,float> splitting = getSplitBasedOnGoodActionAndGini(node);
                        dim = std::get<0>(splitting);
                        splitvalue = std::get<1>(splitting);
                    } 
                    else if (this->splitDecision=="groupingActionGiniExhaustive") {
                        std::tuple<uint,float> splitting = getSplitBasedOnGoodActionAndGiniExhaustive(node);
                        dim = std::get<0>(splitting);
                        splitvalue = std::get<1>(splitting);
                    }
                    else if (this->splitDecision=="brtdp") {
                        std::tuple<uint,float> splitting = getSplitBasedOnBRTDP(node, explorationInformation, bounds);
                        dim = std::get<0>(splitting);
                        splitvalue = std::get<1>(splitting);
                    }
                    else if (this->splitDecision=="brtdpExhaustive") {
                        std::tuple<uint,float> splitting = getSplitBasedOnBRTDPExhaustive(node, explorationInformation, bounds);
                        dim = std::get<0>(splitting);
                        splitvalue = std::get<1>(splitting);
                    }
                    else if (this->splitDecision=="groupQvalue") {
                        std::tuple<uint,float> splitting = getSplitBasedOnGoodActionByQValue(node);
                        dim = std::get<0>(splitting);
                        splitvalue = std::get<1>(splitting); 
                    }
                    else {
                        STORM_LOG_ERROR("Statistic for split " << this->splitDecision << " not implemented!");
                    }

                    bool alreadyExists = false;
                    for (uint i=0; i<node.split_history.size(); i++) {
                        if ((node.split_history[i].first == dim) && (node.split_history[i].second == splitvalue)) {
                            alreadyExists = true;
                            break;
                        }
                    }

                    if (alreadyExists) {
                        std::tuple<int,float> ret = {-1,-1};
                        //logger.open ("DEBUG.OUT", std::ios_base::app);
                        //logger << "Already exists!!! " << dim << " " << splitvalue;//<< "\n";
                        //logger.close();
                        return ret;
                    } else {
                        std::tuple<int, float> ret = {dim, splitvalue};
                        //logger.open ("DEBUG.OUT", std::ios_base::app);
                        //logger << dim << " " << splitvalue;//<< "\n";
                        //logger.close();
                        return ret;
                    }
                }


            // update one node for given rewards and action
            void updateNode(dtstrat::dt::Node<ValueType, float>* changeNode, dtstrat::dt::Node<ValueType, float> nextNode, int reward, int action) const {
                //float alpha = 0.9;
                //float gamma = 0.9;
                int bestaction = nextNode.get_max();
                // Q(s,a) = sum_(s' succ) delta(s,a,s')*Q(s') 
                // Q(s') = max_a Q(s',a)
                float change_value = this->alpha*(reward + this->gamma * nextNode.q_values[bestaction] - changeNode->q_values[action]);
                changeNode->updateQValue(action, changeNode->q_values[action] + change_value);
                //changeNode->q_values.qvalues[action] = changeNode.q_values.qvalues[action] + change_value;
                
                // safety option
                //if (changeNode.history.size() < 1000) {
                    //changeNode.history.push_back(change_value);
                //}
            
                return;
            };

            bool stateEqualsValuation(storm::storage::sparse::StateValuations::StateValuation& value, dtstrat::dt::Observation<ValueType>& observe) const {
                bool ret = true;
                if (value.integerValues.size() + value.booleanValues.size() == 0) {
                    // if the valuation is empty, we don't have if it's equal to the observation
                    ret = false;
                }
                for (uint i=0; i<value.integerValues.size(); i++) {
                    ret = (value.integerValues[i] == observe.observation_values[i]);
                    if (!ret) break;
                }
                if (!ret) {
                    for (int j=0;j<value.booleanValues.size(); j++) {
                        ret = (value.booleanValues[j] == observe.observation_values[j+value.integerValues.size()]);
                        if (!ret) break;
                    }
                }
                return ret;
            }

            int getCurrentStateIdInStorm(dtstrat::dt::Observation<ValueType>& state) const {
                int index = -1;
                for (uint i=0; i<allStateObservations.getSize(); i++) {
                    dt::Observation<ValueType> observation = allStateObservations.getObservation(i);
                    if (observation==state) {
                        index = i;
                        break;
                    }
                }
                return index;
            }

            int possibleAction(uint actionIndex, uint32_t const& currentStateId, storm::modelchecker::exploration_detail::ExplorationInformation<uint32_t, ValueType>& explorationInformation) const {
                auto unexploredIt = explorationInformation.findUnexploredState(currentStateId);
                if (unexploredIt != explorationInformation.unexploredStatesEnd()) {
                    // State is not in explored states -> we don't know anything about it
                    //STORM_PRINT(" --- State " << currentStateId << " is unexplored"<<std::endl);
                    return false;
                }
                auto rowGroup = explorationInformation.getRowGroup(currentStateId);
                //STORM_PRINT(" --- State " << currentStateId << " Storm-Action " << actionIndex << " RowGroup "<< rowGroup << " [" << explorationInformation.getStartRowOfGroup(rowGroup) << "," << explorationInformation.getStartRowOfGroup(rowGroup+1)<< "]");
                // the action was before and is not possible in this state
                if (actionIndex<explorationInformation.getStartRowOfGroup(rowGroup)) {
                    //STORM_PRINT("    --> too small"<<std::endl);
                    return -1;
                } else if (actionIndex>=explorationInformation.getStartRowOfGroup(rowGroup + 1)){
                    //STORM_PRINT("    --> too big"<<std::endl);
                    return 1;
                } else {
                    //STORM_PRINT("   --> IS OK!" <<std::endl);
                    return 0;
                }
            }

                uint32_t getStormActionFromDT(uint dTActionId, uint32_t const& currentStateId, storm::modelchecker::exploration_detail::ExplorationInformation<uint32_t, ValueType>& explorationInformation, std::vector<int> const& stormDtstratActionMapping, bool& worked) const {
                    //STORM_PRINT("// Search " << currentStateId << " and " << dTActionId << std::endl);
                    uint32_t stormIndex = 0;
                    // Get the action index
                    bool found = false;
                    for (uint i = 0; i<stormDtstratActionMapping.size(); i++) {
                        // i is storm index
                        uint const iDTIdx = stormDtstratActionMapping[i];
                        if (iDTIdx == dTActionId) {
                            // Check if this action is possible, unless continue searching
                            int result = possibleAction(i, currentStateId, explorationInformation);
                            //STORM_PRINT("result is "<<result<<std::endl);
                            if (result==0) {
                                found = true;
                                stormIndex = i;
                                break;
                            } else if (result==1){
                                //STORM_PRINT("Not found, because too big"<<std::endl);
                                break;
                            }
                        }
                    }
                    if (!found) {
                        worked = false;
                    }
                    //STORM_LOG_ASSERT(found, "NOT FOUND THE INDEX IN DTSTRAT-STORM MAPPING!");
                    return stormIndex;
                }

                void updateNodeTrueQValue(dtstrat::dt::Node<ValueType, float>* changeNode, uint action, storm::modelchecker::exploration_detail::ExplorationInformation<uint32_t, ValueType>& explorationInformation, dtstrat::dt::Observation<ValueType> state, dtstrat::dt::DecisionTree<ValueType, float>& decisionTree, uint stormStateId, uint count) const {
                    // Q(s,a) = sum_(s' succ) delta(s,a,s')*Q(s') 
                    // Q(s') = max_a Q(s',a) = nextNode.q_values.qvalues[bestaction]
                    int currentStateI = getCurrentStateIdInStorm(state);
                    if (currentStateI == -1) { 
                        STORM_PRINT("ERROR: could not find "<< state.print() << std::endl);
                        return;
                    }
                    uint currentStateId = currentStateI;
                    bool worked = true;
                    auto stormActionId = getStormActionFromDT(action, currentStateId, explorationInformation, stormDtstratActionMapping, worked);
                    if (!worked) {
                        //float change_value = 
                        changeNode->updateQValue(action, changeNode->q_values[action]);
                        return;
                    }
                    std::vector<storm::storage::MatrixEntry<uint32_t, ValueType>> const& row = explorationInformation.getRowOfMatrix(stormActionId);
                    std::vector<std::pair<int, double>> probabilities(row.size());
                    std::transform(row.begin(), row.end(), probabilities.begin(),
                            [] (storm::storage::MatrixEntry<uint32_t, ValueType> const& entry) {return entry.getColumnValuePair();});
                    // probabilities contains tuples of (state, probability)
                    // that is it contains all possible successors of the current state toghether with their proabilities


                    //std::vector<float> qValueOfSuccessors;
                    // delta(s,a,s')?
                    //float oldQvalue = changeNode->q_values.qvalues[action];
                    float newValue = 0.0;
                    //STORM_PRINT("  - state " << state.print() << " action " << action);
                    for (uint i = 0; i<row.size(); i++) {
                        uint stateId = probabilities[i].first;
                        //STORM_PRINT(" succ_" << i+1 << " " << stateId);
                        float qval = 0;
                        if (allStateObservations.existsStateObservation(stateId)) {
                            // Only get the qvalue if we have seen the successor before
                            dtstrat::dt::Observation<ValueType> stateObservation = allStateObservations.getObservation(stateId);
                            // if the successor would be a target state, we don't want to ask the DT about the qvalues
                            //STORM_PRINT("=" << stateObservation.print());
                            for (auto targetState : this->targetStates) {
                                if (stateObservation == targetState) {
                                    qval = 1.0;
                                    break;
                                }
                            }
                            if (qval<1) {
                                auto node = decisionTree.get_leaf_node(stateObservation);
                                qval = node->q_values[node->get_max()];
                            }
                        } else {
                            // otherwise, set the qvalue to zero
                            // we cannot calculate the node to which the state corresponds because we don't know it and its observations
                            qval = 0.0;
                        }
                        //qValueOfSuccessors.push_back(qval);
                        //newValue += probabilities[i].second * qValueOfSuccessors[i];

                        // include penalty for same-state-transitions (loops)
                        if (stateId == currentStateId) {
                            newValue += probabilities[i].second * qval * this->gamma;
                        } else {
                            newValue += probabilities[i].second * qval;
                        }
                        //STORM_PRINT(" newvalue="<< probabilities[i].second << "*" << qval << "=" <<newValue);
                    }
                    // qValueOfSuccessors contains the maximum possible qvalue of the next state

                    //STORM_PRINT(std::endl);
                    newValue *= pow(this->gamma, count);
                    float change_value = changeNode->updateQValue(action, newValue);
                    
                    if (this->splitDecision=="groupingActionGini" | this->splitDecision=="groupingActionGiniExhaustive") {
                        changeNode->updateCompleteHistory(state, action, change_value, newValue);
                    } 
                    else if (this->splitDecision=="brtdp" | this->splitDecision=="brtdpExhaustive") {
                        changeNode->updateBRTDPList(stormStateId);
                    }

                    //STORM_PRINT(" s: " << state.print() << " a: " << action << " n: " << changeNode->num << " c: " << change_value << " q: " << changeNode->q_values.qvalues[action] << std::endl;);
                    //dtstrat::log::logExp logg = dtstrat::log::logExp();
                    //logg.splitDebugLogger(state, action, change_value, changeNode->q_values.qvalues[action], changeNode->num);
                    return;
                }

                // vector of observations (states), actions and rewards
                void performUpdateBatch(dtstrat::dt::DecisionTree<ValueType, float>& decisionTree, storm::modelchecker::exploration_detail::ExplorationInformation<StateType, ValueType>& explorationInformation, storm::modelchecker::exploration_detail::Bounds<StateType, ValueType>& bounds) {
                    // The last element is the target state
                    //STORM_PRINT("update: Size of explorationInfo "<< explorationInformation.getActionCount()<<std::endl);
                    //STORM_PRINT("update: Size of stormdtstrat "<<stormDtstratActionMapping.size()<<std::endl);
                    //STORM_PRINT("UDPDATE"<<std::endl);
                    if (deleteTree) {
                        decisionTree.clear();
                    }

                    std::tuple<dtstrat::dt::Observation<ValueType>,int,int,bool> lastElement = updateBatch.back();
                    // Target State
                    dtstrat::dt::Observation<ValueType> targetState = std::get<0>(lastElement);
                    this->addTargetState(targetState);
                    updateBatch.pop_back();
                    updateBatchStormStates.pop_back();
                    // start from the last element
                    uint count = 1;
                    while (updateBatch.size() > 0) {
                        std::tuple<dtstrat::dt::Observation<ValueType>,int,int,bool> currentElement  = updateBatch.back();
                        updateBatch.pop_back();

                        // get the data from the batch element (state, action, reward)
                        // reward is mostly 0 for reachability queries
                        dtstrat::dt::Observation<ValueType> state = std::get<0>(currentElement);
                        int currentAction = std::get<1>(currentElement);
                        if (currentAction<0) {
                            // This is a target state
                            this->addTargetState(state);
                            continue;
                        }

                        // get reward from last element, because I want the reward 1 from the target state
                        int currentReward = std::get<2>(lastElement);

                        dtstrat::dt::Node<ValueType, float>* const currentNode = decisionTree.get_leaf_node(state);
                        dtstrat::dt::Node<ValueType, float> lastNode = *(decisionTree.get_leaf_node(std::get<0>(lastElement)));
                        
                        if (useTrueQValue) {
                            uint stormState = 0;
                            if (this->splitDecision=="brtdp") {
                                stormState = this->updateBatchStormStates.back();
                                this->updateBatchStormStates.pop_back();
                            }
                            updateNodeTrueQValue(currentNode, currentAction, explorationInformation, state, decisionTree, stormState, count);
                            count++;
                        } else {
                            updateNode(currentNode, lastNode, currentReward, currentAction);
                        }

                        //if (currentNode.observation_history.size() < 10000) {
                            currentNode->updateObservations(state);
                        //}
                        lastElement = currentElement;

                        // check if the DT should be split
                        
                        if (deleteTree) {
                            bool split = decide_split(*currentNode);
                            if (split) {
                                // perform the split
                                //std::cout << "We split!" << std::endl;
                                std::tuple <int, float> splitvals = get_split(*currentNode, explorationInformation, bounds);
                                int dimension = get<0>(splitvals);
                                float value = get<1>(splitvals);
                                if (dimension!=-1 && value>=0) {
                                    decisionTree.split(*currentNode, dimension, value);
                                }
                            }
                        }
                        
                    }
                    // decide if to split outside the function (should make it faster)
                    for (dtstrat::dt::Node<ValueType, float>* leafNode : decisionTree.getLeafNodes()) {
                        //dtstrat::dt::Node<ValueType>* leafNode = *leafNodeAddr;
                        bool split = decide_split(*leafNode);
                        if (split) {
                            // perform the split
                            //std::cout << "We split!" << std::endl;
                            std::tuple <int, float> splitvals = get_split(*leafNode, explorationInformation, bounds);
                            int dimension = get<0>(splitvals);
                            float value = get<1>(splitvals);
                            if (dimension!=-1 && value>=0) {
                                decisionTree.split(*leafNode, dimension, value);
                            }
                        }
                    }
                    return;
                };

                std::string get_last_split(std::string const &str, const char delim)
                {
                    size_t start;
                    size_t end = 0;
                    std::string last;

                    while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
                    {
                        end = str.find(delim, start);
                        last = str.substr(start, end - start);
                    }

                    return last;
                };

                float evaluate_term(nlohmann2::json term, std::vector<float>  state) {
                    float coeff, intercept;
                    std::string var;
                    if (term.contains("coeff") && term.contains("var")) {
                        term["coeff"].get_to(coeff);
                        term["var"].get_to(var);
                        int v = std::stoi(get_last_split(var, '_'));
                        return coeff * state[v];
                    } else {
                        return term["intercept"].get_to(intercept);
                    }
                };

                int evaluate_node(nlohmann2::json split, const std::vector<float>& state) {
                    // return 0 means false edge, 1 means true edge, TODO multi
                    int ret = 0;

                    // split is either an axis-parallel expression or a list of linear terms
                    auto lhs = split["lhs"];
                    auto op = split["op"];

                    if (op == "multi") {
                        // TODO do something else
                        //  and return
                        STORM_PRINT("Not implemented operator 'multi' in the decisionTree from dtcontrol");
                    }

                    float rhs;
                    split["rhs"].get_to(rhs);
                    float lhs_eval;

                    if (lhs.contains("coeff") and lhs.contains("var")) {
                        // Then this is an axis-parallel expression
                        lhs_eval = evaluate_term(lhs, state);
                    }
                    else {
                        // In this case its a linear split
                        lhs_eval = 0;
                        for (auto & term : lhs) {
                            lhs_eval += evaluate_term(term, state);
                        }
                    }

                    if (op == "<=") {
                        if (lhs_eval <= rhs)
                            ret = 1;
                    }
                    else if (op == ">=") {
                        if (lhs_eval >= rhs)
                            ret = 1;
                    }
                    else if (op == "==") {
                        if (lhs_eval == rhs)
                            ret = 1;
                    }
                    else if (op == "multi") {
                        // ???
                    }
                    else {
                        // Control shouldn't come here
                        // TODO: How should exceptions be handled?
                        throw(std::exception());
                    }

                    return ret;
                };

                std::vector<std::string> get_decision(nlohmann2::json tree, const std::vector<float>& state) {
                    bool edge;
                    // If split contains a value, then this is an inner node and needs to be evaluated
                    if (tree.contains("split") && tree["split"] != nullptr) {
                        edge = evaluate_node(tree["split"], state);
                        for (auto & child : tree["children"]) {
                            if (child["edge_label"] == (edge == 1 ? "true" : "false")) {
                                return get_decision(child, state);
                            }
                        }
                        // TODO Handle case where tree has no children
                    } else {
                        // If split doesn't exist, this is a leaf node
                        return tree["actual_label"];
                    }
                };

            void get_lookup(const std::string& filename) {
                std::ifstream inFile(filename);
                int SIZE = std::count(std::istreambuf_iterator<char>(inFile),
                           std::istreambuf_iterator<char>(), '\n');
                lookUpTable = boost::numeric::ublas::mapped_matrix<float>(SIZE, 1000, SIZE * 1000);

                std::ifstream infile(filename);
                std::string line;
                while (std::getline(infile, line))
                {
                    std::istringstream iss(line);
                    std::vector<double> observed;
                    int temp;
                    while (iss >> temp) {
                        observed.push_back(temp);
                    }
                    int a, b;
                    b = observed.back();
                    observed.pop_back();
                    dtstrat::dt::Observation<double> newObservation(observed.size());
                    newObservation.set(observed);
                    int found_index = this->allStateObservations.findObservation(newObservation);
                    if (found_index>=0) {
                        a = found_index;
                    } else {
                        this->allStateObservations.pushback(newObservation);
                        a = allStateObservations.getSize() - 1;
                    }
                    boost::numeric::ublas::matrix_row<boost::numeric::ublas::mapped_matrix<int>> aRow(lookUpTable, a);
                    uint test = 0;
                    for (boost::numeric::ublas::matrix_row<boost::numeric::ublas::mapped_matrix<int> >::iterator it2 = aRow.begin(); it2 != aRow.end(); it2++) test++;
                    lookUpTable(a,test) = b;
                    // process pair (a,b)
                }
            }

            int get_lookup_decision(dtstrat::dt::Observation<double> state) {
                int index = allStateObservations.findObservation(state);
                float returnval = lookUpTable(index,0);
                return returnval;
            }

            std::vector<unsigned int> get_lookup_decisions(dtstrat::dt::Observation<double> state) {
                int index = allStateObservations.findObservation(state);
                std::vector<unsigned int> returnval;
                if (index<0) {
                    return returnval;
                }
                boost::numeric::ublas::matrix_row<boost::numeric::ublas::mapped_matrix<int>> aRow(lookUpTable, index);
                for (boost::numeric::ublas::matrix_row<boost::numeric::ublas::mapped_matrix<int> >::iterator it2 = aRow.begin(); it2 != aRow.end(); it2++) {
                    returnval.push_back((unsigned int) *it2);
                }
                return returnval;
            }

            void setModuleChoice(const std::string& method) {
                if (method=="random") {
                    this->moduleMethodChoice = random;
                } else if (method=="treeSize") {
                    this->moduleMethodChoice = treeSize;
                } else if (method=="smallest") {
                    this->moduleMethodChoice = smallest;
                } else {
                    this->moduleMethodChoice = random;
                }
            }

            uint getModuleToPlay() {
                std::vector<std::pair<uint,int>> availableModules;
                for (uint i=0; i<this->multiDts.size(); i++) {
                    auto const& dt = this->multiDts[i];
                    if (dt.setFromFile()) {
                        availableModules.push_back(std::make_pair(i,dt.size()));
                    }
                }
                if (this->moduleMethodChoice==random) {
                    std::default_random_engine generator;
                    std::uniform_int_distribution<int> distribution(0,availableModules.size());
                    return availableModules[distribution(generator)].first;
                } else if (this->moduleMethodChoice==treeSize) {
                    std::sort(availableModules.begin(), availableModules.end(),
                              [] (std::pair<uint,int> const& a, std::pair<uint,int> const& b) { return a.second>b.second; } );
                    return availableModules[0].first;
                } else if (this->moduleMethodChoice==smallest) {
                    std::sort(availableModules.begin(), availableModules.end(),
                              [] (std::pair<uint,int> const& a, std::pair<uint,int> const& b) { return a.second<b.second; } );
                    return availableModules[0].first;
                } else {
                    STORM_PRINT("Not impelemnted module choice " << this->moduleMethodChoice << std::endl);
                    return availableModules[0].first;
                }
                return 0;
            }

            std::pair<uint, std::vector<unsigned int>> get_multi_module_decisions(dtstrat::dt::Observation<double> state) {
                uint module_to_play = getModuleToPlay();
                vector<uint32_t> decision = this->multiDts[module_to_play].get_decisions(state);
                return std::make_pair(module_to_play, decision);
            }

            std::vector<unsigned int> get_module_decisions(dtstrat::dt::Observation<double> state, std::vector<dtstorm::Action> const &allActions, std::vector<storm::prism::Module>& all_modules) {
                std::vector<unsigned int> returnval;
                dt::Node<ValueType, int> leafNode = *(this->jsonDT.get_leaf_node(state));
                 for (auto module : leafNode.modules) {
                    for (auto action : leafNode.q_values) {
                        for (unsigned int i=0; i<allActions.size(); i++){
                            auto modules = allActions[i].getModules();
                            auto cmdIndices = allActions[i].getCommandIndices();
                            std::string modulename = all_modules[module].getName();
                            //cout << " module: " << modi;
                            if (std::find(modules.begin(), modules.end(), modulename) != modules.end()) {
                                if (std::find(cmdIndices.begin(), cmdIndices.end(), action) != cmdIndices.end()) {
                                    returnval.push_back(i);
                                }
                            }
                        }
                    }
                }
                return returnval;
            }

        };
    }
}

#endif
