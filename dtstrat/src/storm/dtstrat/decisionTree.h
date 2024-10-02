#ifndef DECISION_TREE_H
#define DECISION_TREE_H

/*
This header contains classes to create a fully working decision tree, that contains qvalues in its leaf nodes
*/

#include <optional>
#include <list>
#include <string>
#include <tuple>
#include <limits>
#include <vector>
#include <cstring>
#include <algorithm>
//#include "learning_algorithm.h"
#include "storm/dtstrat/json.hpp"
#include "storm/storage/sparse/StateValuations.h"
#include <random>

using namespace std;

namespace dtstrat {
    namespace dtstorm {
        // needed for DTStrat to remember actions and how to map them between storm and DTStrat
        class Action {
        public:
            Action(std::vector<std::string> mods, std::vector<int> cmdIndices, std::string lab) {
                for (uint i = 0; i<mods.size(); i++) {
                    modules.push_back(mods[i]);
                    commandIndices.push_back(cmdIndices[i]);
                }
                label = lab;
            }

            bool operator==(Action const& one) const {
                bool ret = true;
                if (one.size() != this->size()) {return false;}
                if (one.getLabel() != this->getLabel()) {return false;}
                std::vector<std::string> oneModules = one.getModules();
                std::vector<int> oneCmdIndices = one.getCommandIndices();
                for (uint i = 0; i<modules.size(); i++) {
                    if (modules[i] != oneModules[i]) {
                        ret = false;
                        break;
                    }
                    if (commandIndices[i] != oneCmdIndices[i]) {
                        ret = false;
                        break;
                    }
                }
                return ret;
            }

            int size() const {
                return modules.size();
            }

            std::string getLabel() const {
                return label;
            }

            std::vector<std::string> getModules() const {
                return modules;
            }

            std::vector<int> getCommandIndices() const {
                return commandIndices;
            }

            void appendModule(std::string newModule) {
                modules.push_back(newModule);
                return;
            }

            void appendCommand(int newCommand) {
                commandIndices.push_back(newCommand);
                return;
            }
        private:
            // An action is characterized by all modules that it appears in:
            std::vector<std::string> modules;
            // and the corresponding command indices from storm (basically "lines" in the prism-model)
            std::vector<int> commandIndices;
            // and the label (e.g. 'reset')
            std::string label;
        };
    }


    namespace dt {

		// Class to store an observeration of the environment
		template <typename obValueType>
		class Observation {
			public:
			//... ATTRIBUTES ...
				// integer values for the observation itself
				vector<obValueType> observation_values;

				// number of observation values
				uint size;

			//... FUNCTIONS ...// 
				// creates an observation, default size is 6
				Observation(uint isize=0) {
					// create a new observation and initliaze it with 0 everywhere
					size = isize;
					observation_values.resize(size);
					for (uint i=0; i<size; i++){
						observation_values[i] = 0;
					};
				};

				// set all the qvalues to one value
				void setOne(obValueType value) {
					for (uint i=0; i<size; i++){
						observation_values[i] = value;
					};
				};

				// sets new values to an observation
				void set(const obValueType* values) {
					//memcpy(observation_values, values, size * sizeof(int));
					this->observation_values.clear();
					for(int i=0; i<size; i++) {
						observation_values.push_back(values[i]);
					}
				};

				void set(std::vector<obValueType> values) {
				    this->observation_values.clear();
				    this->observation_values.reserve(values.size());
				    for (int value : values) {
				        this->observation_values.emplace_back(value);
				    }
				};

				// allows to access the values by using [], e.g. values[1]
				int operator[](const uint i) const {
					return observation_values[i];
				};

				bool operator == (const Observation& d) const {
					if (size!=d.size) return false;
					else{
						for(uint i=0; i<size; i++) {
							if(observation_values[i]!=d[i]){
								return false;
							}
						}
						return true;
					}
				};

				// outputs the observation values as a string
				string print(std::vector<int> const& stateVariableModuleAssociation={}) const {
                	string ret = "";
                    bool useSVMA = (!stateVariableModuleAssociation.empty());
					for(uint i=0; i<size; i++) {
                        if (useSVMA) ret += to_string(stateVariableModuleAssociation[i]) + " ";
						ret+= to_string((int)observation_values[i]);
						ret+= " ";
					}
					return ret;
				};

				string printToCSV(std::vector<int> const& stateVariableModuleAssociation={}) const {
					string ret = "";
                    bool useSVMA = (!stateVariableModuleAssociation.empty());
					for(uint i=0; i<size; i++) {
                        if (useSVMA) ret += to_string(stateVariableModuleAssociation[i]) + ",";
						ret+= to_string(observation_values[i]);
						ret+= ",";
					}
					return ret;
				}

                std::vector<double> get_observation_vector() const {
				    std::vector<double> ret;
				    for(uint i=0; i<size; i++) {
				        ret.push_back(observation_values[i]);
				    }
				    return ret;
				}
		};

		// Class to store qvalues
		class QValues {
			public:
			//... ATTRIBUTES ...
				// float values for the qvalues itself
				vector <float> qvalues;

				// number of qvalues
				uint size;

			//... FUNCTIONS ...// 
				// creates Qvalues, default size is 1
				QValues(uint isize=1) {
					// create new qvalues and initialize them with 0 everywhere
					size = isize;
					for(uint i=0; i<size; i++) {
						qvalues.push_back(0.0);
					}
				};

				// returns the index of the maximum value of the stored qvalues
				int get_max() const {
					// run through the qvalues and get the index of the maximum value
					float maxval = -numeric_limits<float>::max();
					int max_indx = -1;
					for(uint i=0; i<size; i++) {
						if (maxval<qvalues[i]) {
							maxval = qvalues[i];
							max_indx = i;
						}
					}
					return max_indx;
				};

				// creates a string
				std::string print() const {
					std::string ret = "[";
					for(uint i=0; i<size; i++) {
						ret += to_string(qvalues[i]) + ",";
					}
					ret += "]";
					return ret;
				};
		};


        // Class to store actions
        class Actions {
        public:
            //... ATTRIBUTES ...
            // vector of actions
            std::vector<uint32_t> acts;

            //... FUNCTIONS ...//
            // creates Actions
            Actions(vector<uint32_t> vec = {}) {
                // create new vector and push the values of the vector
                for(unsigned int i = 0; i < vec.size(); ++i)
                    acts.push_back(vec[i]);
            };

            // creates a string
            std::string print() const {
                std::string ret = "[";
                for(uint i=0; i<acts.size(); i++) {
                    ret += to_string(acts[i]) + ",";
                }
                ret += "]";
                return ret;
            };
        };

		// Class to store a vectyor of observations
        template <typename obValueType>
		class StateObservations {
			public:
				StateObservations() {};

				Observation<obValueType> getObservation(uint index) const {
					if (index<observations.size()){
						return observations[index];
					} else {
						return dtstrat::dt::Observation<obValueType>();
					}
				};

				bool existsStateObservation(uint index) const {
					bool ret = false;
					if (index<observations.size()) {
						ret = (observations[index].size>0);
					}
					return ret;
				};

				void pushback(Observation<obValueType> newObservation) {
				    observations.push_back(newObservation);
				    numDimensions = newObservation.size;
				}

				int findObservation(Observation<obValueType> newObservation) {
				    int returnval = -1;
				    for (uint i=0; i<observations.size(); i++) {
				        auto observed = observations[i];
				        if (observed==newObservation) {
				            returnval = i;
				            break;
				        }
				    }
				    return returnval;
				}

				void addObservation(std::vector<int64_t>& integerValues, std::vector<bool>& booleanValues, std::vector<float>& floatValues, uint state) {
                    if (numDimensions == 0) {
                        numDimensions+=integerValues.size();
                        numDimensions+=booleanValues.size();
                        numDimensions+=floatValues.size();
                    }
					if (state > observations.size()) {
                    	observations.resize(state);
					}
					Observation<obValueType> newObservation(this->numDimensions);
					std::vector<obValueType> values;
					values.resize(numDimensions);
					for (uint i=0; i<integerValues.size(); i++) {
                        values[i] = integerValues[i];
                    }
                    for (uint j=0;j<booleanValues.size(); j++) {
                        values[integerValues.size()+j] = booleanValues[j];
                    }
                    for (uint k=0; k<floatValues.size(); k++) {
                        values[integerValues.size()+booleanValues.size()+k] = floatValues[k];
                    }
					newObservation.set(values);
					if (state == observations.size()) {
						observations.push_back(newObservation);
					} else {
						observations[state] = newObservation;
					}
				}

				uint getSize() const {return observations.size();};

			private:
				std::vector<Observation<obValueType>> observations;
				uint numDimensions=0;
		};

        template <typename obValueType>
		struct HistoryElement{
			dtstrat::dt::Observation<obValueType> state;
			uint action;
			float deltaQ;
			float qValue;
		};

		// Class to store a node of a decision tree
		template <typename obValueType, typename qValueType>
		class Node {
			public:
			//... ATTRIBUTES ...
				// unique value
				uint num;

				// Left and right child of the node
				// default is nullptr
				Node* left_child;
				Node* right_child;
				Node* parent;

				// if the node is a branching node, the dimension gives the dimension to split on
				// the value gives the value to split on
				uint indexToSplit;
				float value;

				// determines whether the node is a leaf node (true) or a branching node (false)
				bool is_leaf = true;
				bool isLeft = false;

				// if the node is a leaf node, the qvalues must be set
				std::vector<qValueType> q_values;
				std::vector<uint> modules;

				// history of updates in the qvalues, needed to determine when t split (Pyeatt&Howe)
				std::vector<float> history;

                // store actions in case of leaf nodes
                Actions actions;

				// history of observations that lead to this node, needed for determining where to split (Pyeatt&Howe)
				std::vector<Observation<obValueType>> observation_history;

				// split history (how do we come to this node)
				std::vector<std::pair<int, int>> split_history;

				// new history list
				std::vector<dtstrat::dt::HistoryElement<obValueType>> complete_history;

				// storm state list
				std::vector<uint> stormStateList;

			//... FUNCTIONS ...
				// creates a leaf node with only zeros in the qvalues, default size of the qvalues is 1
				Node(uint indexToSplitIn=0, float valIn=0.0, uint size=1, uint inum=0, vector<uint32_t> vec = {}) {
					num = inum;
					indexToSplit = indexToSplitIn;
					value = valIn;
					q_values.resize(size);
                    actions = Actions(vec);
					left_child = nullptr;
					right_child = nullptr;
					parent = nullptr;
				};

				// evaluates the node to determine which child to visit, true is the left child, false is the right
				bool evaluate(Observation<obValueType> input_value) const {
					// check how the node would evaluate on some observation
					bool return_value;
					// decisive criterion is whether the value of the specific dimension is larger or equal than the defined value
					if (float(input_value[indexToSplit]) <= value) return_value = true;
					else return_value = false;
					return return_value;
				};

				// returns a pointer to the child node, given by a boolean value (true for left, false for right)
				Node* get_child(bool choice) {
					if (is_leaf) return this;
					if (choice) return left_child;
					else return right_child;
				};

				list <Observation<obValueType>> getObservations(const std::vector<dtstrat::dt::HistoryElement<obValueType>>& history) {
					list<Observation<obValueType>> result;
					for (auto el : history) {
						result.push_back(el.state);
					}
					return result;
				}

				// get and set functions to set the attributes 
				// NOTE: theoretically not needed, because currently they are public
				uint get_dim() const {return indexToSplit;};
				float get_val() const {return value;};
				void set_dim(uint idim) {indexToSplit=idim;};
				void set_val(float ival) {value=ival;};

				void set_qvals(std::vector<float> newvals) {
				    q_values.resize(newvals.size());
					for(uint i=0; i<newvals.size(); i++){
						q_values[i] = newvals[i];
					}
				};

                int get_max() const {
                    // run through the qvalues and get the index of the maximum value
                    float maxval;
                    maxval = -numeric_limits<float>::max();
                    int max_indx = -1;
                    for(uint i=0; i<q_values.size(); i++) {
                        if (maxval<q_values[i]) {
                            maxval = q_values[i];
                            max_indx = i;
                        }
                    }
                    return max_indx;
                }

				std::vector<uint> getSortedActions() const {
					std::vector<std::pair<float,uint>> qvaluesAndIndex;
					for (uint i = 0; i<q_values.size(); i++) {
						qvaluesAndIndex.push_back(std::make_pair(q_values[i],i));
					}
					std::sort(qvaluesAndIndex.begin(), qvaluesAndIndex.end(), greater <>()); 
					std::vector<uint> sortedActionIndices;
					std::vector<uint> actionIndicesForSameQValue;
					float lastQvalue = 10.0;
					unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    				std::default_random_engine e(seed);
					for (uint i = 0; i<q_values.size(); i++) {
						if (lastQvalue != qvaluesAndIndex[i].first) {
							std::shuffle (actionIndicesForSameQValue.begin(), actionIndicesForSameQValue.end(), e);
							for (auto element : actionIndicesForSameQValue){
								sortedActionIndices.push_back(element);
							}
							actionIndicesForSameQValue.clear();
						}
						actionIndicesForSameQValue.push_back(qvaluesAndIndex[i].second);
						lastQvalue = qvaluesAndIndex[i].first;
					}
					std::shuffle(actionIndicesForSameQValue.begin(), actionIndicesForSameQValue.end(), e);
					for (auto element : actionIndicesForSameQValue){
						sortedActionIndices.push_back(element);
					}
					return sortedActionIndices;
				};

				float updateQValue(uint action, float newValue) {
					float change = newValue - this->q_values[action];
					this->history.push_back(change);
					this->q_values[action] = newValue;
					return change;
				}

				void updateObservations(dtstrat::dt::Observation<obValueType> state) {
					this->observation_history.push_back(state);
				}

				void updateCompleteHistory(dtstrat::dt::Observation<obValueType> state, uint action, float change, float qval) {
					HistoryElement<obValueType> element{state, action, change, qval};
					this->complete_history.push_back(element);
				}

				void updateBRTDPList(uint stormStateId) {
					this->stormStateList.push_back(stormStateId);
				}

				std::string getDotLabel(uint k=0, bool dtcontrol=false) const {
				    ostringstream label_stream;
				    if (!is_leaf) {
                        label_stream << this->num << " [label=\"x_" << indexToSplit << " <= " << value << "\"];";
                    }
				    else {
				        std::vector<uint> sorted_actions;
				        if (dtcontrol) {
                            for (auto i : this->actions.acts) sorted_actions.push_back((uint) i);
				            //for (auto i : this->q_values) sorted_actions.push_back((uint) i);
				        }
				        else sorted_actions = getSortedActions();
				        label_stream << this->num << " [label=\"";
				        uint counter = 0;
				        for (uint i=0; i<sorted_actions.size(); i++) {
				            if (counter>=k) break;
				            auto val = sorted_actions.at(i);
				            label_stream << "actions[" << i << "] = " << val << ",";
				            counter++;
				        }
				        label_stream << "\"];";
				    }
				    return label_stream.str();
				}
		};

		// Class to store a decision tree
		template <typename obValueType, typename qValueType>
		class DecisionTree {
			private:
			//... ATTRIBUTES ...
				// stores the root node
				Node<obValueType, qValueType>* root;
				// stores the full size of the decision tree (number of nodes)
				int number_of_nodes=1;
				// for speed-up, store all leaf-nodes
				std::vector<dtstrat::dt::Node<obValueType, qValueType>*> leafNodes;

                // set if DT is created from a file/string (dot/json)
                bool isSetFromFile = false;

				int sizeOfNodes = 0;

            //... FUNCTIONS ...
                bool isLeafNode(nlohmann2::json tree) {
                    return (!tree.contains("split") || tree["split"] == nullptr);
                }

                std::tuple <int, float> getJSplitValues(nlohmann2::json split) {
                    auto op = split["op"];
                    if (op == "<=") {
                        // everything is fine
                    } else {
                        STORM_PRINT("We can't do any other operator than <="<<std::endl);
                    }
                    //std::string rh;
                    float splitval=-1;
                    if (!split.contains("rhs") || split["rhs"] ==nullptr) {
                        STORM_PRINT(" No rhs... " << split << std::endl);
                    }
                    else {
                        split["rhs"].get_to(splitval);
                    }

                    // split is either an axis-parallel expression or a list of linear terms
                    float coeff=1;
                    int varDim = -1;
                    if (!split.contains("lhs") || split["lhs"]==nullptr) {
                        STORM_PRINT(" No lhs... " << split << std::endl);
                    }
                    else {
                        auto lhs = split["lhs"];
                        float coeff;
                        if (!lhs.contains("coeff") || lhs["coeff"]==nullptr) {
                            STORM_PRINT(" No coeff... " << split << std::endl);
                        }
                        else {
                            lhs["coeff"].get_to(coeff);
                        }
                        std::string var;
                        if (!lhs.contains("var") || lhs["var"] == nullptr) {
                            STORM_PRINT(" No var " << split << std::endl);
                        }
                        else {
                            lhs["var"].get_to(var);
                            auto bla = var.substr(var.find("_")+1);
                            varDim = std::stoi(var.substr(var.find("_")+1));
                        }
                    }
                    //STORM_PRINT("Split on " << varDim << " with value " <<splitval/coeff <<" corresp " << split<<std::endl);
                    std::tuple<int,float> ret = {varDim, splitval/coeff};
                    return ret;
                };

                void get_subtree(nlohmann2::json tree, Node<obValueType, qValueType>& currentNode) {
                    if (!isLeafNode(tree)) {
                        auto splitvals = getJSplitValues(tree["split"]);
                        int dimension = get<0>(splitvals);
                        float value = get<1>(splitvals);
                        this->split(currentNode, dimension, value);
                        auto childs =  tree["children"];
                        std::vector<nlohmann2::json> children;
                        for (auto& child : childs) {
                            children.push_back(child);
                        }
                        STORM_LOG_ASSERT((children.size()<3), "More than two children!");
                        get_subtree(children[0], *(currentNode.left_child));
                        get_subtree(children[1], *(currentNode.right_child));
                    }
                    else {
                        int action = 0;
                        if (!tree.contains("actual_label") || tree["actual_label"] ==nullptr) {
                            STORM_PRINT("No actual label " << tree << std::endl);
                        }
                        else {
                            std::vector<nlohmann2::json> labels;
                            for (auto& label : tree["actual_label"]) {
                                labels.push_back(label);
                            }
                            //cout<< "label size: " << labels.size() << endl;
                            //STORM_LOG_ASSERT((labels.size()<2), "Too many labels!");
                            if(labels.size()>=2) {
                                /* if label size (number of actions) is greater or equal to 2, store all the actions */
                                std::random_device rd;
                                std::mt19937 g(rd());
                                // std::shuffle(labels.begin(), labels.end(), g);
                                for (uint i = 0; i < labels.size(); ++i) {
                                    std::string action_string = labels[i].dump();
                                    action_string = action_string.substr(1, action_string.size() - 1);
                                    action = std::stoi(action_string);
                                    currentNode.actions.acts.push_back(action);
                                    //cout << action_string << endl;
                                }
                            }
                            else{
                                std::string action_string = labels[0].dump();
                                action_string = action_string.substr(1, action_string.size() - 1);
                                action = std::stoi(action_string);
                                currentNode.actions.acts.push_back(action);
                                //cout << action_string << endl;
                            }
                        }
                    }
                }

                std::string getDotChildString(Node<obValueType, qValueType>* n, uint k=0) const {
                    std::string result;
                    if (n->is_leaf) {
                        result = n->getDotLabel(k, this->fromDtControl) + "\n";
                    }
                    else {
                        result = n->getDotLabel(k, this->fromDtControl) + "\n";
                        result += std::to_string(n->num) + " -> " + std::to_string(n->left_child->num) +
                                  " [label=\"true\"];\n";
                        result += getDotChildString(n->left_child, k);
                        result += std::to_string(n->num) + " -> " + std::to_string(n->right_child->num) +
                                  " [label=\"false\"];\n";
                        result += getDotChildString(n->right_child, k);
                    }
                    return result;
                }

			public:
                bool fromDtControl = false;
			//... FUNCTIONS ...
				// creates a decision tree, where size determines the number of qvalues 
				DecisionTree(int size) : isSetFromFile(false) {
					this->root = new Node<obValueType, qValueType>(0,0,size);
					this->leafNodes.push_back(root);
				};

				Node<obValueType, qValueType> getRoot(){
					return *(this->root);
				}

				DecisionTree() : isSetFromFile(false){
				    this->number_of_nodes=0;
				};

				explicit DecisionTree(std::string dot_file, bool module=false) : isSetFromFile(true) {
                    std::vector<Node<obValueType, qValueType>*> nodes;
                    std::vector<std::tuple<int,int,bool>> edges;
                    std::ifstream infile(dot_file);
                    std::string line;
                    int linecount = 0;
                    while (std::getline(infile, line)) {
                        if (linecount == 0) {
                            linecount++;
                            continue;
                        }
                        std::istringstream iss(line);
                        uint_fast64_t nodeIndex;
                        iss >> nodeIndex;
                        char symbol;
                        if (!(iss >> symbol)) { break; }
                        if (symbol == '[') {
                            // node definition
                            std::string label;
                            iss >> label;
                            if (label[7] == 'x') {
                                // decision node
                                int dec_dim = std::stoi(label.substr(9));
                                iss >> symbol;
                                std::string rest;
                                iss >> rest;
                                iss >> rest;
                                std::size_t found = rest.find('"');
                                float dec_val = std::stof(rest.substr(0, found));
                                Node<obValueType, qValueType>* newNode = new Node<obValueType, qValueType>(dec_dim, dec_val);
                                if(nodes.size()<=nodeIndex) {
                                    nodes.resize(nodeIndex);
                                    nodes.push_back(newNode);
                                } else {
                                    nodes[nodeIndex] = newNode;
                                }
                            } else if (label[7] == 'a') {
                                // leaf node
                                Node<obValueType, qValueType>* newNode = new Node<obValueType, qValueType>();
                                newNode->q_values.resize(0);
                                while (iss >> symbol && (symbol=='=' || symbol==',')) {
                                    std::string rest;
                                    iss >> rest;
                                    std::size_t found = rest.find(',');
                                    uint action = std::stoi(rest.substr(0, found));
                                    newNode->q_values.push_back(action);
                                    newNode->actions.acts.push_back(action);
                                }
                                if(nodes.size()<=nodeIndex) {
                                    nodes.resize(nodeIndex);
                                    nodes.push_back(newNode);
                                } else {
                                    nodes[nodeIndex] = newNode;
                                }
                            } else if (label[7] == 'm' && module) {
                                // leaf node
                                Node<obValueType, qValueType>* newNode = new Node<obValueType, qValueType>();
                                newNode->q_values.resize(0);
                                newNode->is_leaf=true;
                                // read modules
                                while (iss >> symbol && (symbol=='=' || symbol==',')) {
                                    std::string rest;
                                    iss >> rest;
                                    while(rest.size()>1){
                                        std::size_t found = rest.find(',');
                                        if (found==std::string::npos) {break;}
                                        uint action = std::stoi(rest.substr(0, found));
                                        rest = rest.substr(found+1);
                                        newNode->modules.push_back(action);
                                    }
                                }
                                iss >> label;
                                while (iss >> symbol && (symbol=='=' || symbol==',')) {
                                    std::string rest;
                                    iss >> rest;
                                    while(rest.size()>1){
                                        std::size_t found = rest.find(',');
                                        if (found==std::string::npos) {break;}
                                        uint action = std::stoi(rest.substr(0, found));
                                        rest = rest.substr(found+1);
                                        newNode->q_values.push_back(action);
                                        newNode->actions.acts.push_back(action);
                                    }
                                }
                                if(nodes.size()<=nodeIndex) {
                                    nodes.resize(nodeIndex);
                                    nodes.push_back(newNode);
                                } else {
                                    nodes[nodeIndex] = newNode;
                                }
                            }
                            else { std::cout << "ERROR " << label[7] << std::endl; }
                        }
                        // line starts with -> to define an edge
                        else if (symbol == '-') {
                            // edge definition
                            iss >> symbol;
                            // consumes >
                            int node2;
                            iss >> node2;
                            // defines the "target node"
                            bool left;
                            std::string rest;
                            iss >> rest;
                            // if the edge is annotated with "false", it's the right child, if "true", it's the left child
                            left = (rest[8]=='t');
                            edges.push_back(std::make_tuple(nodeIndex, node2, left));
                        }
                        linecount++;
                    }
                    for(auto edge: edges) {
                        int nodeIndex = std::get<0>(edge);
                        int node2 = std::get<1>(edge);
                        bool left = std::get<2>(edge);
                        if (left) {
                            nodes[nodeIndex]->left_child = nodes[node2];
                            nodes[nodeIndex]->is_leaf = false;
                            nodes[node2]->parent = nodes[nodeIndex];
                            nodes[node2]->isLeft = true;
                        }
                        else {
                            nodes[nodeIndex]->right_child = nodes[node2];
                            nodes[nodeIndex]->is_leaf = false;
                            nodes[node2]->parent = nodes[nodeIndex];
                        }
                    }
                    this->number_of_nodes = nodes.size();
                    this->root = nodes[0];
                }

				explicit DecisionTree(nlohmann2::json tree) : isSetFromFile(true) {
					this->root = new Node<obValueType, qValueType>(0,0,1);
					this->sizeOfNodes = 1;
					this->leafNodes.push_back(root);
					get_subtree(tree, *root);
				}

                int get_number_of_nodes() const {
				    return number_of_nodes;
				}

                bool setFromFile() const {
                    return this->isSetFromFile;
                }

				void clear() {
					if (this->number_of_nodes<2) {
						return;
					}
					//STORM_PRINT("We clear the DT!" << std::endl);
					std::vector<Node<obValueType, qValueType>*> stack;
					stack.push_back(root);
					DFS_clear(this->root, stack);
					this->number_of_nodes = 1;
					this->leafNodes.clear();
					this->root = new Node<obValueType, qValueType>(0,0,this->sizeOfNodes);
				}

				void DFS_clear(Node<obValueType, qValueType>* currentNode, std::vector<Node<obValueType, qValueType>*>& stack) {
					if (!(currentNode->left_child) && !(currentNode->right_child)) {
						// node has no children -> leaf node
						// remove it
						if (!currentNode->parent) {
							delete currentNode;
							stack.pop_back();
							return;
						}
						if (currentNode->isLeft) {
							currentNode->parent->left_child=nullptr;
						}
						else {
							currentNode->parent->right_child=nullptr;
						}
						delete currentNode;
						stack.pop_back();
					} 
					else {
						stack.push_back(currentNode->left_child);
						stack.push_back(currentNode->right_child);
					}
					if (stack.size()>0) {
						auto nextNode = stack.back();
						DFS_clear(nextNode, stack);
					}
				}

				// returns the result of the decision tree for an observation, that is the index of the maximum qvalue
                // of the corresponding leaf node
				// NOTE: giving the index only does not make a lot of sense, does it?
				int get_result(Observation<obValueType> input_value) const {
					bool reached_end = false;
					Node<obValueType, qValueType> current_node = *root;
					// runs through the decision tree
					while (!reached_end) {
						// if it has reached a leaf node, it is done
						if (current_node.is_leaf) reached_end = true;
						// if it is in a branching node, evaluate the node and move to the respective child node
						else {
							current_node = *current_node.get_child(current_node.evaluate(input_value));
						}
					}
					return current_node.get_max();
				};

                int get_decision(Observation<obValueType> input_values) const {
                    bool reached_end = false;
                    Node<obValueType, qValueType> current_node = *root;
                    // runs through the decision tree
                    while (!reached_end) {
                        // if it has reached a leaf node, it is done
                        if (current_node.is_leaf) reached_end = true;
                            // if it is in a branching node, evaluate the node and move to the respective child node
                        else {
                            current_node = *current_node.get_child(current_node.evaluate(input_values));
                        }
                    }
                    //return current_node.q_values.qvalues[0];
                    return current_node.actions.acts[0];
                }

                vector<uint32_t> get_decisions(Observation<obValueType> input_values) const {
                    bool reached_end = false;
                    Node <obValueType, qValueType> current_node = *root;
                    // runs through the decision tree
                    while (!reached_end) {
                        // if it has reached a leaf node, it is done
                        if (current_node.is_leaf) reached_end = true;
                            // if it is in a branching node, evaluate the node and move to the respective child node
                        else {
                            current_node = *current_node.get_child(current_node.evaluate(input_values));
                        }
                    }
                    return current_node.actions.acts;
                }

				// returns a sorted list of indices for the qvalues
				// that is returns the possible actions sorted by their qvalue
				std::vector<uint> get_results(Observation<obValueType> input_value) const {
					bool reached_end = false;
					Node<obValueType, qValueType> current_node = *root;
					// runs through the decision tree
					while (!reached_end) {
						// if it has reached a leaf node, it is done
						if (current_node.is_leaf) reached_end = true;
						// if it is in a branching node, evaluate the node and move to the respective child node
						else {
							current_node = *current_node.get_child(current_node.evaluate(input_value));
						}
					}
					return current_node.getSortedActions();
				};

				// splits a node, by giving it the dimension and the value to split on
				void split(Node<obValueType, qValueType>& node, int dimension, float splitvalue) {
					// Create new child nodes
					Node<obValueType, qValueType>* left_node = new Node<obValueType, qValueType>(0, 0.0, node.q_values.size(), number_of_nodes);
					//cout << "left node number " << left_node->num << endl;
					left_node->split_history = std::vector<std::pair<int, int>>(node.split_history);
					left_node->split_history.push_back(std::make_pair(dimension,int(splitvalue)));
					left_node->parent = &node;
					left_node->isLeft = true;
					number_of_nodes++;
					Node<obValueType, qValueType>* right_node = new Node<obValueType, qValueType>(0, 0.0, node.q_values.size(), number_of_nodes);
					right_node->split_history = std::vector<std::pair<int, int>>(node.split_history);
					right_node->split_history.push_back(std::make_pair(dimension,int(splitvalue)));
					right_node->parent = &node;
					//cout << "right node number " << right_node->num << endl;
					number_of_nodes++;
					this->leafNodes.erase(std::remove(this->leafNodes.begin(), this->leafNodes.end(), &node), this->leafNodes.end());
					this->leafNodes.push_back(left_node);
					this->leafNodes.push_back(right_node);

					// the children get the same qvalues as their parent node
					// they don't. In our application, the qvalues would never decrease, so better learn again!?
					//left_node->set_qvals(node.q_values);
					//right_node->set_qvals(node.q_values);

					// the parent node is no longer a leaf node
					node.is_leaf = false;
					node.observation_history.clear();
					node.history.clear();

					// set dimension and v	alue to split on
					node.set_dim(dimension);
					node.set_val(splitvalue);

					// add the new nodes as the children to the parent node
					node.left_child = left_node;
					node.right_child = right_node;
					//std::ofstream logger;
                    //logger.open ("DEBUG.OUT", std::ios_base::app);
                    //logger << " " << node.q_values.print() <<"\n";
                    //logger.close();
				};

				// similar to get_result, but outputs the corresponding leaf-node
				Node<obValueType, qValueType>* get_leaf_node(Observation<obValueType> input_value) {
					bool reached_end = false;
					Node<obValueType, qValueType>* current_node_adr = root;
					// runs through the decision tree
					while (!reached_end) {
						// if it has reached a leaf node, it is done
						if (current_node_adr->is_leaf) reached_end = true;
						// if it is in a branching node, evaluate the node and move to the respective child node
						else {
							bool evaluation = (*current_node_adr).evaluate(input_value);
							Node<obValueType, qValueType>* const cl = (*current_node_adr).get_child(evaluation);
							current_node_adr = cl;
						}
					}
					return current_node_adr;
				};

				std::vector<dtstrat::dt::Node<obValueType, qValueType>*> getLeafNodes() const {
					return this->leafNodes;
				}

				std::string getDotString(uint k=0) const {
				    std::string result = "digraph {\n";
				    Node<obValueType, qValueType>* n = root;
				    result += getDotChildString(n, k);
				    result += "}\n";
				    return result;
				}

				int size() const {
					return this->number_of_nodes;
				}
		};
	}
}
#endif
