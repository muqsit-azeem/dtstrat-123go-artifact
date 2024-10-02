//
// Created by steffi on 13.12.21.
//

#pragma once
#include <vector>
#include "storm/dtstrat/decisionTree.h"
#include "storm/storage/prism/Module.h"

namespace dtstrat {
    enum moduleChoice { random, treeSize, smallest };

    template <typename obValueType, typename qValueType>
    class DTConfig {
    private:
        std::vector<dtstrat::dt::DecisionTree<obValueType, qValueType>> decisonTrees;
        std::map<uint, uint> moduleDTMapping;
        std::vector<uint_fast64_t> moduleOrdering;
        std::vector<bool> simpleDT;

        std::vector<uint> getDTStratActionFromModuleAction(
                std::map<std::string, std::map<uint, std::vector<uint>>> const& moduleActionToDTStratMapping,
                std::string const &module, std::vector<uint> const& actions) const {
            std::vector<uint> ret;
            for (uint const &action: actions) {
                if (moduleActionToDTStratMapping.find(module) != moduleActionToDTStratMapping.end()) {
                    std::map<uint, std::vector<uint>> cmdIndicesActionMapping = moduleActionToDTStratMapping.at(module);
                    if (cmdIndicesActionMapping.find(action) != cmdIndicesActionMapping.end()) {
                        auto bla = cmdIndicesActionMapping.at(action);
                        ret.insert(ret.end(), bla.begin(), bla.end());
                    } else {
                        STORM_PRINT("Module has no such action!" << std::endl);
                    }
                } else {
                    STORM_PRINT("Module has no actions!" << std::endl);
                }
            }
            return ret;
        }

        std::vector<uint_fast64_t> getModuleOrderToPlay(moduleChoice const& moduleMethodChoice) {
            std::vector<std::pair<uint,int>> availableModules;
            for (uint i=0; i<decisonTrees.size(); i++) {
                auto const& dt = decisonTrees[i];
                if (dt.setFromFile()) {
                    availableModules.push_back(std::make_pair(i,dt.size()));
                }
            }
            if (moduleMethodChoice==random) {
                std::default_random_engine generator;
                std::shuffle(std::begin(availableModules), std::end(availableModules), generator);
            } else if (moduleMethodChoice==treeSize) {
                std::sort(availableModules.begin(), availableModules.end(),
                          [] (std::pair<uint,int> const& a, std::pair<uint,int> const& b) { return a.second>b.second; } );
            } else if (moduleMethodChoice==smallest) {
                std::sort(availableModules.begin(), availableModules.end(),
                          [] (std::pair<uint,int> const& a, std::pair<uint,int> const& b) { return a.second<b.second; } );
            } else {
                STORM_PRINT("Not impelemnted module choice " << moduleMethodChoice << std::endl);
            }
            std::vector<uint_fast64_t> modules;
            modules.reserve(availableModules.size());
            for (auto const& element : availableModules) {
                modules.push_back(element.first);
            }
            return modules;
        }

    public:
        DTConfig() {};

        explicit DTConfig(std::vector<storm::prism::Module> const& allModules) {
            std::iota(moduleOrdering.begin(), moduleOrdering.end(), 0);
        };

        explicit DTConfig(std::map<uint, dtstrat::dt::DecisionTree<obValueType, qValueType>&> const& dts) {
            moduleDTMapping.clear();
        }

        DTConfig(DTConfig<obValueType, qValueType> const& config) {
            setDecisionTrees(config.getDecisionTrees());
            setModuleDTMapping(config.getModuleDTMapping());
            setModuleOrdering(config.getModuleOrdering());
        }

        std::vector<dtstrat::dt::DecisionTree<obValueType, qValueType>> getDecisionTrees() const {
            return decisonTrees;
        }

        std::map<uint, uint> getModuleDTMapping() const {
            return moduleDTMapping;
        }

        std::vector<uint_fast64_t> getModuleOrdering() const {
            return moduleOrdering;
        }

        void allSimple() {
            for (uint i=0; i<decisonTrees.size(); i++) {
                simpleDT[i] = true;
            }
        }

        void setDecisionTrees(std::vector<dtstrat::dt::DecisionTree<obValueType, qValueType>> const& newDecisonTrees) {
            decisonTrees = newDecisonTrees;
            simpleDT = std::vector<bool>(decisonTrees.size(), false);
        }

        void setModuleOrdering(std::vector<uint_fast64_t> const& newOrdering) {
            moduleOrdering = newOrdering;
        }

        void setModuleDTMapping(std::map<uint, uint> const& newMapping) {
            moduleDTMapping = newMapping;
        }

        void setModuleDTMapping(uint const& numberOfModules, uint const& dtIndex) {
            moduleDTMapping.clear();
            for (uint_fast64_t i=0; i<numberOfModules; i++) {
                moduleDTMapping[i] = dtIndex;
            }
        }

        void setModuleOrdering(moduleChoice const& modChoice) {
            std::vector<uint_fast64_t> moduleOrder = getModuleOrderToPlay(modChoice);
            setModuleOrdering(moduleOrder);
        }

        void setOrderOfModule(uint module, uint order) {
            auto currentPlace = std::find(moduleOrdering.begin(),moduleOrdering.end(), module);
            if (currentPlace !=moduleOrdering.end()) {
                if (currentPlace-moduleOrdering.begin() != order) {
                    moduleOrdering.erase(currentPlace);
                } else {
                    return;
                }
            }
            moduleOrdering.insert(moduleOrdering.begin()+order, module);
        }

        std::vector<uint> getDecision(dtstrat::dt::Observation<double> state,
                                      std::vector<storm::prism::Module> const& allModules,
                                      std::map<std::string, std::map<uint, std::vector<uint>>> const& moduleActionToDTStratMapping) const {
            std::vector<uint> possibleActionIds;
            assert(!decisonTrees.empty());
            if (simpleDT[0]) {
                possibleActionIds = {decisonTrees[0].get_decision(state)};
            } else {
                for (auto const &module: moduleOrdering) {
                    vector<uint32_t> decisions = decisonTrees[moduleDTMapping.at(module)].get_decisions(state);
                    std::string moduleName = allModules[module].getName();
                    std::vector<uint> possibleActions = getDTStratActionFromModuleAction(moduleActionToDTStratMapping,
                                                                                         moduleName, decisions);
                    possibleActionIds.insert(possibleActionIds.end(), possibleActions.begin(),
                                             possibleActions.end());
                }
            }
            return possibleActionIds;
        }

        uint getHighestOrderModule() {
            if (moduleOrdering.empty()) return 0;
            else return moduleOrdering[0];
        }

        uint getDTofHighestOrderModule() {
            if (moduleDTMapping.empty()) return 0;
            else return moduleDTMapping[getHighestOrderModule()];
        }

        uint numberOfDTs() const {
            return decisonTrees.size();
        }

        std::vector<float> getSplit() const {
            std::vector<float> res(decisonTrees.size());
            for (auto const& module : moduleOrdering) {
                res[moduleDTMapping.at(module)]+=1;
            }
            for (auto& el : res) {
                el = el/moduleOrdering.size();
            }
            return res;
        }

        friend bool operator== (const DTConfig<obValueType, qValueType>& config1, const DTConfig<obValueType, qValueType>& config2) {
            if (config2.getDecisionTrees().size() != config1.getDecisionTrees().size() ||
                    config2.getModuleDTMapping().size() != config1.getModuleDTMapping().size() ||
                    config2.getModuleOrdering().size() != config1.getModuleOrdering().size()) {
                return false;
            } else {
                auto otherOrdering = config1.getModuleOrdering();
                auto thisOrdering = config2.getModuleOrdering();
                for (uint i=0; i<thisOrdering.size();i++) {
                    if (otherOrdering[i]!=thisOrdering[i]) {
                        return false;
                    }
                }
                auto otherMapping = config1.getModuleDTMapping();
                for (auto const& element :  config2.getModuleDTMapping()) {
                    if (otherMapping.find(element.first)==otherMapping.end()) {
                        return false;
                    } else if (!(otherMapping[element.first] ==element.second)) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool operator== (const DTConfig<obValueType, qValueType>& config) {
            if (decisonTrees.size() != config.getDecisionTrees().size() ||
            moduleDTMapping.size() != config.getModuleDTMapping().size() ||
            moduleOrdering.size() != config.getModuleOrdering().size()) {
                return false;
            } else {
                auto otherOrdering = config.getModuleOrdering();
                for (uint i=0; i<moduleOrdering.size();i++) {
                    if (otherOrdering[i]!=moduleOrdering[i]) {
                        return false;
                    }
                }
                auto otherMapping = config.getModuleDTMapping();
                for (auto const& element : moduleDTMapping) {
                    if (otherMapping.find(element.first)==otherMapping.end()) {
                        return false;
                    } else if (otherMapping[element.first]!=element.second) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool operator< (const DTConfig<obValueType, qValueType>& config) {
            if (decisonTrees.size() < config.getDecisionTrees().size() ||
                moduleDTMapping.size() < config.getModuleDTMapping().size() ||
                moduleOrdering.size() < config.getModuleOrdering().size()) {
                return true;
            } else {
                auto otherOrdering = config.getModuleOrdering();
                for (uint i=0; i<moduleOrdering.size();i++) {
                    if (otherOrdering[i]<moduleOrdering[i]) {
                        return true;
                    }
                }
                auto otherMapping = config.getModuleDTMapping();
                for (auto const& element : moduleDTMapping) {
                    if (otherMapping.find(element.first)<otherMapping.end()) {
                        return false;
                    } else if (otherMapping[element.first]<element.second) {
                        return false;
                    }
                }
            }
            return false;
        }

    };

    template <typename obValueType, typename qValueType>
    class ConfigSet {
    public:
        typedef typename std::vector<DTConfig<obValueType, qValueType>>::iterator iterator;
        typedef typename std::vector<DTConfig<obValueType, qValueType>>::const_iterator const_iterator;

        ConfigSet() = default;

        void insert(DTConfig<obValueType, qValueType> const& config) {
            for (auto const& element : configs) {
                if (element==config) {
                    return;
                }
            }
            configs.push_back(config);
        }

        const_iterator begin() {
            return configs.begin();
        }

        const_iterator end() {
            return configs.end();
        }

        bool empty() {
            return configs.empty();
        }

        void resetDTs(std::vector<dtstrat::dt::DecisionTree<obValueType, qValueType>> const& newDTs) {
            for (auto it = configs.begin(); it<configs.end(); it++) {
                (*it).setDecisionTrees(newDTs);
            }
        }

        DTConfig<obValueType, qValueType> getBestConfig() const {
            assert(currentBest<configs.size());
            return configs[currentBest];
        }


    private:
        std::vector<DTConfig<obValueType, qValueType>> configs;
        uint currentBest = 0;
    };

}