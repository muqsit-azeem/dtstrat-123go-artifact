//
// Created by steffi on 13.12.21.
//

#pragma once
#include <string>
#include "storm/dtstrat/decisionTree.h"

namespace storm{
    namespace modelchecker {
        namespace dtlearning {
            std::string exec(const char *cmd) {
                std::array<char, 128> buffer{};
                std::string result;
                std::unique_ptr < FILE, decltype(&pclose) > pipe(popen(cmd, "r"), pclose);
                if (!pipe) {
                    throw std::runtime_error("popen() failed!");
                }
                while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
                    result += buffer.data();
                }
                return result;
            }

            template<typename ValueType>
            dtstrat::dt::DecisionTree<ValueType, int> learnDTControlTree(std::string const &state_action_csv, std::string jsonFile="DtControl_DT.json") {
                /* for now dtcontrol setup is manual: see storm/dtstrat/dtcontrol.sh for instructions*/
                std::string cmd = ". venv/bin/activate; "
                                  "dtcontrol "
                                  "--input " + state_action_csv +
                                  " --use-preset minnorm "
                                  "--output stdout:json "
                                  "--rerun "
                                  "2>\\dev\\null";
                //cout << "calling dtcontrol : "<<cmd;
                const char *run_cmd = cmd.c_str();
                const std::string &result = exec(run_cmd);
                STORM_LOG_ASSERT(result != "", "dtcontrol not installed? : \nPlease run the following commands"
                                               " to install dtcontrol in the \"Working directory\" of storm:\n"
                                               "#   python -m venv venv\n"
                                               "#   source venv/bin/activate\n"
                                               "#   pip install dtcontrol\n"
                );
                auto start = result.find("START") + 5;
                auto end = result.find("END");
                const std::string jsonStringOfDT = result.substr(start, end - start);
                auto mochs = storm::settings::getModule<storm::settings::modules::ModelCheckerSettings>();
                if (mochs.getPrintDTControlDT()) {
                    std::ofstream logger;
                    logger.open(jsonFile);
                    logger << jsonStringOfDT;
                    logger.close();
                }
                nlohmann2::json json_obj = nlohmann2::json::parse(jsonStringOfDT);
                auto tee = dtstrat::dt::DecisionTree<ValueType, int>(json_obj);
                tee.fromDtControl = true;
                return tee;
            };

            void learnSklearnTree(std::string const& state_action_csv) {
                std::string cmd = ". venv/bin/activate; "
                                 "python decisiontree_sklearn.py "
                                 + state_action_csv;

                cout << "running command: " << cmd << endl;

                const char *run_cmd = cmd.c_str();
                const std::string &result = exec(run_cmd);
            }

            template<typename ValueType>
            std::vector<dtstrat::dt::DecisionTree<ValueType, int>> learnMultiSklearnTrees(std::string const& state_action_csv,
                                        std::vector<storm::prism::Module> const& all_modules, bool const& onlyModuleSVs) {
                std::vector<dtstrat::dt::DecisionTree<ValueType, int>> multiDts;
                for (uint i = 0; i < all_modules.size(); i++) {
                    std::string filename = "DT_module_" + std::to_string(i) + ".dot";
                    ifstream f(filename.c_str());
                    if (f.good()) {
                        remove(filename.c_str());
                    }
                }
                std::string cmd = ". venv/bin/activate; "
                                  "python module_dt.py "
                                  + state_action_csv;
                if (onlyModuleSVs) cmd += " statesPerModule";
                cout << "running command: " << cmd << endl;
                const char *run_cmd = cmd.c_str();
                const std::string &result = exec(run_cmd);
                for (uint i = 0; i < all_modules.size(); i++) {
                    try {
                        std::string filename = "DT_module_" + std::to_string(i) + ".dot";
                        ifstream f(filename.c_str());
                        if (f.good()) {
                            dtstrat::dt::DecisionTree<ValueType, int> first = dtstrat::dt::DecisionTree<ValueType, int>(
                                    filename);
                            multiDts.push_back(first);
                        } else {
                            dtstrat::dt::DecisionTree<ValueType, int> first = dtstrat::dt::DecisionTree<ValueType, int>(
                                    1);
                            multiDts.push_back(first);
                        }
                    } catch (const std::exception &e) {
                        STORM_LOG_THROW(true, storm::exceptions::NotSupportedException, "Problem!!!");
                    }
                }
                return multiDts;
            }
        }
    }
}