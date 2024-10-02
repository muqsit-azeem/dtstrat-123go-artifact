#ifndef DT_LOG_H
#define DT_LOG_H

#include <iostream>
#include <fstream>
#include "storm/dtstrat/decisionTree.h"
#include <ctime>
#include <string>

const std::string filename = "results.txt";
const std::string debugfilename = "debug.csv";
const std::string expResFilename = "exp-res.csv";
const std::string lowerBoundsFilename = "lower-bound.csv";
const std::string upperBoundsFilename = "upper-bound.csv";

// Log functionalities to print debug information or experimental results
namespace dtstrat {
    namespace log {
        // to log the experiments
        class logExp {
            public:
                // Hyperparameter
                // stores if the DTStrat-algorithm has to be used at all
                bool dtstrat;
                // probability to determine how likely we want to do "random exploration" (=BRTDP)
                double epsilon = 0.0;
                // learning rate
                double alpha = 0.0;
                // discount factor, always 1 in reachability setting
                double gamma = 0.0;
                // P&H: how often does a node to be visited before we want to update it
                uint histmin = 0;
                // determines how often the update should be performed
                uint updatesteps = 0;
                // if the true Q-value calculation was used
                bool trueQ = true;
                // batchupdate
                bool batchupdate = false;
                // splitdecision
                string splitdecision = "";

                // Results from experiments
                // How many simulations were run
                uint iterations = 0;
                // how many steps (state to state) in total
                uint total = 0;
                // how long the algorithm took
                double totaltime = 0.0;
                // size of the DT
                uint sizeOfDT = 0;
                // the model that was used
                std::string model = "";
                // the property that was verified
                std::string property = "";
                // lower bound for the initial state
                double lower_b_init;
                double upper_b_init;

                // get the modelname from the complete path
                void setModel(std::string value) {
                    char delimiter = '/';
                    std::string tmp = "";
                    std::vector<std::string> parts;

                    std::istringstream iss(value);
                    while(std::getline(iss, tmp, delimiter)) {
                        parts.push_back(tmp);
                    }
                    this-> model = parts.back();
                }

                // write the experimental results into a file
                void writeToFile(){
                    std::ofstream logger;
                    logger.open (expResFilename, std::ios_base::app);
                    time_t rawtime;
                    struct tm * timeinfo;
                    char buffer[80];

                    time (&rawtime);
                    timeinfo = localtime(&rawtime);

                    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
                    std::string str(buffer);
                    logger << str << ";" << model << ";" << property << ";" << dtstrat << ";" << epsilon  << ";" << alpha << ";" << gamma << ";" << histmin << ";" << updatesteps << ";" << trueQ << ";" << batchupdate << ";" << iterations << ";" << total << ";" <<  totaltime << "\n";
                    logger.close();
                    return;
                }

                void writeTitle() {
                    std::ifstream logfile(expResFilename);
                    if (logfile.is_open()) {
                        std::string line;
                        if (!std::getline(logfile, line)) {
                            logfile.close();
                            writeTitleToFile();
                        }
                    } else {
                        writeTitleToFile();
                    }
                }

                void writeTitleToFile() {
                    std::ofstream logger;
                    logger.open (expResFilename);
                    logger << "Time;Model;Property;Dtstrat;epsilon;alpha;gamma;histmin;updatesteps;trueQ;batchupdate;splitdecision;sizeOfTheDT;iterations;total_iterations;totaltime[ms]\n";
                    logger.close();
                }

                void writeFirstPart() {
                    checkLastLine();
                    std::ofstream logger;
                    logger.open (expResFilename, std::ios_base::app);
                    time_t rawtime;
                    struct tm * timeinfo;
                    char buffer[80];

                    time (&rawtime);
                    timeinfo = localtime(&rawtime);

                    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
                    std::string str(buffer);
                    std::string prop = property;
                    if (property.length() > 30){
                        prop = property.substr(0,30);
                    }
                    std::replace( prop.begin(), prop.end(), ';', '/');
                    logger << str << ";" << model << ";" << prop << ";" << dtstrat << ";" << epsilon  << ";" << alpha << ";" << gamma << ";" << histmin << ";" << updatesteps << ";" << trueQ << ";" << batchupdate << ";" << splitdecision << ";";
                    logger.close();
                    return;
                }

                void printFirstPart() {
                    time_t rawtime;
                    struct tm * timeinfo;
                    char buffer[80];

                    time (&rawtime);
                    timeinfo = localtime(&rawtime);

                    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
                    std::string str(buffer);
                    std::cout << "EXPERIMENTRESULTS" << std::endl;
                    std::string prop = property;
                    if (property.length() > 30){
                        prop = property.substr(0,30);
                    }
                    std::replace( prop.begin(), prop.end(), ';', '/');
                    std::cout << str << ";" << model << ";" << prop << ";" << dtstrat << ";" << epsilon  << ";" << alpha << ";" << gamma << ";" << histmin << ";" << updatesteps << ";" << trueQ << ";" << batchupdate << ";" << splitdecision << ";" << std::endl;
                    return;
                }

                void writeSecondPart() {
                    std::ofstream logger;
                    logger.open(expResFilename, std::ios_base::app);
                    logger << sizeOfDT << ";" << iterations << ";" << total << ";" <<  totaltime << "\n";
                    logger.close();
                }

                void printSecondPart() {
                    std::cout << "EXPERIMENTRESULTS" << std::endl;
                    std::cout << sizeOfDT << ";" << iterations << ";" << total << ";" <<  totaltime << std::endl;
                }

                void  checkLastLine() {
                    std::ifstream logfile(expResFilename);
                    if (logfile.is_open()) {
                        logfile.seekg(-1, std::ios_base::end);
                        if(logfile.peek() == '\n') {
                            logfile.seekg(-1, std::ios_base::cur);
                            int i = logfile.tellg();
                            for(; i>0; i--) {
                                if(logfile.peek() == '\n') {
                                    logfile.get();
                                    break;
                                }
                                logfile.seekg(i, std::ios_base::beg);
                            }
                        }
                        std::string lastline;
                        getline(logfile, lastline);
                        logfile.close();
                        if (lastline.back() == ';') {
                            // last line was incomplete -> timeout
                            // just write empty spaces in there
                            std::ofstream logger;
                            logger.open(expResFilename, std::ios_base::app);
                            logger << "0;0;-1\n";
                            logger.close();
                        }
                    }
                }

                /*void splitDebugLogger(dtstrat::dt::Observation state, uint action, float updatevalue, float qvalue, int node) {
                    std::ofstream logger;
                    logger.open(debugfilename, std::ios_base::app);
                    logger << state.printToCSV() << "," << action << "," << updatevalue << "," << qvalue << "," << node << "\n";
                    logger.close();
                }*/

                void writelowerbounds(){
                    std::ofstream logger;
                    logger.open(lowerBoundsFilename, std::ios_base::app);
                    logger << lower_b_init << "\n";
                    logger.close();
                }

                void writeupperbounds(){
                    std::ofstream logger;
                    logger.open(upperBoundsFilename, std::ios_base::app);
                    logger << upper_b_init << "\n";
                    logger.close();
                }
        };
        /*
        // generate a result-file and the headline
        void generate () {
            std::ofstream logger;
            logger.open ("results.txt");
            time_t rawtime;
            struct tm * timeinfo;
            char buffer[80];

            time (&rawtime);
            timeinfo = localtime(&rawtime);

            strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
            std::string str(buffer);
            logger << str << "\nstate;chosenAction;DTdecicsion;sizeOfDT" << "\n";
            logger.close();
        }

        // write a state in the result file
        void write (dt::Observation state, int chosenAction, bool dt, int sizeOfDT) {
            std::ofstream logger;
            logger.open (filename, std::ios_base::app);
            logger << "[" << state.print() << "];" << chosenAction << ";" << dt << ";" << sizeOfDT << "\n";
            logger.close();
            return;
        }

        // debug information into a the debug-file
        void debuglogger(std::string line) {
            std::ofstream logger;
            logger.open (debugfilename, std::ios_base::app);
            time_t rawtime;
            struct tm * timeinfo;
            char buffer[80];

            time (&rawtime);
            timeinfo = localtime(&rawtime);

            strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
            std::string str(buffer);
            logger << str << ": " << line << "\n";
            logger.close();
            return;
        }

        // print newstart into the results-file (if a new simulation is started)
        void newStart () {
            std::ofstream logger;
            logger.open (filename, std::ios_base::app);
            logger << "Newstart\n";
            logger.close();
            return;
        }*/
    }
}

#endif