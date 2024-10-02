//
// Created by steffi on 03.01.22.
//
#include "storm/storage/DecisionTree.h"
#include "storm/exceptions/OutOfRangeException.h"
#include <iostream>
#include <functional>
#include <random>

namespace storm {
namespace storage {
template<typename DataType>
bool Node<DataType>::evaluate(modelchecker::dtstrat::StateValuationDT inputValues) const {
    return inputValues[indexToSplit]<=value;
}

template<typename DataType>
uint Node<DataType>::size() const {
    return content.size();
}

template<typename DataType>
bool Node<DataType>::isLeafNode() const {
    return isLeaf;
}

template<typename DataType>
uint Node<DataType>::number() const {
    return num;
}

template<typename DataType>
void Node<DataType>::turnToDecisionNode() {
    isLeaf = false;
}

template<typename DataType>
void Node<DataType>::turnToLeafNode(){
    isLeaf = true;
}

template<typename DataType>
void Node<DataType>::setParent(Node<DataType>* parentNode){
    parent = parentNode;
}

template<typename DataType>
void Node<DataType>::setLeft(bool newValue) {
    isLeft = newValue;
}

template<typename DataType>
Node<DataType>* Node<DataType>::getParent() {
    return parent;
}

template<typename DataType>
bool Node<DataType>::isLeftNode() const {
    return isLeft;
}

template<typename DataType>
void Node<DataType>::addAction(uint action) {
    content.push_back(DataType(action));
}

template<typename DataType>
void Node<DataType>::setSplitIndex(uint index){
    indexToSplit = index;
}

template<typename DataType>
void Node<DataType>::setSplitValue(storm::RationalNumber nvalue){
    value = nvalue;
}

template<typename DataType>
void Node<DataType>::setLeftChild(Node* child){
    leftChild = child;
}

template<typename DataType>
void Node<DataType>::setRightChild(Node* child){
    rightChild = child;
}

template<typename DataType>
void Node<DataType>::setChildren(Node* leftChild, Node* rightChild){
    setLeftChild(leftChild);
    setRightChild(rightChild);
}

template<typename DataType>
Node<DataType>* Node<DataType>::getChild(bool choice) {
    if (choice) return leftChild;
    else return rightChild;
}

template<typename DataType>
std::string Node<DataType>::getDotLabel() const {
    std::ostringstream label_stream;
    if (!isLeaf) {
        label_stream << this->num << " [label=\"x_" << indexToSplit << " <= " << value.get_d() << "\"];";
    }
    else {
        std::vector<uint> sorted_actions;
        sorted_actions = getSortedActions();
        label_stream << this->num << " [label=\"";
        uint counter = 0;
        for (uint i=0; i<sorted_actions.size(); i++) {
            auto val = sorted_actions.at(i);
            label_stream << "actions[" << i << "] = " << val << ",";
            counter++;
        }
        label_stream << "\"];";
    }
    return label_stream.str();
}

template<typename DataType>
std::vector<uint> Node<DataType>::getSortedActions() const {
    return content;
}
template<typename DataType>
Node<DataType>::~Node() {
    delete leftChild;
    delete rightChild;
}

template <>
std::vector<uint> Node<float>::getSortedActions() const {
    std::vector<std::pair<float,uint>> qvaluesAndIndex;
    for(uint i=0; i<content.size(); i++) {
        qvaluesAndIndex.push_back(std::make_pair(content[i],i));
    }
    std::sort(qvaluesAndIndex.begin(), qvaluesAndIndex.end(), std::greater <>());
    std::vector<uint> sortedActionIndices;
    std::vector<uint> actionIndicesForSameQValue;
    float lastQvalue = 10.0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine e(seed);
    for (uint i = 0; i<content.size(); i++) {
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
}

template <typename DataType>
bool DecisionTree<DataType>::isLeafNode(nlohmann2::json tree) {
    return (!tree.contains("split") || tree["split"] == nullptr);
}

template <typename DataType>
std::tuple <int, float> DecisionTree<DataType>::getJSplitValues(nlohmann2::json split) {
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
}

template <typename DataType>
void DecisionTree<DataType>::getSubtree(nlohmann2::json tree, Node<DataType>& currentNode) {
    if (!isLeafNode(tree)) {
        auto splitvals = getJSplitValues(tree["split"]);
        int dimension = std::get<0>(splitvals);
        float value = std::get<1>(splitvals);
        this->split(currentNode, dimension, value);
        auto childs =  tree["children"];
        std::vector<nlohmann2::json> children;
        for (auto& child : childs) {
            children.push_back(child);
        }
        getSubtree(children[0], *(currentNode.getChild(true)));
        getSubtree(children[1], *(currentNode.getChild(false)));
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
            if(labels.size()>=2) {
                /* if label size (number of actions) is greater or equal to 2, store all the actions */
                std::random_device rd;
                std::mt19937 g(rd());
                for (uint i = 0; i < labels.size(); ++i) {
                    std::string action_string = labels[i].dump();
                    action_string = action_string.substr(1, action_string.size() - 1);
                    long testAction = std::stol(action_string);
                    STORM_LOG_THROW(testAction < 4294967295,exceptions::OutOfRangeException, "Action is 4294967295 (=-1) for uint maybe.");
                    action = std::stoi(action_string);
                    currentNode.addAction(action);
                }
            }
            else{
                std::string action_string = labels[0].dump();
                action_string = action_string.substr(1, action_string.size() - 1);
                action = std::stoi(action_string);
                currentNode.addAction(action);
            }
        }
    }
}

template<typename DataType>
void DecisionTree<DataType>::split(Node<DataType>& node, uint dimension, storm::RationalNumber splitvalue) {
    // Create new child nodes
    Node<DataType>* leftNode = new Node<DataType>(&node, numberOfNodes, true);
    numberOfNodes++;
    Node<DataType>* rightNode = new Node<DataType>(&node, numberOfNodes, false);
    numberOfNodes++;
    this->leafNodes.erase(std::remove(this->leafNodes.begin(), this->leafNodes.end(), &node), this->leafNodes.end());
    this->leafNodes.push_back(leftNode);
    this->leafNodes.push_back(rightNode);

    // the parent node is no longer a leaf node
    node.turnToDecisionNode();

    // set dimension and value to split on
    node.setSplitIndex(dimension);
    node.setSplitValue(splitvalue);

    // add the new nodes as the children to the parent node
    node.setChildren(leftNode, rightNode);
}

template <typename DataType>
std::string DecisionTree<DataType>::getDotChildString(Node<DataType>* n) const {
    std::string result;
    if (n->isLeafNode()) {
        result = n->getDotLabel() + "\n";
    }
    else {
        result = n->getDotLabel() + "\n";
        result += std::to_string(n->number()) + " -> " + std::to_string(n->getChild(true)->number()) +
                  " [label=\"true\"];\n";
        result += getDotChildString(n->getChild(true));
        result += std::to_string(n->number()) + " -> " + std::to_string(n->getChild(false)->number()) +
                  " [label=\"false\"];\n";
        result += getDotChildString(n->getChild(false));
    }
    return result;
}

template <typename DataType>
DecisionTree<DataType>::DecisionTree(std::string dotFile, bool module) {
    std::vector<Node<DataType>*> nodes;
    std::vector<std::tuple<int,int,bool>> edges;
    std::ifstream infile(dotFile);
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
                Node<DataType>* newNode = new Node<DataType>(dec_dim, storm::RationalNumber(dec_val));
                if(nodes.size()<=nodeIndex) {
                    nodes.resize(nodeIndex);
                    nodes.push_back(newNode);
                } else {
                    nodes[nodeIndex] = newNode;
                }
            } else if (label[7] == 'a') {
                // leaf node
                Node<DataType>* newNode = new Node<DataType>();
                while (iss >> symbol && (symbol=='=' || symbol==',')) {
                    std::string rest;
                    iss >> rest;
                    std::size_t found = rest.find(',');
                    uint action = std::stoi(rest.substr(0, found));
                    newNode->addAction(action);
                }
                if(nodes.size()<=nodeIndex) {
                    nodes.resize(nodeIndex);
                    nodes.push_back(newNode);
                } else {
                    nodes[nodeIndex] = newNode;
                }
            } else if (label[7] == 'm' && module) {
                // leaf node
                Node<DataType>* newNode = new Node<DataType>();
                newNode->turnToLeafNode();
                // read modules
                while (iss >> symbol && (symbol=='=' || symbol==',')) {
                    std::string rest;
                    iss >> rest;
                    while(rest.size()>1){
                        std::size_t found = rest.find(',');
                        if (found==std::string::npos) {break;}
                        uint action = std::stoi(rest.substr(0, found));
                        rest = rest.substr(found+1);
                        // TODO integrate again, if necessary
                        //newNode->modules.push_back(action);
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
                        newNode->addAction(action);
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
            nodes[nodeIndex]->setLeftChild(nodes[node2]);
            nodes[nodeIndex]->turnToDecisionNode();
            nodes[node2]->setParent(nodes[nodeIndex]);
            nodes[node2]->setLeft(true);
        }
        else {
            nodes[nodeIndex]->setRightChild(nodes[node2]);
            nodes[nodeIndex]->turnToDecisionNode();
            nodes[node2]->setParent(nodes[nodeIndex]);
            nodes[node2]->setLeft(false);
        }
    }
    this->numberOfNodes = nodes.size();
    this->root = nodes[0];
}

template <typename DataType>
DecisionTree<DataType>::DecisionTree(nlohmann2::json tree) {
    this->root = new Node<DataType>(0,storm::RationalNumber(0));
    this->leafNodes.push_back(root);
    numberOfNodes++;
    getSubtree(tree, *root);
}

template <typename DataType>
void DecisionTree<DataType>::clear() {
    if (this->numberOfNodes<2) {
        return;
    }
    std::vector<Node<DataType>*> stack;
    stack.push_back(root);
    DFS_clear(this->root, stack);
    this->numberOfNodes = 1;
    this->leafNodes.clear();
    this->root = new Node<DataType>(0,storm::RationalNumber(0));
}

template <typename DataType>
void DecisionTree<DataType>::DFS_clear(Node<DataType>* currentNode, std::vector<Node<DataType>*>& stack) {
    if (!(currentNode->getChild(true)) && !(currentNode->getChild(false))) {
        // node has no children -> leaf node
        // remove it
        if (!currentNode->getParent()) {
            delete currentNode;
            stack.pop_back();
            return;
        }
        if (currentNode->isLeftNode()) {
            Node<DataType>* cug = currentNode->getParent();
            Node<DataType>* parentLeftchild = cug->getChild(true);
            parentLeftchild = nullptr;
        }
        else {
            Node<DataType>* cug = currentNode->getParent();
            Node<DataType>* parentRightChild = cug->getChild(false);
            parentRightChild = nullptr;
        }
        delete currentNode;
        stack.pop_back();
    }
    else {
        stack.push_back(currentNode->getChild(true));
        stack.push_back(currentNode->getChild(false));
    }
    if (stack.size()>0) {
        auto nextNode = stack.back();
        DFS_clear(nextNode, stack);
    }
}

template <typename DataType>
std::string DecisionTree<DataType>::getDotString() const {
    std::string result = "digraph {\n";
    if (this->numberOfNodes>0) {
        Node<DataType>* n = root;
        result += getDotChildString(n);
    }
    result += "}\n";
    return result;
}

template <typename DataType>
uint DecisionTree<DataType>::size() const {
    return this->numberOfNodes;
}

template <typename DataType>
std::vector<uint> DecisionTree<DataType>::getDecisions(modelchecker::dtstrat::StateValuationDT const & inputValues) const {
    if (this->numberOfNodes == 0) return {};
    bool reached_end = false;
    Node<DataType>* currentNode = root;
    // runs through the decision tree
    while (!reached_end) {
        // if it has reached a leaf node, it is done
        if (currentNode->isLeafNode())
            reached_end = true;
        // if it is in a branching node, evaluate the node and move to the respective child node
        else {
            currentNode = currentNode->getChild(currentNode->evaluate(inputValues));
        }
    }
    return currentNode->getSortedActions();
}

template<typename DataType>
DecisionTree<DataType>::~DecisionTree() {
    delete root;
    // leaf nodes will be deleted automatically
    // root deletes its children and so on until all leaf nodes are delted
}

template<typename DataType>
DecisionTree<DataType>& DecisionTree<DataType>::operator=(DecisionTree<DataType> const& copy) {
    if (this == &copy) return *this;
    else {
        delete root;
        root = copy.root;
        numberOfNodes = copy.numberOfNodes;
        leafNodes = copy.leafNodes;
        return *this;
    }
}

template class Node<float>;
template class Node<uint>;

template class DecisionTree<float>;
template class DecisionTree<uint>;
}
}
