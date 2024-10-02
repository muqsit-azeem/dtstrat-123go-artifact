//
// Created by steffi on 03.01.22.
//

#ifndef STORM_DECISIONTREE_H
#define STORM_DECISIONTREE_H

#include <vector>
#include <string>
#include "storm/modelchecker/explorationDT/StateValuationDT.h"
#include "storm/parser/json2.hpp"


namespace storm {
namespace storage {

typedef unsigned int uint;

template<typename DataType>
class Node {
   private:
    // unique value
    uint num;

    // Left and right child of the node
    // default is nullptr
    Node* leftChild;
    Node* rightChild;
    Node* parent;

    // if the node is a branching node, the dimension gives the dimension to split on
    // the value gives the value to split on
    uint indexToSplit;
    storm::RationalNumber value;

    // determines whether the node is a leaf node (true) or a branching node (false)
    bool isLeaf = true;
    bool isLeft = true;

    // stored data
    std::vector<DataType> content;

   public:
    Node(uint indexSplit=0, storm::RationalNumber valIn=0.0, uint num=0, std::vector<DataType> vec= {}) : indexToSplit(indexSplit), value(valIn), num(num), content(vec) {};

    Node(Node* parent, uint num, bool left) : parent(parent), num(num), indexToSplit(0), value(RationalNumber(0)), isLeft(left) {};

    ~Node();

    bool evaluate(modelchecker::dtstrat::StateValuationDT inputValues) const;

    Node* getChild(bool choice);

    std::string getDotLabel() const;
    Node* getParent();

    uint size() const;
    bool isLeafNode() const;
    bool isLeftNode() const;

    uint number() const;

    void turnToDecisionNode();
    void turnToLeafNode();

    void addAction(uint action);

    void setSplitIndex(uint index);
    void setSplitValue(storm::RationalNumber value);
    void setLeftChild(Node* child);
    void setRightChild(Node* child);
    void setChildren(Node* leftChild, Node* rightChild);
    void setParent(Node* parentNode);
    void setLeft(bool newValue);

    std::vector<uint> getSortedActions() const;

};


template<typename DataType>
class DecisionTree {
   private:
    //... ATTRIBUTES ...
    // stores the root node
    Node<DataType>* root;
    // stores the full size of the decision tree (number of nodes)
    int numberOfNodes=0;
    // for speed-up, store all leaf-nodes
    std::vector<Node<DataType>*> leafNodes;

    //... FUNCTIONS ...
    bool isLeafNode(nlohmann2::json tree);
    std::tuple <int, float> getJSplitValues(nlohmann2::json split);
    void getSubtree(nlohmann2::json tree, Node<DataType>& currentNode);

    std::string getDotChildString(Node<DataType>* n) const;
    void split(Node<DataType>& node, uint dimension, storm::RationalNumber splitvalue);

   public:
    DecisionTree() = default;
    DecisionTree(std::string dotFile, bool module=false);
    DecisionTree(nlohmann2::json tree);
    ~DecisionTree();

    DecisionTree& operator=(DecisionTree<DataType> const& copy);

    void clear();
    void DFS_clear(Node<DataType>* currentNode, std::vector<Node<DataType>*>& stack);
    std::string getDotString() const;
    uint size() const;
    std::vector<uint> getDecisions(modelchecker::dtstrat::StateValuationDT const & inputValues) const;
};
}
}

#endif  // STORM_DECISIONTREE_H
