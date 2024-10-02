#!/bin/bash

# Should create the drn file
dtstrat/build/bin/storm $@

# Then we use storm's stable build to evaluate it if full evaluation was chosen. 
if [[ "$@" == *" --evaluationMethod statistical"* ]]
then
    echo "exiting because used statistical evaluation"    
else
    
    dtstrat/storm-stable-build/bin/storm -drn markovchain.drn --prop "P=? [ F\"done\"]"
fi