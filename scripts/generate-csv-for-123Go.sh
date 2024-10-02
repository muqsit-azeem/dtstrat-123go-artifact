#!/usr/bin/sh

DTSTRAT_STORM=../dtstrat/build/bin/storm
NAIVE_MODELS_PATH=../mdp-models/prism/

# For NAIVE
# zeroconf
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/zeroconf_dl.prism --prop $NAIVE_MODELS_PATH/zeroconf_dl.props "deadline_max" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --constants reset=false,N=1000,K=2,deadline=50 --learnParameter 2 --paramName K 

$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/zeroconf_dl.prism --prop $NAIVE_MODELS_PATH/zeroconf_dl.props "deadline_min" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --constants reset=false,N=1000,K=2,deadline=50 --learnParameter 2 --paramName K 

 # Firewire
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/firewire_dl.prism --prop $NAIVE_MODELS_PATH/firewire_dl.props "deadline" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --constants delay=1,deadline=200 --learnParameter 3 --paramName delay

# mer
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/mer.prism --prop $NAIVE_MODELS_PATH/mer.props "p1" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --constants x=0.01,n=1 --learnParameter 1 --paramName n

#### with childtasks
 #pnueli
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/pnueli-zuck.3.prism --prop $NAIVE_MODELS_PATH/pnueli-zuck.props "live" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --childTasks "$NAIVE_MODELS_PATH/pnueli-zuck.3.prism#Pmax=? [F (p1=10)]"

# philosophers
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/philosophers.4.prism --prop $NAIVE_MODELS_PATH/philosophers.4.props "eat" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --childTasks "$NAIVE_MODELS_PATH/philosophers.4.prism#Pmax=? [ F (((p1>=8)&(p1<=9))|((p2>=8)&(p2<=9))|((p3>=8)&(p3<=9))|((p4>=8)&(p4<=9))) ]"

#csma
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/csma.N_2.K_2.prism --prop $NAIVE_MODELS_PATH/csma.props "some_before" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --childTasks $NAIVE_MODELS_PATH'/csma.N_2.K_2.prism#Pmin=? [ F min_backoff_after_success<K ]%'$NAIVE_MODELS_PATH'/csma.N_3.K_2.prism#Pmin=? [ F min_backoff_after_success<K ]'

$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/csma.N_2.K_2.prism --prop $NAIVE_MODELS_PATH/csma.props "all_before_max" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --childTasks $NAIVE_MODELS_PATH'/csma.N_2.K_2.prism#Pmax=? [ !"collision_max_backoff" U "all_delivered" ]%'$NAIVE_MODELS_PATH'/csma.N_3.K_2.prism#Pmax=? [ !"collision_max_backoff" U "all_delivered" ]'
 
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/csma.N_2.K_2.prism --prop $NAIVE_MODELS_PATH/csma.props "all_before_min" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --childTasks $NAIVE_MODELS_PATH'/csma.N_2.K_2.prism#Pmin=? [ !"collision_max_backoff" U "all_delivered" ]'


 # consensus disagree
$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/consensus.2.prism --prop $NAIVE_MODELS_PATH/consensus.props "disagree" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --constants K=2 \
 --childTasks $NAIVE_MODELS_PATH'/consensus.2.prism#Pmax=? [ F "finished" & !"agree" ]%'$NAIVE_MODELS_PATH'/consensus.3.prism#Pmax=? [ F "finished" & !"agree" ]'

$DTSTRAT_STORM --prism $NAIVE_MODELS_PATH/consensus.2.prism --prop $NAIVE_MODELS_PATH/consensus.props "c2" \
 --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable \
 --constants K=8 \
 --childTasks $NAIVE_MODELS_PATH'/consensus.2.prism#Pmin=? [ F "finished"&"all_coins_equal_1" ]%'$NAIVE_MODELS_PATH'/consensus.3.prism#Pmin=? [ F "finished"&"all_coins_equal_1" ]'