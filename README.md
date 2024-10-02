This repository contains the scripts and setup for the evaluation for 1-2-3-Go!

We will run more, possibly even much more, than needed
## Setup
```bash
sudo apt install \
  build-essential \
  git \
  cmake \
  libboost-all-dev \
  libcln-dev libgmp-dev \
  libginac-dev \
  automake \
  libglpk-dev \
  libhwloc-dev \
  libz3-dev \
  libxerces-c-dev \
  libeigen3-dev

# initalize and update the subomdules
git submodule update --init

# build storm
mkdir -p storm/build
pushd storm/build
cmake ..
make -j8
popd

# Manually download and unzip modest: https://www.modestchecker.net/Downloads/
# Now, there should be a directory named Modest

```

## Scripts
The scripts can be divided in 3 parts:
### pre-exec
#### generate data (i.e., the models):
Some models are to be generated because the actual model file is different for different parameters. Same applies for property files.

The ones which do not change can be added to the repo. 
Others shouldn't be added, but generated.

Please follow the following steps to generate the models.
```bash
# create the directory named models
mkdir -p mdp-models

# run the script to generate models
./scripts/generate-models.py

# create jani models
# Different for different models.
 (cd scripts/ ; sh convert-models.sh)


# Pnueli-Zuck
for m in mdp-models/prism/pnueli-zuck*.prism;
do
  ./storm/build/bin/storm-conv --prism $m --tojani mdp-models/jani/$(basename $m).jani --prop mdp-models/prism/pnueli-zuck.props --globalvars
done

```

#### generate models with labeled actions (also called semantic actions):

```bash
# Create models with action labels
./scripts/give_action_names_to_models.py

# Copy the property files
cp mdp-models/prism/*.props mdp-models/action-models-prism/

# Convert these models to JANI 
(cd scripts/ ; sh convert-models.sh)


# create the task-ymls
# first copy
cp benchexec-input/task-yml/* benchexec-input/action-task-yml/

# then edit the path in the yml files
sed -i 's#/mdp-models/prism/#/mdp-models/action-models-prism/#' benchexec-input/action-task-yml/*.yml
sed -i 's#/mdp-models/jani/#/mdp-models/action-models-jani/#' benchexec-input/action-task-yml/*.yml
```


- generate commands to execute:
These scripts are supposed to generate `tool` commands to be executed

### execute run
We will use benchexec to execute the runs. Please go through the file [benchexec-setup.md](benchexec-setup.md).

### post-exec
- post-processing 1: Assemble the data from the output of runs
- generate plots and tables

## High level idea for eval
### Parameter scaling
The models in our benchmark set use 3 kinds of parameters: number of modules, variables, and both.
Moreover, these parameters effect the model size and the time taken to solve it differently.
We do not have time to execute extensive experiments, so we have run a smaller number to get the trend.
N is the number of modules, and K is the variable.

### Setting the baselines


```code
N_Model_Checkers = {Storm with 4 different engines: sparse, dd, dd-to-sparse, hybrid; and Modest in 3: LSS, TODO (finalize the other 2)}

for each model of parameterized MDP
	generate the concrete models
		for each of the N model checkers
			run the model and the property
Plot the data
```

#### Resource limits
We will allocate 2 cores, 8 GB, and 1 hr of CPU time for each run.

### Run 1-2-3-Go
After we have the baselines, we will evaluate different configurations of 1-2-3-Go!

#### One time setup
```bash
# initalize and update the subomdules (in case it isn't there)
git submodule update --init

# build dtstrat
mkdir -p dtstrat/build
pushd dtstrat/build
cmake ..
make -j8 storm-main
popd

# copy the painful script to the dtstrat directory
cp painful-script.sh dtstrat/
cp -r storm/build dtstrat/storm-stable-build
```


```code
#Note: We need to additionally set --overlay-dir . for this to work. 
# Example command
PYTHONPATH=scripts/tool-info-modules benchexec -N 4 --read-only-dir / --overlay-dir . benchexec-input/benchmark-definitions/test-dtstrat.xml --tool-dir dtstrat/
```
