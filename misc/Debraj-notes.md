There are already some scripts to semi-automate things (Steffi has some instructions in a scripts.md). I am writing down what Muqsit told me.

There are three types of examples:

### Examples where parameter type is number of modules module
- philosophers + eat
- rabin + live
- pnueli-zuck + live

There is a script for generating pnueli-zuck models for different number of modules. We may need to create scripts for philosophers and rabin as well.

The following code is for generating data for pnueli-zuck model. 
```bash
build/bin/storm --prism paper_models/pnueli-zuck.3.prism --prop paper_models/pnueli-zuck.props --buildfull --engine learning --no-simplify --dtLearning none --evaluationMethod none --reachable --childTasks "paper_models/pnueli-zuck.3.prism#Pmax=? [F (p1=10)]" > logs/pnueli-zuck-3_dtcontrol_0_3_none.LOG
```
It does not matter what we use as argument of `--prism` or `--prop` because we are not using them. We are using the `--childTasks` option to specify the properties and files. The argument of the `--childTasks` is of the form `prismfile1#property1%prismfile2#property2...`. This is creating a csv file with name `pnueli-zuck.3%live%3%0.csv`.

This command creates the DT and evaluates it:

```bash
build/bin/storm --prism paper_models/pnueli-zuck.6.prism --prop paper_models/pnueli-zuck.props "live" --engine learning --evaluationMethod full --dtLearning dtcontrol --exportscheduler pnueli-zuck-live-NUM-3.dot --datafile pnueli-zuck.3%live%3%0.csv > logs/pnueli-zuck%live%NUM%3%%full-1.LOG 
```
Here the argument of `--prism` is the prism file and the argument of `--prop` is the property. The `--datafile` option is used to specify the dataset that was created by the previous command. There is no command that takes pretrained DT and evaluates it (No need to implement this as this would rather introduce more bugs).
This creates a dot file with name given as the argument of `--exportscheduler`. 

It is also creating two `.drn` files describing the Markov chain and outputs of the dtcontrol (two `.json` file and an `html` file). with `--evaluationMethod full` option, it is not actually evaluate the DT. It is just creating the dot file and the `.drn` files.

To evaluate the markov chain
```bash
storm -drn markovchain.drn --prop "F(done)"
```

### Examples where parameter type is a variable
- firewire + deadline
- mer + p1
- pacman + crash
- zeroconf_dl + deadline_max
- zeroconf_dl + deadline_min

These one does not need a script to generate the new models as we can provide `storm` with the parameter vaues with the `--constants` argument. But we can only train for one parameter value at a time. 

This one is for generating the data:
```bash
build/bin/storm --prism paper_models/zeroconf.prism --prop paper_models/zeroconf.props --buildfull --engine learning --learnParameter 2,3 --paramName K --no-simplify --dtLearning dtcontrol --evaluationMethod exact --constants reset=false,N=20,K=10 --reachable
```
Here the `--prism` and `--prop` files matter. The `--paramName` option is used to specify the name of the parameter we use in learning. We fix the parameter values for which we want to create the training dataset with the `--learnParameter` option. The `--constants` option is used to specify the values of the parameters. For the paramenter given as the argument for `--paramName`, the implimentation would ignore what value we provide as the argument of `--constants`.

Problem if multiple props.

This one is for training the DT and evaluating it:
```bash 
build/bin/storm --prism paper_models/zeroconf.prism --prop paper_models/zeroconf.props "correct_min" --engine learning --evaluationMethod full --constants 'reset=false,N=20,K=10' --dtLearning dtcontrol --exportscheduler zeroconf23.dot --datafile zeroconf%correct_max%K%2_3.csv > logs/zeroconf%correct_max%K%2_3-full-1.LOG
```

### Examples where parameter type are both variables and number of modules

- consensus + c2
- consensus + disagree
- csma + all_before_max
- csma + all_before_min
- csma+some_before

We need to think which one of the above two approaches we should use for these examples.

```
¯\_(ツ)_/¯
```



# TODO

- 0 maps to different action in different things
- remove bscc condition for quicker result