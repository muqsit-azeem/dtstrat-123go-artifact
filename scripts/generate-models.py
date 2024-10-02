#!/usr/bin/env python3
import sys
import os
from itertools import repeat,chain
import yaml

sys.dont_write_bytecode = True  # prevents writing .pyc files

MODELS_DIR="mdp-models/prism"
JANI_MODELS_DIR="mdp-models/jani"

TASK_DIR = "benchexec-input/task-yml"


def generate_models_ij():
    # Model from: https://www.prismmodelchecker.org/casestudies/self-stabilisation.php
    for NUM in range(3,100):
        assert(NUM>=3)

        lines = [
        "// Israeli-Jalfon self stabilising algorithm",
        "// dxp/gxn 10/06/02",
        "mdp",
        "// variables to represent whether a process has a token or not",
        "// note they are global because they can be updated by other processes"]
        for n in range(NUM):
            lines.append(f"global q{n+1} : [0..1] init 1;")
        lines.append(f"global num_tokens_var : [0..{n+1}];")
        lines.append("")
        lines.append("// module of process 1")
        lines.append("module process1")
        lines.append(f"  [] (q1=1) -> 0.5 : (q1'=0) & (q{NUM}'=1) + 0.5 : (q1'=0) & (q2'=1);")
        lines.append("endmodule")
        lines.append("")
        lines.append("// add further processes through renaming")
        for n in range(NUM-1):
            line = f"module process{n+2} = process1 ["
            for k in range(NUM):
                line += f"q{k+1}=q{(n+1+k)%(NUM)+1},"
            line=line[:-1]
            line+="] endmodule"
            lines.append(line)

        with open(f"{MODELS_DIR}/ij.{NUM}.prism", "w+") as f:
            line_sep = repeat("\n")
            line_iter = chain.from_iterable(zip(lines, line_sep))
            f.writelines(line_iter)
            
        with open(f"{MODELS_DIR}/ij.{NUM}.props", "w+") as f:
            line = "\"stable\": Pmax=? [ F  ("
            for n in range(NUM):
                line += "("
                for k in range(NUM):
                    if n==k:
                        line += f"q{k+1}=1 & "
                    else:
                        line += f"q{k+1}=0 & "
                line=line[:-2]
                line += ") | "
            line=line[:-2] 
            line += ") ];"
            f.write(line)


def generate_models_pnueli():
    # The original model starts from 3
    # https://www.prismmodelchecker.org/casestudies/mutual.php
    param_set = {3, 5, 8, 13, 21, 24, 27, 30, 34, 55, 89, 144}

 
    for NUM in param_set:

        lines = []
        lines.append("// mutual exclusion [PZ82]")
        lines.append("// dxp/gxn 19/12/99")
        lines.append("")
        lines.append("mdp")
        lines.append("")
        lines.append("// atomic formula")
        lines.append("// none in low, high, tie")

        line = "formula none_lht = "
        for i in range(1,NUM):
            line+=f"(p{i}<4 | p{i}>13) & "
        line = line[:-3]
        line+=";"

        lines.append(line)

        lines.append("// some in admit")
        line = "formula some_a = "
        for i in range(1,NUM):
            line+=f"(p{i}>=14 & p{i}<=15) | "
        line = line[:-3]
        line+=";"
        lines.append(line)

        lines.append("// some in high, admit")
        line = "formula some_ha	 = "
        for i in range(1,NUM):
            line += f"(p{i}>=4 & p{i}<=5) | (p{i}>=10 & p{i}<=15) | "
        line = line[:-3]
        line+=";"
        lines.append(line)

        lines.append("// none in high, tie, admit")
        line = "formula none_hta = "
        for i in range(1,NUM):
            line+= f"(p{i}>=0 & p{i}<=3) | (p{i}>=7 & p{i}<=8) | "
        line = line[:-3]
        line+=";"
        lines.append(line)

        lines.append("// none in enter")
        line = "formula none_e	 = 	"
        for i in range(1,NUM):
            line+=f"(p{i}<2 | p{i}>3) & "
        line = line[:-3]
        line+=";"
        lines.append(line)

        lines.append("// process 0")
        lines.append("module process0")
        lines.append("")
        lines.append("	p0: [0..15] init 1;")
        lines.append("")	
        lines.append("	[] p0=0 -> (p0'=0);")
        lines.append("	[] p0=0 -> (p0'=1);")
        lines.append("	[] p0=1 -> (p0'=2);")
        lines.append("	[] p0=2 &  (none_lht | some_a) -> (p0'=3);")
        lines.append("	[] p0=2 & !(none_lht | some_a) -> (p0'=2);")
        lines.append("	[] p0=3 -> (p0'=4);")
        lines.append("	[] p0=3 -> (p0'=7);")
        lines.append("	[] p0=4 &  some_ha -> (p0'=5);")
        lines.append("	[] p0=4 & !some_ha -> (p0'=10);")
        lines.append("	[] p0=5 -> (p0'=6);")
        lines.append("	[] p0=6 &  some_ha -> (p0'=6);")
        lines.append("	[] p0=6 & !some_ha -> (p0'=9);")
        lines.append("	[] p0=7 &  none_hta -> (p0'=8);")
        lines.append("	[] p0=7 & !none_hta -> (p0'=7);")
        lines.append("	[] p0=8  -> (p0'=9);")
        lines.append("	[] p0=9  -> 0.5 : (p0'=4) + 0.5 : (p0'=7);")
        lines.append("	[] p0=10 -> (p0'=11);")
        lines.append("	[] p0=11 &  none_lht -> (p0'=13);")
        lines.append("	[] p0=11 & !none_lht -> (p0'=12);")
        lines.append("	[] p0=12 -> (p0'=0);")
        lines.append("	[] p0=13 -> (p0'=14);")
        lines.append("	[] p0=14 &  none_e -> (p0'=15);")
        lines.append("	[] p0=14 & !none_e -> (p0'=14);")
        lines.append("	[] p0=15 -> (p0'=0);")
        lines.append("")	
        lines.append("endmodule")
        lines.append("")
        for i in range(1,NUM):
            line = f"module process{i} = process0 [p0=p{i}, p{i}=p0] endmodule"
            lines.append(line)
            
        with open(f"{MODELS_DIR}/pnueli-zuck.{NUM}.prism","w+") as f:
            for line in lines:
                f.write(line+"\n")
        
        #Generate tasks
        lines = []
        lines.append("format_version: '2.0'")
        lines.append("")
        lines.append("input_files: '../../{MODELS_DIR}/pnueli-zuck.{NUM}.prism'".format(MODELS_DIR=MODELS_DIR,NUM=NUM))
        lines.append("")
        lines.append("options:")
        lines.append("  property_tag: live")
        lines.append("  property_file: pnueli-zuck.props")
        with open(f"{TASK_DIR}/pnueli-zuck.{NUM}.prism.yml","w+") as f:
            for line in lines:
                f.write(line+"\n")
        
        # Task YAML for Modest
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f"../../{JANI_MODELS_DIR}/pnueli-zuck.{NUM}.jani"
        task_yml["options"] = dict(property_tag="live")
        with open(f"{TASK_DIR}/pnueli-zuck.{NUM}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)


def generate_tasks_zeroconf():
    CONSTANTS = 'reset=false,deadline=10,K=1,N='    
    for N in (1000, 2000, 4000, 8000, 16000, 32000):
        task_yml = dict(format_version='2.0')
        options_dict = dict(property_file='zeroconf_dl.props',constants=CONSTANTS+str(N))
        task_yml["options"] = options_dict
        for s in {'deadline_max', 'deadline_min'}:
            task_yml["input_files"] = f'../../{MODELS_DIR}/zeroconf_dl.prism'
            task_yml["options"]["property_tag"] = s
            with open(f"{TASK_DIR}/zeroconf_dl.{s}.N_{N}.prism.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
            
            task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/zeroconf_dl.jani'
            with open(f"{TASK_DIR}/zeroconf_dl.{s}.N_{N}.jani.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
    
    CONSTANTS = 'reset=false,deadline=10,N=1000,K='    
    for K in (1, 2, 3, 4, 8, 16, 32):
        task_yml = dict(format_version='2.0')
        options_dict = dict(property_file='zeroconf_dl.props',constants=CONSTANTS+str(K))
        task_yml["options"] = options_dict
        for s in {'deadline_max', 'deadline_min'}:
            task_yml["input_files"] = f'../../{MODELS_DIR}/zeroconf_dl.prism'
            task_yml["options"]["property_tag"] = s
            with open(f"{TASK_DIR}/zeroconf_dl.{s}.K_{K}.prism.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
            
            task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/zeroconf_dl.jani'
            with open(f"{TASK_DIR}/zeroconf_dl.{s}.K_{K}.jani.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
    
    CONSTANTS = 'reset=false,N=1000,K=1,deadline='    
    for dl in (100, 200, 300, 600, 1000, 1600, 2000, 2400, 2800, 3200):
        task_yml = dict(format_version='2.0')
        options_dict = dict(property_file='zeroconf_dl.props',constants=CONSTANTS+str(dl))
        task_yml["options"] = options_dict
        for s in {'deadline_max', 'deadline_min'}:
            task_yml["input_files"] = f"../../{MODELS_DIR}/zeroconf_dl.prism"
            task_yml["options"]["property_tag"] = s
            with open(f"{TASK_DIR}/zeroconf_dl.{s}.deadline_{dl}.prism.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
            
            task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/zeroconf_dl.jani'
            with open(f"{TASK_DIR}/zeroconf_dl.{s}.deadline_{dl}.jani.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)


def generate_tasks_mer():
    CONSTANTS = 'x=0.01,n='    
    for N in {1, 10, 100, 1000, 3000, 8000, 13000, 21000, 34000, 55000, 89000, 144000}:
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f'../../{MODELS_DIR}/mer.prism'
        options_dict = dict(property_file='mer.props',property_tag="p1",constants=CONSTANTS+str(N))
        task_yml["options"] = options_dict
            
        with open(f"{TASK_DIR}/mer.n_{N}.prism.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
        
        task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/mer.jani'
        with open(f"{TASK_DIR}/mer.n_{N}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
    
    CONSTANTS = 'n=10,x='    
    for x in {0.01, 0.03, 0.05, 0.08, 0.13, 0.21, 0.34, 0.55, 0.89, 1.44}:
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f'../../{MODELS_DIR}/mer.prism'
        options_dict = dict(property_file='mer.props',property_tag="p1",constants=CONSTANTS+str(x))
        task_yml["options"] = options_dict
            
        with open(f"{TASK_DIR}/mer.x_{x}.prism.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
        
        task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/mer.jani'
        with open(f"{TASK_DIR}/mer.x_{x}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)


def generate_tasks_pacman():
    CONSTANTS = 'MAXSTEPS='    
    for MAXSTEPS in {5, 25, 50, 75, 100, 150, 200, 300, 400, 500}:
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f'../../{MODELS_DIR}/pacman.prism'
        options_dict = dict(property_file='pacman.props',property_tag="crash",constants=CONSTANTS+str(MAXSTEPS))
        task_yml["options"] = options_dict
            
        with open(f"{TASK_DIR}/pacman.MAXSTEPS_{MAXSTEPS}.prism.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
        
        task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/pacman.jani'
        with open(f"{TASK_DIR}/pacman.MAXSTEPS_{MAXSTEPS}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)


def generate_tasks_firewire():
    CONSTANTS = 'delay=3,deadline='    
    for deadline in {200, 300, 500, 800, 1300, 2100, 3400, 5500, 8900, 14400}:
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f'../../{MODELS_DIR}/firewire_dl.prism'
        options_dict = dict(property_file='firewire_dl.props',property_tag="deadline",constants=CONSTANTS+str(deadline))
        task_yml["options"] = options_dict
            
        with open(f"{TASK_DIR}/firewire_dl.deadline_{deadline}.prism.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
        
        task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/firewire_dl.jani'
        with open(f"{TASK_DIR}/firewire_dl.deadline_{deadline}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
    
    CONSTANTS = 'deadline=200,delay='    
    for delay in {3, 5, 8, 13, 21, 34, 55, 89, 144}:
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f'../../{MODELS_DIR}/firewire_dl.prism'
        options_dict = dict(property_file='firewire_dl.props',property_tag="deadline",constants=CONSTANTS+str(delay))
        task_yml["options"] = options_dict
            
        with open(f"{TASK_DIR}/firewire_dl.delay_{delay}.prism.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
        
        task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/firewire_dl.jani'
        with open(f"{TASK_DIR}/firewire_dl.delay_{delay}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)

def generate_models_philosophers():
    param_set = {3, 4, 5, 8, 13, 21, 24, 27, 30, 34, 55, 89, 144}

    for NUM_PHILOSOPHERS in param_set:
        # CREATE THE PRISM FILE
        lines = ['// randomized dining philosophers [LR81]\n', '// dxp/gxn 12/12/99\n', '// atomic formulae \n',
                '// left fork free and right fork free resp.\n', 'mdp\n', '\n', 'formula lfree = (p2>=0&p2<=4)|p2=6|p2=10;\n',
                'formula rfree = (p3>=0&p3<=3)|p3=5|p3=7|p3=11;\n', '\n', 'module phil1\n', '\n', '\tp1: [0..11];\n', '\n',
                "\t[] p1=0 -> (p1'=0); // stay thinking\n", "\t[] p1=0 -> (p1'=1); // trying\n",
                "\t[] p1=1 -> 0.5 : (p1'=2) + 0.5 : (p1'=3); // draw randomly\n",
                "\t[] p1=2 &  lfree -> (p1'=4); // pick up left\n", "\t[] p1=2 &  !lfree -> (p1'=2); // left not free\n",
                "\t[] p1=3 &  rfree -> (p1'=5); // pick up right\n", "\t[] p1=3 &  !rfree -> (p1'=3); // right not free\n",
                "\t[] p1=4 &  rfree -> (p1'=8); // pick up right (got left)\n",
                "\t[] p1=4 & !rfree -> (p1'=6); // right not free (got left)\n",
                "\t[] p1=5 &  lfree -> (p1'=8); // pick up left (got right)\n",
                "\t[] p1=5 & !lfree -> (p1'=7); // left not free (got right)\n", "\t[] p1=6  -> (p1'=1); // put down left\n",
                "\t[] p1=7  -> (p1'=1); // put down right\n", "\t[] p1=8  -> (p1'=9); // move to eating (got forks)\n",
                "\t[] p1=9  -> (p1'=10); // finished eating and put down left \n",
                "\t[] p1=9  -> (p1'=11); // finished eating and put down right\n",
                "\t[] p1=10 -> (p1'=0); // put down right and return to think\n",
                "\t[] p1=11 -> (p1'=0); // put down left and return to think\n", '\n', 'endmodule\n', '\n',
                '// construct further modules through renaming\n']

        for n in range(2,NUM_PHILOSOPHERS+1):
            first = n
            second = (n+1)%NUM_PHILOSOPHERS if n+1>NUM_PHILOSOPHERS else n+1
            third = n-1
            line = f"module phil{first} = phil1 [ p1=p{first}, p2=p{second}, p3=p{third} ] endmodule\n"
            lines.append(line)

        with open(f"{MODELS_DIR}/philosophers.{NUM_PHILOSOPHERS}.prism", "w+") as f:
            f.writelines(lines)

        # CREATE THE PROPERTY
        line = '"eat": Pmax=? [ F ('
        for n in range(1,NUM_PHILOSOPHERS+1):
            line+= f"((p{n}>=8)&(p{n}<=9))|"
        line = line[:-1]
        line += ') ];'
        with open(f"{MODELS_DIR}/philosophers.{NUM_PHILOSOPHERS}.props", "w+") as f:
            f.write("// Finally one of the philosophers manages to eat\n")
            f.write(line)
        
        # Task YAML for Storm
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f"../../{MODELS_DIR}/philosophers.{NUM_PHILOSOPHERS}.prism"
        options_dict = dict(property_file=f"philosophers.{NUM_PHILOSOPHERS}.props",property_tag="eat")
        task_yml["options"] = options_dict
                
        with open(f"{TASK_DIR}/philosophers.{NUM_PHILOSOPHERS}.prism.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)
        
        # Task YAML for Modest
        task_yml["input_files"] = f"../../{JANI_MODELS_DIR}/philosophers.{NUM_PHILOSOPHERS}.jani"
        task_yml["options"] = dict(property_tag="eat")
        with open(f"{TASK_DIR}/philosophers.{NUM_PHILOSOPHERS}.jani.yml","w+") as f:
            yaml.dump(task_yml, f, default_flow_style=False)


def generate_models_consensus():
    param_set = {2, 3, 5, 6, 7, 8, 13, 21, 34, 55, 89, 144}

 
    for N in param_set:
        s = f"""// COIN FLIPPING PROTOCOL FOR POLYNOMIAL RANDOMIZED CONSENSUS [AH90] 
    // gxn/dxp 20/11/00

    mdp

    // constants
    const int N={N};
    const int K;
    const int range = 2*(K+1)*N;
    const int counter_init = (K+1)*N;
    const int left = N;
    const int right = 2*(K+1)*N - N;

    // shared coin
    global counter : [0..range] init counter_init;

    module process1
        
        // program counter
        pc1 : [0..3];
        // 0 - flip
        // 1 - write 
        // 2 - check
        // 3 - finished
        
        // local coin
        coin1 : [0..1];	

        // flip coin
        [] (pc1=0)  -> 0.5 : (coin1'=0) & (pc1'=1) + 0.5 : (coin1'=1) & (pc1'=1);
        // write tails -1  (reset coin to add regularity)
        [] (pc1=1) & (coin1=0) & (counter>0) -> (counter'=counter-1) & (pc1'=2) & (coin1'=0);
        // write heads +1 (reset coin to add regularity)
        [] (pc1=1) & (coin1=1) & (counter<range) -> (counter'=counter+1) & (pc1'=2) & (coin1'=0);
        // check
        // decide tails
        [] (pc1=2) & (counter<=left) -> (pc1'=3) & (coin1'=0);
        // decide heads
        [] (pc1=2) & (counter>=right) -> (pc1'=3) & (coin1'=1);
        // flip again
        [] (pc1=2) & (counter>left) & (counter<right) -> (pc1'=0);
        // loop (all loop together when done)
        [done] (pc1=3) -> (pc1'=3);

    endmodule

    // construct remaining processes through renaming"""

        for i in range(1,N):
            s += f"\nmodule process{i+1} = process1[pc1=pc{i+1},coin1=coin{i+1}] endmodule"

        s += """\n
    // labels
    """
        s += f"label \"finished\" = "
        for i in range(N):
            s += f"pc{i+1}=3"
            if i != N-1:
                s += " & "
        s += " ;\n"

        s += f"label \"all_coins_equal_0\" = "
        for i in range(N):
            s += f"coin{i+1}=0"
            if i != N-1:
                s += " & "
        s += " ;\n"

        s += f"label \"all_coins_equal_1\" = "
        for i in range(N):
            s += f"coin{i+1}=1"
            if i != N-1:
                s += " & "
        s += " ;\n"

        s += f"label \"agree\" = "
        for i in range(1,N):
            s += f"coin{i}=coin{i+1}"
            if i != N-1:
                s += " & "
        s += " ;\n"

        s+="""
    // rewards
    rewards "steps"
        true : 1;
    endrewards
    """
        with open(f"{MODELS_DIR}/consensus.{N}.prism", "w+") as f:
            f.writelines(s)
    
    # tasks with N=2 and changing K
    CONSTANTS = 'K='
    N=2
    for K in {2, 3, 5, 8, 13, 21, 34, 55, 89, 144}:        
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f"../../{MODELS_DIR}/consensus.{N}.prism"
        options_dict = dict(property_file='consensus.props',constants=CONSTANTS+str(K))
        task_yml["options"] = options_dict
        for s in {'c2', 'disagree'}:
            task_yml["input_files"] = f"../../{MODELS_DIR}/consensus.{N}.prism"
            task_yml["options"]["property_tag"] = s
            with open(f"{TASK_DIR}/consensus.{s}.K_{K}.prism.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
            
            task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/consensus.{N}.jani'
            with open(f"{TASK_DIR}/consensus.{s}.K_{K}.jani.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
    
    # tasks with K=2 and changing N
    CONSTANTS = 'K=2'
    for N in {2, 3, 5, 6, 7, 8, 13, 21, 34, 55, 89, 144}:        
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f"../../{MODELS_DIR}/consensus.{N}.prism"
        options_dict = dict(property_file='consensus.props',constants=CONSTANTS)
        task_yml["options"] = options_dict
        for s in {'c2', 'disagree'}:
            task_yml["input_files"] = f"../../{MODELS_DIR}/consensus.{N}.prism"
            task_yml["options"]["property_tag"] = s
            with open(f"{TASK_DIR}/consensus.{s}.N_{N}.prism.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
            
            task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/consensus.{N}.jani'
            with open(f"{TASK_DIR}/consensus.{s}.N_{N}.jani.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)


def gen_csma(N, K): 
    s = f"""// CSMA/CD protocol - probabilistic version of kronos model (3 stations)
// gxn/dxp 04/12/01

mdp

// note made changes since cannot have strict inequalities
// in digital clocks approach and suppose a station only sends one message

// simplified parameters scaled
const int sigma=1; // time for messages to propagate along the bus
const int lambda=30; // time to send a message

// actual parameters
const int N = {N}; // number of processes
const int K = {K}; // exponential backoff limit
const int slot = 2*sigma; // length of slot
const int M = floor(pow(2, K))-1 ; // max number of slots to wait
//const int lambda=782;
//const int sigma=26;

//----------------------------------------------------------------------------------------------------------------------------
// the bus
module bus
    
    b : [0..2];
    // b=0 - idle
    // b=1 - active
    // b=2 - collision
    
    // clocks of bus
    y1 : [0..sigma+1]; // time since first send (used find time until channel sensed busy)
    y2 : [0..sigma+1]; // time since second send (used to find time until collision detected)
    
    // a sender sends (ok - no other message being sent)"""

    for i in range(N):
        s += f"\n\t[send{i+1}] (b=0) -> (b'=1);"

    s += "\n\n\t// a sender sends (bus busy - collision)"
    for i in range(N):
        s += f"\n\t[send{i+1}] (b=1|b=2) & (y1<sigma) -> (b'=2);"

    s += "\n\n\t// finish sending"
    for i in range(N):
        s += f"\n\t[end{i+1}] (b=1) -> (b'=0) & (y1'=0);"

    s += "\n\n\t// bus busy"
    for i in range(N):
        s += f"\n\t[busy{i+1}] (b=1|b=2) & (y1>=sigma) -> (b'=b);"

    s += """
    // collision detected
    [cd] (b=2) & (y2<=sigma) -> (b'=0) & (y1'=0) & (y2'=0);
    
    // time passage
    [time] (b=0) -> (y1'=0); // value of y1/y2 does not matter in state 0
    [time] (b=1) -> (y1'=min(y1+1,sigma+1)); // no invariant in state 1
    [time] (b=2) & (y2<sigma) -> (y1'=min(y1+1,sigma+1)) & (y2'=min(y2+1,sigma+1)); // invariant in state 2 (time until collision detected)
    
endmodule

//----------------------------------------------------------------------------------------------------------------------------
// model of first sender
module station1
    
    // LOCAL STATE
    s1 : [0..5];
    // s1=0 - initial state
    // s1=1 - transmit
    // s1=2 - collision (set backoff)
    // s1=3 - wait (bus busy)
    // s1=4 - successfully sent
    
    // LOCAL CLOCK
    x1 : [0..max(lambda,slot)];
    
    // BACKOFF COUNTER (number of slots to wait)
    bc1 : [0..M];
    
    // COLLISION COUNTER
    cd1 : [0..K];
    
    // start sending
    [send1] (s1=0) -> (s1'=1) & (x1'=0); // start sending
    [busy1] (s1=0) -> (s1'=2) & (x1'=0) & (cd1'=min(K,cd1+1)); // detects channel is busy so go into backoff
    
    // transmitting
    [time] (s1=1) & (x1<lambda) -> (x1'=min(x1+1,lambda)); // let time pass
    [end1]  (s1=1) & (x1=lambda) -> (s1'=4) & (x1'=0); // finished
    [cd]   (s1=1) -> (s1'=2) & (x1'=0) & (cd1'=min(K,cd1+1)); // collision detected (increment backoff counter)
    [cd] !(s1=1) -> (s1'=s1); // add loop for collision detection when not important
    
    // set backoff (no time can pass in this state)
    // probability depends on which transmission this is (cd1)"""

    for i in range(K):
        s+= f"\n\t [] s1=2 & cd1={i+1} -> "
        for j in range(2**(i+1)):
            s += f"1/{2**(i+1)} : (s1'=3) & (bc1'={j})"
            if j != 2**(i+1)-1:
                s += " + "
            else:
                s += ";"
    
    s+="""
    // wait until backoff counter reaches 0 then send again
    [time] (s1=3) & (x1<slot) -> (x1'=x1+1); // let time pass (in slot)
    [time] (s1=3) & (x1=slot) & (bc1>0) -> (x1'=1) & (bc1'=bc1-1); // let time pass (move slots)
    [send1] (s1=3) & (x1=slot) & (bc1=0) -> (s1'=1) & (x1'=0); // finished backoff (bus appears free)
    [busy1] (s1=3) & (x1=slot) & (bc1=0) -> (s1'=2) & (x1'=0) & (cd1'=min(K,cd1+1)); // finished backoff (bus busy)
    
    // once finished nothing matters
    [time] (s1>=4) -> (x1'=0);

endmodule

//----------------------------------------------------------------------------------------------------------------------------

// construct further stations through renaming
"""

    for i in range(1,N):
        s += f"module station{i+1}=station1[s1=s{i+1},x1=x{i+1},cd1=cd{i+1},bc1=bc{i+1},send1=send{i+1},busy1=busy{i+1},end1=end{i+1}] endmodule\n"

    s+="""
//----------------------------------------------------------------------------------------------------------------------------

// reward structure for expected time
rewards "time"
    [time] true : 1;
endrewards

//----------------------------------------------------------------------------------------------------------------------------

// labels/formulae
"""

    s+=f"label \"all_delivered\" = "
    for i in range(N):
        s += f"s{i+1}=4"
        if i != N-1:
            s += "&"
    s += ";\n"

    s+=f"label \"one_delivered\" = "
    for i in range(N):
        s += f"s{i+1}=4"
        if i != N-1:
            s += "|"
    s += ";\n"

    s+=f"label \"collision_max_backoff\" = "
    for i in range(N):
        s += f"(cd{i+1}=K & s{i+1}=1 & b=2)"
        if i != N-1:
            s += "|"
    s += ";\n"

    s += "formula min_backoff_after_success = min("
    for i in range(N):
        s += f"s{i+1}=4?cd{i+1}:K+1"
        if i != N-1:
            s += ","
    s += ");\n"

    s += "formula min_collisions = min("
    for i in range(N):
        s += f"cd{i+1}"
        if i != N-1:
            s += ","
    s += ");\n"

    s += "formula max_collisions = max("
    for i in range(N):
        s += f"cd{i+1}"
        if i != N-1:
            s += ","
    s += ");\n"

    with open(f"{MODELS_DIR}/csma.N_{N}.K_{K}.prism", "w+") as f:
        f.writelines(s)

def generate_models_csma():
    #Note: Higher values of second parameter cause too much time in creating models. So we need to be careful
    param_set = {(2,2), (2,3), (2,5), (2,6), (2,7), (2,8), (2,10), (2,11), (2,12), (2,13), (2,17), 
                 (2,2), (3,2), (5,2), (6,2), (7,2), (8,2), (13,2), (21,2), (34,2), (55,2), (89,2), (144,2)}
    
    for (N, K) in param_set:
        print(N,K)
        gen_csma(N, K)
        task_yml = dict(format_version='2.0')
        task_yml["input_files"] = f"../../{MODELS_DIR}/csma.N_{N}.K_{K}.prism"
        options_dict = dict(property_file='csma.props')
        task_yml["options"] = options_dict
        for s in {'all_before_max', 'all_before_min', 'some_before'}:
            task_yml["input_files"] = f"../../{MODELS_DIR}/csma.N_{N}.K_{K}.prism"
            task_yml["options"]["property_tag"] = s
            with open(f"{TASK_DIR}/csma.{s}.N_{N}.K_{K}.prism.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)
            
            task_yml["input_files"] = f'../../{JANI_MODELS_DIR}/csma.N_{N}.K_{K}.jani'
            with open(f"{TASK_DIR}/csma.{s}.N_{N}.K_{K}.jani.yml","w+") as f:
                yaml.dump(task_yml, f, default_flow_style=False)

#generate_models_ij()
generate_models_pnueli()
generate_models_philosophers()
generate_models_consensus()
generate_models_csma()

generate_tasks_zeroconf()
generate_tasks_mer()
generate_tasks_firewire()
generate_tasks_pacman()