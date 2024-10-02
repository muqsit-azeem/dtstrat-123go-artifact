import argparse
import math

parser = argparse.ArgumentParser("Create rabin")
parser.add_argument("range_min", type=int)
parser.add_argument("range_max", type=int)
parser.add_argument("-loops", action='store_true')
args = parser.parse_args()

for NUM_PROCESSES in range(args.range_min, args.range_max + 1):
    # CREATE THE PRISM FILE
    lines = ['// N-processor mutual exclusion [Rab82]\n', '// gxn/dxp 03/12/08\n', '\n',
             '// to remove the need for fairness constraints for this model it is sufficent\n',
             '// to remove the self loops from the model \n', '\n',
             '// the step corresponding to a process making a draw has been split into two steps\n',
             '// to allow us to identify states where a process will draw without knowing the value\n',
             '// randomly drawn\n',
             '// to correctly model the protocol and prevent erroneous behaviour, the two steps are atomic\n',
             '// (i.e. no other process can move one the first step has been made)\n',
             '// as for example otherwise an adversary can prevent the process from actually drawing\n',
             '// in the current round by not scheduling it after it has performed the first step\n', '\n', 'mdp\n',
             '\n', '// size of shared counter\n']
    K = 4+math.ceil(math.log2(NUM_PROCESSES))
    line = f"const int K = {K}; // 4+ceil(log_2 N)\n"
    lines.append(line)
    lines.extend(['\n',
         '// global variables (all modules can read and write)\n',
         'global c : [0..1]; // 0/1 critical section free/taken\n', 'global b : [0..K]; // current highest draw\n',
         'global r : [1..2]; // current round\n', '\n', '// formula for process 1 drawing\n',
         'formula draw = p1=1 & (b<b1 | r!=r1);\n', '\n', '// formula to keep drawing phase atomic\n',
         '// (a process can only move if no other process is in the middle of drawing)\n'])
    addition = ""
    for n in range(2,NUM_PROCESSES+1):
        addition+=f"draw{n}=0&"
    addition = addition[:-1]
    line = f"formula go = ({addition});\n"
    lines.append(line)
    lines.extend(['\n', 'module process1\n', '\n', '\tp1 : [0..2]; // local state\n',
         '\t//  0 remainder\n', '\t//  1 trying\n', '\t//  2 critical section\n',
         '\tb1 : [0..K]; // current draw: bi\n', '\tr1 : [0..2]; // current round: ri\n',
         '\tdraw1 : [0..1]; // performed first step of drawing phase\n', '\n',
                  '\t// remain in remainder\n'])
    if args.loops:
        lines.append("\t[] go & p1=0 -> (p1'=0);\n")
    else:
        lines.append("\t// [] go & p1=0 -> (p1'=0);\n",)
    lines.extend(['\t// enter trying\n', "\t[] go & p1=0 -> (p1'=1);\n",
         '\t// make a draw\n', "\t[] go & draw & draw1=0 -> (draw1'=1);\n",
         "\t[] draw1=1 -> 0.5 : (b1'=1) & (r1'=r) & (b'=max(b,1)) & (draw1'=0)\n",
         "\t         + 0.25 : (b1'=2) & (r1'=r) & (b'=max(b,2)) & (draw1'=0)\n",
         "\t         + 0.125 : (b1'=3) & (r1'=r) & (b'=max(b,3)) & (draw1'=0)\n"])

    value = 16
    for i in range(4,K):
        line = f"\t         + 1/{value} : (b1'={i}) & (r1'=r) & (b'=max(b,{i})) & (draw1'=0)\n"
        value *=2
        lines.append(line)
    value = int(value/2)
    line = f"\t         + 1/{value} : (b1'={K}) & (r1'=r) & (b'=max(b,{K})) & (draw1'=0);\n"
    lines.append(line)
    lines.extend(['\t// enter critical section and randomly set r to 1 or 2\n',
         "\t[] go & p1=1 & b=b1 & r=r1 & c=0 -> 0.5 : (r'=1) & (c'=1) & (b'=0) & (b1'=0) & (r1'=0) & (p1'=2)\n",
         "\t                                  + 0.5 : (r'=2) & (c'=1) & (b'=0) & (b1'=0) & (r1'=0) & (p1'=2);\n",
         '\t// loop when trying and cannot make a draw or enter critical section\n'])
    if args.loops:
        lines.append("\t[] go & p1=1 & r1=r & b1<=b & ((c=0 & b1!=b) | c=1) -> (p1'=p1);\n")
    else:
        lines.append("\t//[] go & p1=1 & r1=r & b1<=b & ((c=0 & b1!=b) | c=1) -> (p1'=p1);\n")
    lines.extend(['\t// leave crictical section\n',
         "\t[] go & p1=2 -> (p1'=0) & (c'=0);\n", '\t// stay in critical section\n'])
    if args.loops:
        lines.append("\t[] go & p1=2 -> (p1'=2);\n")
    else:
        lines.append("\t//[] go & p1=2 -> (p1'=2);\n")
    lines.extend(['\t\n', 'endmodule\n', '\n', '// construct further modules through renaming\n'])

    for n in range(2, NUM_PROCESSES + 1):
        line = f"module process{n} = process1 [p1=p{n}, b1=b{n}, r1=r{n}, draw1=draw{n}, draw{n}=draw1 ] endmodule\n"
        lines.append(line)

    with open(f"rabin.{NUM_PROCESSES}.prism", "w+") as f:
        f.writelines(lines)

    # CREATE THE PROPERTY
    line = '"live": Pmax=? [ F ('
    for n in range(1, NUM_PROCESSES+1):
        line += f"p{n}=2|"
    line = line[:-1]
    line += ') ] ;'

    with open(f"rabin.{NUM_PROCESSES}.props", "w+") as f:
        f.write("// Some process will eventually enter the critical section\n")
        f.write(line)
