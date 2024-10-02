import argparse

parser = argparse.ArgumentParser("Create philosophers")
parser.add_argument("range_min", type=int)
parser.add_argument("range_max", type=int)
args = parser.parse_args()

for NUM_PHILOSOPHERS in range(args.range_min, args.range_max+1):
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

    with open(f"philosophers.{NUM_PHILOSOPHERS}.prism", "w+") as f:
        f.writelines(lines)

    # CREATE THE PROPERTY
    line = '"eat": Pmax=? [ F ('
    for n in range(1,NUM_PHILOSOPHERS+1):
        line+= f"((p{n}>=8)&(p{n}<=9))|"
    line = line[:-1]
    line += ') ];'
    with open(f"philosophers.{NUM_PHILOSOPHERS}.props", "w+") as f:
        f.write("// Finally one of the philosophers manages to eat\n")
        f.write(line)
