from itertools import repeat,chain

for NUM in range(25,32):#range(17,25):
    #NUM = 4
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

    with open(f"paper_models/ij.{NUM}.prism", "w+") as f:
        line_sep = repeat("\n")
        line_iter = chain.from_iterable(zip(lines, line_sep))
        f.writelines(line_iter)
        
    with open(f"paper_models/ij.{NUM}.props", "w+") as f:
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
