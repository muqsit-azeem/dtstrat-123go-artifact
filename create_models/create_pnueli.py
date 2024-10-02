for NUM in range(16,65):#[4,5,6,7,8,9,10,11,12,13,14,15]:
#NUM = 3

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
        
    with open(f"paper_models/pnueli-zuck.{NUM}.prism","w+") as f:
        for line in lines:
            f.write(line+"\n")
