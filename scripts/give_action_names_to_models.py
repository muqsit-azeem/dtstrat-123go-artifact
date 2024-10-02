#!/usr/bin/env python3
import re
from pathlib import Path

def transform(input_model, output_model):
    with open(input_model,"r") as f:
        getFormula = False
        if "philosophers" in input_model or "rabin" in input_model or "pnueli" in input_model:
            getFormula = True
        formulas = []
        formulasNameToValue = {}
        count = 0
        writelines = []
        line = f.readline()
        inmodule = False
        currentmodulename = ""
        keep = []
        formulaKeep = []
        modules = {}
        while line!="":
            if line.startswith("formula") and getFormula:
                formulas.append(line)
                name = line.split(" ")[1].replace('\t','')
                value = line[line.find('=')+1:].replace(';','').replace('\n','').replace('\t','')
                if not value.startswith('('):
                    value = '(' + value + ')'
                formulasNameToValue[name] = value
#            print(line)
            if inmodule:
                if getFormula:
                    for name in formulasNameToValue.keys():
                        if name in line:
                            something = [m.start() for m in re.finditer(name, line)]
                            something.reverse()
                            for el in something:
                                if el+len(name)+1 < len(line) and not line[el+len(name)].isnumeric():
                                    #print(line[:-1]," #### ",line[el:-1])
                                    firstpart = line[:el]
                                    secondpart = line[el+len(name):]
                                    if firstpart.replace(' ','').endswith('('):
                                        rep = formulasNameToValue[name][1:-1]
                                    else:
                                        rep = formulasNameToValue[name]
                                    line = firstpart+rep+secondpart
                if line.replace(" ","").replace("\t","").startswith("["):
                    actionname = line.split('[')[1].split(']')[0]
                    if actionname == "":
                        nl = line
                        newactionname = "STEFFI_"+currentmodulename+"_STEFFI_"+str(count)
                        nl = nl.replace("[]",f"[{newactionname}]")
                        keep.append(nl)
#                        print(newactionname)
                    else:
                        keep.append(line)
                        #if '_' in actionname:
                        #    print(input_model, actionname)
                else:
                    keep.append(line)
            if line.replace(" ","").startswith("module"):
                #print("MODULE",line)
                inmodule = True
                currentmodulename = line.split("=")[0].replace(" ","")[6:].replace("\n","")
                #print("modulename: ", currentmodulename)
                keep.append(line)
                writelines.append("%%%MODULE%$"+currentmodulename)
            if line.replace(" ","").replace("\t","").endswith("endmodule\n") or line.replace(" ","").replace("\t","").endswith("endmodule"):
                if '=' in line:
                    copymodule = line.split('=')[1].split('[')[0].replace(" ","")
                    renamings = line.split('[')[1].split(']')[0].split(',')
                    leftR = [v.split('=')[0].replace(" ","") for v in renamings]
                    trueRightR = [v.split('=')[1].replace(" ","") for v in renamings]
                    rightR = [f"%temp{i}%" for i in range(len(leftR))]
                    renamings = dict(zip(leftR,rightR))
                    renamings[copymodule] = currentmodulename
                    #print(renamings)
                    keep = []
                    for l in modules[copymodule]:
                        newline = l
                        for lk in renamings.keys():
                            if lk in newline:
                                something = [m.start() for m in re.finditer(lk, newline)]
                                robll = False
                                for el in something:
                                    if el+len(lk)<len(newline) and newline[el+len(lk)].isnumeric():
                                        robll = True
                                if not robll:
                                    newline = newline.replace(lk,renamings[lk])
                                else:
                                    something.reverse()
                                    for el in something:
                                        if el + len(lk) < len(newline) and newline[el + len(lk)].isnumeric():
                                            continue
                                        else:
                                            firstpart = newline[:el]
                                            secondpart = newline[el + len(lk):]
                                            newline = firstpart + renamings[lk] + secondpart
                                    print(input_model, "PROBLEM : ",line[:-1])
                        keep.append(newline)
                    for i in range(len(keep)):
                        l = keep[i]
                        for j,el in enumerate(rightR):
                            if el in l:
                                l = l.replace(el,trueRightR[j])
                        keep[i] = l
                modules[currentmodulename] = keep
                inmodule = False
                currentmodulename = ""
                # do something with keep
                keep = []
                #print(modules)
            if not line.endswith("endmodule\n") and not line.replace(" ","").startswith("module") and not inmodule:
                writelines.append(line)
            line = f.readline()
            count += 1
    #for module in modules.keys():
        #print(module)
        ##for l in modules[module]:
        #    print(l[:-1])
    with open(output_model,"w") as g:
        for line in writelines:
            if line.startswith("%%%MODULE"):
                modulename = line.split("%$")[1]
                for l in modules[modulename]:
                    g.write(l)
            else:
                g.write(line)
            
            
if __name__=="__main__":
#   transform("paper_models/philosophers.3.prism","TEST_philo3.prism")
#    transform("paper_models/consensus.4.prism","TEST_consensus.prism")
#    transform("paper_models/csma.2-2.prism","TEST_csma.prism")
    all_files = ['consensus.10.prism',
 'consensus.2.prism',
 'consensus.3.prism',
 'consensus.4.prism',
 'consensus.5.prism',
 'consensus.6.prism',
 'consensus.7.prism',
 'consensus.8.prism',
 'csma.2-2.prism',
 'csma.2-4.prism',
 'csma.2-6.prism',
 'csma.3-2.prism',
 'csma.3-4.prism',
 'csma.3-6.prism',
 'csma.4-2.prism',
 'csma.4-4.prism',
 'csma.4-6.prism',
 'csma.5-2.prism',
 'csma.8-2.prism',
 'firewire.false.prism',
 'firewire.true.prism',
 'firewire_abst.prism',
 'firewire_dl.prism',
 'mer.prism',
 'pacman.prism',
 'philosophers.10.prism',
 'philosophers.20.prism',
 'philosophers.3.prism',
 'philosophers.30.prism',
 'philosophers.4.prism',
 'philosophers.5.prism',
 'philosophers.6.prism',
 'philosophers.7.prism',
 'philosophers.80.prism',
 'pnueli-zuck.3.prism',
 'pnueli-zuck.4.prism',
 'pnueli-zuck.5.prism',
 'pnueli-zuck.6.prism',
 'pnueli-zuck.7.prism',
 'pnueli-zuck.8.prism',
 'pnueli-zuck.9.prism',
 'rabin.10.prism',
 'rabin.3.prism',
 'rabin.4.prism',
 'rabin.5.prism',
 'rabin.6.prism',
 'rabin.7.prism',
 'rabin.8.prism',
 #'running_example.prism',
 'wlan.0.prism',
 'wlan.1.prism',
 'wlan.2.prism',
 'wlan.3.prism',
 'wlan.4.prism',
 'wlan.5.prism',
 'wlan.6.prism',
 #'wlan.prism',
 'wlan_dl.0.prism',
 'wlan_dl.1.prism',
 'wlan_dl.2.prism',
 'wlan_dl.3.prism',
 'wlan_dl.4.prism',
 'wlan_dl.5.prism',
 'wlan_dl.6.prism',
 #'wlan_dl.prism',
 #'zeroconf.prism',
 'zeroconf_dl.prism']
    all_files = ['rabin.50.prism','rabin.80.prism']#,'pnueli-zuck.30.prism',#'philosophers.50.prism',,'pnueli-zuck.60.prism','philosophers.100.prism']
    dtstrat_eva_dir = Path(__file__).resolve().parent.parent
    action_models_dir = dtstrat_eva_dir / "mdp-models" / "action-models-prism"
    for file in (dtstrat_eva_dir / "mdp-models" / "prism").glob("*.prism"):
        transform(str(file),str(action_models_dir / file.name))