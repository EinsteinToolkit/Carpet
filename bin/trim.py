import re
import os

vars = {}
for var in re.findall(r'\S+',os.environ["RDWR_DEBUG_VARS"]):
  vars[var] = 1
pat = r"^RDWR:\s*(" +("|".join(vars.keys()))+r")\s+:="

def trimf(fname,wname):
  with open(fname,"r") as fd:
    with open(wname,"w") as fw:
      for line in fd.readlines():
        line = line.strip()
        if re.match(pat,line, flags=re.IGNORECASE):
          print(line,file=fw)
        elif re.match(r'^[/\\]==',line):
          print(line,file=fw)
        elif re.match(r'^(Iteration|SymBC:)',line):
          print(line,file=fw)
        elif re.match(r'^\s*update boundary',line):
          print(line,file=fw)

trimf("CCTK_Proc0.yes","x.yes")
trimf("CCTK_Proc0.no","x.no")
