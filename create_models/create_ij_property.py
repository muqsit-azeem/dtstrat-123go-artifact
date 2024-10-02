import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument(
    "num",
    type=int,
)

args = parser.parse_args()

with open(f"ij.{args.num}.props","w+") as f:
    line = '"stable": Pmax=? [ F ('
    for i in range(args.num):
        part = '('
        for j in range(args.num):
            if i==j:
                part+= f"q{i+1}=1 &"
            else:
                part+= f"q{j+1}=0 &"
        part = part[:-1]
        part += ") | "
        line += part
    line = line[:-2]
    line += ")];"
    f.write(line)
