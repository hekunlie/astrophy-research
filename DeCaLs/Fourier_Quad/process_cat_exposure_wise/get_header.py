from sys import argv
import os

files = os.listdir(argv[1])

for fn in files:
    try:
        with open(argv[1] + "/%s"%fn,"r") as f:
            heads = f.readline().split()
        break
    except:
        pass

print(heads)
if heads[0] == "#":
    heads.__delitem__(0)
print(heads)
header = []
for tag, h in enumerate(heads):
    header.append("%d\t%s\n"%(tag, h))
with open('cat_header.dat',"w") as f:
    f.writelines(header)