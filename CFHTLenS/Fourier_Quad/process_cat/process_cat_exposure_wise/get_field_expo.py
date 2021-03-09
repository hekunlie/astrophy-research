import os

data_path = "/lustre/home/acct-phyzj/phyzj/CFHT/i"
fns = os.listdir(data_path)

fields = []
expo_name = []
for fn in fns:
    if "w" in fn:
        fields.append(fn)
        expos = os.listdir(data_path + "/%s/science"%fn)
        temp = [fn]
        for expo in expos:
            if ".fits" in expo:
                nm = expo.split("_")[0]
                if nm not in temp:
                    temp.append(nm)
        strs = "\t".join(temp) +"\n"
        expo_name.append(strs)

logs = []
for i in range(len(fields)):
    nms = expo_name[i].split("\n")[0].split("\t")
    for j in range(1,len(nms)):
        expo_path = data_path + "/%s/science/%s_1.fits"%(nms[0],nms[j])
        if not os.path.exists(expo_path):
            logs.append(expo_path+"\n")
        else:
            logs.append("Found\n")

with open("./logs.dat", "w") as f:
    f.writelines(logs)

with open("./field_expo.dat","w") as f:
    f.writelines(expo_name)


