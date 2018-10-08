from subprocess import Popen

source = "dimmerm3"

filter_names = ["sex2_1.5", "sex3_1.5", "sex4_1.5"]
gauss_filters = ["gauss_2.0_5x5","gauss_3.0_5x5","gauss_4.0_7x7"]

for ii, filter_name in enumerate(filter_names):
    # change the filter
    with open("./default.sex_ori","r") as f:
        contents = f.readlines()
    contents[23] = "FILTER_NAME      %s.conv   # name of the file containing the filter"%gauss_filters[ii]
    with open("./default.sex", "w") as f:
        f.writelines(contents)

    # run the estimation code
    cmd = "mpirun -np 28 python snr_est.py snr %s %s"%(source, filter_name)
    a = Popen(cmd, shell=True)
    a.wait()

    # add the data in .cat files to the final catalogs
    cmd = "mpirun -np 14 python snr_est.py add %s %s"%(source, filter_name)
    a = Popen(cmd, shell=True)
    a.wait()


