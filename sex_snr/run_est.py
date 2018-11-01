from subprocess import Popen

sources = ["dimmerm3", "dimmer", "pts", "ptsb"]
for source in sources:

    if "pts" in source:
        max_radius = 8
    else:
        max_radius = 5.5

    filter_names = ["sex2_2", "sex3_2", "sex4_2"]
    gauss_filters = ["gauss_2.0_5x5", "gauss_3.0_5x5", "gauss_4.0_7x7"]

    for ii, filter_name in enumerate(filter_names):
        # change the filter
        with open("./default.sex_ori", "r") as f:
            contents = f.readlines()
        contents[23] = "FILTER_NAME      %s.conv   # name of the file containing the filter"%gauss_filters[ii]
        with open("./default.sex", "w") as f:
            f.writelines(contents)

        # run the estimation code
        cmd = "mpirun -np 14 python snr_est.py snr %s %s %.1f"%(source, filter_name, max_radius)
        a = Popen(cmd, shell=True)
        a.wait()

        # add the data in .cat files to the final catalogs
        cmd = "mpirun -np 14 python snr_est.py add %s %s %.1f"%(source, filter_name, max_radius)
        a = Popen(cmd, shell=True)
        a.wait()

        # check
        # cmd = "mpirun -np 14 python snr_est.py check %s %s %.1f"%(source, filter_name, max_radius)
        # a = Popen(cmd, shell=True)
        # a.wait()

