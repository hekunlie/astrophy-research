import os
from sys import argv
import time

file_num = int(argv[1])

for i in range(file_num):
    with open("./log/cal_%d.dat"%i,"r") as f:
        # contents = f.readlines()

        first_line = f.readline()
        off = -50
        while True:
            f.seek(off, 2)
            lines = f.readlines()
            if len(lines) >= 2:
                last_line = lines[-1]
                break
            off *= 2

    log_time = last_line.split()[1]
    log_time_ = log_time.split(":")
    # print(log_time_, log_time)
    log_hour, log_min, log_sec = int(log_time_[0]), int(log_time_[1]),int(log_time_[2])

    sys_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
    curr_time = sys_time.split()[1].split(":")
    curr_hour, curr_min, curr_sec = int(curr_time[0]), int(curr_time[1]),int(curr_time[2])

    delta_hour, delta_min, delta_sec = curr_hour - log_hour, curr_min - log_min, curr_sec - log_sec

    if delta_hour > 0:
        print("cal_%d.dat %d hours %s %s"%(i, delta_hour, log_time, sys_time))
        continue
    if delta_min > 2:
        print("cal_%d.dat %d min %s %s"%(i, delta_min,log_time, sys_time))
        continue