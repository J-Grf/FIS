#! /usr/local/bin/python3.9

import subprocess as sp
import re

m_min = 232
m_max = 233
time = []
iterations = []
count = 3

for m in range(m_min, m_max):
    print("Running m = " + str(m))
    with open(str(m) + ".txt", "w") as f:
        sp.run(["./main.exe", "GMRES", "../../material/gmres_test_msr.txt", str(m), "NONE"], stdout=f)
    
    with open(str(m) + ".txt", "r") as f:
        data = f.readlines()
        for line in data:
            
            # check if we enter restarted GMRES
            count += 1
            if count == 3 and not re.match(r"restart", line):
                break

            match = re.match(r".* = ([+-]?([0-9]*[.])?[0-9]+)s CPU .*", line)
            if match:
                time.append(float(match.group(1)))
                print(time)
            
            match2 = re.match(r"m = (\d+) iterations = (\d+).*", line)
            if match2:
                iterations.append(int(match2.group(2)))

MIN = min(time)
IDX = time.index(MIN) + m_min
IT = iterations[time.index(MIN)]

print("Minimal time is: " + str(MIN) + " with restartParameter m = " + str(IDX) + " and " + str(IT) + " restartiterations")

