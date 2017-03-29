import numpy as np
import pylab as plt

vals_2_3 = []
vals_2_1 = []

with open("roche.e17875017") as f:
    lines = f.readlines()
    for line in lines:
        if "-inf" not in line:
            val = float(line.split(" ")[1].rstrip())
            vals_2_3.append(val)


with open("roche.e17875024") as f:
    lines = f.readlines()
    for line in lines:
        if "-inf" not in line:
            val = float(line.split(" ")[1].rstrip())
            vals_2_1.append(val)

samp_2_1 = np.random.choice(vals_2_1, 1000)
samp_2_3 = np.random.choice(vals_2_3, 1000)

plt.plot(samp_2_1,'o')
plt.plot(samp_2_3,'o')
plt.ylim(-440,-425)
plt.ylabel("logprob")
plt.xlabel("Sample Number")
plt.savefig("logprob.png", dpi=150)
