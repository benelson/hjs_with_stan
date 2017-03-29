import numpy as np
import matplotlib.pyplot as plt

f1 = 0.25
f2 = 1 - f1

total = 1000

def truncated_power_law(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return list((ag + (bg - ag)*r)**(1./g))



data = truncated_power_law(1., 20., g=1., size=f1*total) + truncated_power_law(2., 20., g=-1, size=f2*total)

file = open('planets_synthetic.in', 'w')

for item in data:
    file.write("%s\n" % item)

