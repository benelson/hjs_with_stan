import numpy as np

class Observation:
    t = None
    rv = None
    err = None
    off = None
    nPoints = 0
    error = 0

class get_obs(Observation):
    def __init__(self, filename):
        """
            Generates fake observations. 
        """
        with open(filename) as f:
            f.readline()
            lines = f.read().splitlines()

        for i, line in enumerate(lines):
            lines[i] = list(map(float, line.split('\t')))

        data = np.transpose(lines)
            
        self.t = data[0]
        self.rv = data[1]
        self.err = data[2]
        try:
            self.off = list(map(int, data[3]))
        except IndexError:
            self.off = [0 for i in self.t]
        self.nPoints = len(data[0])

