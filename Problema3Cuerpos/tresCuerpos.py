import numpy as np
import matplotlib.pylab as plt

data=np.genfromtxt('tresCuerpos.dat');
plt.scatter(data[:,0],data[:,1],s=0.1)

plt.show()

