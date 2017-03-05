import numpy as np
import matplotlib.pyplot as plt

inp=np.genfromtxt('output.txt')

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.xlim((0,1024))
plt.plot(inp)
plt.savefig('Potential.pdf', format='pdf')
plt.close()
