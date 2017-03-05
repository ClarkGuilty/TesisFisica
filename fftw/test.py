import numpy as np
import matplotlib.pyplot as plt

inp=np.genfromtxt('output.txt')

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.xlim((-10,1030))
plt.plot(inp[:,0], label='Input')
plt.plot(inp[:,1], label='OutReal')
#plt.plot(inp[:,2], label='OutImag')
plt.legend()
plt.savefig('Coefs.pdf', format='pdf')
plt.close()
