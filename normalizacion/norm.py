import numpy as np
import matplotlib.pyplot as plt


inp=np.genfromtxt('norm.txt')
fig=plt.figure()
ax=plt.axes()
plt.plot(inp[:,0], inp[:,1], 'k', label='Densidad')
plt.plot(inp[:,0], inp[:,2], 'c', label='Relajacion')
plt.plot(inp[:,0], inp[:,3], 'r', label='Fourier')
plt.plot(inp[:,0], inp[:,4], 'm--', label='Analitico')
plt.plot(inp[:,0], inp[:,5], 'b', label='Fourier Real')
plt.legend(framealpha=0.5, loc=2)
plt.grid()
plt.savefig('norm.pdf', format='pdf')
plt.close()
