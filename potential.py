import numpy as np
import matplotlib.pyplot as plt
Nx=1024
inp=np.genfromtxt('test.txt')

fig=plt.figure()
ax=plt.axes()
plt.xlim((0,Nx))
plt.plot(inp[:,0])
plt.savefig('densidad.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.xlim((0,Nx))
plt.plot(inp[:,1])
plt.savefig('potential.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.xlim((0,Nx))
plt.plot(inp[:,2])
plt.savefig('aceleracion.pdf', format='pdf')
plt.close()
