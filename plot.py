import numpy as np
import matplotlib.pyplot as plt

inp=np.genfromtxt('output.txt', delimiter='  ')

fig=plt.figure()
ax=plt.axes()
plt.imshow(inp, cmap='RdPu', extent=[-1,1,-1,1]) #extent[x,x,y,y]
plt.colorbar()
plt.xlabel(r'Posicion ($x$)')
plt.ylabel(r'Velocidad ($v$)')
plt.title('Espacio de fase condicion inicial')
plt.savefig('inizio.pdf', format='pdf')
plt.close()
