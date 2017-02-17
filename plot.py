import numpy as np
import matplotlib.pyplot as plt
import sys #Para tomar los datos del archivo input
import imageio #Para hacer el gif
import os #Para crear directorio temporal
import shutil #Para eliminar directorio temporal


Nx=1024
inp=np.genfromtxt('output.txt', delimiter='  ')
size=len(inp[:,0])/Nx
"""fig=plt.figure()
ax=plt.axes()
plt.imshow(inp, cmap='RdPu', extent=[-1,1,-1,1]) #extent[x,x,y,y]
plt.colorbar()
plt.xlabel(r'Posicion ($x$)')
plt.ylabel(r'Velocidad ($v$)')
plt.title(r'Espacio de fase $\Delta t$ despues del inicio')
plt.savefig('dt.pdf', format='pdf')
plt.close()"""

os.mkdir('temp')

with imageio.get_writer('./movimiento.gif', mode='I') as writer:
    for i in range(size):
        fig=plt.figure()
        ax=plt.axes()
        plt.xlabel(r'Posicion($x$)')
        plt.ylabel(r'Velocidad($v$)')
        plt.title('Espacio de fase')
        plt.imshow(inp[i*Nx:(i+1)*Nx,:], cmap='RdPu', extent=[-1,1,-1,1])
        plt.colorbar()
        plt.savefig('./temp/'+str(i)+'phase.png', format='png')
        plt.close()

        image=imageio.imread('./temp/'+str(i)+'phase.png')
        writer.append_data(image)
#borra directorio temporal
#shutil.rmtree('temp')
