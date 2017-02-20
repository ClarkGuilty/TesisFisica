import numpy as np
import matplotlib.pyplot as plt
import sys #Para tomar los datos del archivo input
import imageio #Para hacer el gif
import os #Para crear directorio temporal
import shutil #Para eliminar directorio temporal

"""fig=plt.figure()
ax=plt.axes()
plt.imshow(inp, cmap='RdPu', extent=[-1,1,-1,1]) #extent[x,x,y,y]
plt.colorbar()
plt.xlabel(r'Posicion ($x$)')
plt.ylabel(r'Velocidad ($v$)')
plt.title(r'Espacio de fase $\Delta t$ despues del inicio')
plt.savefig('dt.pdf', format='pdf')
plt.close()"""

Nx=1024
inp=np.genfromtxt('output.txt', delimiter='  ')
size=len(inp[:,0])/Nx

os.mkdir('temp')

with imageio.get_writer('./movimiento.gif', mode='I') as writer:
    for i in range(size):
        fig=plt.figure(figsize=(8,8))
        ax=plt.axes()
        plt.xlabel(r'Posicion($x$)')
        plt.ylabel(r'Velocidad($v$)')
        plt.title('Espacio de fase')
        plt.imshow(inp[i*Nx:(i+1)*Nx,:], extent=[-1,1,-1,1], cmap='BuPu') #cmap='RdPu')
        plt.colorbar()
        plt.savefig('./temp/'+str(i)+'phase.png', format='png')
        plt.close()

        image=imageio.imread('./temp/'+str(i)+'phase.png')
        writer.append_data(image)
#borra directorio temporal
#shutil.rmtree('temp')
