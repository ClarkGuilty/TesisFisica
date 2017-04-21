import numpy as np
import matplotlib.pyplot as plt
import sys #Para tomar los datos del archivo input
import imageio #Para hacer el gif
import os #Para crear directorio temporal
import shutil #Para eliminar directorio temporal
from matplotlib import gridspec
import plotly.plotly as py

cons=np.genfromtxt("Constantes.txt")[:,1]
Nx=int(cons[0])
Nv=int(cons[1])
L=cons[2]
L_min=cons[3]
V=cons[4]
V_min=cons[5]
T=int(cons[6])
skip=int(cons[7])
deltat=cons[8]
metodo=np.genfromtxt("Constantes.txt", dtype="string")[9,1]
delx=L/(Nx-1)
x=np.arange(L_min, L_min+L+delx, delx)

dens=np.genfromtxt("dens_dat.txt")
acc=np.genfromtxt("acc_dat.txt")
pot=np.genfromtxt("pot_dat.txt")

phase=np.genfromtxt("phase_dat.txt")
os.mkdir('temp'+metodo)

with imageio.get_writer('./'+metodo+'.gif', mode='I') as writer:
    for i in range(int(T/skip)):
        fig=plt.figure(figsize=(12,8))
        plt.suptitle(metodo+" Method", fontsize=20)
        gs=gridspec.GridSpec(3,4)

        ax1=fig.add_subplot(gs[0,0])
        ax1.plot(x, dens[i*Nx:(i+1)*Nx])
        plt.ylabel('Density (a.u.)')
        plt.ylim((np.min(dens),np.max(dens)))

        ax2=fig.add_subplot(gs[1,0])
        ax2.plot(x, pot[i*Nx:(i+1)*Nx])
        plt.ylim((np.min(pot),np.max(pot)))
        plt.ylabel('Potential (a.u.)')

        ax3=fig.add_subplot(gs[2,0])
        ax3.plot(x, acc[i*Nx:(i+1)*Nx])
        plt.ylim((np.min(acc),np.max(acc)))
        plt.ylabel('Acceleration (a.u.)')
        plt.xlabel(r'Position($x$)')

        ax5=fig.add_subplot(gs[:,1:])
        im=ax5.imshow(np.fliplr(phase[i*Nx:(i+1)*Nx,:]), extent=[L_min,L_min+L,V_min,V_min+V], cmap='BuPu', aspect='auto')
        fig.colorbar(im)
        plt.xlabel(r'Position($x$)')
        plt.ylabel(r'Velocity($v$)')
        plt.title('Phase Space: Time=' + str(i*skip*deltat))

        gs.update(wspace=0.5, hspace=0.5)
        fig = plt.gcf()
        plt.savefig('./temp'+metodo+'/'+str(i)+'phase.png', format='png')
        plt.close()

        image=imageio.imread('./temp'+metodo+'/'+str(i)+'phase.png')
        writer.append_data(image)

#borra directorio temporal
#shutil.rmtree('temp')
