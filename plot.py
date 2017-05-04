import numpy as np
import matplotlib.pyplot as plt
import sys #Para tomar los datos del archivo input
import imageio #Para hacer el gif
import os #Para crear directorio temporal
import shutil #Para eliminar directorio temporal
from matplotlib import gridspec
import plotly.plotly as py

extra=0
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
delx=L/(Nx)
delv=L/(Nx)
x=np.arange(L_min, L_min+L, delx)
v=np.arange(V_min, V_min+V, delv)

dens=np.genfromtxt("dens_dat.txt")
acc=np.genfromtxt("acc_dat.txt")
pot=np.genfromtxt("pot_dat.txt")
vels=np.genfromtxt("vels_dat.txt")

phase=np.genfromtxt("phase_dat.txt")
os.mkdir('temp'+metodo)

print len(phase[:,0])/1024

f_min=np.min(phase[0:Nx,:])
f_max=np.max(phase[0:Nx,:])

for i in range(int(T/skip)):
    print str(i)+ ". El centro de masa esta en x="+str(np.dot(dens[i*Nx:(i+1)*Nx],x)*delx)
    print str(i)+ ". El centro de velocidades esta en v="+str(np.dot(vels[i*Nv:(i+1)*Nv],v)*delv)
    print str(i)+ ". La aceleracion neta es a="+str(np.trapz(acc[i*Nx:(i+1)*Nx],x))

if extra==1:
    with imageio.get_writer('./'+metodo+'.gif', mode='I') as writer:
        for i in range(int(T/skip)):
            fig=plt.figure(figsize=(18,12))
            plt.suptitle("Metodo de "+metodo, fontsize=20)
            gs=gridspec.GridSpec(4,5)

            ax1=fig.add_subplot(gs[0,0])
            ax1.plot(x, dens[i*Nx:(i+1)*Nx])
            plt.ylabel(r'Densidad ($u.m./u.l$)')
            plt.ylim((np.min(dens),np.max(dens)))

            ax2=fig.add_subplot(gs[1,0])
            ax2.plot(x, pot[i*Nx:(i+1)*Nx])
            plt.ylim((np.min(pot),np.max(pot)))
            plt.ylabel(r'Potencial ($u.l.^2/u.t^2$)')

            ax3=fig.add_subplot(gs[2,0])
            ax3.plot(x, acc[i*Nx:(i+1)*Nx])
            plt.ylim((np.min(acc),np.max(acc)))
            plt.ylabel(r'Aceleracion ($u.l./u.t.^2$)')
            plt.xlabel(r'Posicion($u.l.$)')

            ax4=fig.add_subplot(gs[3,0])
            ax4.plot(x, vels[i*Nv:(i+1)*Nv])
            plt.ylim((np.min(vels),np.max(vels)))
            plt.ylabel(r'Densidad de velocidad ($u.m u.t./u.l$)')
            plt.xlabel(r'Velocidad ($u.l/u.t.$)')

            ax5=fig.add_subplot(gs[:,1:])
            im=ax5.imshow((phase[i*Nx:(i+1)*Nx,:]), extent=[L_min,L_min+L,V_min,V_min+V], cmap='BuPu', aspect='auto', vmin=f_min, vmax=f_max)
            fig.colorbar(im)
            plt.xlabel(r'Posicion($u.l.$)')
            plt.ylabel(r'Velocidad ($u.l/u.t.$)')
            plt.title('Espacio de fase: Tiempo=' + str(i*skip*deltat)+r"$u.t.$")

            gs.update(wspace=0.5, hspace=0.5)
            fig = plt.gcf()
            plt.savefig('./temp'+metodo+'/'+str(i)+'phase.png', format='png')
            plt.close()

            image=imageio.imread('./temp'+metodo+'/'+str(i)+'phase.png')
            writer.append_data(image)
if extra==0:
    with imageio.get_writer('./'+metodo+'.gif', mode='I') as writer:
        for i in range(int(T/skip)):
            fig=plt.figure(figsize=(12,8))
            plt.suptitle("Metodo de "+metodo, fontsize=20)
            gs=gridspec.GridSpec(3,4)

            ax1=fig.add_subplot(gs[0,0])
            ax1.plot(x, dens[i*Nx:(i+1)*Nx])
            plt.ylabel(r'Densidad ($u.m./u.l$)')
            plt.ylim((np.min(dens),np.max(dens)))

            ax2=fig.add_subplot(gs[1,0])
            ax2.plot(x, pot[i*Nx:(i+1)*Nx])
            plt.ylim((np.min(pot),np.max(pot)))
            plt.ylabel(r'Potencial ($u.l.^2/u.t^2$)')

            ax3=fig.add_subplot(gs[2,0])
            ax3.plot(x, acc[i*Nx:(i+1)*Nx])
            plt.ylim((np.min(acc),np.max(acc)))
            plt.ylabel(r'Aceleracion ($u.l./u.t.^2$)')
            plt.xlabel(r'Posicion($u.l.$)')

            ax5=fig.add_subplot(gs[:,1:])
            im=ax5.imshow((phase[i*Nx:(i+1)*Nx,:]), extent=[L_min,L_min+L,V_min,V_min+V], cmap='BuPu', aspect='auto', vmin=f_min, vmax=f_max)
            fig.colorbar(im)
            plt.xlabel(r'Posicion($u.l.$)')
            plt.ylabel(r'Velocidad ($u.l/u.t.$)')
            plt.title('Espacio de fase: Tiempo=' + str(i*skip*deltat)+r"$u.t.$")

            gs.update(wspace=0.5, hspace=0.5)
            fig = plt.gcf()
            plt.savefig('./temp'+metodo+'/'+str(i)+'phase.png', format='png')
            plt.close()

            image=imageio.imread('./temp'+metodo+'/'+str(i)+'phase.png')
            writer.append_data(image)


#borra directorio temporal
#shutil.rmtree('temp')
