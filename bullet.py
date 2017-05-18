import sys #Para tomar los datos del archivo input
import imageio #Para hacer el gif
import os #Para crear directorio temporal
import shutil #Para eliminar directorio temporal
import numpy as np
import matplotlib.pyplot as plt

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
delv=V/(Nv)
x=np.arange(L_min, L_min+L, delx)
v=np.arange(V_min, V_min+V, delv)

if os.path.isfile("phase_four_dat.txt"):
    bullet=np.genfromtxt("phase_four_dat.txt")
if os.path.isfile("phase_rela_dat.txt"):
    bullet=np.genfromtxt("phase_rela_dat.txt")

delxs=[]
delvs=[]
times=np.linspace(0, T*deltat, int(T/skip))

for i in range(int(T/skip)):
    bull_act=bullet[i*Nx:(i+1)*Nx,:]
    donde=np.where(bull_act==bull_act[512:-512,512:-512].max())
    if len(donde[0]) == 1:
        donde1=[donde[0][0], donde[1][0]]
        diffx=np.abs(2048-2*donde1[1])*delx
        diffv=np.abs(2048-2*donde1[0])*delv
        print i, donde1
    else:
        donde1=[donde[0][0], donde[1][0]]
        donde2=[donde[0][1], donde[1][1]]
        diffx=np.abs(donde1[1]-donde2[1])*delx
        diffv=np.abs(donde1[0]-donde2[0])*delv
        print i, donde1, donde2
    delxs.append(diffx)
    delvs.append(diffv)

fig=plt.figure()
plt.title(r"Distancia entre maximos del $bullet$ $cluster$")
plt.xlabel(r"Tiempo ($u.t.$)")
plt.ylabel(r"Distancia entre maximos ($u.l.$, $u.l./u.t.$)")
plt.plot(times, delxs, 'bo', label=r"$\Delta x$")
plt.plot(times, delxs, 'b-')
plt.plot(times, delvs, 'ro', label=r"$\Delta v$")
plt.plot(times, delvs, 'r-')
plt.legend(framealpha=0.5)
plt.grid()
plt.savefig("bullet_max.pdf", format='pdf')
plt.close()
