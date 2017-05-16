import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

god=np.genfromtxt("GOD_fourier.txt")
com=god[:,0]
cov=god[:,1]
acc_neta=god[:,2]
deltat=0.5
T=64
itera=16
times=np.linspace(0, 30, itera)

fig=plt.figure(figsize=(12,4))
gs=gridspec.GridSpec(1,3)

ax1=fig.add_subplot(gs[0,0])
ax1.plot(times, com[:itera], 'go', label="N=256")
ax1.plot(times, com[:itera], 'g-')
ax1.plot(times, com[itera:2*itera], 'bo', label="N=512")
ax1.plot(times, com[itera:2*itera], 'b-')
ax1.plot(times, com[2*itera:3*itera], 'ro', label="N=1024")
ax1.plot(times, com[2*itera:3*itera], 'r-')
ax1.plot(times, com[3*itera:4*itera], 'mo', label="N=2048")
ax1.plot(times, com[3*itera:4*itera], 'm-')
ax1.grid()
ax1.legend(framealpha=0.5, loc=2)
plt.xlabel(r"Tiempo ($u.t.$)")
plt.ylabel(r"Centro de masa ($u.l.$)")
plt.title("Evolucion del centro de masa")
plt.ylim((-0.5,0.5))

ax2=fig.add_subplot(gs[0,1])
ax2.plot(times, cov[:itera], 'go', label="N=256")
ax2.plot(times, cov[:itera], 'g-')
ax2.plot(times, cov[itera:2*itera], 'bo', label="N=512")
ax2.plot(times, cov[itera:2*itera], 'b-')
ax2.plot(times, cov[2*itera:3*itera], 'ro', label="N=1024")
ax2.plot(times, cov[2*itera:3*itera], 'r-')
ax2.plot(times, cov[3*itera:4*itera], 'mo', label="N=2048")
ax2.plot(times, cov[3*itera:4*itera], 'm-')
ax2.grid()
ax2.legend(framealpha=0.5, loc=3)
ax2.set_yticks([-0.0, -0.03, -0.06, -0.09, -0.12])
plt.ylim((-0.13,0.01))
plt.xlabel(r"Tiempo ($u.t.$)")
plt.ylabel(r"Centro de velocidad ($u.l./u.t$)")
plt.title("Evolucion del centro de velocidad")

ax3=fig.add_subplot(gs[0,2])
ax3.plot(times, acc_neta[:itera], 'go', label="N=256")
ax3.plot(times, acc_neta[:itera], 'g-')
ax3.plot(times, acc_neta[itera:2*itera], 'bo', label="N=512")
ax3.plot(times, acc_neta[itera:2*itera], 'b-')
ax3.plot(times, acc_neta[2*itera:3*itera], 'ro', label="N=1024")
ax3.plot(times, acc_neta[2*itera:3*itera], 'r-')
ax3.plot(times, acc_neta[3*itera:4*itera], 'mo', label="N=2048")
ax3.plot(times, acc_neta[3*itera:4*itera], 'm-')
ax3.grid()
ax3.legend(framealpha=0.5, loc=2)
ax3.set_yticks([-0.0014, -0.0007, -0.0, 0.0007, 0.0014])
plt.ylim((-0.002,0.002))
plt.xlabel(r"Tiempo ($u.t.$)")
plt.ylabel(r"Aceleracion neta ($u.l./u.t^2$)")
plt.title("Evolucion de la aceleracion neta")

plt.tight_layout()
plt.savefig("GOD.pdf", format='pdf')
plt.close()
