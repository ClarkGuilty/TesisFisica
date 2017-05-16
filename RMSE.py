import numpy as np
import matplotlib.pyplot as plt

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

cons=np.genfromtxt("Constantes.txt")[:,1]
T=int(cons[6])
skip=int(cons[7])
deltat=cons[8]

four_256=np.genfromtxt('phase_four_dat_256.txt')
rela_256=np.genfromtxt('phase_rela_dat_256.txt')

four_512=np.genfromtxt('phase_four_dat_512.txt')
rela_512=np.genfromtxt('phase_rela_dat_512.txt')

four_1024=np.genfromtxt('phase_four_dat_1024.txt')
rela_1024=np.genfromtxt('phase_rela_dat_1024.txt')

four_2048=np.genfromtxt('phase_four_dat_2048.txt')
rela_2048=np.genfromtxt('phase_rela_dat_2048.txt')

RMSE_256=[]
RMSE_512=[]
RMSE_1024=[]
RMSE_2048=[]
times=[]

for i in range(int(T/skip)):
    rmse_temp=rmse((four_256[i*256:(i+1)*256,:]), (rela_256[i*256:(i+1)*256,:]))
    RMSE_256.append(rmse_temp)
    rmse_temp=rmse((four_512[i*512:(i+1)*512,:]), (rela_512[i*512:(i+1)*512,:]))
    RMSE_512.append(rmse_temp)
    rmse_temp=rmse((four_1024[i*1024:(i+1)*1024,:]), (rela_1024[i*1024:(i+1)*1024,:]))
    RMSE_1024.append(rmse_temp)
    rmse_temp=rmse((four_2048[i*2048:(i+1)*2048,:]), (rela_2048[i*2048:(i+1)*2048,:]))
    RMSE_2048.append(rmse_temp)

    times.append(i*deltat*skip)

fig=plt.figure()
plt.xlabel(r"Tiempo ($u.t.$)")
plt.ylabel(r"$RMSD$")
plt.title("Discrepancia entre metodo de relajacion y de Fourier")
plt.grid()
plt.plot(times,RMSE_256, 'go', alpha=0.5, label="N=256")
plt.plot(times,RMSE_256, 'g-')
plt.plot(times,RMSE_512, 'bo', alpha=0.5, label="N=512")
plt.plot(times,RMSE_512, 'b-')
plt.plot(times,RMSE_1024, 'ro', alpha=0.5, label="N=1024")
plt.plot(times,RMSE_1024, 'r-')
plt.plot(times,RMSE_1024, 'mo', alpha=0.5, label="N=2048")
plt.plot(times,RMSE_1024, 'm-')
plt.legend(framealpha=0.5, loc=2)
plt.savefig('RMSE.pdf')
plt.close()
