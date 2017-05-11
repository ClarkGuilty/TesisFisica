movimiento.gif : Constantes.txt acc_dat.txt dens_dat.txt pot_dat.txt vels_dat.txt plot.py
	time python plot.py

Constantes.txt acc_dat.txt dens_dat.txt pot_dat.txt vels_dat.txt: a.out
	time ./a.out

a.out : LB1D.c
	time gcc -lm -lfftw3 LB1D.c

clean:
	rm -fr a.out Constantes.txt acc_dat.txt dens_dat.txt pot_dat.txt phase_four_dat.txt phase_rela_dat.txt vels_dat.txt Fourier.gif Relaxation.gif tempFourier/ tempRelajacion/
