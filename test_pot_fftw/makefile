Potential.pdf : test.py output.txt
	python test.py

output.txt : a.out
	./a.out > output.txt

a.out : pot_fourier.c
	gcc -lm -lfftw3 pot_fourier.c
