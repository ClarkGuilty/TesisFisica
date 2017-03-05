Coefs.pdf : test.py output.txt
	python test.py

output.txt : a.out
	./a.out > output.txt

a.out : fft-test.c
	gcc -lm -lfftw3 fft-test.c
