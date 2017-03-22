norm.pdf : norm.txt norm.py
	python norm.py
norm.txt : a.out
	./a.out > norm.txt

a.out : norm.c
	gcc -lfftw3 -lm -Wall norm.c

clean:
	rm -fr norm.txt a.out norm.pdf
