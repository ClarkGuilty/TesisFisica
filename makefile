inizio.pdf : plot.py output.txt
	python plot.py

output.txt : a.out
	./a.out > output.txt

a.out : LB1D.c
	gcc -lm LB1D.c

clean:
	rm -f a.out output.txt inizio.pdf
