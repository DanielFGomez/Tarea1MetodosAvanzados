tresCuerpos.pdf: tresCuerpos.py tresCuerpos.dat tresCuerposRK.dat tresCuerposA.dat energiaRK.dat energia.dat
	python tresCuerpos.py

tresCuerpos.dat : a.out
	./a.out

tresCuerposRK.dat : a.out
	./a.out

tresCuerposA.dat : a.out
	./a.out

energia.dat : a.out
	./a.out

energiaRK.dat : a.out
	./a.out

a.out : tresCuerpos.c
	cc tresCuerpos.c -lm
