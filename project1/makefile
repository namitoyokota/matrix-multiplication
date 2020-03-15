all:
	make compile
	make submit

compile:
	gcc ./project1.c -o ./project1

submit:
	sbatch ./project1.sbatch

open:
	vi project1.out

clean:
	rm -f project1	
	rm -f *.out
