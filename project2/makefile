swap-5.4:
	module swap gnu7/7.3.0 gnu/5.4.0

swap-7.3:
	module swap gnu/5.4.0 gnu7/7.3.0

submit:
	gcc -O0 ./project2.c -o ./O0
	sbatch ./O0.sbatch
	gcc -O1 ./project2.c -o ./O1
	sbatch ./O1.sbatch
	gcc -O2 ./project2.c -o ./O2
	sbatch ./O2.sbatch
	gcc -O3 ./project2.c -o ./O3
	sbatch ./O3.sbatch

clean:
	rm -f O0
	rm -f O1
	rm -f O2
	rm -f O3