all: clean fast

clean:
	rm -f main sol

main: HadamardMCSim.hpp main.cpp
	g++ -Wall --std=c++0x main.cpp -o main

fast: HadamardMCSim.hpp main.cpp
	g++ -D NDEBUG -Ofast main.cpp -o main

run:
	./main > out.log

sol: sol.c
	gcc -Ofast -o sol sol.c -lm && ./sol 12 5.5 sol_out.log

fast-grading:
	g++ -D NDEBUG\
		-D NOSTDEVS\
		-D DIM=20\
		-D BETA=6\
		-D NEPOCHS=1\
		-D NUMSWEEPS=10000\
		-D INITEPOCHS=0\
		-Ofast main.cpp -o yubo

timing:
	./main > /dev/null 2> &1
