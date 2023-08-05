all: ground-state-annealer

ground-state-annealer: ground-state-annealer.c
	gcc -g -O2 -Wall ground-state-annealer.c -o ground-state-annealer -lm
