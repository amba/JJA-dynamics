all: ground-state-annealer

ground-state-annealer: ground-state-annealer.c
	gcc -g -O2  ground-state-annealer.c -o ground-state-annealer -lm
