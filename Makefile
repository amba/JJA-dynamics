all: ground-state-annealer benchmark

CFLAGS = -g -Ofast -march=native -ffast-math -Wall
LDFLAGS = -lm
ground-state-annealer: ground-state-annealer.c
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)

benchmark: benchmark.c
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o ground-state-annealer benchmark

.PHONY: clean
#ground-state-annealer: ground-state-annealer.c
#	gcc -g -O2 -Wall ground-state-annealer.c -o ground-state-annealer -lm
