all: ground-state-annealer time-evolution long-range-model ground-state-annealer-periodic

CFLAGS = -g -Ofast -march=native -ffast-math -Wall
LDFLAGS = -lm
ground-state-annealer: ground-state-annealer.c
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)

ground-state-annealer-periodic: ground-state-annealer-periodic.c
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)

time-evolution: time-evolution.c
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)

long-range-model: long-range-model.c
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)


clean:
	rm -f *.o ground-state-annealer benchmark time-evolution

.PHONY: clean
#ground-state-annealer: ground-state-annealer.c
#	gcc -g -O2 -Wall ground-state-annealer.c -o ground-state-annealer -lm
