CXX = g++
CC = gcc
CPPFLAGS = -std=c++17
CCFLAGS = -std=gnu11 -Wall -O3

partrace: src/partrace.c
	$(CC) $(CCFLAGS) src/partrace.c -o partrace -lm
	$(CC) $(CCFLAGS) -shared -o src/libpartrace.so -fPIC src/partrace.c

velocities: src/find_grain_velocities.c
	$(CC) $(CCFLAGS) src/find_grain_velocities.c -o find_grain_velocities -lm
	$(CC) $(CCFLAGS) -shared -o src/libgrainvelocities.so -fPIC src/find_grain_velocities.c

temperatures: interpolate_temps.c
	$(CC) $(CCFLAGS) interpolate_temps.c -o interpolate_temps -lm

uv: interpolate_uv.c
	$(CC) $(CCFLAGS) interpolate_uv.c -o interpolate_uv -lm

timescales: interpolate_uv.c
	$(CC) $(CCFLAGS) interpolate_timescales.c -o interpolate_timescales -lm

positions: gather_positions.c
	$(CC) $(CCFLAGS) gather_positions.c -o gather_positions -lm

all: partrace velocities temperatures positions

clean:
	rm -rf partrace find_grain_velocities
	rm -rf src/*.so
