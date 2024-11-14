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

all: src/partrace.c
	$(CC) $(CCFLAGS) src/partrace.c -o partrace -lm
	$(CC) $(CCFLAGS) src/find_grain_velocities.c -o find_grain_velocities -lm
	$(CC) $(CCFLAGS) -shared -o src/libpartrace.so -fPIC src/partrace.c
	$(CC) $(CCFLAGS) -shared -o src/libgrainvelocities.so -fPIC src/find_grain_velocities.c

clean:
	rm -rf partrace find_grain_velocities
	rm -rf src/*.so