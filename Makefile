CXX = g++
CC = gcc
CPPFLAGS = -std=c++17
CCFLAGS = -Wall -O3

compile: src/partrace.c
	$(CC) $(CCFLAGS) src/partrace.c -o partrace -lm

clean:
	rm -rf partrace