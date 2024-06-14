CXX = g++
CC = clang
CPPFLAGS = -std=c++17
CCFLAGS = -Wall -O3

recompile: compile run

run: partrace
	./partrace

compile: src/partrace.c
	$(CC) $(CCFLAGS) src/partrace.c -o partrace

clean:
	rm -rf partrace