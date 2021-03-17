# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
CFLAGS  = -Wall -lgomp

#Give them when compile 
#File name of the matrix
FILENAME = ""
#Is this thing weighted?
WEIGHTED = ""

target: compile run

compile: dijkstra_betweenness.c
	$(CC) dijkstra_betweenness.c mmio.c $(CFLAGS) -o dijkstra_betweenness

run: 
	./dijkstra_betweenness $(FILENAME) $(WEIGHTED)

clean:
	rm dijkstra_betweenness