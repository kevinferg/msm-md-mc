# Makefile

CC = gcc
RM = del
srcs = $(wildcard *.c)
#$(info srcs is [$(srcs)])
objs = $(srcs:.c=.o)
deps = $(objs:.o=.d)
incs = $(wildcard *.h)



all: mdsim.exe

mdsim.exe: $(objs)
	$(CC) $^ -o $@

%.o: %.c
	$(CC) -MMD -MP -c $< -o $@

.PHONY: clean

clean:
	$(RM) $(objs) $(deps)
