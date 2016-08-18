EXTLIBS := -L/usr/lib64 -Ldmn_i -lg2c -lm -lc -ldmn_i

all: example

clean:
	rm -rf *.o example

example: example.o TDoA.o
	g++ -o example example.o TDoA.o $(EXTLIBS)

example.o: example.cxx
	g++ -c example.cxx

TDoA.o: TDoA.h TDoA.cxx
	g++ -c TDoA.cxx
