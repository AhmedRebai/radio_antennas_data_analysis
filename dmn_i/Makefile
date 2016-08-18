FOPTS      := -fno-second-underscore -O -Wpadded -Wpacked
EXTLIBS    :=  -L/usr/lib64 -lg2c -lm -lc
FORTPATH   := fortran

all: libdmn_i.a

clean:
	rm -rf *.o *.a $(FORTPATH)/dmnfgh.o $(FORTPATH)/dmnfghb.o

example: example.o dmn_i.o $(FORTPATH)/dmnfgh.o $(FORTPATH)/dmnfghb.o
	g++ -o example example.o dmn_i.o $(FORTPATH)/dmnfgh.o $(FORTPATH)/dmnfghb.o ${EXTLIBS}

example.o: example.cxx
	g++ -c example.cxx

dmn_i.o: dmn_i.h dmn_i.cxx
	g++ -c dmn_i.cxx

$(FORTPATH)/dmnfgh.o: $(FORTPATH)/dmnfgh.f
	g++ $(FOPTS) -c -o $(FORTPATH)/dmnfgh.o $(FORTPATH)/dmnfgh.f

$(FORTPATH)/dmnfghb.o: $(FORTPATH)/dmnfghb.f
	g++ $(FOPTS) -c -o $(FORTPATH)/dmnfghb.o -c $(FORTPATH)/dmnfghb.f

libdmn_i.a: dmn_i.o $(FORTPATH)/dmnfgh.o $(FORTPATH)/dmnfghb.o
	ar rv libdmn_i.a dmn_i.o $(FORTPATH)/dmnfgh.o $(FORTPATH)/dmnfghb.o
	ranlib libdmn_i.a
