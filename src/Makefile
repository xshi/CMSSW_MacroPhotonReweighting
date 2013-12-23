OBJS = \
	fits.o \
	main.o \
	reweight.o \
	samples.o
HEADERS = \
	fits.h \
	reweight.h \
	samples.h

CC = g++
CFLAGS = -g -O3 -ansi -Wall `root-config --cflags` -I../../MacroLibrary -I$(ROOFITSYS)/include -std=c++0x
LFLAGS = -Wl,-Bstatic -L../../MacroLibrary -lHZZ2l2nu -Wl,-Bdynamic `root-config --libs` -L$(ROOFITSYS)/lib -lRooFitCore -lRooFit -lFoam -lMathMore -lMinuit -ljsoncpp -O3 -ansi --std=c++0x
NAME = photons

$(NAME) : $(OBJS) ../../MacroLibrary/libHZZ2l2nu.a
	$(CC) $(OBJS) -o $(NAME) $(LFLAGS)

fits.o : fits.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c fits.cpp

samples.o : samples.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c samples.cpp

reweight.o : reweight.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c reweight.cpp

main.o : main.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c main.cpp

clean :
	rm -f *.o $(NAME)
