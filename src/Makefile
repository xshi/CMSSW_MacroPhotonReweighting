OBJS = fits.o photonpreselection.o reweight.o RooZPtPdf.o main.o
HEADERS = fits.h photonpreselection.h reweight.h RooZPtPdf.h
CC = g++
CFLAGS = -O3 -ansi -Wall `root-config --cflags` -I../../MacroLibrary -I$(ROOFITSYS)/include -std=c++0x
LFLAGS = -Wl,-Bstatic -L../../MacroLibrary -lHZZ2l2nu -Wl,-Bdynamic `root-config --libs` -L$(ROOFITSYS)/lib -lRooFitCore -lRooFit -lFoam -lMathMore -lMinuit -O3 -ansi --std=c++0x
NAME = photons

$(NAME) : $(OBJS)
	$(CC) $(OBJS) -o $(NAME) $(LFLAGS)

fits.o : fits.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c fits.cpp

photonpreselection.o : photonpreselection.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c photonpreselection.cpp

reweight.o : reweight.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c reweight.cpp

RooZPtPdf.o : RooZPtPdf.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c RooZPtPdf.cpp

main.o : main.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c main.cpp

clean :
	rm -f *.o $(NAME)
