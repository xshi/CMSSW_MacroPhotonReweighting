SRCS = src/*.cpp
HEADERS = src/*.h

NAME = photons

$(NAME) : src/$(NAME)
	cp src/$(NAME) $(NAME)

src/$(NAME) : $(SRCS) $(HEADERS) ../MacroLibrary/libHZZ2l2nu.a
	cd src; make

clean :
	rm $(NAME)
	cd src; make clean
