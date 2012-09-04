SRCS = src/*.cpp
HEADERS = src/*.h

NAME = photons

$(NAME) : $(SRCS) $(HEADERS)
	cd src; make
	cp src/$(NAME) $(NAME)

clean :
	rm $(NAME)
	cd src; make clean
