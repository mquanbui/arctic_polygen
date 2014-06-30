OBJS = offset_main_new.o offset_poly.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -std=c++0x -c $(DEBUG) -I../clipper/cpp/
LFLAGS = $(DEBUG) -L../clipper/cpp/ -lpolyclipping_static

Offset_poly: $(OBJS)
				$(CC) -o Offset_poly $(OBJS) $(LFLAGS)

offset_main_new.o: offset_main_new.cpp offset_poly.h
				$(CC) $(CFLAGS) offset_main_new.cpp

offset_poly.o: offset_poly.cpp offset_poly.h
				$(CC) $(CFLAGS) offset_poly.cpp

clean:
		\rm *.o *~ Offset_poly
