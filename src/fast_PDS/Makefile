OBJS = fpd_main.o fpd_sampling.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -std=c++0x -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

FPD_Sampling: $(OBJS)
				$(CC) $(LFLAGS) $(OBJS) -o FPD_Sampling
				
fpd_main.o: fpd_main.cpp fpd.h
			$(CC) $(CFLAGS) $(LFLAGS) fpd_main.cpp

fpd_sampling.o: fpd_sampling.cpp fpd.h
			$(CC) $(CFLAGS) $(LFLAGS) fpd_sampling.cpp
				
clean:
		\rm *.o *~ FPD_Sampling
