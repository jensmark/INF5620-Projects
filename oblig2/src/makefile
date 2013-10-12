#------------------------------------------------------------------------------
CC            = gcc
CFLAGS        = -g -c
FILE		  = wave2d_compiled

#------------------------------------------------------------------------------

all:
	$(CC) $(CFLAGS) -fPIC $(FILE).c
	$(CC) $(FILE).o -shared -o $(FILE).so

clean:
	rm -fr *.o
	rm -fr *.so
