CFLAGS = -Wall -g -I../boole/src -L../boole/src

ifeq ($(shell hostname),serveur-imath02.univ-tln.fr)
    CFLAGS = -O2 -I../src -L../src
endif


#CFLAGS = -O2 -I../src -L../src

all : degree.exe inv.exe mmf.exe pi.exe

#mkpi.exe pi.exe mkstab.exe 

pi.exe :     pi.c
	gcc $(CFLAGS)   $^ -o $@  -lgmp -lboole  -lm

degree.exe :     degree.c
	gcc $(CFLAGS)   $^ -o $@  -lgmp -lboole  -lm


mkstab.exe :     mkstab.c
	gcc $(CFLAGS)   $^ -o $@  -lgmp -lboole  -lm

mkpi.exe :    mkpi.c
	gcc $(CFLAGS)   $^ -o $@  -lgmp -lboole  -lm

mmf.exe :    mmf.c
	gcc $(CFLAGS)   $^ -o $@  -lgmp -lboole  -lm
inv.exe :    inv.c
	gcc $(CFLAGS)   $^ -o $@  -lgmp -lboole  -lm
clean :
	rm -f *.exe
