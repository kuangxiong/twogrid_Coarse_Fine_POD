all:twogrid_APOD BuildFEM.o

OBJS1=FEMheat.o src/BuildFEM.o
OBJS2=PODheat.o src/BuildFEM.o
OBJS3=twogrid_APOD.o src/BuildFEM.o

include ${PHG_MAKEFILE_INC}

twogrid_APOD.o: /opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so twogrid_APOD.c fun.h
FEMheat.o: /opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so FEMheat.c fun.h
PODheat.o: /opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so PODheat.c fun.h

FEMheat:/opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so ${OBJS1}
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ ${OBJS1}${USER_LIB} ${LIBS}

PODheat:/opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so ${OBJS2}
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ ${OBJS2}${USER_LIB} ${LIBS}

twogrid_APOD:/opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so ${OBJS3}
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ ${OBJS3}${USER_LIB} ${LIBS}
clean:
	rm -f *.o *.text FEMheat twogrid_APOD 

