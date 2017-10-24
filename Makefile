#
FC=pgf90

SRCS =  module_mp_suliaharrington.f90 \
	netcdf_plugin.f90 \
	wrapper.f90

FFLAGS = -r8 -C #-O3 -ffree-form -ffree-line-length-none
INCS =  -I/opt/local/include
LIBS =  -L/opt/local/lib -lnetcdff -lnetcdf

OBJS= ${SRCS:.f90=.o}

agg.exe: ${SRCS}
	$(FC) -c ${FFLAGS} ${INCS} ${SRCS} 
	$(FC) -o $@ ${OBJS} ${LIBS}

clean:
	rm -f *.mod *.o core *~ 

clobber:
	make clean
	rm -f agg.exe
	rm -f OUTPUT.NC
