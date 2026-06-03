F77  = gfortran -fdefault-integer-8 -O1 -std=legacy

EX_OBJS = extrapol.o ewald.o setbox.o shells.o sk0.o vewald_lr_k.o \
          spline.o ujasee.o optpot_Ujas.o

FITPN_OBJS = fillk.o fitpnnew.o basis.o splint3D.o

LINPACK_OBJS = dpoco.o dposl.o blas_min.o

OBJS = $(EX_OBJS) $(FITPN_OBJS) $(LINPACK_OBJS)

.SUFFIXES: .f

.f.o:
	$(F77) -c $<

extrapol: $(OBJS)
	$(F77) -o fse_extrapol $(OBJS)

clean:
	rm -f *.o fse_extrapol
