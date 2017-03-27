
#LIBS := /home/zpw/luochong/intel/mkl/lib/32/libguide*\
	/home/zpw/luochong/intel/mkl/lib/32/libmkl_lapack* \
	/home/zpw/luochong/intel/mkl/lib/32/libmkl_p4* \
	/home/zpw/luochong/intel/mkl/lib/32/libmkl.so\
	/home/zpw/luochong/intel/mkl/lib/32/libvml.so\
 	/home/zpw/luochong/intel/mkl/lib/32/libmkl_vml_p4.so 
LIBS := -lm

F90C := /usr/local/mpich-1.2.4..8a/bin/mpif90
F90FLAGS := -c -r8 -w

#DEBUG := -g

OBJS := sparse_matrix.o \
  amg_solver.o \
  triangulation.o \
  h1_fem_space.o \
  sphere.o \
  fftconv.o \
  main.o

all : main

main : $(OBJS) 
	$(F90C) -o main $(OBJS) $(LIBS)

%.o : %.f90
	$(F90C) $(F90FLAGS) $(DEBUG) -o $@ $<

%.o : %.for
	$(F90C) $(F90FLAGS) $(DEBUG) -o $@ $<
