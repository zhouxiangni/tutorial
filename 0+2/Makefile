
LIBS := /usr/local/intel/mkl/lib/32/libguide* \
        /usr/local/intel/mkl/lib/32/libmkl_lapack* \
        /usr/local/intel/mkl/lib/32/libmkl_p3* \
        /usr/local/intel/mkl/lib/32/libmkl.so \
        /usr/local/intel/mkl/lib/32/libvml.so \
        /usr/local/intel/mkl/lib/32/libmkl_vml_p3*

F90C := ifc
F90FLAGS := -c -r8 -w

#DEBUG := -g

OBJS := sparse_matrix.o \
  amg_solver.o \
  triangulation.o \
  h1_fem_space.o \
  sphere.o \
  main.o

all : main

main : $(OBJS) 
	$(F90C) -o main $(OBJS) $(LIBS)

%.o : %.f90
	$(F90C) $(F90FLAGS) $(DEBUG) -o $@ $<
