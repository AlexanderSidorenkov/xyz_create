gfortran.exe -fdefault-real-8 -o crystal_lattice.o -c crystal_lattice.f90
gfortran.exe -fdefault-real-8 -o graphene graphene.f90 crystal_lattice.o
graphene.exe -o D:\science2\xyz_cells\
pause