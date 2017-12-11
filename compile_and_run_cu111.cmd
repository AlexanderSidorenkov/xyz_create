gfortran.exe -fdefault-real-8 -o crystal_lattice.o -c crystal_lattice.f90
gfortran.exe -fdefault-real-8 -o cu111 cu111.f90 crystal_lattice.o
cu111.exe -o D:\science2\xyz_cells\
pause