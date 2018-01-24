gfortran.exe -fdefault-real-8 -o crystal_lattice.o -c crystal_lattice.f90
gfortran.exe -fdefault-real-8 -o cu100 cu100.f90 crystal_lattice.o
cu100.exe -o D:\science2\xyz_cells\
pause