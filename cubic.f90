program cubic
use crystal_lattice
implicit none

character(len=32):: str

call create(reshape((/ 0., 0., 0./),shape=(/3,1/)),(/ 1/),4*reshape((/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/)))
call replicate((/16,16,16/))
open(out_id,file=trim(output_path())//'cubic.xyz')
write(str,"(A)") 'A'
call write_xyz(out_id,1.,(/str/))
close(out_id)

end program cubic