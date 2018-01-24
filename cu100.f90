program cu100
use crystal_lattice
implicit none

real:: a=3.6147,m=63.5463,angle,mat(3,3)
character(len=32):: cu_atom_names(2)=(/'CU     ','CUfixed'/)
real,parameter:: &
face_centered_100(3,4)=reshape((/ &
0.,		0., 	0.,		&
0.5, 	0.5, 	0.,		&
0.5, 	0., 	0.5, 	&
0.,		0.5, 	0.5 	&
/),shape=(/3,4/))
integer,parameter:: face_centered_ids(4)=(/ 1, 1, 1, 1/)
integer:: i

angle = (0.)*pi/180.
mat = reshape((/ cos(angle), -sin(angle), 0., sin(angle), cos(angle), 0., 0., 0., 1. /),shape=(/3,3/))

call create(face_centered_100,face_centered_ids,e)
call cut_transformed_cell(mat,e,0)
call stretch(a*(e(:,1)+e(:,2)+e(:,3)))
call replicate((/17,20,4/))
do i=1,size(ATOM_IDS)
	if(POSITIONS(3,i)<3.) ATOM_IDS(i) = 2
enddo
open(out_id,file=trim(output_path())//'cu100_rectangle_8L.xyz')
call write_xyz(out_id,m,cu_atom_names)
close(out_id)

end program cu100