program create_graphene
use crystal_lattice
implicit none

real:: a=2.45858,b=3.354,m=12.0107,angle,mat(3,3)
character(len=32):: gr_atom_names(2)=(/'C_a','C_b'/)
real,parameter:: &
hexagonal(3,2)=reshape((/ 0., 0., 0., 1./sqrt(3.)/2, 0.5, 0. /),shape=(/3,2/)),&
hexagonal_lv(3,3)=reshape((/ sqrt(3.)/2, 0.5, 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/)),&
hexagonal_rectangle_lv(3,3)=reshape((/ sqrt(3.), 0., 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/))
integer,parameter:: &
hexagonal_ids(2)=(/ 1, 2/),&
hexagonal_rectangle_ids(4)=(/ 1, 2, 1, 2/)

angle = (90.)*pi/180.
mat = reshape((/ cos(angle), -sin(angle), 0., sin(angle), cos(angle), 0., 0., 0., 1. /),shape=(/3,3/))

!call create_turn_cut(hexagonal,hexagonal_ids,hexagonal_lv,mat,hexagonal_lv,-1)
!call stretch(a*(e(:,1)+e(:,2))+b*e(:,3))
!call replicate((/1,1,1/))
!call shift(20.*e(:,3))
!open(out_id,file='graphene_parallelogram.xyz')
!call write_xyz(out_id,m,gr_atom_names)
!close(out_id)

!call reset()

call create_turn_cut(hexagonal,hexagonal_ids,hexagonal_lv,mat,hexagonal_rectangle_lv,-1)
call stretch(a*(e(:,1)+e(:,2))+b*e(:,3))
call stretch(58.787514/59.005920*(e(:,1)+e(:,2)))
call replicate((/24,24,1/))
call shift(20.*e(:,3))
open(out_id,file='graphene_rectangle.xyz')
call write_xyz(out_id,m,gr_atom_names)
close(out_id)

end program create_graphene