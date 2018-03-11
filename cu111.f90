!> \file
!> \brief Создает .xyz файл c кристаллом меди.
!> \details Кристалл имеет поверхность (111) в плоскости XY. 
!  \warning 
!> \param a постоянная решетки меди.
!> \param m масса атома меди.
!> \param face_centered_111 координаты атомов меди в элементарной ячейке.
!> \param cubic_lv_111 вектора элементарной ячейки для получения прямоугольного кристалла.
!> \param hexagonal_lv вектора ромбической элементарной ячейки.
program cu111
use crystal_lattice
implicit none

real:: a=3.6147,m=63.5463,angle,mat(3,3)
real,parameter:: hexagonal_lv(3,3)=reshape((/ sqrt(3.)/2, 0.5, 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/))
character(len=32):: cu_atom_names(2)=(/'CU     ','CUfixed'/)
real,parameter:: &
face_centered_111(3,4)=reshape((/ &
0., 0., 0., &
-0.5/sqrt(6.), -0.5/sqrt(2.), 1./sqrt(3.), &
1./sqrt(6.), 0. , 1./sqrt(3.), &
-0.5/sqrt(6.), 0.5/sqrt(2.), 1./sqrt(3.) &
/),shape=(/3,4/)),&
cubic_lv_111(3,3)=reshape((/ &
1./sqrt(6.), 1./sqrt(2.), 1./sqrt(3.), &
-2./sqrt(6.), 0. , 1./sqrt(3.), &
1./sqrt(6.), -1./sqrt(2.), 1./sqrt(3.) &
/),shape=(/3,3/))
integer,parameter:: face_centered_ids(4)=(/ 1, 1, 1, 1/)
integer:: i

angle = (90.)*pi/180.
mat = reshape((/ cos(angle), -sin(angle), 0., sin(angle), cos(angle), 0., 0., 0., 1. /),shape=(/3,3/))

!call create_turn_cut(face_centered_111,face_centered_ids,cubic_lv_111,mat,hexagonal_lv,0)
!call stretch(a*(e(:,1)+e(:,2)+e(:,3)))
!call replicate((/1,1,3/))
!do i=1,size(atomid)
!	if(X(3,i)<3.) atomid(i) = 2
!enddo
!call shift(e(:,3)*1.5*a/sqrt(3.))
!call cut()
!call shift(-e(:,3)*1.5*a/sqrt(3.))
!open(out_id,file='cu111_parallelogram_8L.xyz')
!call write_xyz(out_id,m,cu_atom_names)
!close(out_id)
!call shift(e(:,3)*5.5*a/sqrt(3.))
!call cut()
!call shift(-e(:,3)*5.5*a/sqrt(3.))
!open(out_id,file='cu111_parallelogram_4L.xyz')
!call write_xyz(out_id,m,cu_atom_names)
!close(out_id)

!call reset()

call create(face_centered_111,face_centered_ids,cubic_lv_111)
call cut_transformed_cell(mat,e,0)
call stretch(a*(e(:,1)+e(:,2)+e(:,3)))
call replicate((/23,23,3/))
call shift(e(:,3)*1.5*a/sqrt(3.))
call cut()
call shift(-e(:,3)*1.5*a/sqrt(3.))
do i=1,size(ATOM_IDS)
	if(POSITIONS(3,i)<3.) ATOM_IDS(i) = 2
enddo
open(out_id,file=trim(output_path())//'cu111_rectangle_8L.xyz')
call write_xyz(out_id,m,cu_atom_names)
close(out_id)
!call shift(e(:,3)*5.5*a/sqrt(3.))
!call cut()
!call shift(-e(:,3)*5.5*a/sqrt(3.))
!open(out_id,file='cu111_rectangle_4L.xyz')
!call write_xyz(out_id,m,cu_atom_names)
!close(out_id)

end program cu111