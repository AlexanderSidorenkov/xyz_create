!> \file
!> \brief Создает .xyz файл с атомами графена.
!> \details Создается прямоугольный лист графена в плоскости XY повернутый на некоторый угол.
!> \warning Для графита учесть его модификацию
!> \param a постоянная решетки графена.
!> \param b расстояние между листами графена в графите.
!> \param m масса атома углерода.
!> \param hexagonal координаты атомов графена в элементарной ячейке.
!> \param hexagonal_lv вектора элементарной ячейки.
!> \param hexagonal_rectangle_lv вектора прямоугольной элементарной ячейки.
program graphene
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

call create(hexagonal,hexagonal_ids,hexagonal_lv)
call cut_transformed_cell(mat,hexagonal_rectangle_lv,-1)
call stretch(a*(e(:,1)+e(:,2))+b*e(:,3))
!call stretch(61.449900/61.464500*e(:,1)+e(:,2)+e(:,3))
call replicate((/19,11,1/))
call shift(20.*e(:,3))
open(out_id,file=trim(output_path())//'graphene_rectangle_47x47.xyz')
call write_xyz(out_id,m,gr_atom_names)
close(out_id)

end program graphene