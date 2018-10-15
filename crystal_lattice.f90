!> \brief Модуль для создания файлов, содержащих координаты атомов кристаллов.
!> \details Все переменные доступны при подключении модуля.
!> Переменные POSITIONS, ATOM_IDS и LATTICE_VECTORS могут быть входными неявно для некоторых подпрограмм.
module crystal_lattice
implicit none

public
integer:: out_id=2017 !< Идентификатор вывода
integer:: t0=101 !< Сколько раз копировать ячейку для поиска элементарной ячейки.
real,allocatable:: POSITIONS(:,:) !< Координаты атомов в трехмерной декартовой системе.
integer,allocatable:: ATOM_IDS(:) !< Типы атомов.
real:: LATTICE_VECTORS(3,3) !< Вектора текущей ячейки. 
real,parameter:: smallvalue=10.**(-5) !< Маленькое число. 
real,parameter:: e(3,3)=reshape((/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/)) !< Единичная матрица 3 на 3. 
real,parameter:: pi=3.14159265358979323846 !< Число пи.
contains

!> Копирует текущую ячейку t(1) раз по первому вектору решетки, t(2) раз по второму, t(3) раз по третьему.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
!> \param[in,out] LATTICE_VECTORS
subroutine replicate(t)
	integer:: t(3) !< входной массив
	integer:: k,i,j,n 
	n = size(POSITIONS,dim=2)
	call allocate_more_X(n*t(1)*t(2)*t(3))
	do k=1,3
		do i=1,t(k)-1
			do j=1,n
				POSITIONS(:,n*i+j) = POSITIONS(:,j)+i*LATTICE_VECTORS(:,k)
				ATOM_IDS(n*i+j) = ATOM_IDS(j)
			enddo
		enddo
		n = n*t(k)
	enddo
	do k=1,3
		LATTICE_VECTORS(:,k) = LATTICE_VECTORS(:,k)*t
	enddo
end subroutine

!> Прибавляет ко всем координатам вектор.
!> \param[in,out] POSITIONS
subroutine shift(vec)
	real:: vec(3) !< Вектор сдвига.
	integer:: i
	do i=1,size(POSITIONS,dim=2)
		POSITIONS(:,i) = POSITIONS(:,i)+vec
	enddo
end subroutine

!> Растягивает текущую ячейку по X, Y, Z.
!> \param[in,out] POSITIONS
!> \param[in,out] LATTICE_VECTORS
subroutine stretch(vec)
	real:: vec(3) !< По X, Y, Z ячейка растянется в vec(1), vec(2), vec(3) раз соответственно.
	integer:: i,k
	do i=1,size(POSITIONS,dim=2)
		POSITIONS(:,i) = POSITIONS(:,i)*vec
	enddo
	do k=1,3
		LATTICE_VECTORS(:,k) = LATTICE_VECTORS(:,k)*vec
	enddo
end subroutine

!> Делает преобразование координат.
!> \param[in,out] POSITIONS
!> \param[in,out] LATTICE_VECTORS
subroutine transform(mat)
	real:: mat(3,3) !< Матрица перобразования 3 на 3.
	real:: vec(3)
	integer:: i,k,kk
	do i=1,size(POSITIONS,dim=2)
		do k=1,3
			vec(k) = sum(mat(:,k)*POSITIONS(:,i))
		enddo
		POSITIONS(:,i) = vec
	enddo
	do k=1,3
		do kk=1,3
			vec(kk) = sum(mat(:,kk)*LATTICE_VECTORS(:,k))
		enddo
		LATTICE_VECTORS(:,k) = vec
	enddo
end subroutine

!> Делает преобразование координат и находит элементарную ячейку для новых векторов решетки.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
!> \param[in,out] LATTICE_VECTORS
subroutine cut_transformed_cell(mat,directions,add_j)
	real:: directions(3,3) !< Новые вектора решетки - матрица 3 на 3.
	real:: mat(3,3) !< Матрица перобразования 3 на 3.
	integer:: t(3),j
	integer:: add_j !< Сдвиг номера атома в начале координат.
	t = (/t0,t0,t0/)
	j = (t0*t0*int(real(t0)/2)+t0*int(real(t0)/2)+int(real(t0)/2))*size(POSITIONS,dim=2)+1
	call transform(mat)
	call replicate(t)
	call cut_cell(directions,j+add_j)
end subroutine

!> Убирает атомы вне ячейки.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
subroutine cut()
	real,allocatable:: Xt(:,:)
	integer,allocatable:: atomidt(:)
	real:: dps(3),surf_norm(3,3),dlvs(3)
	integer:: i,k,n,kk
	logical:: in_cell(3)
	allocate(Xt(3,size(POSITIONS,dim=2)))
	allocate(atomidt(size(POSITIONS,dim=2)))
	n=0
	do k=1,3
		do kk=1,3
			surf_norm(kk,k) = LATTICE_VECTORS(mod(kk,3)+1,mod(k,3)+1)*LATTICE_VECTORS(mod(kk+1,3)+1,mod(k+1,3)+1)-&
			LATTICE_VECTORS(mod(kk+1,3)+1,mod(k,3)+1)*LATTICE_VECTORS(mod(kk,3)+1,mod(k+1,3)+1)
		enddo
		dlvs(k) = sum(surf_norm(:,k)*LATTICE_VECTORS(:,k))/sqrt(sum(surf_norm(:,k)**2))
	enddo
	do i=1,size(POSITIONS,dim=2)
		do k=1,3
			dps(k) = sum(surf_norm(:,k)*POSITIONS(:,i))/sqrt(sum(surf_norm(:,k)**2))
			in_cell(k) = dps(k)>-smallvalue .and. dps(k)<dlvs(k)-smallvalue
		enddo
		if(in_cell(1) .and. in_cell(2) .and. in_cell(3)) then
			n = n+1
			Xt(:,n) = POSITIONS(:,i)
			atomidt(n) = ATOM_IDS(i)
		endif
	enddo
	print*,n
	deallocate(POSITIONS)
	deallocate(ATOM_IDS)
	allocate(POSITIONS(3,n))
	allocate(ATOM_IDS(n))
	POSITIONS(:,:) = Xt(:,:n)
	ATOM_IDS(:) = atomidt(:n)
	deallocate(Xt)
	deallocate(atomidt)
end subroutine

!> Находит элементарную ячейку для новых векторов решетки.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
!> \param[in,out] LATTICE_VECTORS
subroutine cut_cell(directions,j)
	real:: directions(3,3) !< Новые вектора решетки - матрица 3 на 3.
	real:: dpl2,dpp
	integer:: i,k
	integer:: j !< Номер атома в начале координат.
	call shift(-POSITIONS(:,j))
	LATTICE_VECTORS = 1./smallvalue
	do i=1,size(POSITIONS,dim=2)
		do k=1,3
			dpl2 = ( sum(directions(:,k)**2)*sum(POSITIONS(:,i)**2) - sum(directions(:,k)*POSITIONS(:,i))**2 )/(sum(directions(:,k)**2))
			if(dpl2<0. .or. sqrt(dpl2)<smallvalue) then
				dpp = sqrt(sum(POSITIONS(:,i)**2))
				if(dpp<sqrt(sum(LATTICE_VECTORS(:,k)**2)) .and. dpp>smallvalue .and. ATOM_IDS(i)==ATOM_IDS(j)) then
					LATTICE_VECTORS(:,k) = POSITIONS(:,i)
				endif
			endif
		enddo
	enddo
	do k=1,3
		if(LATTICE_VECTORS(k,k)<0.) LATTICE_VECTORS(:,k) =  -LATTICE_VECTORS(:,k)
	enddo
	call cut()
end subroutine

!> Инициализирует массивы POSITIONS, ATOM_IDS и LATTICE_VECTORS значениями входных переменных.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
!> \param[in,out] LATTICE_VECTORS
!> \todo Добавить проверку размеров входных массивов.
subroutine create(x1,ids,lv)
	real:: x1(:,:) !< Трехмерные координаты атомов
	real:: lv(:,:) !< Три трехмерных вектора текущей ячейки
	integer:: ids(:) !< Типы атомов
	call allocate_more_X(size(x1,dim=2))
	POSITIONS = x1
	LATTICE_VECTORS = lv
	ATOM_IDS = ids
end subroutine

!> Сбрасывает текущие положения атомов и вектора решетки.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
!> \param[in,out] LATTICE_VECTORS
subroutine reset()
	LATTICE_VECTORS = 0.
	if(allocated(POSITIONS)) deallocate(POSITIONS,ATOM_IDS)
	allocate(POSITIONS(3,0),ATOM_IDS(0))
	print*,'reset ',size(POSITIONS,dim=2)
end subroutine

!> Выводит информацию об атомах в формате .xyz. 
!> \param[in] POSITIONS
!> \param[in] ATOM_IDS
!> \param[in] LATTICE_VECTORS
subroutine write_xyz(out_id,m,names)
	integer:: out_id !< Идентификатор вывода
	integer:: i
	real:: m(:) !< Массы атомов
	character(len=32):: names(:) !< Названия атомов
	write(out_id,*) size(POSITIONS,dim=2)
	write(out_id,'(A,9f16.6,A)') 'Lattice="',LATTICE_VECTORS,' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do i=1,size(POSITIONS,dim=2)
		write(out_id,'(7f17.7,A,A)') POSITIONS(:,i),0.,0.,0.,m(ATOM_IDS(i)),'    ',names(ATOM_IDS(i))
	enddo
	print*,size(POSITIONS,dim=2),' atoms written'
end subroutine

!> Выводит информацию об атомах в формате .xyz. Информация выводится в порядке следования имен атомов.
!> \param[in] POSITIONS
!> \param[in] ATOM_IDS
!> \param[in] LATTICE_VECTORS
subroutine write_xyz_by_id(out_id,m,names)
	integer:: out_id
	integer:: i,j
	real:: m(:) !< Массы атомов
	character(len=32):: names(:) !< Названия атомов
	write(out_id,*) size(POSITIONS,dim=2)
	write(out_id,'(A,9f16.6,A)') 'Lattice="',LATTICE_VECTORS,' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do j=1,size(names)
		do i=1,size(POSITIONS,dim=2)
			if (ATOM_IDS(i)==j) then
				write(out_id,'(7f17.7,A,A)') POSITIONS(:,i),0.,0.,0.,m(ATOM_IDS(i)),'    ',names(ATOM_IDS(i))
			endif
		enddo
	enddo
	print*,size(POSITIONS,dim=2),' atoms written'
end subroutine

!> Объединеняет .xyz файлы
subroutine concatenate_xyz_files(files_ids,out_id)
	integer:: files_ids(:)
	integer:: out_id
	integer:: i,j,n(size(files_ids))
	character(len=128):: str1,str2,pl,pl_i
	character(len=512):: line
	real:: lv(9),lv_i(9)
	lv = 0
	do i=1,size(files_ids)
		read(files_ids(i),*) n(i)
		read(files_ids(i),'(A9,9f16.6,A14,A)') str1,lv_i,str2,pl_i
		do j=1,9
			if (lv(j)<lv_i(j)) lv(j) = lv_i(j)
		enddo
		if (i>1) then
			if (pl/=trim(pl_i)) print*, 'properties do not match'
		else
			pl = trim(pl_i)
		endif
	enddo
	write(out_id,*) sum(n)
	write(out_id,'(A,9f16.6,A,A)') 'Lattice="',lv,' " Properties=',pl
	do i=1,size(files_ids)
		do j=1,n(i)
			read(files_ids(i),'(A)') line
			write(out_id,*) trim(line)
		enddo
	enddo
	print*,sum(n),' atoms written'
end subroutine

!> Увеличивает массивы координат атомов и их типов.
!> \param[in,out] POSITIONS
!> \param[in,out] ATOM_IDS
subroutine allocate_more_X(N)
	real,allocatable:: Xt(:,:)
	integer,allocatable:: atomidt(:)
	integer:: N !< Новое количество атомов
	if(N>=size(POSITIONS,dim=2)) then
		if(allocated(POSITIONS)) then
			allocate(Xt(3,size(POSITIONS,dim=2)))
			allocate(atomidt(size(POSITIONS,dim=2)))
			Xt = POSITIONS
			atomidt = ATOM_IDS
			deallocate(POSITIONS)
			deallocate(ATOM_IDS)
			allocate(POSITIONS(3,N))
			allocate(ATOM_IDS(N))
			POSITIONS(:,:size(Xt,dim=2)) = Xt
			ATOM_IDS(:size(Xt,dim=2)) = atomidt
			deallocate(Xt)
			deallocate(atomidt)
		else
			allocate(POSITIONS(3,N))
			allocate(atomidt(N))
		endif
	else
		print*,'N<size(POSITIONS,dim=2)=',size(POSITIONS,dim=2)
		stop
	endif
	print*,size(POSITIONS,dim=2)
end subroutine 

!> Возвращает путь вывода введенный при запуске программы.
function output_path()
	integer:: i
	character(len=128):: arg
	character(len=128):: output_path !< Путь вывода.
	i = 0
	do while (i<=command_argument_count())
		i = i+1
		call get_command_argument(i,arg)
		select case (arg)
		case('-o')
			i = i+1
			call get_command_argument(i,output_path);
		end select
	enddo
	print*,output_path
	return
end function

end module crystal_lattice