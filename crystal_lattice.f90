module crystal_lattice
implicit none

public
integer:: out_id=2017,t0=101
real,allocatable:: POSITIONS(:,:)
integer,allocatable:: ATOM_IDS(:)
real:: LATTICE_VECTORS(3,3),smallvalue=10.**(-5)
real,parameter:: e(3,3)=reshape((/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/))
real,parameter:: pi=3.14159265358979323846
contains

subroutine replicate(t)
	integer:: t(3),k,i,j,n
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

subroutine shift(vec)
	real:: vec(3)
	integer:: i
	do i=1,size(POSITIONS,dim=2)
		POSITIONS(:,i) = POSITIONS(:,i)+vec
	enddo
end subroutine

subroutine stretch(vec)
	real:: vec(3)
	integer:: i,k
	do i=1,size(POSITIONS,dim=2)
		POSITIONS(:,i) = POSITIONS(:,i)*vec
	enddo
	do k=1,3
		LATTICE_VECTORS(:,k) = LATTICE_VECTORS(:,k)*vec
	enddo
end subroutine

subroutine transform(mat)
	real:: mat(3,3),vec(3)
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

subroutine cut_transformed_cell(mat,directions,add_j)
	real:: directions(3,3),mat(3,3)
	integer:: t(3),j,add_j
	t = (/t0,t0,t0/)
	j = (t0*t0*int(real(t0)/2)+t0*int(real(t0)/2)+int(real(t0)/2))*size(x1,dim=2)+1
	call transform(mat)
	call replicate(t)
	call cut_cell(directions,j+add_j)
end subroutine
	
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
	
subroutine cut_cell(directions,j)
	real:: directions(3,3),dpl2,dpp
	integer:: i,j,k
	call shift(-POSITIONS(:,j))
	LATTICE_VECTORS = 1./smallvalue
	do i=1,size(POSITIONS,dim=2)
		do k=1,3
			dpl2 = ( sum(directions(:,k)**2)*sum(POSITIONS(:,i)**2) - sum(directions(:,k)*POSITIONS(:,i))**2 )/(sum(directions(:,k)**2))
			if(dpl2<0. .or. sqrt(dpl2)<smallvalue) then
				dpp = sqrt(sum(POSITIONS(:,i)**2))
				if(dpp<sqrt(sum(LATTICE_VECTORS(:,k)**2)) .and. dpp>smallvalue .and. ATOM_IDS(i)==ATOM_IDS(j)) LATTICE_VECTORS(:,k) = POSITIONS(:,i)
			endif
		enddo
	enddo
	do k=1,3
		if(LATTICE_VECTORS(k,k)<0.) LATTICE_VECTORS(:,k) =  -LATTICE_VECTORS(:,k)
	enddo
	call cut()
end subroutine

subroutine create(x1,ids,lv)
	real:: x1(:,:),lv(:,:)
	integer:: ids(:)
	call allocate_more_X(size(x1,dim=2))
	POSITIONS = x1
	LATTICE_VECTORS = lv
	ATOM_IDS = ids
end subroutine

subroutine reset()
	LATTICE_VECTORS = 0.
	if(allocated(POSITIONS)) deallocate(POSITIONS,ATOM_IDS)
	allocate(POSITIONS(3,0),ATOM_IDS(0))
	print*,'reset ',size(POSITIONS,dim=2)
end subroutine

subroutine write_xyz(out_id,m,str)
	integer:: out_id,i
	real:: m
	character(len=32):: str(:)
	write(out_id,*) size(POSITIONS,dim=2)
	write(out_id,'(A,9f16.6,A)') 'Lattice="',LATTICE_VECTORS,' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do i=1,size(POSITIONS,dim=2)
		write(out_id,'(7f17.7,A,A)') POSITIONS(:,i),0.,0.,0.,m,'    ',str(ATOM_IDS(i))
	enddo
	print*,size(POSITIONS,dim=2),' atoms written'
end subroutine

subroutine allocate_more_X(N)
	real,allocatable:: Xt(:,:)
	integer,allocatable:: atomidt(:)
	integer:: N
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
	
end module crystal_lattice