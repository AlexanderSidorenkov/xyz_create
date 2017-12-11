module crystal_lattice
implicit none

public
integer:: out_id=2017,t0=101
real,allocatable:: X(:,:)
integer,allocatable:: atomid(:)
real:: lattice_vectors(3,3),smallvalue=10.**(-5)
real,parameter:: e(3,3)=reshape((/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/))
real,parameter:: pi=3.14159265358979323846
contains

subroutine replicate(t)
	integer:: t(3),k,i,j,n
	n = size(X,dim=2)
	call allocate_more_X(n*t(1)*t(2)*t(3))
	do k=1,3
		do i=1,t(k)-1
			do j=1,n
				X(:,n*i+j) = X(:,j)+i*lattice_vectors(:,k)
				atomid(n*i+j) = atomid(j)
			enddo
		enddo
		n = n*t(k)
	enddo
	do k=1,3
		lattice_vectors(:,k) = lattice_vectors(:,k)*t
	enddo
end subroutine

subroutine shift(vec)
	real:: vec(3)
	integer:: i
	do i=1,size(X,dim=2)
		X(:,i) = X(:,i)+vec
	enddo
end subroutine

subroutine stretch(vec)
	real:: vec(3)
	integer:: i,k
	do i=1,size(X,dim=2)
		X(:,i) = X(:,i)*vec
	enddo
	do k=1,3
		lattice_vectors(:,k) = lattice_vectors(:,k)*vec
	enddo
end subroutine

subroutine turn(mat)
	real:: mat(3,3),vec(3)
	integer:: i,k,kk
	do i=1,size(X,dim=2)
		do k=1,3
			vec(k) = sum(mat(:,k)*X(:,i))
		enddo
		X(:,i) = vec
	enddo
	do k=1,3
		do kk=1,3
			vec(kk) = sum(mat(:,kk)*lattice_vectors(:,k))
		enddo
		lattice_vectors(:,k) = vec
	enddo
end subroutine

subroutine create_turn_cut(x1,ids,lv,mat,directions,add_j)
	real:: x1(:,:),lv(:,:),directions(3,3),mat(3,3)
	integer:: ids(:)
	integer:: t(3),j,add_j
	t = (/t0,t0,t0/)
	j = (t0*t0*int(real(t0)/2)+t0*int(real(t0)/2)+int(real(t0)/2))*size(x1,dim=2)+1
	call create(x1,ids,lv)
	call turn(mat)
	call replicate(t)
	call cut_cell(directions,j+add_j)
end subroutine
	
subroutine cut()
	real,allocatable:: Xt(:,:)
	integer,allocatable:: atomidt(:)
	real:: dps(3),surf_norm(3,3),dlvs(3)
	integer:: i,k,n,kk
	logical:: in_cell(3)
	allocate(Xt(3,size(X,dim=2)))
	allocate(atomidt(size(X,dim=2)))
	n=0
	do k=1,3
		do kk=1,3
			surf_norm(kk,k) = lattice_vectors(mod(kk,3)+1,mod(k,3)+1)*lattice_vectors(mod(kk+1,3)+1,mod(k+1,3)+1)-&
			lattice_vectors(mod(kk+1,3)+1,mod(k,3)+1)*lattice_vectors(mod(kk,3)+1,mod(k+1,3)+1)
		enddo
		dlvs(k) = sum(surf_norm(:,k)*lattice_vectors(:,k))/sqrt(sum(surf_norm(:,k)**2))
	enddo
	do i=1,size(X,dim=2)
		do k=1,3
			dps(k) = sum(surf_norm(:,k)*X(:,i))/sqrt(sum(surf_norm(:,k)**2))
			in_cell(k) = dps(k)>-smallvalue .and. dps(k)<dlvs(k)-smallvalue
		enddo
		if(in_cell(1) .and. in_cell(2) .and. in_cell(3)) then
			n = n+1
			Xt(:,n) = X(:,i)
			atomidt(n) = atomid(i)
		endif
	enddo
	print*,n
	deallocate(X)
	deallocate(atomid)
	allocate(X(3,n))
	allocate(atomid(n))
	X(:,:) = Xt(:,:n)
	atomid(:) = atomidt(:n)
	deallocate(Xt)
	deallocate(atomidt)
end subroutine
	
subroutine cut_cell(directions,j)
	real:: directions(3,3),dpl2,dpp
	integer:: i,j,k
	call shift(-X(:,j))
	lattice_vectors = 1./smallvalue
	do i=1,size(X,dim=2)
		do k=1,3
			dpl2 = ( sum(directions(:,k)**2)*sum(X(:,i)**2) - sum(directions(:,k)*X(:,i))**2 )/(sum(directions(:,k)**2))
			if(dpl2<0. .or. sqrt(dpl2)<smallvalue) then
				dpp = sqrt(sum(X(:,i)**2))
				if(dpp<sqrt(sum(lattice_vectors(:,k)**2)) .and. dpp>smallvalue .and. atomid(i)==atomid(j)) lattice_vectors(:,k) = X(:,i)
			endif
		enddo
	enddo
	do k=1,3
		if(lattice_vectors(k,k)<0.) lattice_vectors(:,k) =  -lattice_vectors(:,k)
	enddo
	call cut()
end subroutine

subroutine create(x1,ids,lv)
	real:: x1(:,:),lv(:,:)
	integer:: ids(:)
	call allocate_more_X(size(x1,dim=2))
	X = x1
	lattice_vectors = lv
	atomid = ids
end subroutine

subroutine reset()
	lattice_vectors = 0.
	deallocate(X,atomid)
	allocate(X(3,0),atomid(0))
	print*,'reset ',size(X,dim=2)
end subroutine

subroutine write_xyz(out_id,m,str)
	integer:: out_id,i
	real:: m
	character(len=32):: str(:)
	write(out_id,*) size(X,dim=2)
	write(out_id,'(A,9f16.6,A)') 'Lattice="',lattice_vectors,' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do i=1,size(X,dim=2)
		write(out_id,'(7f17.7,A,A)') X(:,i),0.,0.,0.,m,'    ',str(atomid(i))
	enddo
	print*,size(X,dim=2),' atoms written'
end subroutine

subroutine allocate_more_X(N)
	real,allocatable:: Xt(:,:)
	integer,allocatable:: atomidt(:)
	integer:: N
	if(N>=size(X,dim=2)) then
		if(allocated(X)) then
			allocate(Xt(3,size(X,dim=2)))
			allocate(atomidt(size(X,dim=2)))
			Xt = X
			atomidt = atomid
			deallocate(X)
			deallocate(atomid)
			allocate(X(3,N))
			allocate(atomid(N))
			X(:,:size(Xt,dim=2)) = Xt
			atomid(:size(Xt,dim=2)) = atomidt
			deallocate(Xt)
			deallocate(atomidt)
		else
			allocate(X(3,N))
			allocate(atomidt(N))
		endif
	else
		print*,'N<size(X,dim=2)=',size(X,dim=2)
		stop
	endif
	print*,size(X,dim=2)
end subroutine 
	
end module crystal_lattice

!simple(3,1)=reshape((/ 0.0, 0.0, 0.0 /),shape=(/3,1/)),&
!volume_centered(3,2)=reshape((/ 0.0, 0.0, 0.0, 0.5, 0.5, 0.5 /),shape=(/3,2/)),&
!face_centered(3,4)=reshape((/ 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5 /),shape=(/3,4/)),&
!cubic_lv(3,3)=reshape((/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /),shape=(/3,3/))
