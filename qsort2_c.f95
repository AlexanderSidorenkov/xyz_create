! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! and makes the same rearrangement in the second integer array

! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
! Made F conformant by Walt Brainerd
! The second integer array added by Alexander Sidorenkov

module qsort2_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A,B)
  real, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A)==size(B)) then
     if(size(A) > 1) then
        call Partition(A, B, iq)
        call QsortC(A(:iq-1), B(:iq-1))
        call QsortC(A(iq:), B(iq:))
     endif
  else
     print*, 'ERROR: sizes of arrays are not equal.'
	 stop
  endif
end subroutine QsortC

subroutine Partition(A, B, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  real :: tempA
  integer :: tempB
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        ! exchange B(i) and B(j)
        tempA = A(i)
        A(i) = A(j)
        A(j) = tempA
        tempB = B(i)
        B(i) = B(j)
        B(j) = tempB
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort2_c_module