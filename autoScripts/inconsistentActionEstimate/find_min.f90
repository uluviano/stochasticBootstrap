program fit
   implicit none
   integer ND
   real*8, allocatable :: xx(:,:),act(:),dist(:)
   integer, allocatable :: s_order(:)
   integer nOp,i,j
  

   ND=0
   do 
      read(5,*,end=111)
      ND=ND+1
   enddo
111 continue
   rewind(5)

   nOp=nn+1
   allocate(xx(nOp,ND),act(ND),dist(ND),s_order(ND))
   do i=1,ND
      read(5,*)xx(:,i),act(i)
   enddo

   do i=1,ND
      do j=1,ND
         dist(j)=sum((xx(:,i)-xx(:,j))**2)
      enddo
      dist(i)=1.d10
      call HPSORT(ND,dist,s_order)
      do j=1,ND
         ! find the closest neighbour with lower action. Write its order (j)
         if(act(s_order(j))<act(i))then
            write(6,'(i10,100f12.7)')j,xx(:,i),sqrt(dist(j))
            goto 112
         endif
      enddo
            write(6,'(i10,100f12.7)')ND,xx(:,i),sqrt(dist(j))
112   continue
   enddo
     

contains
SUBROUTINE HPSORT(N,RA,s_order)
  implicit none
  integer N,s_order(N)
  real*8 RA(N)

  integer L,IR,I,J,sord
  real*8 RRA
  do i=1,n
     s_order(i)=i
  enddo
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
    sord=s_order(L)
  else
    RRA=RA(IR)
    sord=s_order(IR)
    RA(IR)=RA(1)
    s_order(ir)=s_order(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      s_order(1)=sord
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
    if(J < IR)then
      if(RA(J) < RA(J+1))  J=J+1
    end if
    if(RRA < RA(J))then
      RA(I)=RA(J)
      s_order(i)=s_order(j)
      I=J 
      J=J+J
    else
      J=IR+1
    end if
    goto 20
  end if
  RA(I)=RRA
  s_order(I)=sord
  goto 10
END subroutine HPSORT
!===========================================================
!===========================================================
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
end program fit
