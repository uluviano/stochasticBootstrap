function relu (x)
real*8 relu, x
relu = (x + abs(x))/2
end

function RHOf (x,y)
real*8 RHOf, x,y
complex*16 z
z=CMPLX(x,y)
RHOf = ABS(z/(1 + SQRT(1-z))**2)
end

program MC
   implicit none
   integer NT,nZFile,nZ,nSpin, nPoints,seedInput
   integer, allocatable :: nZList(:)
   integer nops
   integer externalInOPE
   integer, allocatable :: spinh(:)
   real*8, allocatable :: op_bounds(:,:)
   real*8, allocatable :: DDini(:),DDn(:),opstep(:)
   real*8, allocatable :: zz(:,:),zzt(:,:)
   real*8, allocatable :: gg1(:,:,:,:),gg2(:,:,:,:)
   real*8  Dstar,Dext,Temp
   real*8, allocatable :: lambda(:),eps(:),z1(:),z2(:),err(:),rho(:),rhon(:),DD(:),DDi(:,:),step(:)
   real*8, allocatable :: gt1(:),gt2(:),b(:),c(:),d(:)
   real*8, allocatable :: MM(:,:),AA(:,:),res(:)
   integer, allocatable :: seed(:)
   character*200   filepar



   integer i,j,jj,k,ii,it,acc,positiveOPEs,areOPEspos,invSampleDensity,nGoodZ,validPoint
   real*8 g1,g2,DDmin,DDmax,SS,dx,actt,actn,MCstep,r,coleps,wall,adaptStep,lambda0,zPointDensity
   !functions
   real*8 relu, RHOf
   adaptStep=1.01d0

   read(5,*)filepar
   read(5,*)nZ
   read(5,*)lambda0
   if(lambda0<0.344d0)stop "lambda_0 is too small"
   read(5,*)Dstar
   read(5,*)NT
   read(5,*)positiveOPEs
   read(5,*)externalInOPE
   read(5,*)Temp
   read(5,*)MCstep
   read(5,*)coleps
   read(5,*)wall
   read(5,*)invSampleDensity
   read(5,*)seedInput
   read(5,*)nops !should not include the external operator
   if(nops>nZ+2)stop "Need more z Points in order to have a well defined quadratic form."
   allocate(spinh(nops),op_bounds(2,0:nops),opstep(0:nops),DDini(0:nops))
   read(5,*)i,op_bounds(:,0),DDini(0),opstep(0)
   if(i.ne.-1)stop "the first operator should be external"
   if(DDini(0)<op_bounds(1,0).or. DDini(0)>op_bounds(2,0)) stop "initial condition outside boudaries"
   do i=1,nops
      read(5,*)spinh(i),op_bounds(:,i),DDini(i),opstep(i)
      spinh(i)=spinh(i)/2+1 !
      if(DDini(i)<op_bounds(1,i).or. DDini(i)>op_bounds(2,i))then
         write(6,*)spinh(i),op_bounds(:,i),DDini(i)
         stop "initial condition outside boudaries"
      endif
   enddo
   !now read the file with the values of the CBs
   open(50,file=filepar,status="old")
   read(50,*)nSpin,nZFile,nPoints
   do i=1,nops
      if(spinh(i)>nSpin)stop "spin spectrum not present"
   enddo


   allocate(zz(2,nZ),zzt(2,nZFile),gg1(4,nPoints,nZ,nSpin),gg2(4,nPoints,nZ,nSpin),DDi(nPoints,nSpin),step(nSpin))
   allocate(err(nZ),nZList(nZ))
   read(50,*)zzt(:,:)

   !Select actual zPoint Sample based on lambda0
   allocate(lambda(nZFile))
   
   nGoodZ=-1
   do j=1,nZFile
      lambda(j)= RHOf(zzt(1,j),zzt(2,j)) + RHOf(1-zzt(1,j),-zzt(2,j))
      if(lambda(j)>lambda0)then
         nGoodZ=j-1
         exit
      endif
   enddo
   

   !initialize random Gen
   call random_seed(size=k)
   allocate(seed(1:k))
   seed(:) = seedInput
   call random_seed(put=seed)

   !create nZList
   if(nGoodz.eq.-1)stop "Your \lambda_0 is too large or invalid."
   if(nGoodz<nZ +1)stop "Not Enough Points at required Lambda. Try a larger lambda or a smaller nZ"
   zPointDensity=nGoodz/dble(nZ)
   do i=1,nZ
      do
         call random_number(r)
         nZList(i)=floor((i-1)*zPointDensity + r*zPointDensity) + 1
         validPoint=0
         if(i.eq.1) then
            validPoint=validPoint+1
         else
            if(nZList(i).ne.nZList(i-1))validPoint=validPoint+1
         endif
         if(validPoint.eq.1)exit
      enddo
   enddo

   do j=1,nSpin
      read(50,*)DDmin,DDmax
      step(j)=(DDmax-DDmin)/dble(nPoints-1)
      do i=1,nPoints
         DDi(i,j)=DDmin+dble(i-1)*step(j)
      enddo
   enddo

   allocate(gt1(nPoints),gt2(nPoints),b(nPoints),c(nPoints),d(nPoints))
   do i=1,nSpin
      jj=1
      do j=1,nZFile
         read(50,*)gt1(:),gt2(:)
         if(jj>nZ)cycle
         if(j.eq.nZList(jj))then
            call spline (DDi(:,i),gt1,b,c,d,nPoints)
            gg1(1,:,jj,i)=gt1(:)
            gg1(2,:,jj,i)=b(:)
            gg1(3,:,jj,i)=c(:)
            gg1(4,:,jj,i)=d(:)
            call spline (DDi(:,i),gt2,b,c,d,nPoints)
            gg2(1,:,jj,i)=gt2(:)
            gg2(2,:,jj,i)=b(:)
            gg2(3,:,jj,i)=c(:)
            gg2(4,:,jj,i)=d(:)
            jj=jj+1
         endif
      enddo
   enddo
   close(50)
      
   allocate(eps(nZ),z1(nZ),z2(nZ))

   !estimate the error
   ii=int((Dstar-DDi(1,1))/step(1))+1
   dx = Dstar - DDi(ii,1)
   jj=1
   do j=1,nZFile
      if(jj>nZ)cycle
      if(j.eq.nZList(jj))then
         z1(jj)=((1.d0-zzt(1,j))**2+zzt(2,j)**2)**0.5d0
         z2(jj)=(zzt(1,j)**2+zzt(2,j)**2)**0.5d0
         !and fill good zz array
         zz(:,jj)=zzt(:,j)
         eps(jj)=z1(jj)-z2(jj)
         g1= gg1(1,ii,jj,1) + dx*(gg1(2,ii,jj,1) + dx*(gg1(3,ii,jj,1) + dx*gg1(4,ii,jj,1)))
         g2= gg2(1,ii,jj,1) + dx*(gg2(2,ii,jj,1) + dx*(gg2(3,ii,jj,1) + dx*gg2(4,ii,jj,1)))
         err(jj)=(z1(jj)*g1-z2(jj)*g2)**2
         jj=jj+1
      endif
   enddo


   allocate(DD(0:nops),DDn(0:nops)) ! this includes the external operator (DD(1))
   allocate(AA(nops,nops),rho(nops),rhon(nops))
   allocate(res(nZ),MM(nZ,nops))



   !MCstep = 0.01d0

   DD(:)=DDini(:)
   if(externalInOPE.eq.1)DD(1)=DD(0)
      call compute_SS(actt,DD,rho)
   do it=1,NT
      do i=0,nops
        call random_number(r)
        DDn(i)=DD(i)+opstep(i)*MCstep*(r-0.5)
      enddo

      Dext=DDn(0)
      if(externalInOPE.eq.1)DDn(1)=DDn(0)
      call compute_SS(actn,DDn,rhon)
      acc=0
      if(actn<actt)then
         acc=1
      else
         call random_number(r)
         if(r<exp((actt-actn)/Temp))acc=1
      endif
      if(acc.eq.1)then
         actt=actn
         DD(:)=DDn(:)
         rho(:)=rhon(:)
            MCstep=MCstep*adaptStep
      else
            MCstep=MCstep/adaptStep
      endif

      if(MOD(it,invSampleDensity).eq.0)write(6,'(t1,2000f29.16)')DD(:),actt,rho(:)
   enddo

contains
!
   subroutine compute_SS(SS,DD,rho)
   implicit none
   real*8 DD(0:nops),AAb(nops,nops),rho(nops),SS,oob_op,relu

   integer newnegs,ii,i,j,k,itt,jtt,sph, nopst,negOPEs(nops),redOPs(nops),rank,lwork,info
   real*8 rcond,g1,g2,dx,Dext
   real*8, allocatable :: work(:), MMb(:,:),MMren(:,:),epsRenb(:),epsRen(:),jpvt(:),AAt(:,:)
   !Real (Kind=dp) :: rcond

    lwork = 3*nOps + 64*(nOps+1)
    allocate(jpvt(nOps),work(lwork))


   
   AA(:,:)=0
!  stop

   Dext=DD(0)
   if(externalInOPE.eq.1)DD(1)=DD(0)
      oob_op=0
      do i=0,nops
         oob_op=oob_op+ (relu(-DD(i)+op_bounds(1,i)) + relu(-op_bounds(2,i) +DD(i)))**2
      enddo
   do j=1,nZ
      z1(j)=((1.d0-zz(1,j))**2+zz(2,j)**2)**Dext
      z2(j)=(zz(1,j)**2+zz(2,j)**2)**Dext
      eps(j)=z1(j)-z2(j)
   enddo

   !estimate the error
   ii=int((Dstar-DDi(1,1))/step(1))+1
   dx = Dstar - DDi(ii,1)
   do j=1,nZ
      g1= gg1(1,ii,j,1) + dx*(gg1(2,ii,j,1) + dx*(gg1(3,ii,j,1) + dx*gg1(4,ii,j,1)))
      g2= gg2(1,ii,j,1) + dx*(gg2(2,ii,j,1) + dx*(gg2(3,ii,j,1) + dx*gg2(4,ii,j,1)))
      err(j)=(z1(j)*g1-z2(j)*g2)**2
   enddo

   nopst=nOps
   redOPs(:)=0
   do i=1,nops
      sph=spinh(i)
      !check collision
      do j=i+1,nops
         if(sph.eq.spinh(j) .and. abs(DD(i)-DD(j))<coleps)then
            !write(6,*)i,j
            redOPs(i)=1
            nopst=nopst-1
            exit
         else
            redOPs(i)=0
         endif
      enddo

      if(redOPs(i).eq.0)then
         ii=int((DD(i)-DDi(1,sph))/step(sph))+1
         dx = DD(i) - DDi(ii,sph)
         do j=1,nZ
            g1= gg1(1,ii,j,sph) + dx*(gg1(2,ii,j,sph) + dx*(gg1(3,ii,j,sph) + dx*gg1(4,ii,j,sph)))
            g2= gg2(1,ii,j,sph) + dx*(gg2(2,ii,j,sph) + dx*(gg2(3,ii,j,sph) + dx*gg2(4,ii,j,sph)))
            MM(j,i)=z1(j)*g1-z2(j)*g2
         enddo
      else
         do j=1,nZ
            MM(j,i)=0
         enddo
      endif

   enddo

   !add first time through

   allocate(epsRenb(nZ),epsRen(nZ),MMren(nZ,nops))
   allocate(MMb(nZ,nops))
   !initialize internal indices
   do k=1,nZ
        epsRen(k)=-eps(k)/sqrt(err(k))
        epsRenb(k)=epsRen(k)
        do i=1,nOps
            MMren(k,i) = MM(k,i)/sqrt(err(k))
            MMb(k,i) = MMren(k,i)
        enddo
   enddo
    jpvt(1:nOps) = 0
      rcond = 0.0000000001

Call dgelsy(nZ, nOps, 1, MMren, nZ, epsRen, nZ, jpvt, rcond, rank, work, lwork, info)

rho(:) = epsRen(1:nOps)
deallocate(MMren,jpvt,work)



   !positive OPE enforcing
   if(positiveOPEs.eq.1)then
      !start checking
      do i=1,nops
         negOPEs(i)=0
      enddo
      newnegs=0
      do i=1,nops
         if(rho(i)<0)then
             nopst=nopst-1
             newnegs=newnegs+1
             negOPEs(i)=1
         endif
      enddo
      if(newnegs.eq.0)then
         areOPEspos=1
      else
         areOPEspos=0
      endif
      if(nopst<1)then
         areOPEspos=1
         rho(:)=0
      endif


      !initialize temps
      do while(areOPEspos.eq.0)
         allocate(MMren(nZ,nOpst))
         !initialize internal indices
         itt=0
         do i=1,nops
            if(negOPEs(i).eq.0 .and. redOPs(i).eq.0)then
               itt=itt+1
               do j=1,nZ
                     MMren(j,itt)= MMb(j,i)
                     epsRen(j)=epsRenb(j)
               enddo
            endif
         enddo
    lwork = 3*nOpst + 64*(nOpst+1)
    allocate(jpvt(nOpst),work(lwork))

      Call dgelsy(nZ, nOpst, 1, MMren, nZ, epsRen, nZ, jpvt, rcond, rank, work, lwork, info)
         !write(6,*)epsRen(1:nOpst)
      deallocate(MMren,jpvt,work)

            itt=0
         do i=1,nops
            rho(i)=0
            if(negOPEs(i).eq.0 .and. redOPs(i).eq.0)then
               itt=itt+1
               rho(i)=epsRen(itt)
            endif
         enddo

         !check if negative OPEs are present
         newnegs=0
         !nopst=nops
         do i=1,nops
            if(rho(i)<0)then
                nopst=nopst-1
                newnegs=newnegs+1
                negOPEs(i)=1
            endif
         enddo
         if(newnegs.eq.0)then
            areOPEspos=1
         else
            areOPEspos=0
         endif
         !deallocate(AAt)
            if(nopst<1)then
               areOPEspos=1
               rho(:)=0
            endif
      enddo
   endif

   do k=1,nZ
      res(k)=(sum(MM(k,:)*rho(:))+eps(k))/sqrt(err(k))
   enddo
   SS=log(sum(res(:)**2)/dble(nZ))
           !write(6,*)SS,rho(:)
              SS=SS+wall*oob_op

   end subroutine compute_SS

   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
   implicit none
   integer n
   real*8 x(n), y(n), b(n), c(n), d(n)
   integer i, j, gap
   real*8 h

   gap = n-1
   if(n<4)stop 'err 1'
   !
   ! step 1: preparation
   !
   d(1) = x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i = 2, gap
     d(i) = x(i+1) - x(i)
     b(i) = 2.0*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
   end do
   !
   ! step 2: end conditions 
   !
   b(1) = -d(1)
   b(n) = -d(n-1)
   c(1) = 0.0
   c(n) = 0.0
   if(n /= 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
   end if
   !
   ! step 3: forward elimination 
   !
   do i = 2, n
     h = d(i-1)/b(i-1)
     b(i) = b(i) - h*d(i-1)
     c(i) = c(i) - h*c(i-1)
   end do
   !
   ! step 4: back substitution
   !
   c(n) = c(n)/b(n)
   do j = 1, gap
     i = n-j
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
   end do
   !
   ! step 5: compute spline coefficients
   !
   b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
   do i = 1, gap
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
   end do
   c(n) = 3.0*c(n)
   d(n) = d(n-1)
   end subroutine spline
   subroutine invert(a,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A (ALE: c=a in output
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
  a(:,:)=c(:,:) !overwrite the input matrix
  end subroutine invert
end program MC
