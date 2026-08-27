   program tabpol
   implicit none

 integer:: i,ii,init,j,k,kp,lconf,nat,ncel,ncount,lvec(29)   
 integer,allocatable::lwei(:),ncen(:)
 real(8)::rcut2,rsq,x,xref,y,yref,z,zref,vol,ap(3,3), rcell(3,3)
 real(8),allocatable::xxx(:),yyy(:),zzz(:)  
 character(80):: linea  
 character(10),allocatable:: atm(:) 

 open(1,file='CONFIG',form='formatted',status='old')
 open(2,file='tabpol.dat',form='formatted')

 read(1,'(a80)') linea
 write(*,'(a80)') linea
 read(1,*) lconf 
 write(*,*) lconf
!c        write(2,'(a80)') linea

  write(*,*) ' Ingrese # de centros, # numero de atomos '
  read (*,*) ncel, nat
  allocate(lwei(nat),ncen(ncel))
  allocate(atm(nat))
  allocate(xxx(nat),yyy(nat),zzz(nat))
  write(2,*) ncel, ' Centros'
   do j = 1, 3
     read(1,'(3f20.10)') (ap(i,j),i=1,3)
     write(*,'(3f20.10)') (ap(i,j),i=1,3)
  enddo

     do k = 1, nat
        read(1,'(a10)') atm(k)
!c        write(*,'(a10)') atm(k) 
        read(1,'(3g20.10)') xxx(k),yyy(k),zzz(k)
!c        write(*,'(3g20.10)') xxx(k),yyy(k),zzz(k)
       if(lconf.gt.0) then
        read(1,'(a80)') linea
        read(1,'(a80)') linea
        write(*,'(a80)') linea
       endif
     enddo

     i = 0
  do k = 1, nat
    if (atm(k).eq.'Mg    '.or.atm(k).eq.'Zr   '.or. atm(k).eq.'Ti   ') then
     i = i + 1
     if(i.gt.ncel) stop ' i muy grande'
     ncen(i) = k
!c    write (2,'(i5)') ncen(i)
   endif
  enddo

  call invert2(ap,rcell,vol)
   ii = 0
   init = 1

  do i = 1, ncel 
    ncount = 0
    k = ncen(i)
    xref = xxx(k)
    yref = yyy(k)
    zref = zzz(k)

    do kp = 1, nat
      if (k.ne.kp) then
        if(atm(kp).eq.'MgA  '.or.atm(kp).eq.'Pb-sh') then
           rcut2 = 19.765
        elseif(atm(kp).eq.'MgA-sh '.or.atm(kp).eq.'Pbs-sh') then
          rcut2 = 19.765
        elseif(atm(kp).eq.'Ca  '.or.atm(kp).eq.'Ca-sh ') then
          rcut2 = 19.765
        elseif(atm(kp).eq.'Sr  '.or.atm(kp).eq.'Sr-sh ') then
          rcut2 = 19.765
        elseif(atm(kp).eq.'Ba  '.or.atm(kp).eq.'Ba-sh') then
          rcut2 = 19.765
        elseif(atm(kp).eq.'Pb  '.or.atm(kp).eq.'Pb-sh') then
         rcut2 = 19.765
        elseif(atm(kp).eq.'O  '.or.atm(kp).eq.'O-sh') then
          rcut2 = 11.0
        elseif(atm(kp).eq.'Os '.or.atm(kp).eq.'Os-sh') then
          rcut2 = 12.0
        elseif(atm(kp).eq.'Mg-sh'.or.atm(kp).eq.'Zr-sh'.or.atm(kp).eq.'Ti-sh') then
          rcut2 = 2.0 
        elseif(atm(kp).eq.'Fe-sh') then
          rcut2 = 2.0 
        else
          rcut2 = 0.0
        endif

       x = xxx(kp) - xref
       y = yyy(kp) - yref
       z = zzz(kp) - zref
       call image(x,y,z,ap,rcell) 
       rsq = x*x + y*y + z*z
       if (rsq.lt.rcut2) then
          ncount = ncount + 1
          lvec(ncount) = kp

       endif

     endif
   enddo  ! kp

   if (ncount.ne.29) write(*,*) 'error vecinos', k, ncount


   write (2,'(i8,i5,3x,a10)') k, ncount,atm(k)
   write (2,'(29i7)') (lvec(j), j=1,29)

  enddo  ! k

!c   ---   weight  ----
   do kp = 1, nat
      if(atm(kp).eq.'Ba '.or.atm(kp).eq.'Ba-sh') then
         lwei(kp) = 8
      elseif(atm(kp).eq.'Pbs '.or.atm(kp).eq.'Pbs-sh') then
         lwei(kp) = 8
      elseif(atm(kp).eq.'Pb  '.or.atm(kp).eq.'Pb-sh ') then
         lwei(kp) = 8
      elseif(atm(kp).eq.'Bi '.or.atm(kp).eq.'Bi-sh') then
         lwei(kp) = 8
      elseif(atm(kp).eq.'Ca '.or.atm(kp).eq.'Ca-sh') then
         lwei(kp) = 8
      elseif(atm(kp).eq.'O  '.or.atm(kp).eq.'O-sh ') then
         lwei(kp) = 2 
      elseif(atm(kp).eq.'Os '.or.atm(kp).eq.'Os-sh ') then
         lwei(kp) = 2
      elseif(atm(kp).eq.'Mg-sh') then
         lwei(kp) = 1
      elseif(atm(kp).eq.'Fe-sh') then
         lwei(kp) = 1
      elseif(atm(kp).eq.'Ti-sh') then
         lwei(kp) = 1
      else
         lwei(kp) = 100
      endif
   enddo

        write (2,'(20i5)') (lwei(j), j=1,nat)

    stop
    end

!************************************************************************
!*                                                                      *
!************************************************************************

    subroutine image (x,y,z, ap, rcell)

      implicit none

      real(8):: ap(3,3), rcell(3,3),ux,uy,uz,x,y,z


      ux=(rcell(1,1)*x+rcell(1,2)*y+rcell(1,3)*z)
      uy=(rcell(2,1)*x+rcell(2,2)*y+rcell(2,3)*z)
      uz=(rcell(3,1)*x+rcell(3,2)*y+rcell(3,3)*z)

        ux = ux - nint(ux)
        uy = uy - nint(uy)
        uz = uz - nint(uz)

!c -----   calculo nuevas posiciones  ----

      x=(ap(1,1)*ux+ap(1,2)*uy+ap(1,3)*uz)
      y=(ap(2,1)*ux+ap(2,2)*uy+ap(2,3)*uz)
      z=(ap(3,1)*ux+ap(3,2)*uy+ap(3,3)*uz)

      return
      end

!************************************************************************
!*                                                                      *
!************************************************************************

      subroutine invert2(a,b,d)
        implicit none

        real(8):: a(3,3), b(3,3), d,r
!c
!c***********************************************************************
!c
!c     subroutine to invert a 3 * 3 matrix using cofactors
!cc
!c***********************************************************************
!c

!c
!c     calculate adjoint matrix
      b(1,1)=a(2,2)*a(3,3)-a(3,2)*a(2,3)
      b(2,1)=a(3,1)*a(2,3)-a(2,1)*a(3,3)
      b(3,1)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
      b(1,2)=a(3,2)*a(1,3)-a(1,2)*a(3,3)
      b(2,2)=a(1,1)*a(3,3)-a(3,1)*a(1,3)
      b(3,2)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
      b(1,3)=a(1,2)*a(2,3)-a(2,2)*a(1,3)
      b(2,3)=a(2,1)*a(1,3)-a(1,1)*a(2,3)
      b(3,3)=a(1,1)*a(2,2)-a(2,1)*a(1,2)
!c
!c     calculate determinant
      d=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
      r=0.d0
      if(abs(d)>0.0) r=1.0/d
!c
!c     complete inverse matrix
      b(1,1)=r*b(1,1)
      b(2,1)=r*b(2,1)
      b(3,1)=r*b(3,1)
      b(1,2)=r*b(1,2)
      b(2,2)=r*b(2,2)
      b(3,2)=r*b(3,2)
      b(1,3)=r*b(1,3)
      b(2,3)=r*b(2,3)
      b(3,3)=r*b(3,3)

      return
      end

!************************************************************************
!*                                                                      *
!************************************************************************

