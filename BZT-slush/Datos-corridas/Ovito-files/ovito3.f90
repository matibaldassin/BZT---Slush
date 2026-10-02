 program ovito
   implicit none

 ! serie con mpol_xx.dat   
  integer:: nat,ncel 

 ncel = 21952 ; nat = ncel*10

 write(*,*) ' voy a gen ovito'
 call gen_ovito(ncel,nat)


    stop
    end

!*************************************************************************
!
!*************************************************************************


subroutine gen_ovito(ncel,nat)
 implicit none
 integer::ncel,nat

 integer::i,ii,jj,k,ncen(ncel),nptyp(ncel)
 integer :: pos_inicio, pos_final
 real(8)::ap(3,3),xxx(nat),yyy(nat),zzz(nat)
 real(8)::px(ncel),py(ncel),pz(ncel),pmod(ncel)
 character(10)::atm(nat)
 character(130)::texto, nueva_linea, parte1, parte2, nuevo_lattice

  write(*,*) 'voy a config' 

  call read_config(nat,ap,xxx,yyy,zzz,atm)
  write(*,*) 'Lei config' 
 !    ---  find cell centers -------
     i = 0
    do k = 1, nat
      if (atm(k).eq.'Zr   '.or. atm(k).eq.'Ti   ') then
       i = i + 1
       ncen(i) = k
      endif
    enddo
!  armo lista ordenada
   do i = 1, ncel
    k = ncen(i)
    if(atm(k)== 'Zr    ') then
      nptyp(i) = 2
    else
       nptyp(i) = 1
    endif
   enddo

  write(*,*) 'ordeno' 

   open (11,file='P-local')
!     read(11,'(a80)') linea
     do i = 1, ncel
       read(11,*) ii,px(i),py(i),pz(i)
       pmod(i) = sqrt(px(i)**2 + py(i)**2 + pz(i)**2)
     enddo
     close(11)
  write(*,*) 'Lei mpola' 

! modifico texto
     texto=  'Lattice="112.00 0.0 0.0 0.0  112.00 0.0 0.0 0.0  112.93" Properties=species:I:1:pos:R:3:dipoles:R:3:dipoles_mag:R:1'
    pos_inicio = index(texto, 'Lattice="') + len('Lattice="') - 1
    pos_final = index(texto(pos_inicio+1:), '"') + pos_inicio - 1

    ! Extraer las partes antes y después de "Lattice"
    parte1 = texto(1:pos_inicio-1)
    parte2 = texto(pos_final+1:)

    print *, parte1
   print *, parte2
    ! Convertir los valores reales a texto con formato
    nuevo_lattice = ""

    i = 0
    do jj = 1, 3
      do ii = 1, 3
        i = i + 1
        if (i > 1) nuevo_lattice = trim(nuevo_lattice) // " "  ! Agregar espacio
        if (ii==jj) then 
         write(nuevo_lattice(len_trim(nuevo_lattice)+1:), '(F7.2)') ap(ii,jj)
        else 
         write(nuevo_lattice(len_trim(nuevo_lattice)+1:), '(F5.1)') ap(ii,jj)
        endif
!         print *, trim(nuevo_lattice) 
      end do
    enddo

!    print *, trim(parte1)
!    print *, trim(nuevo_lattice)
!    print *, trim(parte2)
    
    ! Reconstruir la línea con las partes correctamente ajustadas
    nueva_linea = trim(parte1) // '"' // trim(nuevo_lattice) // trim(parte2)
!      print *, 
     open(9,file='ovito.xyz')
     write(9,*) ncel
!     write(9,*) trim(adjustl(nuevo_linea)) 
     write(9,'(a130)') nueva_linea 

     do i = 1, ncel
       k = ncen(i)
       write(9,'(i4,7f10.4)') nptyp(i),xxx(k),yyy(k),zzz(k),px(i),py(i),pz(i),pmod(i)
     enddo  ! k

     close(9)

  return
  end

!*************************************************************************
!
!*************************************************************************
 subroutine read_config(nat,ap,x,y,z,atm)
 implicit none
 integer:: nat
 real(8)::ap(3,3),x(nat),y(nat),z(nat)
 character(10)::atm(nat)

 integer:: lconf,i,j,k
 real(8)::ux,uy,uz,vol,rcell(3,3)
 character(80)::linea

   open(1,file='CONFIG',form='formatted',status='old')

   read(1,'(a80)') linea
   read(1,*) lconf

!   ---------------------------------------
!   read CONFIG
!
   do j = 1, 3
     read(1,'(3f20.10)') (ap(i,j),i=1,3)
! pongo en 0 elementos 
     do i = 1, 3
       if(i/=j) ap(i,j) = 0.0d0
     enddo
     write(*,'(3f20.10)') (ap(i,j),i=1,3)
    enddo

   do k = 1, nat
       read(1,'(a10)') atm(k)
        read(1,'(3g20.10)') x(k),y(k),z(k)
       if(lconf.gt.0) then
        read(1,'(a80)') linea
        read(1,'(a80)') linea
       endif
     enddo

  close(1)

! llevo a celdas entre [0,1)

 call invert2(ap,rcell,vol)

 do k = 1, nat
   
    ux=(rcell(1,1)*x(k)+rcell(1,2)*y(k)+rcell(1,3)*z(k))
    uy=(rcell(2,1)*x(k)+rcell(2,2)*y(k)+rcell(2,3)*z(k))
    uz=(rcell(3,1)*x(k)+rcell(3,2)*y(k)+rcell(3,3)*z(k))

    if (ux < 0.0) ux = ux + 1.0
    if (uy < 0.0) uy = uy + 1.0
    if (uz < 0.0) uz = uz + 1.0

!c -----   calculo nuevas posiciones  ----
    x(k)=(ap(1,1)*ux+ap(1,2)*uy+ap(1,3)*uz)
    y(k)=(ap(2,1)*ux+ap(2,2)*uy+ap(2,3)*uz)
    z(k)=(ap(3,1)*ux+ap(3,2)*uy+ap(3,3)*uz)

 enddo

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


! subroutine modifica_texto(texto,
