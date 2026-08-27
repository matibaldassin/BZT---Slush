 program genconf

 implicit none

 integer:: nat,ncel,k,kk,ix,iy,iz,j,i,nc(3)
 real(8):: celx,cely,celz,x,y,z,  &
     ap(3,3), ap2(3,3), xxx(200000), yyy(200000), zzz(200000)

character(80):: linea
character(10):: atm(200000)


  open(1,file='config',form='formatted',status='old')
  open(3,file='datos',form='formatted')
  open(2,file='CONFIG',form='formatted')



  read(1,*) nat
  read(1,*) nc(:)
  write(2,*) 'BULK nx,ny,nz =',nc(:)
  read(1,'(a80)') linea
  write(2,'(a80)') linea

  do j = 1, 3
    read(1,*) (ap(i,j),i=1,3)
  enddo

  do j = 1, 3
   ap2(:,j) = ap(:,j) * nc(j)
  enddo

!  escribo numero de celdas
   write(3,*) nc(1),nc(2),nc(3)

! escribo numero de centros y numero de atomos
  write(3,*) nc(1)*nc(2)*nc(3), nat*nc(1)*nc(2)*nc(3)
  

do j = 1, 3
  write(2,'(3f20.10)') (ap2(i,j),i=1,3)
enddo


do k = 1, nat
  read(1,'(a10)') atm(k)
!  write(2,'(a10,i10)') atm(k), k
!	  write(2,'(a80)') linea
          read(1,*) xxx(k),yyy(k),zzz(k)
 enddo


! c   ---  construyo mi macrocelda  ----
kk = 0
do iz = 1, nc(3)
 do iy = 1, nc(2)
  do ix = 1, nc(1)
   celx = (ix-1)*ap(1,1)+(iy-1)*ap(1,2)+(iz-1)*ap(1,3)
   cely = (ix-1)*ap(2,1)+(iy-1)*ap(2,2)+(iz-1)*ap(2,3)
   celz = (ix-1)*ap(3,1)+(iy-1)*ap(3,2)+(iz-1)*ap(3,3)
    do k = 1, nat
      kk = kk + 1
      x = celx + xxx(k)
      y = cely + yyy(k)
      z = celz + zzz(k)
      write(2,'(a10,i10)') atm(k), kk
      write(2,'(3f20.10)') x, y, z
    enddo
  enddo
 enddo
enddo

stop

end







