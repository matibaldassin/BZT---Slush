
!parameter(vac=12)
integer, parameter :: centros=2744

integer celda, pasos, n !, vacancias(vac), l, centro(vac)
real(8) px(centros), py(centros), pz(centros), pox, poy, poz, p(centros), pol
real(8) pvacx, pvacy, pvacz
character bla

pasos = 51   ! cantidad de pasos a leer     

open(1,file="entrada.in")
open(2,file="salida.out")
!open(3,file="pol-vacancias.out")
!open(4, file="centros.in")

!do i=1, vac

!   read(4,*) l, centro(i)

!   n = (l-1)/2

!   vacancias(i) = n
   

!enddo

   pvacx=0.0
   pvacy=0.0
   pvacz=0.0
 

   do j=1, centros

   px(j) = 0.0
   py(j) = 0.0
   pz(j) = 0.0
   
   p(j) = 0.0
   
   enddo

   pol = 0.0

do i=1, pasos

read(1,*) bla

   do j=1, centros

   read(1,*) n, pox, poy, poz

   px(j) = px(j) + pox
   py(j) = py(j) + poy
   pz(j) = pz(j) + poz
   
   pol = (pox**2 + poy**2 + poz**2)**0.5
   p(j) = p(j) + pol
   
   enddo
  

enddo



do i=1, centros

px(i) = px(i)/pasos 
py(i) = py(i)/pasos
pz(i) = pz(i)/pasos
pol = (px(i)**2 + py(i)**2 + pz(i)**2)**0.5

write(2,2000) i, px(i), py(i), pz(i), pol

!if (i.lt.vac/2+1) then 
!                     write(3,1000) i, centro(i), px(i)/pasos, py(i)/pasos, pz(i)/pasos
!                     pvacx = pvacx + px(i)/pasos
!                     pvacy = pvacy + py(i)/pasos
!                     pvacz = pvacz + pz(i)/pasos
!                  endif

!     do j=vac/2+1, vac
!            if (i.eq.(vacancias(j))) then 
!                                         write(3,1000) i, centro(j), px(i)/pasos, py(i)/pasos, pz(i)/pasos
!                                         pvacx = pvacx + px(i)/pasos
!                                         pvacy = pvacy + py(i)/pasos
!                                         pvacz = pvacz + pz(i)/pasos
!                                     endif 
     !enddo

enddo

!write(3,*) 
!write(3,3000) pvacx/vac, pvacy/vac, pvacz/vac
!write(3,4000) pvacx/centros, pvacy/centros, pvacz/centros


1000 format(I5, 2x, I10, 2x, 3(F7.2, 2x))
2000 format(I5,  2x, 4(F7.2, 2x))
3000 format('promedio vacancias: ',3(F7.2, 2x))
4000 format('promedio vacancias sobre total del sistema: ',3(F7.2, 2x))

end
