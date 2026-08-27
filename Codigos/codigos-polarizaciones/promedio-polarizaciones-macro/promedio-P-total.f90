integer pasos,n
real(8) px,py,pz,pxt,pyt,pzt, pola
character bla


open(1, file="entrada.dat")
open(2, file="pola-total-prom")

pasos=51


pxt=0
pyt=0
pzt=0

read(1,*) bla
read(1,*) bla


do i=1, pasos


  read(1,*) n,px,py,pz
  pxt=pxt+px
  pyt=pyt+py
  pzt=pzt+pz


enddo

pxt=pxt/pasos
pyt=pyt/pasos
pzt=pzt/pasos


pola=(pxt**2+pyt**2+pzt**2)**0.5


write(2,100) pxt, pyt, pzt, pola 


100 format (F6.2, "  ", F6.2, "  ", F6.2, "  ",F6.2)



end
