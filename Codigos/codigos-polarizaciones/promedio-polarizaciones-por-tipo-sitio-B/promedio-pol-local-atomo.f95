
integer, parameter :: celdasZr = 123, celdasTi = 2621
real(8) pxsum, pysum, pzsum, pxsumZr, pysumZr, pzsumZr, pxsumTi, pysumTi, pzsumTi

open(1, file="entrada.dat")

open(11, file="salida.dat")


pxsum = 0
pysum = 0
pzsum = 0


do i=1, celdasZr

   read(1,*) num, px, py, pz

   pxsum = px + pxsum
   pysum = py + pysum
   pzsum = pz + pzsum
   
enddo





pxsumZr = pxsum/celdasZr
pysumZr = pysum/celdasZr
pzsumZr = pzsum/celdasZr



pxsum = 0
pysum = 0
pzsum = 0


do i=1, celdasTi

   read(1,*) num, px, py, pz

   pxsum = px + pxsum
   pysum = py + pysum
   pzsum = pz + pzsum
   
enddo





pxsumTi = pxsum/celdasTi
pysumTi = pysum/celdasTi
pzsumTi = pzsum/celdasTi

write(11,*) " Polarización esfera Zr                Polarización matriz BaTiO3 "
write(11,"(3(f6.2, 3x), 6x, 3(f6.2, 3x), 6x )") pxsumZr, pysumZr, pzsumZr, pxsumTi, pysumTi, pzsumTi




end
