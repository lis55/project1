#! /bin/bash -f
#
# Scriptfile
#
# rho    = density
# lmax   = number of steps
# nequil = number of equilibration steps
# dr     = maximum displacement
# npart  = number of particles

for rho in 0.1 0.2; do

cat > input  <<endofdata
  ibeg  nequil   lmax  nsamp  
     0  10000    10000      1        
    dr
  0.09 
ndispl
    50  
 npart    temp    rho
   500     2.0  ${rho} 
 epsilon sigma  mass  cutoff
 1.0     1.0    1.0   3.5
endofdata
time  ../Source/mc_nvt
rm input
done
exit
