program main

use cosmo_module
use omp_lib
use mpi
use mcmc
use nr

implicit none

integer (kind=4), parameter :: Rdim=500, kdim=1000
integer (kind=4) :: i
real (kind=8), dimension(Rdim) :: radius, density,redshift,growthfactor
real (kind=8), dimension(Rdim) :: rho_DM_AC,rho_DM,radius2
real (kind=4), dimension(kdim) :: k1,pk
real (kind=8) :: Mstar_center,x,mass,AF
logical :: yesorno

!$ print *, "compiled with -fopenmp"
	idum = 2214524
	call NFW(1d14,0d0,Rdim,radius,rho_DM) 
	call DensityProfileStars(1d14,0d0,0d0,1d0,Rdim,radius2,density,Mstar_center)
	call AdiabaticContraction(radius,rho_DM,1d14,0d0,Rdim,Mstar_center,rho_DM_AC)
	call mcmc_chain(AF)

end

