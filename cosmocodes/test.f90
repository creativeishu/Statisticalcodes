program main

use cosmo_module
use omp_lib
use mpi
use mcmc
use nr

implicit none

integer (kind=4), parameter :: Rdim=100, kdim=500
integer (kind=4) :: i
real (kind=8), dimension(Rdim) :: radius,rhoNFW, rhoGAS,rhoBCG
real (kind=8), dimension(Rdim) :: rho_DM_AC,rho_DM,radius2
real (kind=4), dimension(kdim) :: k1,pk
real (kind=8), dimension(6) :: params
real (kind=8) :: Mstar_center,x,mass,AF,redshift
logical :: yesorno
real (kind=8), dimension(kdim_linear) :: Pk_nonlinear,k_nonlinear
real (kind=8), dimension(kdim_linear) :: Pk1h_nonlinear,Pk2h_nonlinear


mass = 2d14
redshift = 0d0
call NFW(mass,0d0,Rdim,radius,rhoNFW)
call DensityProfileGas(mass,redshift,0.15d0,Rdim,radius,rhoGAS)
call DensityProfileStars(mass,redshift,1d0,3d0,Rdim,&
				radius,rhoBCG,Mstar_center)
call AdiabaticContraction(radius,rhoNFW,mass,redshift,Rdim,Mstar_center,rho_DM_AC)

!do i=1,Rdim
!	write(*,*) radius(i),rhoNFW(i),rhoGAS(i),rhoBCG(i),rho_DM_AC(i)
!end do

baryons=.True.
AC=.True.
Mcrit=1d13
call pknonlinear(redshift,k_nonlinear,Pk_nonlinear,pk1h_nonlinear,pk2h_nonlinear)
do i=1, kdim_linear
	write(12,*) k_nonlinear(i),pk_nonlinear(i),pk1h_nonlinear(i),pk2h_nonlinear(i)
end do



end

