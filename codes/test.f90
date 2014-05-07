program main

use cosmo_module
use omp_lib

implicit none

integer (kind=4), parameter :: Rdim=1000
integer (kind=4) :: i
real (kind=8), dimension(Rdim) :: radius, density,redshift,growthfactor
real (kind=4), dimension(100) :: k1,pk
!$ print *, "compiled with -fopenmp"


	call NFW(1d15,0d0,Rdim,radius,density) 
	call growthfactorDE(redshift,growthfactor)
!do i=1,Rdim
!	write(*,*) redshift(i),growthfactor(i)
!end do

	call linearPk(0.01,10.0,100,k1,pk)
	do i=1,100
		write(*,*) k1(i),pk(i)
	end do

end
