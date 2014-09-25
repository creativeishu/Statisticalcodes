program main

use cosmo_module
use omp_lib
use mcmc
use nr
use fisher

implicit none

integer (kind=4), parameter :: Rdim=100, kdim=500,ldim=1000,nopd=7,nopb=8
integer (kind=4) :: i,j
real (kind=8), dimension(Rdim) :: radius,rhoNFW, rhoGAS,rhoBCG
real (kind=8), dimension(Rdim) :: rho_DM_AC,rho_DM,radius2
real (kind=4), dimension(kdim) :: k1,pk
real (kind=8), dimension(6) :: params
real (kind=8) :: Mstar_center,x,mass,AF,redshift,norm,z1,z2,fgas
logical :: yesorno
real (kind=8), dimension(kdim_linear) :: Pk_nonlinear,k_nonlinear
real (kind=8), dimension(kdim_linear) :: Pk1h_nonlinear,Pk2h_nonlinear
real (kind=8), dimension(100) :: zz, nn, kii,gg
real (kind=8) :: lmin,lmax,xxx,delta_p,acceptedfraction
real (kind=8), dimension(ldim) :: ell, C_ell
real (kind=8), dimension(ldim,redshiftbins*(redshiftbins+1)/2,redshiftbins*(redshiftbins+1)/2) :: cov,inv
real (kind=8), dimension(nopd) :: parametersd
real (kind=8), dimension(nopb) :: parametersb

parametersd = (/0.279d0,0.817d0,0.701d0,0.96d0,0.0462d0,-1.0d0,0d0/)
parametersb = (/0.279d0,0.817d0,0.701d0,0.96d0,0.0462d0,-1.0d0,0d0,13.d0/)

!deltaparameters = (/0.01,0.01,0.04,0.02,0.006,0.05,0.2/)

	call set_fiducial_parameters()
	ell_min = 1d1
	ell_max = 1d5



	redshift = 0d0
	fgas = 0.15d0
do j=1,6
	redshift = 0.0d0 + float(j-1)*0.4
	write(*,*) redshift
	do i =1,100
		mass = 7d0 + float(i-1)/99*(11d0)
		mass = 10**mass
		
		write(60+j,*) mass,concentration(mass,redshift)
	end do
end do

stop
	redshift=4.0d0
	call pknonlinear(redshift,k_nonlinear,Pk_nonlinear,pk1h_nonlinear,pk2h_nonlinear)
	do j=1,kdim_linear
		write(50,*) k_nonlinear(j),pk_nonlinear(j)
		end do

	Baryons = .True.
	AC = .False.
	do i=1,8
		Mcrit = 9d0 + float(i-1)
		Mcrit = 10**Mcrit
		write(*,*) Mcrit
		call pknonlinear(redshift,k_nonlinear,Pk_nonlinear,pk1h_nonlinear,pk2h_nonlinear)
	do j=1,kdim_linear
		write(60+i,*) k_nonlinear(j),pk_nonlinear(j)
		end do
	end do

	Baryons = .True.
	AC = .True.
	do i=1,8
		Mcrit = 9d0 + float(i-1)
		Mcrit = 10**Mcrit
		write(*,*) Mcrit
		call pknonlinear(redshift,k_nonlinear,Pk_nonlinear,pk1h_nonlinear,pk2h_nonlinear)
	do j=1,kdim_linear
		write(50+i,*) k_nonlinear(j),pk_nonlinear(j)
		end do
	end do
stop

!	CALL write_cl_data(51,52)
!	BARYONS = .tRUE.
!	ac = .tRUE.
!	mCRIT = 1D12
!	CALL WRITE_CL_DATA(53,54)

!STOP
!	call read_cl_data_DMO()
!	call calculatefishermatrix()
!	call write_fishermatrix(91)
	idum = 545786589
	call read_cl_data_BAR()
	baryons = .True.
	AC = .True.
	write(*,*) chi2(parametersb,nopb)
	baryons = .False.
	write(*,*) chi2(parametersd,nopd)
stop
	call read_cl_data_DMO()
	baryons = .True.
	AC = .True.
	write(*,*) chi2(parametersb,nopb)
	baryons = .False.
	write(*,*) chi2(parametersd,nopd)


!	call mcmc_chain(2222)

end

