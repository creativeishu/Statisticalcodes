module cosmo_module

use nr

implicit none

real (kind=8), Parameter :: c_light = 3.0d5		! in km/sec^2
real (kind=8), Parameter :: pi = 4d0*atan(1d0)
real (kind=8), Parameter :: G_constant = 6.6738480d-11	! in m^3/kg/sec^2

! Cosmology
real (kind=8), Parameter :: Omega_m = 0.279d0	
real (kind=8), Parameter :: Omega_b = 0.0462d0	
real (kind=8), Parameter :: Omega_r = 0d0	
real (kind=8), Parameter :: Omega_k = 0d0	
real (kind=8), Parameter :: Omega_l = 1d0 - Omega_m - Omega_r - Omega_k
real (kind=8), Parameter :: HubbleConstant = 0.701d0	
real (kind=8), Parameter :: n_s = 0.960d0	
real (kind=8), Parameter :: sigma_8 = 0.817d0	
real (kind=8), Parameter :: w0 = -1.0d0	
real (kind=8), Parameter :: wa = 0.0d0	

real (kind=8), Parameter :: T_cmb = 2.72d0
real (kind=8), Parameter :: RhoC = 2.778d11
real (kind=8), Parameter :: RhoM = RhoC * Omega_m

real (kind=8), parameter :: overdensityCR = 2.0d2
real (kind=8), Parameter :: deltaC = 1.686d0
real (kind=8) :: overdensity_matter

!baryon parameters
logical :: baryons
logical :: AC
real (kind=8) :: Mcrit 
real (kind=8), Parameter :: beta_Mcrit = 2d0

!linear Power spectrum
integer (kind=4), Parameter :: kdim_linear=500
real (kind=8), dimension(kdim_linear) :: k_linear, Pk_linear, y2pk_linear,dk_linear

!Growthfactor
integer (kind=4), Parameter :: zdim_gf=500
real (kind=8), dimension(zdim_gf) :: redshift_gf, gf_gf, y2gf_gf

contains

!-------------------------------------------------------------------------------

real (kind=8) function likelihood(chi2)
	real (kind=8), intent(in) :: chi2
	likelihood = dexp(-chi2/2d0)
	
	end function likelihood

!-------------------------------------------------------------------------------

subroutine Initialise_cosmo_module(redshift)
	integer (kind=4) :: i
	real (kind=8), intent(in) :: redshift
	real (kind=4), dimension(kdim_linear) :: kk, pkk
	real (kind=4) :: kmin,kmax
	
	kmin = 1.e-3
	kmax = 1.e2
	
	call linearPk(kmin,kmax,kdim_linear,kk,pkk)
	k_linear = dble(kk)
	pk_linear = dble(pkk)
	call spline(k_linear,pk_linear,kdim_linear,0.0d0,0.0d0,y2pk_linear)
	do i=1,kdim_linear
		if (i.eq.1) then
			dk_linear(i) = k_linear(i+1) - k_linear(i)
		else
			dk_linear(i) = k_linear(i) - k_linear(i-1)
		endif
		end do
		
!	overdensity_matter = overdensityCR2M(overdensityCR,redshift)
	overdensity_matter = 2d2
	call growthfactorDE(redshift_gf,gf_gf)
	call spline(redshift_gf,gf_gf,zdim_gf,0.0d0,0.0d0,y2gf_gf)
	
	end subroutine Initialise_cosmo_module


!-------------------------------------------------------------------------------

real (kind=8) function Ez(redshift)
	real (kind=8), intent(in) :: redshift
	real (kind=8) :: a	

	a = 1d0/(1d0+redshift)
	Ez  = 	(Omega_m/a**3 + &
		Omega_k/a**2 + &
		Omega_r/a**4 + &
		Omega_l/a**(3.0*(1d0+w0+wa))/exp(3.0*wa*(1d0-a)))**0.5
	end function Ez


!-------------------------------------------------------------------------------
! Comoving distance

real (kind=8) function ComovingDistance(z)
	real (kind=8), intent(in) :: z
	real (kind=8) :: ss

	call qgaus(func,0.d0,z,ss)
	ComovingDistance=c_light/1d2*ss
	end function ComovingDistance
!-----------------------------------
real (kind=8) function func(z)
	real (kind=8), intent(in) :: z
	func=1.d0/Ez(z)
	return
end function

!-------------------------------------------------------------------------------

real (kind=8) function concentration(mass,redshift)
	real (kind=8), intent(in) :: mass,redshift
	real (kind=8), parameter :: w = 0.029d0
	real (kind=8), parameter :: m = 0.097d0
	real (kind=8), parameter :: alpha = -110d0
	real (kind=8), parameter :: beta = 2469d0
	real (kind=8), parameter :: gammaa = 16.89d0
	real (kind=8) :: a, b

	a = w*redshift - m
	b = alpha/(redshift + gammaa) + beta/(redshift + gammaa)**2
	concentration = a * log10(mass) + b
	concentration = 10**concentration
	if (concentration.lt.4d0) then
		concentration = 4d0
		end if
	end function concentration

!-------------------------------------------------------------------------------

real (kind=8) function VirialRadius(mass,redshift)
	real (kind=8), intent(in) :: mass, redshift

	VirialRadius =	1.63d-2 * &
			mass**(1./3.) * &
			Ez(redshift)**(-2./3.) * &
			1d-3
	end function VirialRadius

!-------------------------------------------------------------------------------

real (kind=8) function Rs(mass,redshift)
	real (kind=8), intent(in) :: mass, redshift

	Rs = VirialRadius(mass,redshift)/concentration(mass,redshift)
	end function Rs

!-------------------------------------------------------------------------------

subroutine NFW(mass,redshift,Rdim,radius,density) ! density in Msun/mpc^3
	integer (kind=4), intent(in) :: Rdim
	real (kind=8), intent(in) :: mass, redshift
	real (kind=8), intent(out), dimension(Rdim) :: radius, density
	real (kind=8) :: Rvir, Rs, conc
	Real (kind=8) :: Rho_s
	integer (kind=4) :: i

	Rvir = VirialRadius(mass,redshift)
	conc = concentration(mass,redshift)
	Rs = Rvir/conc
	Rho_s = mass / 4.0 / pi / Rs**3 &
		 / (log(1 + conc) - conc / (1 + conc))
	
	forall (i=1:Rdim)
		radius(i) = log10(Rvir/Rdim) + float(i-1)/(Rdim-1)*(log10(Rvir)-log10(Rvir/Rdim))
		radius(i) = 10**radius(i)
		endforall
	density = Rho_s / (radius/Rs) / (1d0+radius/Rs)**2
	end subroutine NFW

!-------------------------------------------------------------------------------

subroutine DensityProfileGas(mass,redshift,fgas,Rdim,radius,density)
	integer (kind=4), intent(in) :: Rdim
	integer (kind=4) :: i
	real (kind=8), intent(in) :: mass, redshift,fgas
	real (kind=8), intent(out), dimension(Rdim) :: radius, density
	real (kind=8) :: xeq,GammaGas,power,norm,Rvir,conc,R_s

	Rvir = VirialRadius(mass,redshift)
	conc = concentration(mass,redshift)
	R_s = Rvir/conc
	xeq=conc/sqrt(5d0)
    GammaGas = 1d0+ (((1d0+xeq)*(log(1d0+xeq))-(xeq))/&
    			((1d0+3d0*xeq)*(log(1d0+xeq))))
    power = 1d0/(GammaGas-1d0)
	
	norm = 0d0
	do i=1,Rdim
		radius(i) = log10(Rvir/Rdim) + float(i-1)/(Rdim-1)*(log10(Rvir)-log10(Rvir/Rdim))
		radius(i) = 10**radius(i)
	    density(i) = (log(1d0 + radius(i)/R_s) / (radius(i)/R_s) )**power
	    if (i>1) then
		norm = norm + density(i) * radius(i)**2 * (radius(i)-radius(i-1))
		endif
		enddo
	norm = mass*fgas/4d0/pi/norm
	density = density*norm
	end subroutine DensityProfileGas

!-------------------------------------------------------------------------------

subroutine DensityProfileStars(mass,redshift,rmean,rsigma,Rdim,&
				radius,density,Mstar_center)
	integer (kind=4), intent(in) :: Rdim	
	integer (kind=4) :: i

	real (kind=8), intent(in) :: mass,redshift,rmean,rsigma
	real (kind=8), intent(out), dimension(Rdim) :: radius,density
	real (kind=8), intent(out) :: Mstar_center
 	real (kind=8), Parameter :: m10=11.59
    real (kind=8), Parameter :: m11=1.195
    real (kind=8), Parameter :: n10=0.035
    real (kind=8), Parameter :: n11=-0.0247
    real (kind=8), Parameter :: beta10=1.376
    real (kind=8), Parameter :: beta11=-0.826
    real (kind=8), Parameter :: gamma10=0.608
    real (kind=8), Parameter :: gamma11=0.329
    real (kind=8) :: m,n,beta,gm,norm,Rmin,Rvir
    
    m=1d1**(m10+m11*(redshift/(redshift+1)))
    n=n10+n11*(redshift/(redshift+1))
    beta=beta10+beta11*(redshift/(redshift+1))
    gm=gamma10+gamma11*(redshift/(redshift+1))
    Mstar_center=mass*(2d0*n*((mass/m)**(-beta)+(mass/m)**gm)**(-1.))

    Rmin = 1d-3
    density = exp(-((radius-Rmean*Rmin)/(Rsigma*Rmin))**2)
    Rvir = VirialRadius(mass,redshift)

	norm = 0d0
	do i=1,Rdim
		radius(i) = log10(Rvir/Rdim) + float(i-1)/(Rdim-1)*(log10(Rvir)-log10(Rvir/Rdim))
		radius(i) = 10**radius(i)
    	density(i) = exp(-((radius(i)-Rmean*Rmin)/(Rsigma*Rmin))**2)
	    if (i>1) then
			norm = norm + density(i) * radius(i)**2 * &
				(radius(i)-radius(i-1))
			endif
		enddo
    norm = Mstar_center / 4.0 / pi /norm
    density = norm*density
	end subroutine DensityProfileStars

!-------------------------------------------------------------------------------

subroutine AdiabaticContraction(radius,rho_DM,mass,redshift,Rdim,M_BCG,rho_DM_AC)
	integer (kind=4), intent(in) :: Rdim
	integer (kind=4) :: i
	real (kind=8), intent(in) :: mass, redshift,M_BCG
	real (kind=8), dimension(Rdim) :: radius, rho_DM
	real (kind=8), dimension(Rdim) :: rho_DM_AC
	real (kind=8) :: alpha,fd,Rvir,conc,R_s
	real (kind=8), dimension(Rdim) :: Mi,ratio,x,xx,y2,rho_DM_AC_temp
	
   	alpha = 0.68d0
	Rvir = VirialRadius(mass,redshift)
	conc = concentration(mass,redshift)
	R_s = Rvir/conc    
	x = radius/R_s
	Mi = mass * (log(1+x) - x/(1+x)) / (log(1+conc) - conc/(1+conc))
	fd = 1.d0-M_BCG/mass
	ratio = 1.d0 + alpha * (Mi/(fd*Mi + M_BCG) - 1.d0)
	rho_DM_AC_temp = rho_DM * ratio**(-2.d0)
	xx = x*ratio
	call spline(xx,rho_DM_AC_temp,Rdim,0.0d0,0.0d0,y2)
	do i=1,Rdim
	    call splint(xx,rho_DM_AC_temp,y2,Rdim,x(i),rho_DM_AC(i))
	    end do
	end subroutine AdiabaticContraction

!-------------------------------------------------------------------------------

real (kind=8) function Fgas_Mhalo(Mhalo)
	real (kind=8), intent(in) :: Mhalo
	Fgas_Mhalo = Omega_b/Omega_m / (1d0 + (Mcrit/Mhalo)**beta_Mcrit)
	return
	end function Fgas_Mhalo

!-------------------------------------------------------------------------------

subroutine UK_single(radius,density,Rdim,k,U_K)
	integer (kind=4) :: i
	integer (kind=4), intent(in) :: Rdim
	real (kind=8), intent(in), dimension(Rdim) :: radius, density
	real (kind=8), dimension(Rdim) :: dr,kr
	real (kind=8), intent(in) :: k
	real (kind=8), intent(out) :: U_K
	real (kind=8) :: norm
	
	do i = 1,Rdim
		if (i.eq.1) then
			dr(i) = radius(i+1) - radius(i)
		else
			dr(i) = radius(i) - radius(i-1)
		endif
		end do
	kr = radius * k		
	U_K = sum(density*radius**2*dr*sin(kr)/kr)
	norm = sum(density*radius**2*dr)
	U_K = U_K/norm

	end subroutine UK_single

!-------------------------------------------------------------------------------

subroutine UK_array(radius,density,Rdim,k,U_K,kdim)
	integer (kind=4) :: i
	integer (kind=4), intent(in) :: Rdim, kdim
	real (kind=8), intent(in), dimension(Rdim) :: radius, density
	real (kind=8), dimension(Rdim) :: dr,kr
	real (kind=8), intent(in), dimension(kdim) :: k
	real (kind=8), intent(out), dimension(kdim) :: U_K
	real (kind=8) :: norm
	
	do i = 1,Rdim
		if (i.eq.1) then
			dr(i) = radius(i+1) - radius(i)
		else
			dr(i) = radius(i) - radius(i-1)
		endif
		end do
	
	do i=1,kdim
		kr = radius * k(i)		
		U_K(i) = sum(density*radius**2*dr*sin(kr)/kr)
		norm = sum(density*radius**2*dr)
		U_K(i) = U_K(i)/norm
		end do
	end subroutine UK_array

!-------------------------------------------------------------------------------

real (kind=8) function overdensityCR2M(overdensity_crit,redshift)
	real (kind=8), intent(in) :: overdensity_crit, redshift
	real (kind=8) :: Om
	
	Om=Omega_m*(1.d0+redshift)**3.d0/(Ez(redshift)**2.d0)
	overdensityCR2M  = overdensity_crit/Om
	return
	end function overdensityCR2M

!-------------------------------------------------------------------------------

subroutine TinkerFparameters(overdensityM,redshift,params)
	integer (kind=4), Parameter :: dim=9
	real (kind=8), dimension(dim) :: overdensity,A,aa,b,c
	real (kind=8), dimension(dim) :: y2_A,y2_aa,y2_b,y2_c
	real (kind=8), intent(in) :: overdensityM,redshift
	real (kind=8), intent(out), dimension(4) :: params
	real (kind=8) :: alpha
	
	!table 2 from Tinker et al 2008
	overdensity = (/200d0,300d0,400d0,600d0,800d0,1200d0,1600d0,2400d0,3200d0/)
	A = (/0.186d0,0.200d0,0.212d0,0.218d0,0.248d0,0.255d0,0.260d0,0.260d0,0.260d0/)
	aa = (/1.47d0,1.52d0,1.56d0,1.61d0,1.87d0,2.13d0,2.30d0,2.53d0,2.60d0/)
	b = (/2.57d0,2.25d0,2.05d0,1.87d0,1.59d0,1.51d0,1.46d0,1.44d0,1.41d0/)
	c = (/1.19d0,1.27d0,1.34d0,1.45d0,1.58d0,1.80d0,1.97d0,2.24d0,2.44d0/)
	
	 call spline(overdensity,A,dim,0.d0,0.d0,y2_A)
	 call spline(overdensity,aa,dim,0.d0,0.d0,y2_aa)
	 call spline(overdensity,b,dim,0.d0,0.0d0,y2_b)
	 call spline(overdensity,c,dim,0.d0,0.d0,y2_c)	
	
	call splint(overdensity,A,y2_A,dim,overdensityM,params(1))	 
	call splint(overdensity,aa,y2_aa,dim,overdensityM,params(2))
	call splint(overdensity,b,y2_b,dim,overdensityM,params(3))
	call splint(overdensity,c,y2_c,dim,overdensityM,params(4))
	
	alpha = 10d0**(-(0.75d0/(dlog10(overdensityM/75.0d0)))**1.2d0)
	params(1) = params(1) * (1 + redshift)**(-0.14d0)
	params(2) = params(2) * (1 + redshift)**(-0.06d0)
	params(3) = params(3) * (1 + redshift)**alpha
	
	end subroutine TinkerFparameters

subroutine Fsigma(params,nu,f)
	integer (kind=4), Parameter :: dim=4
	real (kind=8), intent(in), dimension(dim) :: params
	real (kind=8), intent(in) :: nu
	real (kind=8), intent(out) :: f
	real (kind=8) :: sigma
	
	sigma = (deltaC / nu)
	f = params(1) * ((sigma/params(3))**(-params(2)) + 1) * exp(-params(4)/sigma**2)	
	end subroutine Fsigma

!-------------------------------------------------------------------------------

subroutine TinkerBparameters(overdensityM,params)
	integer (kind=4), Parameter :: dim=6
	real (kind=8), intent(out), dimension(dim) :: params
	real (kind=8), intent(in) :: overdensityM
	real (kind=8) :: y
	
	y = dlog10(overdensityM)
	params(1) = 1.0d0 + 0.24d0 * y * dexp(-(4.0d0/y)**4d0)
	params(2) = 0.44d0 * y - 0.88d0
	params(3) = 1.83d0
	params(4) = 1.5d0
	params(5) = 0.019d0 + 0.107d0 * y + 0.19d0 * dexp(-(4.0d0/y)**4d0)
	params(6) = 2.4d0
	end subroutine TinkerBparameters

subroutine bias(params,nu,b)
	integer (kind=4), Parameter :: dim=6
	real (kind=8), intent(in), dimension(dim) :: params
	real (kind=8), intent(in) :: nu
	real (kind=8), intent(out) :: b
	
	b = 1.0d0 - (params(1) * sqrt(nu)**params(2) / (sqrt(nu)**params(2) + &
			deltaC**params(2))) + &
			params(3) * sqrt(nu)**params(4) + params(5) * sqrt(nu)**params(6)
	end subroutine bias

!-------------------------------------------------------------------------------

subroutine growthfactorDE(redshift,growthfactor)
	integer (kind=4) :: i,j
	real (kind=8), dimension(zdim_gf), intent(out) :: redshift, growthfactor
	real (kind=8), dimension(zdim_gf) :: g, zz, a, z, E, dedz
	real (kind=8) :: a0, h, k1y, k1z


	a0=1.0d-4
	g(1)=a0
	zz(1)=1.0
	do i=1,zdim_gf
		a(i) = a0 + float(i-1)/(zdim_gf-1)*(1.0-a0)
		z(i) = (1.-a(i))/a(i)
		E(i) = Ez(z(i))
		dedz(i) = (1d0/2d0/E(i)) * ( -3d0*omega_m/a(i)**4  + &
			  omega_l * (a(i)**(-3d0*(1.+w0+wa))) * &
			  exp(-3d0*wa*(1d0-a(i)))  * &
			  (3d0*wa-3d0*(1d0+w0+wa)/a(i)))
		if (i>1) then
			h=a(i) - a(i-1)
			k1y = zz(i-1)*h
			k1z = (-(3d0/a(i-1) + dedz(i-1)/E(i-1))*zz(i-1) + &
				3d0*omega_m*g(i-1)/2d0/(a(i-1)**5)/E(i-1)**2)*h
			g(i) = g(i-1) + k1y
			zz(i) = zz(i-1) + k1z
			endif
			end do
		
	g = g/g(zdim_gf)
		
	do i=1,zdim_gf
		redshift(i) = z(zdim_gf-i+1)
		growthfactor(i) = g(zdim_gf-i+1)
		end do
		
	end subroutine growthfactorDE

!-------------------------------------------------------------------------------

real (kind=8) function sigmasquare_single(mass)
	real (kind=8), intent(in) :: mass
	real (kind=8) :: ss
	real (kind=8) :: kmin,kmax
	real (kind=8), dimension(kdim_linear) :: dk,x,W

	real (kind=8) :: radius

	
	radius = (mass/((4.d0/3.d0)*RhoM*pi))**(1.d0/3.d0)

	x = radius*k_linear
	W = 3d0 * (sin(x) - x*cos(x)) / x**3
	
	sigmasquare_single = sum(pk_linear * k_linear**2 * dk_linear * W**2)/2d0/pi**2

	return
	end function sigmasquare_single


!-------------------------------------------------------------------------------

real (kind=8) function mass2nu(mass)
	real (kind=8), intent(in) :: mass
	real (kind=8) :: sigma_M
	
	sigma_M = sigmasquare_single(mass)
	mass2nu = (deltaC/sigma_M)
	return
	end function mass2nu


!-------------------------------------------------------------------------------
subroutine pknonlinear(redshift,k_nonlinear,Pk_nonlinear,pk1h_nonlinear,pk2h_nonlinear)
	real (kind=8), intent(out), dimension(kdim_linear) :: Pk_nonlinear,k_nonlinear
	real (kind=8), intent(out), dimension(kdim_linear) :: Pk1h_nonlinear,Pk2h_nonlinear	
	real (kind=8), intent(in) :: redshift
	integer (kind=4), Parameter :: Mdim=500, Rdim=400
	integer (kind=4) :: i
	real (kind=8) :: Mmin, Mmax,growthfactor,fdm
	real (kind=8), dimension(Mdim) :: Mass_array,nu_array,fgas,dnu_array
	real (kind=8), dimension(Mdim) :: fnu, bnu
	real (kind=8), dimension(Rdim,Mdim) :: radius_matrix, density_matrix
	real (kind=8), dimension(Rdim,Mdim) :: densitygas_matrix, densityMstar_matrix
	real (kind=8), dimension(kdim_linear,Mdim) :: U_K
	real (kind=8), dimension(4) :: params_fnu
	real (kind=8), dimension(6) :: params_bias
	real (kind=8) :: f0, b0, Mstar_center,b0_norm
	
	call Initialise_cosmo_module(redshift)
	
	call splint(redshift_gf,gf_gf,y2gf_gf,zdim_gf,redshift,growthfactor)	
	fdm = 1d0 - Omega_b/Omega_m
	call TinkerFparameters(overdensity_matter,redshift,params_fnu)
	call TinkerBparameters(overdensity_matter,params_bias)
	
	Mmin = 1.d10
	Mmax = 1.d16
	
	do i = 1,Mdim
		Mass_array(i) = dlog10(Mmin)+float(i-1)/(Mdim-1)*(dlog10(Mmax)-dlog10(Mmin))
		Mass_array(i) = 10**Mass_array(i)
		nu_array(i) = mass2nu(Mass_array(i))/growthfactor
		if (baryons) then
			fgas(i) = Fgas_Mhalo(Mass_array(i))
		else
			fgas(i) = 1d0 - fdm
			endif
			
		call Fsigma(params_fnu,nu_array(i),fnu(i))
		call bias(params_bias,nu_array(i),bnu(i))
		
		call NFW(Mass_array(i),redshift,Rdim,radius_matrix(:,i),density_matrix(:,i))
		if (baryons) then
			call DensityProfileStars(Mass_array(i),redshift,0d0,3d0,Rdim,&
				radius_matrix(:,i),densityMstar_matrix(:,i),Mstar_center)
			call DensityProfileGas(Mass_array(i),redshift,fgas(i),Rdim,&
									radius_matrix(:,i),densitygas_matrix(:,i))
			if (AC) then
				call AdiabaticContraction(radius_matrix(:,i),density_matrix(:,i),&
						Mass_array(i),redshift,Rdim,Mstar_center,density_matrix(:,i))
				endif
			density_matrix(:,i) = fdm*density_matrix(:,i) + &
									densityMstar_matrix(:,i) + &
									densitygas_matrix(:,i)	
			endif
		
		call UK_array(radius_matrix(:,i),density_matrix(:,i),Rdim,&
						k_linear,U_K(:,i),kdim_linear)
		
		if (i.gt.1) then
			dnu_array(i) = nu_array(i)-nu_array(i-1)
			endif
		end do
		dnu_array(1) = dnu_array(2)
	
	
	b0_norm = sum(fnu*bnu*dnu_array)
	bnu = bnu/b0_norm
	f0 = 1d0 - sum((fdm + fgas)*fnu*dnu_array)
	b0 = (1d0 - sum((fdm + fgas)*fnu*bnu*dnu_array))/f0
	
	do i =1,kdim_linear
		if (k_linear(i).lt.0.01) then
			pk1h_nonlinear(i) = 0d0
		else
		pk1h_nonlinear(i) = sum(dnu_array * (fdm+fgas) * fnu * &
								Mass_array * U_K(i,:)**2)
		endif
		pk2h_nonlinear(i) = sum(dnu_array * (fdm+fgas) * fnu * bnu * &
								U_K(i,:))
		end do
	Pk1h_nonlinear = Pk1h_nonlinear/RhoM
	Pk2h_nonlinear = (f0*b0 + Pk2h_nonlinear)**2 * Pk_linear

	Pk_nonlinear = Pk1h_nonlinear + Pk2h_nonlinear
	k_nonlinear = k_linear
	end subroutine



!-------------------------------------------------------------------------------

!real (kind=8) function ComovingDistance(redshift)
!	real (kind=8), intent(in) :: redshift
!	real (kind=8) :: z
!	
!	call qgauss(1d0/Ez(z), 0d0, redshift, ComovingDistance)
!	ComovingDistance = ComovingDistance * c_light / hubble / 1d2
!	end function ComovingDistance

!-------------------------------------------------------------------------------
!=======================================================================
! Eisenstein and Hu linear power spectrum/Transfer function

subroutine linearPk(kmin,kmax,numk,k1,pk)


    real    omega0,f_baryon,hubble,Tcmb,kmax,kmin,omegab,omhh
	real	k,tf_full,tf_baryon,tf_cdm,tf_nowiggles,tf_zerobaryon,k_peak
	integer numk, numk1,i
	real (kind=4), dimension(numk) :: k1,tk,pk
	real xx,integral,integrand,&
			normalisation,sigma8,ns,pi


	omega0 = Omega_m
	omegab = Omega_b
	ns = n_s
	sigma8 = sigma_8
	Tcmb = T_cmb
	hubble = HubbleConstant

	f_baryon=omegab/omega0
	omhh = omega0*hubble*hubble
	pi=4.0*atan(1.0)

    call TFset_parameters(omhh, f_baryon, Tcmb)

	if (kmax.le.0) kmax=10.
	if (numk.le.0) numk=50

	integral=0.
	do i=1,numk
	 	k=10.**(i*(log10(kmax/kmin)/numk))*kmin
        call TFtransfer_function(k*hubble,omhh,f_baryon,tf_full,tf_baryon,tf_cdm)
		k1(i)=k
		tk(i)=tf_full
	end do
	
	do i=1,numk-1
		xx=8.*k1(i)
		integrand=(tk(i)**2.)*(k1(i)**(ns+2.))*&
			((3.*((sin(xx)-xx*cos(xx))/xx**3.))**2.)*(k1(i+1)-k1(i))
		integral=integral+integrand
	end do
	
	normalisation=(sigma8**2)*(2.*pi**2)/(integral)
	do i=1,numk
		pk(i)=normalisation*((k1(i))**ns)*(tk(i)**2)
	end do

end subroutine linearPk
!-------------------------------------

!c
!c
!c PART I:------------------- FITTING FORMULAE ROUTINES ----------------- 
!c
!c There are two routines and a set of functions.  
!c   TFset_parameters() sets all the scalar parameters, while 
!c   TFtransfer_function() calculates various transfer functions 
!c
!c Global variables -- We've left many of the intermediate results as
!c global variables in case you wish to access them, e.g. by declaring
!c them as a common block in your main program. 
!c
!c Note that all internal scales are in Mpc, without any Hubble constants! 
!c
!-------------------------------------

subroutine TFset_parameters(omhh,f_baryon,Tcmb)

	real y,omhh,obhh,Tcmb,hubble
	real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,sound_horizon,k_silk,&
	alpha_c,beta_c,alpha_b,beta_b,f_baryon,beta_node
	common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,&
	sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,beta_node

!c Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
!c Input omhh -- The density of CDM and baryons, in units of critical dens,
!c                multiplied by the square of the Hubble constant, in units
!c                of 100 km/s/Mpc */
!c       f_baryon -- The fraction of baryons to CDM */
!c       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
!c		the default reached by inputing Tcmb=0 -- reset on output. */
!c Output nothing, but set many global variables in common block 
!c       GLOBALVARIABLES. You can access them yourself, if you want:
!c
!c	theta_cmb,	/* Tcmb in units of 2.7 K */ 
!c	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
!c	k_equality,	/* Scale of equality, in Mpc^-1 */
!c	z_drag,		/* Redshift of drag epoch */
!c	R_drag,		/* Photon-baryon ratio at drag epoch */
!c	R_equality,	/* Photon-baryon ratio at equality epoch */
!c	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
!c	k_silk,		/* Silk damping scale, in Mpc^-1 */
!c	alpha_c,	/* CDM suppression */
!c	beta_c,		/* CDM log shift */
!c	alpha_b,	/* Baryon suppression */
!c	beta_b,		/* Baryon envelope shift */


	if (f_baryon.le.0) f_baryon=1.e-5
	if (Tcmb.le.0) Tcmb=2.728
        if (omhh.le.0.0) then
	   write(6,*) 'TFset_parameters(): Illegal input'  
	  ! pause
	end if

        obhh = omhh*f_baryon
        theta_cmb = Tcmb/2.7

        z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*omhh*theta_cmb**(-2.) 

          z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
          z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
        z_drag = 1291e0 * omhh**(0.251)/(1e0 + 0.659*omhh**(0.828)) * z_drag
 
        R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag) 
        R_equality = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_equality) 

        sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )/(1.+sqrt(R_equality)))

        k_silk = 1.6*obhh**(0.52)*omhh**(0.73)*(1e0 + (10.4*omhh)**(-0.95))

          alpha_c = ((46.9*omhh)**(0.670)*(1e0+(32.1*omhh)**(-0.532)))
          alpha_c = alpha_c**(-f_baryon) 
	alpha_c = alpha_c*((12.0*omhh)**(0.424)*(1e0 + (45.0*omhh)**(-0.582)))**(-f_baryon**3.)

    
          beta_c = 0.944/(1+(458.*omhh)**(-0.708))
          beta_c = 1.+beta_c*((1.-f_baryon)**((0.395*omhh)**(-0.0266))- 1e0)
	  beta_c = 1./beta_c

          y = (1e0+z_equality)/(1e0+z_drag)
          alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)))
        alpha_b = 2.07*k_equality*sound_horizon*(1.+R_drag)**(-0.75)*alpha_b


        beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt((17.2*omhh)**2.+1e0)

        beta_node = 8.41*omhh**(0.435)

        return
end subroutine TFset_parameters

!-------------------------------------

subroutine TFtransfer_function(k,omhh,f_baryon,tf_full,tf_baryon,tf_cdm)

!c  Calculate transfer function from the fitting parameters stored in
!c  GLOBALVARIABLES.
!c
!c  Input: 
!c	 k -- wavenumber in Mpc^{-1}  
!c        omhh -- The density of CDM and baryons, in units of critical dens,
!c                multiplied by the square of the Hubble constant, in units
!c                of 100 km/s/Mpc */
!c        f_baryon -- The fraction of baryons to CDM */
!c	
!c  Output:
!c	 tf_full -- The full fitting formula, eq. (16), for the matter
!c	            transfer function. 
!c	 tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
!c	 tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
!c


	real k,tf_full,tf_baryon,tf_cdm,q,ks,omhh,s_tilde
	real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,sound_horizon,&
	k_silk,alpha_c,beta_c,alpha_b,beta_b, f_baryon,beta_node
	common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag, R_drag,R_equality,&
	sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,beta_node

        if (k.le.0) then
           write(6,*) 'TFtransfer_function(): Illegal k'
        !   pause
        end if 

	    q = k/13.41/k_equality
	    ks = k*sound_horizon


        tf_cdm = 1./(1.+(ks/5.4)**4.)
	    tf_cdm = tf_cdm*TF_pressureless(q,1.,beta_c) +(1.-tf_cdm)*TF_pressureless(q,alpha_c,beta_c)


	      s_tilde = sound_horizon/(1.+(beta_node/ks)**3.)**(1./3.) 
	      tf_baryon = TF_pressureless(q,1.,1.)/(1.+(ks/5.2)**2.)
	      tf_baryon = tf_baryon + alpha_b/(1.+(beta_b/ks)**3)*exp(-(k/k_silk)**(1.4))
	      tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
	    tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

         return

end subroutine TFtransfer_function


real function TF_pressureless(q,a,b)

	  real q,a,b
	
	  TF_pressureless = Log(exp(1.)+1.8*b*q)
	  TF_pressureless = TF_pressureless/(TF_pressureless + (14.2/a + 386/(1.+69.9*q**1.08))*q**2)

	  return

end	function TF_pressureless

!-------------------------------------

!c
!c
!c
!c
!c
!c
!c PART II:------------------- Scaling Functions ROUTINES ----------------- 
!c
!c       omhh -- The density of CDM and baryons, in units of critical dens,
!c                multiplied by the square of the Hubble constant, in units
!c                of 100 km/s/Mpc */
!c       f_baryon -- The fraction of baryons to CDM */
!c
!c
!c	TF_zerobaryon:     
!c	  Input:  q = k/omhh * (Tcmb/2.7)**2    (k in Mpc^{-1})
!c	  Output: zero baryon TF Eq(29)
!c	TF_nowiggles:      
!c	  Input:  k = wavenumber in Mpc^{-1}, omhh, f_baryon, Tcmb
!c	  Output: shape approximation TF  Eq(30-31)
!c	  Calls: TF_zerobaryon,sound_horizon_fit,alpha_gamma
!c 	sound_horizon_fit: 
!c         Input:  omhh,f_baryon	
!c	  Output: approximate sound horizon in Mpc	
!c	kpeak:		   
!c	  Input:  omhh,f_baryon
!c         Output: first peak location in Mpc
!c	  Calls:  sound_horizon_fit
!c	alpha_gamma:	   
!c	  Input: omhh,f_baryon
!c	  Output: effective small scale suppression

!-------------------------------------

real function TF_zerobaryon(q)

	  real q
	  TF_zerobaryon = log(2.0*exp(1.)+1.8*q)
	  TF_zerobaryon = TF_zerobaryon/(TF_zerobaryon+(14.2 + 731.0/(1+62.5*q))*q**2)

	  return
end function TF_zerobaryon

!-------------------------------------

real function TF_nowiggles(k,omhh,f_baryon,Tcmb)
	  real k,omhh,f_baryon,q_eff,a,Tcmb
	
	    if (Tcmb.le.0) Tcmb=2.728
	    a = alpha_gamma(omhh,f_baryon)
	    q_eff = k/omhh*(Tcmb/2.7)**2
	    q_eff = q_eff/(a+(1.-a)/(1.+(0.43*k*sound_horizon_fit(omhh,f_baryon))**4))

	    TF_nowiggles = TF_zerobaryon(q_eff)
	return
end function TF_nowiggles

!-------------------------------------

real function sound_horizon_fit(omhh,f_baryon)
	 real omhh,obhh,f_baryon
	 obhh = f_baryon*omhh
         sound_horizon_fit = 44.5*log(9.83/omhh) /sqrt(1.+10.0*obhh**(0.75))
	return
end function sound_horizon_fit

!-------------------------------------

real function k_peak(omhh,f_baryon)
	 real omhh,obhh,f_baryon
	 obhh = f_baryon*omhh
         k_peak = 5.*3.14159/2.*(1.+0.217*omhh)/ sound_horizon_fit(omhh,f_baryon)

	 return
end function k_peak

!-------------------------------------

real function alpha_gamma(omhh,f_baryon)
	 real omhh,f_baryon
         alpha_gamma = 1.-0.328*log(431.0*omhh)*f_baryon  + 0.38*log(22.3*omhh)*(f_baryon)**2
	return
end function alpha_gamma


!=======================================================================


end module cosmo_module
