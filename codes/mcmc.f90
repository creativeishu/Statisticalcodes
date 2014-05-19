module mcmc
use nr

implicit none

integer (kind=4), Parameter :: mcmc_steps = 100000, NOP = 2
real (kind=8), Parameter :: step_multiplier = 1.d0
character (len=20) :: outfilename='mcmc_output.dat'
real (kind=8), dimension(NOP) :: MinParameter, MaxParameter,Std_dev

contains

!-------------------------------------------------------------------------------
subroutine Initialise()
	MinParameter = (/0d0,0d0/)
	MaxParameter = (/1d2,1d2/)
	Std_dev = (/1d-1,1d-1/)
	
	end subroutine Initialise

!-------------------------------------------------------------------------------

subroutine FirstPosition(parameters)
	integer (kind=4) :: i
	real (kind=8), intent(out), dimension(NOP) :: parameters
	do i = 1,NOP
		parameters(i) = MinParameter(i) + &
					ran2(idum) * &
					(MaxParameter(i) - MinParameter(i))
		end do
	end subroutine FirstPosition

!-------------------------------------------------------------------------------

logical function MetHas(P_old,P_new)
	real (kind=8), intent(in) :: P_old, P_new
	real (kind=8) :: random_num
	if (P_new.gt.P_old) then
		MetHas = .True.
	else
		random_num = ran2(idum)
		if (P_new/P_old.ge.random_num) then
			MetHas = .True.
		else
			MetHas = .False.
		end if
	end if
	return
	end function MetHas

!-------------------------------------------------------------------------------

subroutine Step(PrevStep, NextStep)
	integer (kind=4) :: i
	real (kind=8), intent(in), dimension(NOP) :: PrevStep
	real (kind=8), intent(out), dimension(NOP) :: NextStep
	
	do i=1,NOP
		Nextstep(i) = PrevStep(i) + gasdev(idum) * &
						step_multiplier * &
						Std_dev(i)
		do while (Nextstep(i).lt.MinParameter(i).or.Nextstep(i).gt.MaxParameter(i))
			Nextstep(i) = PrevStep(i) + gasdev(idum) * &
							step_multiplier * &
							Std_dev(i)
			end do			
			end do
	end subroutine Step

!-------------------------------------------------------------------------------

subroutine mcmc_chain(AcceptedFraction)
	real (kind=8), intent(out) :: AcceptedFraction
	real (kind=8), dimension(NOP) :: Params_array_old, Params_array_new	
	real (kind=8) :: likelihood_old, likelihood_new
	real (kind=8) :: chi2
	integer (kind=4) :: i, AcceptedPoints, weight
	logical :: Accept_point
	
	open(11,file=outfilename)
	
	Call Initialise()
	Call FirstPosition(Params_array_old)
	likelihood_old = ran2(idum)
	
	AcceptedPoints = 0
	weight = 0
	do i = 1,mcmc_steps
		weight = weight+1
		Call Step(Params_array_old,Params_array_new)
		likelihood_new = ran2(idum)
		Accept_point = MetHas(likelihood_old,likelihood_new)
		
		if (Accept_point) then
			likelihood_old = likelihood_new
			Params_array_old = Params_array_new
			AcceptedPoints = AcceptedPoints + 1
			write(11,*) likelihood_new, weight, Params_array_new
			weight = 0
			endif

		if (mod(i,5000)==0) then
			AcceptedFraction = float(AcceptedPoints)/(i)
			write(*,*) "Accepted Fraction after",i,'steps is: ',AcceptedFraction
			end if
		end do
		close(11)
		
		AcceptedFraction = float(AcceptedPoints)/(i-1)
		write(*,*) "Final accepted fraction after ",i-1,"steps is:", AcceptedFraction
	end subroutine mcmc_chain

!-------------------------------------------------------------------------------

end module mcmc