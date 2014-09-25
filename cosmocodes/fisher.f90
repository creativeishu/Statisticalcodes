module fisher
use cosmo_module
use nr

implicit none

integer (kind=4), parameter :: numberofparameters=7
real (kind=8), dimension(numberofparameters) :: fiducialparameters
real (kind=8), dimension(numberofparameters) :: deltaparameters
real (kind=8), dimension(numberofparameters) :: errorparameters

real (kind=8), dimension(numberofparameters,numberofparameters) :: fishermatrix
real (kind=8), dimension(numberofparameters,numberofparameters) :: covariance
real (kind=8), dimension(numberofparameters,numberofparameters) :: correlation

real (kind=8), dimension(numberofparameters,numberofparameters) :: a_ellipse
real (kind=8), dimension(numberofparameters,numberofparameters) :: b_ellipse
real (kind=8), dimension(numberofparameters,numberofparameters) :: th_ellipse

real (kind=8), dimension(3) :: ConfidenceLevelMultiplier

contains

!======================================================================================
real (kind=8) function chisquare(parameters,nop)
	integer (kind=4) :: nop
	real (kind=8), intent(in), dimension(numberofparameters) :: parameters
	chisquare = chi2(parameters,nop)
	end function chisquare

!======================================================================================

subroutine write_fishermatrix(number)
	integer (kind=4) :: i,number
	do i=1,numberofparameters
		write(number,*) fishermatrix(i,:)
		end do
	end subroutine write_fishermatrix

!======================================================================================
subroutine calculatefishermatrix()
	integer (kind=4) :: i,j
	real (kind=8), dimension(numberofparameters) :: parametersplus
	real (kind=8), dimension(numberofparameters) :: parametersminus
	real (kind=8), dimension(numberofparameters) :: parametersplusplus
	real (kind=8), dimension(numberofparameters) :: parametersplusminus
	real (kind=8), dimension(numberofparameters) :: parametersminusplus
	real (kind=8), dimension(numberofparameters) :: parametersminusminus
	real (kind=8) :: fiducial_chi2
	
	fiducial_chi2 = chisquare(fiducialparameters,numberofparameters)

	do i=1,numberofparameters
		do j=i,numberofparameters
		
			if (i.eq.j) then
				parametersplus = fiducialparameters
				parametersplus(i) = parametersplus(i) + deltaparameters(i)	
			
				parametersminus = fiducialparameters
				parametersminus(i) = parametersminus(i) + deltaparameters(i)	
			
				Fishermatrix(i,j) = (chisquare(parametersplus,numberofparameters) - &
									2*fiducial_chi2 + chisquare(parametersminus,numberofparameters)) / &
									deltaparameters(i)**2			
			else
				
				parametersplusplus = fiducialparameters
				parametersplusminus = fiducialparameters
				parametersminusplus = fiducialparameters
				parametersminusminus = fiducialparameters

				parametersplusplus(i) = parametersplusplus(i) + deltaparameters(i)
				parametersplusplus(j) = parametersplusplus(j) + deltaparameters(j)

				parametersminusminus(i) = parametersminusminus(i) - deltaparameters(i)
				parametersminusminus(j) = parametersminusminus(j) - deltaparameters(j)

				parametersplusminus(i) = parametersplusminus(i) + deltaparameters(i)
				parametersplusminus(j) = parametersplusminus(j) - deltaparameters(j)

				parametersminusplus(i) = parametersminusplus(i) - deltaparameters(i)
				parametersminusplus(j) = parametersminusplus(j) + deltaparameters(j)
				
				Fishermatrix(i,j) = (chisquare(parametersplusplus,numberofparameters) + &
									chisquare(parametersminusminus,numberofparameters) - &
									chisquare(parametersplusminus,numberofparameters) - &
									chisquare(parametersminusplus,numberofparameters)) / &
									4d0/deltaparameters(i)/deltaparameters(j)
				Fishermatrix(j,i) = Fishermatrix(i,j)
			end if
			
			write(*,*) "Fisher matrix element",i,j,Fishermatrix(i,j)
			end do
			end do

	Fishermatrix = Fishermatrix/2d0
	end subroutine calculatefishermatrix

!======================================================================================
subroutine fisher2covariance()
	integer (kind=4) :: i,i1,i2
	integer (kind=4), dimension(numberofparameters) :: indx
	real (kind=8) :: d
	real (kind=8), dimension(numberofparameters,numberofparameters) :: fishertemp

	fishertemp = Fishermatrix
	call ludcmp(fishertemp,numberofparameters,numberofparameters,indx,d)

	do i1=1,numberofparameters
		do i2=1,numberofparameters
				covariance(i1,i2)=0.
				end do
				covariance(i1,i1)=1.
				end do

	do i=1,numberofparameters
		call lubksb(fishertemp,numberofparameters,numberofparameters,indx,covariance(1,i))
		end do
	end subroutine fisher2covariance

!======================================================================================
subroutine ellipseparameters()
	integer (kind=4) :: i,j
	
	do i=1,numberofparameters
		do j=i,numberofparameters
			a_ellipse(i,j) = (covariance(i,i)+covariance(j,j))/2 + &
					sqrt((covariance(i,i)-covariance(j,j))**2/4 + covariance(i,j)**2)
					
			b_ellipse(i,j) = (covariance(i,i)+covariance(j,j))/2 - &
					sqrt((covariance(i,i)-covariance(j,j))**2/4 + covariance(i,j)**2)
					
			th_ellipse(i,j) = 0.5*atan(2d0*covariance(i,j)/ &
							(covariance(i,i)-covariance(j,j)))
			
			if (i.ne.j) then
				a_ellipse(j,i) = a_ellipse(i,j)
				b_ellipse(j,i) = b_ellipse(i,j)
				th_ellipse(j,i) = th_ellipse(i,j)
				end if
			end do
			end do
	end subroutine ellipseparameters

!======================================================================================
subroutine getsigma()
	integer (kind=4) :: i
	do i=1,numberofparameters
		errorparameters(i) = sqrt(covariance(i,i))
		end do
	end subroutine getsigma
	
!======================================================================================
subroutine FigureOfMerit(i,j,CL,fom)
	integer (kind=4), intent(in) :: i,j,CL
	integer (kind=8), intent(out) :: fom
	
	ConfidenceLevelMultiplier = (/1.52d0,2.48d0,3.44d0/)
	fom = pi * ConfidenceLevelMultiplier(CL)**2 * a_ellipse(i,j) * b_ellipse(i,j) 

	end subroutine FigureOfMerit

!======================================================================================

end module fisher
