module fisher
use cosmo_module

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
real (kind=8) function chisquare(parameters)
	real (kind=8), intent(in), dimension(numberofparameters) :: parameters
	chisquare=0d0
	end function chisquare

!======================================================================================
subroutine calculatefishermatrix()
	integer (kind=4) :: i,j
	real (kind=8), dimension(numberofparameters) :: parametersplus
	real (kind=8), dimension(numberofparameters) :: parametersminus
	real (kind=8), dimension(numberofparameters) :: parametersplusplus
	real (kind=8), dimension(numberofparameters) :: parametersplusminus
	real (kind=8), dimension(numberofparameters) :: parametersminusplus
	real (kind=8), dimension(numberofparameters) :: parametersminusminus
	real (kind=8) :: fiducial_chi2,chi2
	
	fiducial_chi2 = chisquare(fiducialparameters)

	do i=1,numberofparameters
		do j=i,numberofparameters
		
			if (i.eq.j) then
				parametersplus = fiducialparameters
				parametersplus(i) = parametersplus(i) + deltaparameters(i)	
			
				parametersminus = fiducialparameters
				parametersminus(i) = parametersminus(i) + deltaparameters(i)	
			
				Fishermatrix(i,j) = (chisquare(parametersplus) - &
									2*fiducial_chi2 + chisquare(parametersminus)) / &
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
				
				Fishermatrix(i,j) = (chisquare(parametersplusplus) + &
									chisquare(parametersminusminus) - &
									chisquare(parametersplusminus) - &
									chisquare(parametersminusplus)) / &
									4d0/deltaparameters(i)/deltaparameters(j)
				Fishermatrix(j,i) = Fishermatrix(i,j)
			end if
			
			
			end do
			end do

	Fishermatrix = Fishermatrix/2d0
	end subroutine calculatefishermatrix

!======================================================================================
subroutine fisher2covariance()
	integer (kind=4) :: i
	integer (kind=4), dimension(numberofparameters) :: indx
	real (kind=8) :: d
	real (kind=8), dimension(numberofparameters,numberofparameters) :: fishertemp

	fishertemp = Fishermatrix
	call ludcmp(fishertemp,numberofparameters,numberofparameters,indx,d)

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

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
!        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

!======================================================================================
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END


!======================================================================================

end module fisher