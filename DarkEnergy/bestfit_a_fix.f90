! 0911130 by hongbo
! //397SNa[Constitution] + R + A / //
! 307Union + R + A / SuQP's Model, 2 Bins with z2 fixed to 1 & a fixed to -1

program bestfit
implicit none

!  type snadata
!    real(kind=4) :: z,mu,sigma
!  end type
!  type(snadata) :: sd(sn_num)

  integer,parameter :: sn_num = 307, n=100
  real(kind=8) :: sd(3,sn_num)
  common sd
  character(80) :: filename = "sn_z_mu_dmu.txt"
  logical :: alive
  character(6) :: sn_name
  real(kind=8),external :: chi2
  real(kind=8),parameter :: zero = 1.0D-1, step=1.0D-4
  integer :: i,j,countcycle =0
  integer,parameter :: dimen = 6
  real(kind=8) :: beginx(dimen), gradient(dimen),stepvector(dimen,dimen)
! beginx = [omegam0,z1,a0,a1,a2,a]


!-------------------!
!---read sna data---!
!-------------------!
  write(*,*) 'Reading Supernovea Data, Uion Set'

  inquire (file=filename, exist=alive)
    if (.not.alive) then
       write(*,*) filename,"doesn't exist"
       stop
    else 
       write(*,*) "data file exists"
    end if

  open (unit=44, file=filename)
  do i=1,sn_num,1
     read(44,*) sn_name, sd(:,i)
!     write(*,*) sn_name, sd(i)
  end do
  close(44)

!--------------!
!---best fit---!
!--------------!

  beginx(1) = 0.2
  beginx(2) = 0.4
  beginx(3) = -0.6
  beginx(4) = -1.5
  beginx(5) = 0.5
  beginx(6) = -1.0




!-Unit Matrix-!

  do i=1,dimen,1
     do j=1,dimen,1
        if (i==j) then
           stepvector(i,j) = step
        else
           stepvector(i,j) = 0.0
        end if
     end do
  end do

!---get gradient---!
	do i=1,dimen-1,1
	  gradient(i) = (chi2(beginx+stepvector(i,:))-chi2(beginx))/step
	end do

	gradient(dimen)=0

  open (unit=88,file = "bestfit.txt")

  do while (dot_product(gradient,gradient)>zero)
	beginx=beginx-step*gradient
	do i=1,dimen-1,1
	  gradient(i) = (chi2(beginx+stepvector(i,:))-chi2(beginx))/step
	end do
        countcycle = countcycle + 1
        write(*,*) countcycle, beginx(:), dot_product(gradient,gradient),chi2(beginx)
        write(88,*) countcycle, beginx(:), dot_product(gradient,gradient),chi2(beginx)
  end do

  write(*,*) "The minimal chi2 is", chi2(beginx)
  write(*,*) "best fit : "
  write(*,*) "omega = ",beginx(1),"z1=",beginx(2),"a0=",beginx(3),"a1=",beginx(4),"a2=",beginx(5),"a=",beginx(6)
  

  write(88,*) "The minimal chi2 is", chi2(beginx)
  write(88,*) "best fit : "
  write(88,*) "omega = ",beginx(1),"z1=",beginx(2),"a0=",beginx(3),"a1=",beginx(4),"a2=",beginx(5),"a=",beginx(6)
  close(88)

end program

!-------------------------------------------------------------------------
function chi2(param)  ! -2*LnLikelihood of Supernovea + CMB[R] +BAO[A]
implicit none

  integer,parameter :: sn_num = 307
  real(kind=8) :: chi2, chi2sn,chi2R,chi2A,R,A
  real(kind=8) :: param(6),abar,bbar,cbar
  real(kind=8),external :: muth
  real(kind=8),external :: Dl
  real(kind=8),external :: hubble
  real(kind=8) :: sd(3,sn_num)
  common sd
  integer :: i

    abar=0.0
    bbar=0.0
    cbar=0.0
    do i=1,sn_num,1
       cbar=cbar+1.0/(sd(3,i)**2)
    end do       
! not good. actually we can calculate cbar only once 
!instead of calculating it everytime when call the function chi2sn
  
    do i=1,sn_num,1
       bbar=bbar+(sd(2,i)-muth(param,sd(1,i)))/(sd(3,i)**2)
    end do

    do i=1,sn_num,1
       abar=abar+((sd(2,i)-muth(param,sd(1,i)))**2)/(sd(3,i)**2)
    end do
  
  chi2sn = abar - (bbar**2)/cbar

  R=sqrt(param(1))*Dl(param,1.089D3)

  chi2R = ((R-1.715)**2)/(0.021**2)
  
  A=sqrt(param(1))*((Dl(param,3.5D-1)/0.35)**(2.0/3.0))/(hubble(param,3.5D-1)**(1.0/3.0))

  chi2A = ((A-0.472)**2)/(0.017**2)

  chi2 = chi2sn + chi2R +chi2A

return
end function
!--------------------------------------------------------------------------

function Dl(param,z)   ! /int_{0}^{z} dt/H(t) lum-distance=(1+z)Dl
! this part is to do trape method of integration
! read notes/../integration/trape.f90 for detail of this method
implicit none

  real(kind=8) :: Dl,z,param(6)
  real(kind=8),external :: hubble
  real(kind=8) :: xmin=0.0,xmax 
  integer,parameter :: seg=10000  
  real(kind=8) :: datalist(seg+1)   
  real(kind=8) :: width  
  real(kind=8) :: s
  integer :: i

  s=0.0
  xmax=z
  width=(xmax-xmin)/seg
  do i=0,seg,1
     datalist(i+1)=1.0/hubble(param,xmin+i*width)
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  Dl=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  

return
end function

!---------------------------------------------------------------------

function muth(param,z) ! /mu_th = 5*lg(luminositydistance)+0[not /mu_0]
implicit none

  real(kind=8),external :: Dl
  real(kind=8) :: muth
  real(kind=8) :: z,param(6)

  muth=5*log10(1+z) + 5*log10(Dl(param,z))

return
end function

!-----------------------------------------------------------------------!
!-----------------------------Hubble Parameter--------------------------!
!-----------------------------------------------------------------------!

function hubble(param,z)
implicit none

real(kind=8) :: hubble, z, param(6)
real(kind=8),external :: Intergal1,Intergal2,Intergal3

! param = [omegam0,z1,a0,a1,a2,a]

if ( z < param(2) ) then

  hubble = sqrt(param(1)*((1+z)**3) + (1- param(1))*Intergal1(param(2),param(3),param(4),z))

else if ( z>=param(2) .and. z<=1.0 ) then

  hubble = sqrt(param(1)*((1+z)**3) + (1- param(1))*Intergal2(param(2),param(3),param(4),param(5),z))

else if ( z>1.0 ) then
  
  hubble = sqrt(param(1)*((1+z)**3) + (1- param(1))*Intergal3(param(2),param(3),param(4),param(5),param(6),z))

end if

return
end function


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal1(z1,a0,a1,z)
implicit none

real(kind=8) :: Intergal1,z1,a0,a1,z

Intergal1 = ((1+z)**(3*(1+a0-(a1-a0)/z1)))*exp(3*z*(a1-a0)/z1)

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal2(z1,a0,a1,a2,z)
implicit none

real(kind=8) :: Intergal2,z1,a0,a1,a2,z

Intergal2 = ((1+z1)**(3*(1+a0-(a1-a0)/z1)))*(((1+z)/(1+z1))**(3*(1+a1-(a2-a1)/(1-z1)-z1*(a2-a1)/(1-z1))))&
&*exp(3*(a1-a0) + 3*(z-z1)*(a2-a1)/(1-z1))

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal3(z1,a0,a1,a2,a,z)
implicit none

real(kind=8) :: Intergal3,z1,a0,a1,a2,a,z

Intergal3=((1+z1)**(3*(1+a0-(a1-a0)/z1)))*(((1+1.0)/(1+z1))**(3*(1+a1-(a2-a1)/(1.0-z1)-z1*(a2-a1)/(1.0-z1))))&
&*exp(3*(a1-a0) + 3*(a2-a1)) !!!*(((1+z)/(1+1.0))**(3*(1+a)))

return
end function
