! 0911127 by hongbo
! //397SNa[Constitution] + R + A / OmegaCDM model//
! 307Union + R + A / OmegaCDM model

program wCDM
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
  integer :: i
  real(kind=8),external :: chi2
  real(kind=8) :: beginomega,beginw,gradient(2)
  real(kind=8),parameter :: zero = 1.0D0, step=1.0D-4
  integer :: countcycle =0



!---read sna data---!
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

!---best fit---!

  beginomega = 0.2
  beginw = -0.8
  gradient(1) = (chi2(beginomega+step,beginw)-chi2(beginomega,beginw))/step
  gradient(2) = (chi2(beginomega,beginw+step)-chi2(beginomega,beginw))/step

  do while (dot_product(gradient,gradient)>zero)
	beginomega=beginomega-step*gradient(1)
        beginw=beginw-step*gradient(2)
	gradient(1) = (chi2(beginomega+step,beginw)-chi2(beginomega,beginw))/step
	gradient(2) = (chi2(beginomega,beginw+step)-chi2(beginomega,beginw))/step
        countcycle = countcycle + 1
        write(*,*) countcycle, beginomega, beginw, dot_product(gradient,gradient)
  end do

  write(*,*) "The minimal chi2 is", chi2(beginomega,beginw)
  write(*,*) "best fit : "
  write(*,*) "omega = ",beginomega,"w=",beginw
!!!  write(*,*) "cycle",countcycle
  
  open (unit=88,file = "bestfit.txt")
  write(88,*) "The minimal chi2 is", chi2(beginomega,beginw)
  write(88,*) "best fit : "
  write(88,*) "omega = ",beginomega,"w=",beginw
  close(88)

end program

!-------------------------------------------------------------------------
function chi2(omegam0,w)  ! -2*LnLikelihood of Supernovea + CMB[R] +BAO[A]
implicit none

  integer,parameter :: sn_num = 307
  real(kind=8) :: chi2, chi2sn,chi2R,chi2A,R,A
  real(kind=8) :: omegam0,w,abar,bbar,cbar
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
       bbar=bbar+(sd(2,i)-muth(omegam0,w,sd(1,i)))/(sd(3,i)**2)
    end do

    do i=1,sn_num,1
       abar=abar+((sd(2,i)-muth(omegam0,w,sd(1,i)))**2)/(sd(3,i)**2)
    end do
  
  chi2sn = abar - (bbar**2)/cbar

  R=sqrt(omegam0)*Dl(omegam0,w,1.089D3)

  chi2R = ((R-1.715)**2)/(0.021**2)
  
  A=sqrt(omegam0)*((Dl(omegam0,w,3.5D-1)/0.35)**(2.0/3.0))/(hubble(omegam0,w,3.5D-1)**(1.0/3.0))

  chi2A = ((A-0.472)**2)/(0.017**2)

  chi2 = chi2sn + chi2R +chi2A

return
end function
!--------------------------------------------------------------------------

function Dl(omegam0,w,z)   ! /int_{0}^{z} dt/H(t) lum-distance=(1+z)Dl
! this part is to do trape method of integration
! read notes/../integration/trape.f90 for detail of this method
implicit none

  real(kind=8) :: Dl,omegam0,w,z
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
     datalist(i+1)=1.0/hubble(omegam0,w,xmin+i*width)
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  Dl=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  

return
end function

!------------------------------------------------------------------------

function hubble(omegam0,w,z) ! hubble parameter in wCDM model
implicit none

  real(kind=8) :: hubble, omegam0,w,z

  hubble=sqrt(omegam0*((1+z)**3) + (1- omegam0)*((1+z)**(3*(1+w))))

return
end function

!---------------------------------------------------------------------

function muth(omegam0,w,z) ! /mu_th = 5*lg(luminositydistance)+0[not /mu_0]
implicit none

  real(kind=8),external :: Dl
  real(kind=8) :: muth
  real(kind=8) :: omegam0,w,z

  muth=5*log10(1+z) + 5*log10(Dl(omegam0,w,z))

return
end function

!-----------------------------------------------------------------------
