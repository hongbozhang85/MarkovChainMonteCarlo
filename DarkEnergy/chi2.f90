  program main
  implicit none

  integer,parameter :: sn_num = 397
  real(kind=4) :: sd(3,sn_num+12)
  common sd
  character(80) :: filename = "sn397hubble12.dat"
  logical :: alive
  character(6) :: sn_name
  real(kind=4),external :: chi2,muth,Dl
  integer :: i
  real :: aaa(3)=(/ 0.274,0.72,0.046 /)

  write(*,*) muth(0.274,0.088), Dl(0.274,0.088)


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
  do i=1,sn_num+12,1
     read(44,*) sn_name, sd(:,i)
!     write(*,*) sn_name, sd(i)
  end do
  close(44)

  write(*,*) chi2(aaa)
  

  end program

!---------------------

function chi2(params)  ! -2*LnLikelihood of Supernovea + CMB[R] +BAO[A]
implicit none

  integer,parameter :: sn_num = 397
  real(kind=4) :: chi2,chi2sn,chi2h
  real(kind=4) :: params(3),abar,bbar,cbar
  real(kind=4),external :: muth,Dl,hubble,RDl,chi2cmb,chi2bao
  real(kind=4) :: sd(3,sn_num+12)
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
       bbar=bbar+(sd(2,i)-muth(params(1),sd(1,i)))/(sd(3,i)**2)
    end do

    do i=1,sn_num,1
       abar=abar+((sd(2,i)-muth(params(1),sd(1,i)))**2)/(sd(3,i)**2)
    end do
  
  chi2sn = abar - (bbar**2)/cbar

  chi2h=0.0

  do i=1,12,1
    chi2h= chi2h+((params(2)*hubble(params(1),sd(1,i+sn_num))-sd(2,i+sn_num))**2)/(sd(3,i+sn_num)**2)
  end do

  chi2= chi2sn + chi2h + chi2cmb(params) + chi2bao(params)
  
  write(*,*) chi2cmb(params),chi2sn, chi2h, chi2bao(params)

return
end function
!--------------------------------------------------------------------------

  function chi2cmb(params)
  implicit none

  integer,parameter :: num_params=3 ! params=(omegam0,h100,omegab0)
  real :: chi2cmb,params(num_params),R,la,zs,cmbarray(3),g1,g2
  real :: cov(3,3)
  integer :: rr,cc
  real,external :: RDl,rs
  
  data((cov(rr,cc),rr=1,3),cc=1,3)/ 1.8,27.968,-1.103,27.968,5667.577,-92.263,-1.103,-92.263,2.923 /

  g1= 0.0783*((params(3)*(params(2)**2))**(-0.238))/(1+39.5*((params(3)*(params(2)**2))**0.763))

  g2= 0.560/(1+21.1*((params(3)*(params(2)**2))**1.81))

  zs= 1048*(1+0.00124*((params(3)*(params(2)**2))**(-0.738)))*(1+g1*((params(1)*(params(2)**2))**g2))

  R=sqrt(params(1))*RDl(params(1),zs)

  la=3.1415926*RDl(params(1),zs)/rs(params,zs)

  cmbarray(1)= la-302.1
  cmbarray(2)= R-1.710
  cmbarray(3)= zs-1090.04

  write(*,*) rs(params,zs),la,cmbarray(:)

  chi2cmb= dot_product(cmbarray,matmul(cov,cmbarray))

  return
  end function

!---------------------

  function chi2bao(params)
  implicit none

  integer,parameter :: num_params=3  ! params=(omegam0,h100,omegab0)
  real :: chi2bao, params(num_params), bao1,bao2,zd,b1,b2 !bao1=Dv(0.35)/Dv(0.2); bao2=rs(zd)/Dv(0.275)
  real,external :: Dl,rs,hubble

  b1= 0.313*((params(1)*params(2)*params(2))**(-0.419))*(1+0.607*((params(1)*params(2)*params(2))**0.674))

  b2= 0.238*((params(1)*params(2)*params(2))**0.223)

zd=1291*((params(1)*(params(2)**2))**0.251)*(1+b1*((params(3)*(params(2)**2))**b2))/(1+0.659*((params(1)*(params(2)**2))**0.828))

  bao1= (((0.35/hubble(params(1),0.35))**(1.0/3.0))*(Dl(params(1),0.35)**(2.0/3.0)))/&
&(((0.2/hubble(params(1),0.2))**(1.0/3.0))*(Dl(params(1),0.2)**(2.0/3.0)))

  bao2= rs(params,zd)/(((0.275/hubble(params(1),0.275))**(1.0/3.0))*(Dl(params(1),0.275)**(2.0/3.0)))

  chi2bao= ((bao1-1.736)**2)/(0.065**2)+((bao2-0.139)**2)/(0.0037**2)

!  write(*,*) bao1,bao2,rs(params,zd)

  return
  end function


!---------------------

  function rs(params,z)  ! comoving sound horizon * H0
  implicit none

  real :: rs, params(3),z, cs,a
  real,external :: hubble
  real :: zmin=1.0E-10,zmax
  integer,parameter :: seg=100 ! seg=1E4,rs=6.57E-2,la=159; seg=1E5,rs=6.15E-2,la=170
  real :: datalist(seg+1)
  real :: width,s
  integer :: i

!  cs=sqrt(3+9*a*params(3)*params(2)*params(2)/(4*2.469E-5))

! use trape method, a bad method in this case--delete this note

  s=0.0
  zmax=1.0/(1.0+z)
  width=(zmax-zmin)/seg
  do i=0,seg,1
     a=zmin+i*width
     datalist(i+1)=1.0/(a*a*hubble(params(1),1/a-1)*sqrt(3+9*a*params(3)*params(2)*params(2)/(4*2.469E-5)))
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  rs=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)

  return
  end function

!---------------------


function Dl(omegam0,z)   ! /int_{0}^{z} dt/H(t) lum-distance=(1+z)Dl
! this part is to do trape method of integration
! read notes/../integration/trape.f90 for detail of this method
implicit none

  real(kind=4) :: Dl,omegam0,z
  real(kind=4),external :: hubble
  real(kind=4) :: xmin=0.0,xmax 
  integer,parameter :: seg=100  
  real(kind=4) :: datalist(seg+1)   
  real(kind=4) :: width  
  real(kind=4) :: s
  integer :: i

  s=0.0
  xmax=z
  width=(xmax-xmin)/seg
  do i=0,seg,1
     datalist(i+1)=1.0/hubble(omegam0,xmin+i*width)
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  Dl=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  

return
end function

!------------------------------------------------------------------------

function RDl(omegam0,z)   ! /int_{0}^{z} dt/H(t) lum-distance=(1+z)Dl
! dived integrate interval into 2 parts, and integrate them seperately

implicit none

  real(kind=4) :: RDl,omegam0,z
  real(kind=4),external :: hubble
  real(kind=4) :: xmin=0.0,xmax 
  integer,parameter :: xmid=5.0,seg= 10000, seg1=500 ! seg1--(0,xmid), seg--(xmid,z) 
  real(kind=4) :: datalist(seg+1), datalist1(seg1+1)   
  real(kind=4) :: width, width1
  real(kind=4) :: s,s1
  integer :: i

  s=0.0
  s1=0.0
  xmax=z
  width1=(xmid-xmin)/seg1
  width=(xmax-xmid)/seg

  do i=0,seg1,1
     datalist1(i+1)=1.0/hubble(omegam0,xmin+i*width1)
  end do

  do i=0,seg,1
     datalist(i+1)=1.0/hubble(omegam0,xmid+i*width)
  end do

  do i=0,seg1,1
     s1=s1+datalist1(i+1)
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  RDl=width1*(s1-datalist1(1)/2.0-datalist1(seg1+1)/2.0)+width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  

return
end function

!------------------------------------------------------------------------

function hubble(omegam0,z) ! hubble parameter in LCDM model
implicit none

  real(kind=4) :: hubble, omegam0,z

  hubble=sqrt(omegam0*((1+z)**3) + 1.929*4.15E-5*((1+z)**4) + 1- omegam0)

return
end function

!---------------------------------------------------------------------

function muth(omegam0,z) ! /mu_th = 5*lg(luminositydistance)+0[not /mu_0]
implicit none

  real(kind=4),external :: Dl
  real(kind=4) :: muth
  real(kind=4) :: omegam0,z

  muth=5*log10(1+z) + 5*log10(Dl(omegam0,z))


return
end function

!-----------------------------------------------------------------------
