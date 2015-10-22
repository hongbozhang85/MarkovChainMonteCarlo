! /3rd/final/task_3/dl

module CalcLike
 use CMB_Cls
 use cmbtypes
 use cmbdata
 use mpk
 use Random
 use settings
 use ParamDef
 use snovae
 use WeakLen
 use Lya
 implicit none

 logical :: Use_HST = .true.
 logical :: Use_Age_Tophat_Prior = .true.
 logical :: Use_CMB = .true.
 logical :: Use_BBN = .false.
 logical :: Use_Clusters = .false.
 
 integer :: H0_min = 40, H0_max = 100
 real :: Omk_min = -0.3, Omk_max = 0.3
 real :: Use_min_zre = 0
 integer :: Age_min = 10, Age_max = 20
 real :: Temperature  = 1

contains

  function GenericLikelihoodFunction(Params) 
    type(ParamSet)  Params 
    real :: GenericLikelihoodFunction
    real(kind=8) :: PY(num_params)

   !Used when you want to plug in your own CMB-independent likelihood function:
   !set generic_mcmc=.true. in settings.f90, then write function here returning -Ln(Likelihood)
   !Parameter array is Params%P

!----------------------------------------------------- 
    PY=Params%P   
    GenericLikelihoodFunction =  chi2(PY)/2.0  !LogZero 
!----------------------------------------------------
!    stop 'GenericLikelihoodFunction: need to write this function!'

  end function GenericLikelihoodFunction

!-------------------------------------------------------------------------!
!-------------------------------added code--------------------------------!
!-------------------------------------------------------------------------!
subroutine read_sn_data

  integer,parameter :: sn_num = 397
  real(kind=8) :: sd(3,sn_num+12), mytest(num_params)=(/ 0.27,0.72,0.05,-1,0.0 /)
  common sd
  character(80) :: filename = "sn397hubble12.dat"
  logical :: alive
  character(6) :: sn_name
  integer :: i

!---read sna data---!
  write(*,*) 'Reading Supernovea Data, Constitution Set & 12 Hubble Data'

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
!     write(*,*) sn_name, sd(:,i)
  end do
  close(44)
!----read data ends---!
   
    write(*,*) mytest(:) !!!chi2(mytest)

end subroutine

!---------------------

function chi2(P)  ! -2*LnLikelihood of Supernovea + CMB[3] + BAO[2] + Hubble[12]
implicit none

  integer,parameter :: sn_num = 397
  real(kind=8) :: chi2,chi2sn,chi2h
  real(kind=8) :: P(num_params),abar,bbar,cbar
!  real(kind=4),external :: muth,Dl,hubble,RDl,chi2cmb,chi2bao
  real(kind=8) :: sd(3,sn_num+12)
  common sd
  integer :: i
  logical,save :: do_read_data = .true.

  if (do_read_data) then 
     call read_sn_data
     do_read_data = .false.
  end if


    abar=0.0
    bbar=0.0
    cbar=0.0
    do i=1,sn_num,1
       cbar=cbar+1.0/(sd(3,i)**2)
    end do       
! not good. actually we can calculate cbar only once 
!instead of calculating it everytime when call the function chi2sn
  
    do i=1,sn_num,1
       bbar=bbar+(sd(2,i)-mymuth(P,sd(1,i)))/(sd(3,i)**2)
    end do

    do i=1,sn_num,1
       abar=abar+((sd(2,i)-mymuth(P,sd(1,i)))**2)/(sd(3,i)**2)
    end do
  
  chi2sn = abar - (bbar**2)/cbar

  chi2h=0.0

  do i=1,12,1
    chi2h= chi2h+((P(2)*myhubble(P,sd(1,i+sn_num))-sd(2,i+sn_num))**2)/(sd(3,i+sn_num)**2)
  end do

  chi2= chi2sn + chi2h + chi2cmb(P) + chi2bao(P)
  
!  write(*,*) chi2cmb(params),chi2sn, chi2h, chi2bao(params)

return
end function
!--------------------------------------------------------------------------

  function chi2cmb(P)
  implicit none

!  integer,parameter :: num_params=9 !P = [omegam,h100,omegab,z1,z2,a0,a1,a2,a3]
  real(kind=8) :: chi2cmb,P(num_params),R,la,zs,cmbarray(3),g1,g2
  real(kind=8) :: cov(3,3)
  integer :: rr,cc
!  real,external :: RDl,rs
  
  data((cov(rr,cc),rr=1,3),cc=1,3)/ 1.8,27.968,-1.103,27.968,5667.577,-92.263,-1.103,-92.263,2.923 /

  g1= 0.0783*((P(3)*(P(2)**2))**(-0.238))/(1+39.5*((P(3)*(P(2)**2))**0.763))

  g2= 0.560/(1+21.1*((P(3)*(P(2)**2))**1.81))

  zs= 1048*(1+0.00124*((P(3)*(P(2)**2))**(-0.738)))*(1+g1*((P(1)*(P(2)**2))**g2))

  R=sqrt(P(1))*myRDl(P,zs)

  la=3.1415926*myRDl(P,zs)/rs(P,zs)

  cmbarray(1)= la-302.1
  cmbarray(2)= R-1.710
  cmbarray(3)= zs-1090.04

!  write(*,*) rs(params,zs),la,cmbarray(:)

  chi2cmb= dot_product(cmbarray,matmul(cov,cmbarray))

  return
  end function

!---------------------

  function chi2bao(P)
  implicit none

!  integer,parameter :: num_params=9 !P = [omegam,h100,omegab,z1,z2,a0,a1,a2,a3]
  real(kind=8) :: chi2bao, P(num_params), bao1,bao2,zd,b1,b2 !bao1=Dv(0.35)/Dv(0.2); bao2=rs(zd)/Dv(0.275)
!  real,external :: Dl,rs,hubble

  b1= 0.313*((P(1)*P(2)*P(2))**(-0.419))*(1+0.607*((P(1)*P(2)*P(2))**0.674))

  b2= 0.238*((P(1)*P(2)*P(2))**0.223)

  zd=1291*((P(1)*(P(2)**2))**0.251)*(1+b1*((P(3)*(P(2)**2))**b2))/(1+0.659*((P(1)*(P(2)**2))**0.828))

  bao1= (((0.35/myhubble(P,3.5D-1))**(1.0/3.0))*(myDl(P,3.5D-1)**(2.0/3.0)))/&
&(((0.2/myhubble(P,2.0D-1))**(1.0/3.0))*(myDl(P,2.0D-1)**(2.0/3.0)))

  bao2= rs(P,zd)/(((0.275/myhubble(P,2.75D-1))**(1.0/3.0))*(myDl(P,2.75D-1)**(2.0/3.0)))

  chi2bao= ((bao1-1.736)**2)/(0.065**2)+((bao2-0.139)**2)/(0.0037**2)

!  write(*,*) bao1,bao2,rs(params,zd)

  return
  end function


!---------------------

  function rs(P,z)  ! comoving sound horizon * H0
  implicit none

  real(kind=8) :: rs, P(num_params),z, cs,a
!  real,external :: hubble
  real(kind=8) :: zmin=1.0D-10,zmax
  integer,parameter :: seg=100 ! seg=1E4,rs=6.57E-2,la=159; seg=1E5,rs=6.15E-2,la=170
  real(kind=8) :: datalist(seg+1)
  real(kind=8) :: width,s
  integer :: i

!  cs=sqrt(3+9*a*params(3)*params(2)*params(2)/(4*2.469E-5))

! use trape method, a bad method in this case--delete this note

  s=0.0
  zmax=1.0/(1.0+z)
  width=(zmax-zmin)/seg
  do i=0,seg,1
     a=zmin+i*width
     datalist(i+1)=1.0/(a*a*myhubble(P,1/a-1)*sqrt(3+9*a*P(3)*P(2)*P(2)/(4*2.469D-5)))
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  rs=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)

  return
  end function


!---------------------


function myDl(P,z)   ! /int_{0}^{z} dt/H(t) lum-distance=(1+z)Dl
! this part is to do trape method of integration
! read notes/../integration/trape.f90 for detail of this method
implicit none

  real(kind=8) :: myDl,z
    real(kind=8) :: P(num_params)
!!!  real(kind=8),external :: hubble
  real(kind=8) :: xmin=0.0,xmax 
  integer,parameter :: seg=100  
  real(kind=8) :: datalist(seg+1)   
  real(kind=8) :: width  
  real(kind=8) :: s
  integer :: i

  s=0.0
  xmax=z
  width=(xmax-xmin)/seg
  do i=0,seg,1
     datalist(i+1)=1.0/myhubble(P,xmin+i*width)
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  myDl=width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  

!return
end function

!------------------------------------------------------------------------

function myRDl(P,z)   ! /int_{0}^{z} dt/H(t) lum-distance=(1+z)Dl
! this part is to do trape method of integration
! read notes/../integration/trape.f90 for detail of this method
implicit none

  real(kind=8) :: myRDl,z
    real(kind=8) :: P(num_params)
!!!  real(kind=8),external :: hubble
  real(kind=8) :: xmin=0.0,xmax 
  integer,parameter ::seg= 10000, seg1=500
  real(kind=8),parameter ::  xmid=5.0
  real(kind=8) :: datalist(seg+1), datalist1(seg1+1)   
  real(kind=8) :: width, width1
  real(kind=8) :: s,s1
  integer :: i


  s=0.0
  s1=0.0
  xmax=z
  width1=(xmid-xmin)/seg1
  width=(xmax-xmid)/seg

  do i=0,seg1,1
     datalist1(i+1)=1.0/myhubble(P,xmin+i*width1)
  end do

  do i=0,seg,1
     datalist(i+1)=1.0/myhubble(P,xmid+i*width)
  end do

  do i=0,seg1,1
     s1=s1+datalist1(i+1)
  end do

  do i=0,seg,1
     s=s+datalist(i+1)
  end do

  myRDl=width1*(s1-datalist1(1)/2.0-datalist1(seg1+1)/2.0)+width*(s-datalist(1)/2.0-datalist(seg+1)/2.0)
  

!return
end function


!------------------------------------------------------------------------


!------------------------------------------!
!------------Hubble Parameter--------------!
!------------------------------------------!

function myhubble(P,z)
implicit none

real(kind=8) :: myhubble, z, P(num_params)
!real,external :: Intergal1,Intergal2,Intergal3,Intergal4
real(kind=8),parameter :: z3=1.8, a=-1.0

! P = [omegam,h100,omegab,z1,z2,a0,a1,a2,a3]

if ( z < P(4) ) then

  myhubble = sqrt(P(1)*((1+z)**3) + 4.15E-5*((1+z)**4)/(P(2)*P(2)) +(1- P(1))*Intergal1(P(4),P(6),P(7),z))

else if ( z>=P(4) .and. z<=P(5) ) then

  myhubble = sqrt(P(1)*((1+z)**3) + 4.15E-5*((1+z)**4)/(P(2)*P(2)) + (1- P(1))*Intergal2(P(4),P(5),P(6),P(7),P(8),z))

else if ( z>P(5) .and. z<=z3 ) then
  
  myhubble = sqrt(P(1)*((1+z)**3)+ 4.15E-5*((1+z)**4)/(P(2)*P(2)) + (1- P(1))*Intergal3(P(4),P(5),P(6),P(7),P(8),P(9),z3,z))

else if ( z>z3 ) then
  
  myhubble =sqrt(P(1)*((1+z)**3)+ 4.15E-5*((1+z)**4)/(P(2)*P(2)) +(1- P(1))*Intergal4(P(4),P(5),P(6),P(7),P(8),P(9),z3,a,z))

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

function Intergal2(z1,z2,a0,a1,a2,z)
implicit none

real(kind=8) :: Intergal2,z1,z2,a0,a1,a2,z
!real,external :: Intergal1

Intergal2 = Intergal1(z1,a0,a1,z1)*(((1+z)/(1+z1))**(3*(1+a1-(1+z1)*(a2-a1)/(z2-z1))))*exp(3*(z-z1)*(a2-a1)/(z2-z1))

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal3(z1,z2,a0,a1,a2,a3,z3,z)
implicit none

real(kind=8) :: Intergal3,z1,z2,a0,a1,a2,a3,z3,z
!real,external :: Intergal2

Intergal3=Intergal2(z1,z2,a0,a1,a2,z2)*(((1+z)/(1+z2))**(3*(1+a2-(1+z2)*(a3-a2)/(z3-z2))))*exp(3*(z-z2)*(a3-a2)/(z3-z2))

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal4(z1,z2,a0,a1,a2,a3,z3,a,z)
implicit none

real(kind=8) :: Intergal4,z1,z2,a0,a1,a2,a3,z3,a,z
!real,external :: Intergral3

Intergal4 = Intergal3(z1,z2,a0,a1,a2,a3,z3,z3)*(((1+z)/(1+z3))**(3*(1+a)))

return
end function


!-----------------------------------------------------------------------

function mymuth(P,z) ! /mu_th = 5*lg(luminositydistance)+0[not /mu_0]
implicit none

!!!  real(kind=8),external :: Dl
  real(kind=8) :: mymuth
    real(kind=8) :: P(num_params) 
  real(kind=8) :: z

  mymuth=5*log10(1+z) + 5*log10(myDl(P,z))

!return
end function

!-----------------------------------------------------------------------!
!---------------------------added code ends-----------------------------!
!-----------------------------------------------------------------------!
  
  function GetLogPrior(CMB, Info) !Get -Ln(Prior)
    real GetLogPrior
    real Age
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info

    GetLogPrior = logZero
 
    if (.not. generic_mcmc) then
     if (CMB%H0 < H0_min .or. CMB%H0 > H0_max) return
     if (CMB%omk < Omk_min .or. CMB%omk > Omk_max .or. CMB%Omv < 0) return
     if (CMB%zre < Use_min_zre) return

     Age = GetAge(CMB, Info)
      !This sets up parameters in CAMB, so do not delete unless you are careful!
 
     if (Use_Age_Tophat_Prior .and. (Age < Age_min .or. Age > Age_max) .or. Age < 0) return
    
    end if
    GetLogPrior = 0
 
  end function GetLogPrior

  function GetLogLike(Params) !Get -Ln(Likelihood)
    type(ParamSet)  Params 
    Type (CMBParams) CMB
    real GetLogLike
    real dum(1,1)
 
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
        return
    end if

    if (generic_mcmc) then
        GetLogLike = GenericLikelihoodFunction(Params) 
        if (GetLogLike /= LogZero) GetLogLike = GetLogLike/Temperature

    else

     call ParamsToCMBParams(Params%P,CMB)

     GetLogLike = GetLogLikePost(CMB, Params%Info,dum,.false.)
    end if 
   end function GetLogLike

    
  function GetLogLikePost(CMB, Info, inCls, HasCls) 
    real GetLogLikePost
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info
    real, intent(in):: inCls(:,:)
    logical, intent(in) :: HasCls
    real acl(lmax,num_cls)
    integer error

    if (generic_mcmc) stop 'GetLogLikePost: not supported for generic'


    GetLogLikePost  = GetLogPrior(CMB, Info)
    if ( GetLogLikePost >= logZero) then
       GetLogLikePost = logZero
       
    else 
          
       GetLogLikePost = GetLogLikePost + sum(CMB%nuisance(1:nuisance_params_used)**2)/2
          !Unit Gaussian prior on all nuisance parameters
       if (Use_HST) GetLogLikePost = GetLogLikePost + (CMB%H0 - 72)**2/(2*8**2)  !HST 
       if (Use_BBN) GetLogLikePost = GetLogLikePost + (CMB%ombh2 - 0.022)**2/(2*0.002**2) 
          !I'm using increased error bars here
   
       if (Use_CMB .or. Use_LSS) then
          if (HasCls) then
           acl = inCls
           error =0
          else
           call GetCls(CMB, Info, acl, error)
          end if
         if (error /= 0) then
          GetLogLikePost = logZero 
         else
          if (Use_CMB) GetLogLikePost = CMBLnLike(acl, CMB%norm(norm_SZ),CMB%nuisance) + GetLogLikePost
          if (Use_mpk) GetLogLikePost = GetLogLikePost + LSSLnLike(CMB, Info%theory)
          if (Use_WeakLen) GetLogLikePost = GetLogLikePost + WeakLenLnLike(CMB, Info%theory)     
          if (Use_Lya) GetLogLikePost = GetLogLikePost +  LSS_Lyalike(CMB, Info%Theory)
          if ( GetLogLikePost >= logZero) then
            GetLogLikePost = logZero
          end if
         end if
         if (Use_SN .and. GetLogLikePost /= logZero ) then
            if (Info%Theory%SN_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%SN_loglike
            else
             GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
  
  
         end if
   
       else
         if (Use_SN) GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
       end if

      if (Use_Clusters .and. GetLogLikePost /= LogZero) then
          GetLogLikePost = GetLogLikePost + &
                 (Info%Theory%Sigma_8-0.9)**2/(2*0.05**2)
          stop 'Write your cluster prior in calclike.f90 first!'
      end if

     if (GetLogLikePost /= LogZero) GetLogLikePost = GetLogLikePost/Temperature
   

    end if

  end function GetLogLikePost

end module CalcLike
