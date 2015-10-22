 program main
 implicit none

 real,external :: chi2h,hubble
 real :: a(2)

 a(1)=0.274
 a(2)=0.72
 
 write(*,*) hubble(0.0,0.0),hubble(1.0,3.0)
 write(*,*) chi2h(a)

 end program

!-----------------------------

  function chi2h(params)
  implicit none

  integer,parameter :: num_params=2
  real :: chi2h, params(num_params)  !params=(omegam0,h100)
  character(80) :: hubble12="hubble12.txt"
  real :: hub12(3,12)
  integer :: i
  real,external :: hubble
  logical :: alive

  write(*,*) "reading hubble data"

  inquire (file=hubble12, exist=alive)
  if (.not.alive) then
       write(*,*) hubble12,"doesn't exist"
       stop
    else 
       write(*,*) "data file exists"
    end if

  open (unit=77, file=hubble12)
  do i=1,12,1
     read(77,*) hub12(:,i)
     write(*,*) hub12(:,i)
  end do
  close(77)
  
  chi2h=0.0

  do i=1,12,1
    chi2h= chi2h+((params(2)*hubble(params,hub12(1,i))-hub12(2,i))**2)/(hub12(3,i)**2)
  end do

  end function

!-------------------------

function hubble(params,z) ! hubble parameter in LCDM model
implicit none

  real(kind=4) :: hubble,params(2),z

  hubble=sqrt(params(1)*((1+z)**3) + 4.15E-5*((1+z)**4)+ 1- params(1))

return
end function

