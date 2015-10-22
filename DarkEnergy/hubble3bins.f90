!-----------------------------------------------------------------------!
!-----------------------------Hubble Parameter--------------------------!
!-----------------------------------------------------------------------!

function hubble(P,z)
implicit none

real :: hubble, z, P(9) !P(num_params)
real,external :: Intergal1,Intergal2,Intergal3，Intergal4
real,parameter :: z3=1.8, a=-1.0

! P = [omegam,h100,omegab,z1,z2,a0,a1,a2,a3]

if ( z < P(4) ) then

  hubble = sqrt(P(1)*((1+z)**3) + 4.15E-5*((1+z)**4)/(P(2)*P(2)) +(1- P(1))*Intergal1(P(4),P(6),P(7),z))

else if ( z>=P(4) .and. z<=P(5) ) then

  hubble = sqrt(P(1)*((1+z)**3) + 4.15E-5*((1+z)**4)/(P(2)*P(2)) + (1- P(1))*Intergal2(P(4),P(5),P(6),P(7),P(8),z))

else if ( z>=P(5) .and. z<=z3 ) then
  
  hubble = sqrt(P(1)*((1+z)**3)+ 4.15E-5*((1+z)**4)/(P(2)*P(2)) + (1- P(1))*Intergal3(P(4),P(5),P(6),P(7),P(8),P(9),z3,z))

else if ( z>z3 ) then
  
  hubble =sqrt(P(1)*((1+z)**3)+ 4.15E-5*((1+z)**4)/(P(2)*P(2)) +(1- P(1))*Intergal4(P(4),P(5),P(6),P(7),P(8),P(9),z3,a,z))

end if

return
end function


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal1(z1,a0,a1,z)
implicit none

real :: Intergal1,z1,a0,a1,z

Intergal1 = ((1+z)**(3*(1+a0-(a1-a0)/z1)))*exp(3*z*(a1-a0)/z1)

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal2(z1,z2,a0,a1,a2,z)
implicit none

real :: Intergal2,z1,z2,a0,a1,a2,z
real,external :: Intergal1

Intergal2 = Intergal1(z1,a0,a1,z1)*(((1+z)/(1+z1))**(3*(1+a1-(1+z1)*(a2-a1)/(z2-z1))))*exp(3*(z-z1)*(a2-a1)/(z2-z1))

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergal3(z1,z2,a0,a1,a2,a3,z3,z)
implicit none

real :: Intergal3,z1,z2,a0,a1,a2,a3,z3,z
real,external :: Intergal2

Intergal3=Intergal2(z1,z2,a0,a1,a2,z2)*(((1+z)/(1+z2))**(3*(1+a2-(1+z2)*(a3-a2)/(z3-z2))))*exp(3*(z-z2)*(a3-a2)/(z3-z2))

return
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function Intergral4(z1,z2,a0,a1,a2,a3,z3,a,z)
implicit none

real :: Intergral4,z1,z2,a0,a1,a2,a3,z3,a,z
real,external :: Intergral3

Intergral4 = Intergral3(z1,z2,a0,a1,a2,a3,z3,z3)*(((1+z)/(1+z3))**(3*(1+a)))

return
end function
