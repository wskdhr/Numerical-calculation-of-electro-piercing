title

 'Optimization of Piercing Device'
COORDINATES
YCYLINDER('r','z')

SELECT
ERRLIM=1e-3
TERRLIM=1e-3
SMOOTHINIT=ON
stages=7

variables

 Phi1

definitions
   R1= 600e-6{the outer radius of the egg}
   R2= R1-15e-6{pvs radius}
   R3= R2-50e-6{yolk radius}
!  Rhai=(R2+R3)/2
   Rn= 18.5e-6{needle radius}
   theta=20*PI/180 {needle angle}
    dn= 10e-6! 600e-6{distance between needle and chorion}
   dp=10e-6!  600e-6{distance between plate electrode and chorion}
   length=stage*0.25e-3
   Ztop=R1+length{upper electrode location}
   Zbtm= -R1-dp {lower electrode location}
   EPS0=8.854e-12    {vacuum permittivity}
    eps=78*EPS0  {external relative permittivity}
   Rmax=2.5e-3{cylinder radius}
!Rmax=3e-3 {cylinder radius}
 !R_isulation=0.25e-3*stage

R_isulation=1e-3
Phi0=150
equations

 div(eps*grad(Phi1)) = 0       

boundaries

 REGION 0 'box'

    START(0,R1+dn)
    LINE TO(Rn,R1+dn+Rn/tan(theta))
    LINE TO(Rn,Ztop)
    LINE TO(Rmax,Ztop)
    LINE TO(Rmax,Zbtm)
    LINE TO(0,Zbtm)
    LINE TO CLOSE

 REGION 1 'extracell'


    START(0,R1+dn)
    VALUE(Phi1)=150    LINE TO(Rn,R1+dn+Rn/tan(theta))
    VALUE(Phi1)=150    LINE TO(Rn,Ztop)
    NATURAL(Phi1)=0    LINE TO(Rn+R_isulation,Ztop)
    Value(Phi1)=-150    LINE TO(Rmax,Ztop)
    NATURAL(Phi1)=0    LINE TO(Rmax,Zbtm)
    Natural(Phi1)=0  LINE TO(0,Zbtm)
    NATURAL(Phi1)=0        LINE TO(0,-R1)
    ARC(CENTER=0,0)  TO (R1,0)
    ARC(CENTER=0,0)  TO (0,R1)
    NATURAL(Phi1)=0
    LINE TO CLOSE

 REGION 2 'chorion'
eps=30*EPS0
    START(0,-R1)
    ARC(CENTER=0,0)  TO (R1,0)
    ARC(CENTER=0,0)  TO (0,R1)
    NATURAL(Phi1)=0         LINE TO CLOSE



     LINE TO CLOSE
Feature 'insu'
 start(Rn,Ztop) line to (Rn+R_isulation,Ztop)
Feature 'negative'
 start(Rn+R_isulation,Ztop) line to (Rmax,Ztop)
!
!FEATURE 'degree_0_30'
!    START(0,Rhai)  ARC(CENTER=0,0)  to (0.5*Rhai,sqrt(0.75)*Rhai)
!FEATURE 'degree_30_60'
!    START(0.5*Rhai,sqrt(0.75)*Rhai)  ARC(CENTER=0,0)  to (sqrt(0.75)*(R2+R3)/2,0.25*(R2+R3))
!FEATURE 'degree_60_90'
 !   START(sqrt(0.75)*(R2+R3)/2,0.25*(R2+R3))  ARC(CENTER=0,0)  to (0.5*(R2+R3),0) 
!FEATURE 'degree_90_120'
!    START(0.5*(R2+R3),0)  ARC(CENTER=0,0)  to (sqrt(0.75)*(R2+R3)/2,-0.25*(R2+R3))
!FEATURE 'degree_120_150'
!    START(sqrt(0.75)*(R2+R3)/2,-0.25*(R2+R3))  ARC(CENTER=0,0)  to (0.25*(R2+R3),-sqrt(0.75)*(R2+R3)/2)
!FEATURE 'degree_150_180'
!    START(0.25*(R2+R3),-sqrt(0.75)*(R2+R3)/2)  ARC(CENTER=0,0)  to (0,-0.5*(R2+R3))



FEATURE 'membrane c'
    START(0,-R1)  ARC(CENTER=0,0)  TO (R1,0)
    ARC(CENTER=0,0)  TO (0,R1)

monitors

contour(Phi1) as 'Potential'

plots


surface(Phi1) as 'potential'
contour(Phi1) as 'potential'
contour(log10(sqrt(dz(Phi1))*dz(Phi1)+dr(Phi1)*dr(Phi1))) as 'Electric field strength'

contour(dz(Phi1)) zoom (0,R1-100e-6,105e-6,105e-6) as 'Electric field strength'

vector(-dz(Phi1),-dr(Phi1)) as 'Electric Field'

histories
!history(Line_INTEGRAL(dz(Phi1), 'degree_0_30'),Line_INTEGRAL(dz(Phi1), 'degree_30_60'),Line_INTEGRAL(dz(Phi1), 'degree_60_90'),Line_INTEGRAL(dz(Phi1), 'degree_90_120'),Line_INTEGRAL(dz(Phi1), 'degree_120_150'),Line_INTEGRAL(dz(Phi1), 'degree_150_180'))
history(dz(Phi1)) at (0,R1-7.5e-6) vs length as "Electric Field Strength" export format "#1" file="length of needle with E.txt"
end

