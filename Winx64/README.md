T-matrix-abs-host
FORTRAN T-matrix code for particles in absorbing media, dielectric hosts are available
Program       qCQ-3  (NGauss<=192)

This code was developed by Nikolai G. Khlebtsov during 1987-2022.
The formulation for particles in absorbing media is given in JQRT paper

Nikolai G. Khlebtsov,  Extinction and scattering of light by nonspherical particles in absorbing media
JQSRT, 2022, https://doi.org/10.1016/j.jqsrt.2022.108069

 I would highly appreciate the proper citation and informing me of any problems encountered 
 with this code.  Please send your message to the following         
 e-mail address:  khlebtsov@ibppm.ru                              

 WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
 IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
 INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
 VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
 THE AUTHOR AND MY ORGANIZATION DISCLAIM ALL LIABILITY FOR
 ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM.  

This code is an extended precision version of CQ-3.for 
with double precision T-matrix calculations
Calculation of the light scattering cross sections Qext,sca,abs
for randomly oriented MONODISPERSE particles (rods and disks)
In addition, the orientation averaged Scattering matrix Fij(beta_sca)
and derivative parameters Ivv, Itot, Ivh/Ivv are calculated at a fixed wavelength
for details see Khlebtsov, JQSRT, 2022.	
In addition to orientation-averaged cross sections,
Qext,sca,abs(theta) , and the forward scattering amplitudes
A11(0)=A1, A22(0)=A2   are calculated for TE and TM
polarizations of the incident wave with respect to the plane
(a_sym;k) where a_sym is the symmetry axis, k is the inc. wave vector
for different orientation angles theta : cos (theta)=(a*k)/|a||k|
The Mueller matrix is calculated in (S) and (L) representations
The Stokes parameters are defined by     
(Il, Ir, U, V)   (U, V)   Bohren and Huffman book  
Il=El(El*), Ir=Er(Er*), Q=ElEl* + ErEr*;
U=Re[El(Er*)+Er(El*)  V=Im[El(Er*)-Er(El*)];
where indexes r,l correspond to perpendicular
and parallel directions with respect to the reference plane (Hulst)
 Vectors :(k;s)=Er x El is defined as in Bohren and Huffman book.
The signs of final amplitudes Aij and A1,A2 are changed according to 
Mishchenko (2000) book definition
complex AL = (2*pi/lambda(in vacuum)*Rev*refmed     |
Four particle shapes are considered:
 Nshape=1  spheroid: length aL=2*a, thickness d=2*b
semiaxis (a,b,b), a is rotation axis
eqivolume radius = b*(a/b)**(1/3)
axis ratio e=a/b (<1 or >1), parameter epsiln=b/a (>1 or <1)

Nshape=2  right circular cylinder or disk
length aL= 2a, and thickness d=2*b (rotation axis), 
axis ratio e=a/b (<1 or >1), epsiln = (b/a >1 or <1
eqivolume radius = b*[(3a/2b)]**(1/3)

Nshape=3  PROLATE s-cylinder with SEMI-SPHERICAL ends,
total length aL= 2a (rotation axis),total thickness d=2*b
Equivolume radius = b*(1+3(a-b)/2b)**1/3,
the aspect ratio e=a/b>1, the parameter epsiln = b/a<1

OBLATE s-disk with thickness L=2a (rotation axis)
radius of rotation is  b, the disk diameter is d=2*b
the equivolume radius = a[1+(3/2)*((b-a)/a)**2+(3*pi/4)*((b-a)/a)]**1/3
the aspect ratio e=a/b<1, the parameter epsiln=b/a>1
Nshape=4  Dog-bone rods with elliptical caps
Four parameters have to be defined:
The total length L=L_aver
the minimal central diameter d=2*b
the maximal end-diameter d1=d*(1+par_hi)  
the end-cap height b_c=(d1/2)*par_c
The averaged aspect ratio equals e_aver = L_aver/d_aver,
where d_aver=(d+d1)/2  is the AVERAGED diameter
For calculations we have used also
epsiln=1/e e=q/b=L/d d=2b is the CENTRAL diameter
For details, see Khlebtsov et al.,J. Phys. Chem. C 2011, 115, 6317–6323
NOTES
NOTE1: par_hi must be less than 0.5
NOTE2: if Nshape=3 par_c and par_hi will be used incorrectly
Use Nshape=4 (!) for par_c<=1 and par_hi<=0.5 !!!!!!!!
FNOTE 3: in this version, Nshape=4 works only for dog-bone rods
