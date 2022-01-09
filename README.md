# T-matrix-abs-host
FORTRAN T-matrix code for particles in absorbing media, dielectric hosts are available
C|                 Program       qCQ-3  (NGauss<=192)            |
C|                                                                  |
C|  This code is an extended precision version of CQ-3.for 
c   with double precision T-matrix calculations
c
C   Calculation of the light scattering cross sections Qext,sca,abs
c   for randomly oriented MONODISPERSE particles (rods and disks)
c   In addition, the orientation averaged Scattering matrix Fij(beta_sca)
c   and derivative parameters Ivv, Itot, Ivh/Ivv are calculated at a wixed wavelength
c   for Deatails see Khlebtsov, JQSRT, 2022.	
c   In addition to orientation-averaged cross sections,
c   Qext,sca,abs(theta) , and the forward scattering amplitudes
c    A11(0)=A1, A22(0)=A2   are calculated for TE and TM
c   polarizations of the incident wave with respect to the plane
c   (a_sym;k) where a_sym is the symmetry axis, k is the inc. wave vector
c   for different orientation angles theta : cos (theta)=(a*k)/|a||k|
c   
c   The Mueller matrix is calculated in (S) and (L) representations
c
c   The Stokes parameters are defined by     
C|  (Il,Ir,U,V)   (U,V)   Bohren and Huffman book  
C|  Il=El(El*),Ir=Er(Er*),Q=ElEl* + ErEr*;U=Re[El(Er*)+Er(El*)  |
C|  V=Im[El(Er*)-Er(El*)]; where indexes r,l correspond to perpendicular
c   and parallel directions with respect to the reference palne (Hulst)
c   Vectors :(k;s)=ErxEl is defined as in Bohren and Huffman book.
c 
C  The signs of final amplitudes Aij and A1,A2 are changed according to 
c  Mishchenko (2000) book definition
C|==================================================================|
C|  complex AL = (2*pi/lambda(in vacuum)*Rev*refmed                 |
c***********************************************************
c      Four particle shapes are considered:
c  Nshape=1  spheroid: length aL=2*a, thickness d=2*b
c            semiaxis (a,b,b), a is rotation axis
c		   eqivolume radius = b*(a/b)**(1/3)
c            axis ratio e=a/b (<1 or >1), parameter epsiln=b/a (>1 or <1)
c
c  Nshape=2  right circular cylinder or disk
c            length aL= 2a, and thickness d=2*b (rotation axis), 
c            axis ratio e=a/b (<1 or >1), epsiln = (b/a >1 or <1
c		   eqivolume radius = b*[(3a/2b)]**(1/3)
c
c  Nshape=3  PROLATE s-cylinder with SEMI-SPHERICAL ends,
c            total length aL= 2a (rotation axis),total thickness d=2*b
c            Equivolume radius = b*(1+3(a-b)/2b)**1/3,
c            the aspect ratio ratio e=a/b>1, the parameter epsiln = b/a<1
c
c            OBLATE s-disk with thicknes L=2a (rotation axis)
c            radius of rotation is  b, the disk diamater is d=2*b
c            the equivolume radius = a[1+(3/2)*((b-a)/a)**2+(3*pi/4)*((b-a)/a)]**1/3
c            the aspect ratio e=a/b<1, the parameter epsiln=b/a>1
c  Nshape=4  Dog-bone rods with elliptical caps
c		   Four parameters have to be defeined:
c		   The total length L=L_aver
c		   the minimal central diameter d=2*b
c		   the maximal end-diameter d1=d*(1+par_hi)  
c		   the end-cap height b_c=(d1/2)*par_c
c			The averaged aspect ratio equals e_aver = L_aver/d_aver,
c             where d_aver=(d+d1)/2  is the AVERAGED diameter
c			For calculations we have used also
c             epsiln=1/e e=q/b=L/d d=2b is the CENTRAL diameter
c	        For details, see Khlebtsov et al.,J. Phys. Chem. C 2011, 115, 6317–6323
c	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c             NOTE1: par_hi must be less than 0.5
c			NOTE2: if Nshape=3 par_c and par_hi will be used incorrectly
c			Use Nshape=4 (!) for par_c<=1 and par_hi<=0.5 !!!!!!!!	
c
c			FNOTE 3: in this version, Nshape=4 works only for dog-bone rods	
