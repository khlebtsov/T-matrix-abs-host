T-matrix-abs-host
C                                                                  
C                 Program       qCQ-4A.for   (NGauss<=192)           
c****
C	In this version, the complex (constant) refractive indexes of medium (n1) and particles (n2)
c	can be defined in the parameter file QC-4A.par                                                            
c*****
C  This code is an extended precision version of CQ-4A.for 
c   with double precision T-matrix calculations
c  This code was developed by Nikolai G. Khlebtsov during 1987-2022.
c  The formulation for particles in absorbing media is given in JQSRT paper
c  Nikolai G. Khlebtsov,  Extinction and scattering of light by nonspherical particles in absorbing media
c  JQSRT, 2022, https://doi.org/10.1016/j.jqsrt.2022.108069
c
c   I would highly appreciate the proper citation and informing me of any problems encountered 
c   with this code.  Please send your message to the following         
c   e-mail address:  khlebtsov@ibppm.ru                              
c
c  WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
c  IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
c  INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
c  VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
c  THE AUTHOR AND MY ORGANIZATION DISCLAIM ALL LIABILITY FOR
c  ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM.  
c
C   Calculation of the light scattering cross sections Qext,sca,abs
c   for randomly oriented MONODISPERSE particles (rods and disks)
c   In addition, the orientation averaged Scattering matrix Fij(beta_sca)
c   and derivative parameters Ivv, Itot, Ivh/Ivv are calculated at a wixed wavelength
c  
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
C|===================================================|
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
c            OBLATE s-disks with circular edges . 
c		   The disk thicknes is L=2a (the rotation axis)
c            the radius of rotation is  b, the disk diamater is d=2*b
c            the equivolume radius = a[1+(3/2)*((b-a)/a)**2+(3*pi/4)*((b-a)/a)]**1/3
c            the aspect ratio e=a/b<1, the parameter epsiln=b/a>1
c  Nshape=4  Dog-bone rods with elliptical caps and disks with elliptical edges
c			
c		Dog-bone rods: e=Lav/d_av>1
c			Four parameters have to be defeined:
c		   The total length L=L_aver
c		   the minimal central diameter d=2*b
c		   the maximal end-diameter d1=d*(1+par_hi)  
c		   the end-cap height b_c=(d1/2)*par_c
c		   The averaged aspect ratio equals e_aver = L_aver/d_aver,
c		   where d_aver=(d+d1)/2  is the AVERAGED diameter
c		   For calculations we have used also
c		   epsiln=1/e e=q/b=L/d d=2b is the CENTRAL diameter
c		   For details, see Khlebtsov et al.,J. Phys. Chem. C 2011, 115, 6317–6323
c		!!!!!!!!
c             NOTE1: par_hi must be less than 0.5
c			NOTE2: if Nshape=3 par_c and par_hi will be used incorrectly
c			Use Nshape=4 (!) for par_c<=1 and par_hi<=0.5 !!!!!!!!	
c
c			Disks with elliptical edges
c		ec=par_c. If par_c=1 we have a disk with the circular edge, Nshape=3, L/d<1.
c				  If par_c=0.001, we have the right circular disk,  Nshape=2, L/d<1 
c		The equivoloume radius = a*[ec**2+(3/2)(e-ec)**2+(3*pi/4)*ec*(e-ec)]**(1./3.)
c		ec=b_c/a, e=b/a, b=d/2 is the rotation radius, 2a=L is the disk thickness
c		NOTE!!!: below the symbol 'e' stands for the aspect ratio e=eps=b/a>1
c		Y=X*RT X=k1*a  FW1=d(RT)/d(theta)/RT/sin(theta)
c
c		if z=cos(theta)>=z1=a/sqrt(a**2+b**2) then
c		RT=1/z FW1=T/sz=RT, sz=sin(theta)
c
c		if z<z1 RT=p*[1+det] FW1=(1/sz/p)*dp/d(theta)-[1/2/sz/(1-q-det)]*dq/d(theta)
c		det=dsqrt(1.-q), q=q1*[1+ec**2*z*z/sz/sz], q1=(e*e-2*e*ec)/(e-ec)**2
c		dp/d(theta)=dpdt=(e-ec)*z/p1*[1-2*sz*sz(1-ec**2)/p1]
c		p1=1+(ec**1-1)*z*z
c 		
c		The input parameters block:
      read(14,*) aLav,dav
c	The average length and diameter of partciles
c	aLav>dav for rods (prolate particles) and aLav<dav for disks (oblate particles)
	read(14,*) refm1,refm2,NRefm
c	Real and Imag part of refmed if Nrefm=0 (the refractive index of an external medium) 
c	If Nrefm=1,2,3 the external medium is  P3HT, PMMA, and water, respectively 
c     For water we use the following dispersion formula (Khlebtsov, Nanoscale 2020)
c     W2=WAVEL*WAVEL  REFMED=1.32334d0+3.479d+03/W2-5.111d+07/W2/W2
c	For a specific external medium, replace P3HT or PMMA subroutine with your subroutine MYrefm
c	for example, you can define NRefm=2 and replace 
C     		call PMMA(wavel,refmed,epsmed)
c	with
c			call MYrefm(wavel,refmed)
c    here refmed is the complex*32 variable; epsmed=refmed**2 and it is not used here
c
	read(14,*) crm2,cim2,Nm2
c	The  complex absolute refractive index of the particle m2=(crm2,cim2);
c	To use this optio, set Nm2=1. Otherwise, set Nm2=0
c
	read(14,*) KD, KS, As
c	KD and As are parameters for size-coorection of the particle refreactive index
c	if KS=0 the TS set by Olmon et al. is used
c	if KS=1 the SC set by Olmon et al. is used
c	for details: see Khlebtsov et al. Nanoscale, 2020
	read(14,*) wavelbeg,wavelend,wavelstep
c	wavelength range and step
	read(14,*) NeedFij
c	if NeedFif=0 Fij is not calculated
      read(14,*) Nshape,Nmetal,NAu
c	Nshape=1,2,3,4 dedine four particle shapes
c	Nmetal=0 for Ag and 1 for Au
c	if NAu=1 or 2 the Olmon or g34JCH3 splines are used, respectively 
	read(14,*) par_c,par_hi
c	The shape parameters for Nshape=4 case
c     For details, see Khlebtsov et al.,J. Phys. Chem. C 2011, 115, 6317–6323
      read(14,*) Norient,Ntheta
c	Norient is the number of orientation points for averaging
c     FSrandom(i,j) FLrandom(i,j) 4x4 matrixes over ZU=cos(Theta_u) and HI
c     see Khlebtsov, JQSRT, 2022, https://doi.org/10.1016/j.jqsrt.2022.108069
c	Ntheta is the number of angles between the incident wave vector k and the symmetry axis
c	for calculating the integral cross-sections for TE and TM polarizations of the incident wave
	read(14,*) Bet_beg,Bet_end,Bet_step 
c	The range of scattering angles for calculating FSrandom(i,j) FLrandom(i,j) 4x4 matrixes
      read(14,*) N,Mmax
c	The T-matrix Tmn,m'n' orders n_max=N and m_max=Mmax for calculating the light scattering parameters
	read(14,*) ERR_Q, ERR_Aij
c	the tolerance error parameters for the efficiencies and the amplitude matrix
      read(14,*) NGauss
c	The number of Gauss points in the T-matrix quadratures
	NI=1
c	The number of subintervals in the T-matrix integrals. Here NI is set to be 1
C--------------------------------------------------------
