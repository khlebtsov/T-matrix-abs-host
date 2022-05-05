C                                                                  
C                 Program       qCQ-4A.for   (NGauss<=192)           
c****
C	A new subroutine qRTFW1 has been implemented in a previous qCQ-3A version                                                                 
c	A new cq-4A.par file is used instead of previous cq-4.par
c*****
C  This code is an extended precision version of CQ-4.for 
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

c************************************************************************	
c	Main program
c
	IMPLICIT real*16(A-H,O-Z)

	real*16 F(4,4),FS(4,4),FSran(4,4), FLran(4,4),ZUW(48),ZUX1(48) 
	real*16 MG(4,4),MH(4,4)
      complex*32 AL, refm

	COMPLEX*32 AI1,n2,eps2,refmed, epsmed,
     +           A1, A2, A11, A12, A21, A22


	COMPLEX*32 TE_0,TM_0,TE_m,TM_m
	
      Character f_line*120,outfile*120
c
 	common/T/ TE_0(1000),TM_0(1000),TE_m(74000),TM_m(74000)
c	TM dimension for M=0 (N/2)x(N/2)
c	Here for Nmax=60 Dim_TM0=30*30=900 is set 1000
c	TM dimension (M>0): 
c	N=100 Dim_TM=338350 N=60 Dim_TM=73810 N=50 Dim_TM=42925
c	Here Nmax=60 and Dim_TM is set 74000
c	NGauss (max) is set 192

c	for q programs
      character (9) ::today
      write (*,*)'Enter Output File name:'
      read (*,'(a)') outfile
      OPEN (UNIT=16, FILE=outfile)
C     Write program-name and date
      call date_and_time(today)
      write (16,*) today
      write(16,*) 'qCQ-4A.for'
C--------------------------------------------------------------------
C     SAVE Input Data file "***.par" in "outfile"
c     for control and arhive purposes
      OPEN(unit=14,file='CQ-4A.par',status='old')
 111   read(14,'(a)',end=222,err=222)f_line
c       write(16,'(a80)')f_line
       write(16,'(a)')f_line
       go to 111
 222  rewind(14)
C--------------------------------------------------------------------

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
	read(14,*) crn2,cin2,Nm2
c	The  complex absolute refractive index of the particle n2=(crn2,cin2);
c	If this option is used, set Nm2=1. Otherwise, set Nm2=0
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
c	Nshape=1,2,3,4 define four particle shapes
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
c	the toleranse error parameters for the efficiencies and the amplitude matrix
      read(14,*) NGauss
c	The number of Gauss points in the T-matrix quadratures
	NI=1
c	The number of subintervals in the T-matrix integrals. Here NI is set to be 1
C--------------------------------------------------------

	if (NRefm.gt.3) then
		write (16,*) 'Check NRefm', NRefm
		STOP
	end if

C--------------------------------------------------------
c	the average size parameters of rods
C     aLav and dav(central diameter)
	a=aLav/2.
	arod=a
	b_av=dav/2.
c---------------------------------------------------------
c	For Nshape=1-3  b=b_av=dav/2 
c	For Nshape=4 b=b_av/(1+par_hi/2)
c      and the aspect ratio e =a/b DEPENDS on the par_hi parameter!!! 
	b=b_av/(1.q0+par_hi/2.q0)
	e=a/b
C--------------------------------------------------------------   end
c	Parameter epsiln=1/e 
	epsiln=1.q0/e
C--------------------------------------------------------------   end
	AI1=(0.q0,1.q0)
      pi=qacos(-1.q0)
	CALL qGAU_short(Norient,ZUW,ZUX1)
C--------------------------------------------------------------   beg
		if (Nshape.eq.1) then
			a0=(a*b*b)**(1./3.)
			go to 1123
		end if

		if (Nshape.eq.2) then
			a0=b*(3.q0*a/2.q0/b)**(1./3.) 
			go to 1123
		end if

		if (Nshape.eq.3) then
		    if (e.ge.1.q0) then
c           Equivolume radius a0= b*(1+3(a-b)/2b)**1/3,
			   a0=b*(1.q0+3.q0*(a-b)/2.q0/b)**(1./3.)
	        else
c           Equivolume radius a0= a[1+(3/2)*((b-a)/a)**2+(3*pi/4)*((b-a))/a)]**1/3
	           a0=a*(1.q0+1.5d0*((b-a)/a)**2.q0+
     +	          (0.75d0*pi)*((b-a)/a))**(1.q0/3.q0)
			end if
			go to 1123
		 end if
c
c	For Nshape=4 
		if (Nshape.eq.4) then
	if (e.lt.1.) goto 1120
c		Case 1: dog-bone rods	
c			end-cap diameter d1=2*b1 b1=b*(1+par_hi)			
			b1=b*(1.q0+par_hi)
c			end-cap second size bc=par_c*b1
			bc=b1*par_c
c			rods length without end-caps: 2a=aL-2*bc		
			ac=arod-bc
c	Calcualtion of the equivolume radius:
c	Radius of a sphere for side rod curve
		p=ac/b
		zz=2.q0*p*par_hi/(p*p+par_hi*par_hi)
		xx=zz/p
		zp=zz*p
			if(zp.lt.0.001) then			
			ff=1.q0+(zp/3.q0)*(1.q0+(3./20.)*zz*zz)+zp*zp/20.q0
			else
				fi=1.q0-( qasin(zz)/zz+qsqrt((1.q0-zz)*(1.q0+zz)) )/2.q0
				ff=1.q0+(2.q0*(1.q0+xx)/xx/xx)*fi-p**2/3.q0
			end if

c           The equivolume radius a0= b[(b1/b)**3*c+(3/2)*p*f(x)]**(1/3)
c		f(x)=1 + 2*(1+x)/[x**2]*fi(z) - p**2/3
c		fi(z)=1 - (1/2)*[arcsin(z)/z+sqrt(1-z**2)]
c		x=b/R  p=a/b  z=p*x  R=[a**2+(b1-b)**2]/[2*(b1-b)]

		a0=b*((b1/b)**3.q0*par_c + 1.5d0*p*ff )**(1./3.)
		end if
c---------------------------------------------------------
c		The disk with elliptical edge
c		This shape incldudes the right OBLATE cylinders (e.g.,par_c=0.001)
c		and OBLATE s-cylinders if par_c=1.
1120	ec=par_c
	ee=1./e
c	e=a/b<1 eps=1./e=b/a
	eec2=(ee-ec)*(ee-ec)
	ec2=ec*ec
c		The equivoloume radius Rev= a0=a*[ec**2+(3/2)(e-ec)**2+(3*pi/4)*ec*(e-ec)]**(1./3.)
c		ec=b_c/a=par_c, e=b/a, b=d/2 is the rotation radius, 2a=L is the disk thickness
	sc_a=(ec2+1.5*eec2+0.75*pi*ec*(ee-ec))**(1./3.)
c	
	a0=a*sc_a
c---------------------------------------------------------

1123		continue
c---------------------------------------------------------
c****************
C	The effective electron scattering length is defined with a0=Rev and As parameter:
	aLeff=a0/As
c*****************	
C--------------------------------------------------------------
	write (16,*)
	WRITE (16,191) a0,e
191      FORMAT ('a0=Rev=',E11.4,' e=a/b=',E11.4) 
c-----------------------------------------------------------------------
C--------------------------------------------------------
            write(16,*)
            WRITE (16,6)
6           FORMAT (1H ,'Wavel',6x,'Qext',10X,'Qsca',9x,'Qabs')

	write(16,192)
192	format(1h ,'Wavel',6x,'Ivh/Ivv[90]',3x,'Ivv_norm',
     +       5x,'FS11_norm',2x,'sca_ang')

C--------------------------------------------------------
c  NOTE!!! c	If wavel_sca.ne.0 then only one wavel_sca will be used

      do 999 wavel=wavelbeg,wavelend,wavelstep
	write (*,765) wavel
765	format(1h , f7.2)
c-----------------------------------------------------------------------
C--------------------------------------------------------------------

	if(NRefm.eq.0) then
		refmed=dcmplx(refm1,refm2)
	end if

c	Complex refmed and epsmed for P3HT polymer
	if(NRefm.eq.1) then
		call P3HT2020DF(wavel,refmed,epsmed)
	end if 

c	Real refmed and epsmed for PMMA polymer
	if(NRefm.eq.2) then
		call PMMA(wavel, refmed, epsmed)
	end if

c	Real refmed and epsmed for water

	if(NRefm.eq.3) then
         W2=wavel*wavel
         refmed=1.32334d0+3.479d+03/W2-5.111d+07/W2/W2
	end if

C--------------------------------------------------------------------
c	Calculate the absolute and relative complex refractive index of the particle
c	If Nm2=1 Otherwise, skip these calcuations
	if(Nm2.eq.1) then
		n2=dcmplx(crn2,cin2)
		refm=n2/refmed
		goto 1111
	end if
C--------------------------------------------------------------------
C     Calculate the absolute refractive index n2 of Au or Ag particle            
C     with account for size limiting effect                         
c	Nmetal=1 Au
c	Nmetal=0 Ag
	if (Nmetal.eq.1) then 
		 if (NAu.eq.1) then
		    call Olmon(wavel,aLeff,n2,eps2,KD,KS)
			refm=n2/refmed
		 else
  			call g34JCh3(wavel,aLeff,n2,KD)
			refm=n2/refmed
		 end if
	else        
		 call Ag2021(wavel,aLeff,n2,eps2)
			refm=n2/refmed
	end if

C--------------------------------------------------------------------
1111	continue
C--------------------------------------------------------------------
	  AL0=2.q0*pi*a0/wavel	  
	  AL=AL0*REFMED
C--------------------------------------------------------------   end
c	Calculation of ALL T0 and Tm matrices for given N and Mmax
c	The particle parameters are specified by Al, EPSILN, Nshape, par_c,par_hi, and Refm
c	The otuput matrices are saved in TE_0,TM_0, TE_m and TM_m one-dimensional arrays
c
        call qT(AL,EPSILN,Nshape,par_c,par_hi,REFM,Mmax,N,NGauss)

c	T-matrix has been calculated
C--------------------------------------------------------------------
c	T-matrix has been calculated
c     PRINT first T-matrix elements  
c
c	write(16,*) 'TE_0 TM_0 TE_m TM_m i=1-4'
c 	do 333 i=1,2
c		write(16,544)i
c		write(16,543) TE_0(i)
c		write(16,543) TM_0(i)
c		write(16,543) TE_m(i)
c		write(16,543) TM_m(i)
c543	format(1h ,2(2x,e11.4))
c544	FORMAT (1h ,i2)
c333   end do
C--------------------------------------------------------------------

c	Calculation of Qext,sca,abs for random orientations
c     and for TE,TM polarizations with respect
c      to the plane (symmetry axis, k). Theta is the angle between
c	the symmetry axis and the plane wave vector k.
c	Define Theta angles for given Ntheta
       i_th=0
1122         if(Ntheta.eq.0) then
                 Theta=0.q0
              else
                 Theta=i_th*90.q0/Ntheta
              end if   

	Z=qcos(theta*pi/180.q0)
c	Calculate A1(0), A2(0) light scattering amplitudes for 
c	TE [E perp. (a,k)] and TM [E parallel (a,k)] for different
c	orientation angles theta : cos (theta)=(ak)/|a||k|
c	NOTE: M_max is actual M that ensures ERR_Q accuracy of calculations
	call qQesa(AL,Mmax,N,Z,A1,A2,Qext,Qsca,Qabs,
     +           QexE,QexM,QscE,QscM,QabE,QabM,ERR_Q,M_max)


c******************************************************************
c     TEST for the scattering amplitudes
c	Here: the test is excluded
C     Optical Theorem for the scattering amplitudes A1(0) and A2(0)
c	reads as follows:

c      Qext1=4.*Imag(A1/AL)/Real(AL)
c      Qext2=4.*Imag(A2/AL)/Real(AL)
c	write(16,876)wavel, Qext1, Qext2,A1,A2
c876	format (1x,f6.1,2(1x,e11.5),2(1x,E11.4),2(1x,E11.4))
c******************************************************************
c	Print the output quantities
      if (Ntheta.eq.0) then
	    WRITE (16,991) WAVEL,Qext, Qsca, Qabs
		goto 994
	end if
	if (theta.eq.0) WRITE (16,991) WAVEL,Qext,Qsca,Qabs             


          WRITE (16,992) WAVEL,QexE, QscE, QabE,Theta
          WRITE (16,993) WAVEL,QexM, QscM, QabM,Theta

991    format(/1x,f6.1,3(2x,e11.6),2x,'Random')
992    format(1x,f6.1,3(2x,e11.6),2x,f4.1,2x,'TE')
993    format(1x,f6.1,3(2x,e11.6),2x,f4.1,2x,'TM')

994		continue
	if (i_th.lt.Ntheta) then
	i_th=i_th+1
	goto 1122
	end if

c	End loop over Theta. All cross sections have been calculated
C--------------------------------------------------------------------
C--------------------------------------------------------------------
c	Calculate the average Mueller matrix only for wavel=wavel_sca
C	Output Ivh/Ivv[90], Itot_norm, Ivv_norm
c	The normalization coefficient anor=(8.q0/3.q0)/x**2
c	Itot_nor=Itot*anor  Ivv_nor=Ivv*anor

	if (NeedFij.eq.0) goto 999
C--------------------------------------------------------------------

c	write(16,192)
c192	format(2x,'Wavel',4x,'Ivh/Ivv[90]',3x,'Ivv_norm',
c     +       5x,'phaseF',6x,'sca_ang')

c-----------------------------------------------------------------------
c	WRITE (16,91) wavel,A0,E,REFM,Norient,N,MMAX,NGauss,NI
C|==================================================================|
C|  BETTA - the scattering angle, ZU is the cosine of the polar orienation angle
C|  of the symmetry axis,  HI is the azimuth in the incident wave refernce frame                             |
C|==================================================================|
c	start loop ober beta (scattering angle)
	DO 99 Bet=Bet_beg,Bet_end,Bet_step
		BETTA=Bet*PI/180.q0
		CBET=qcos(BETTA)
		SBET=qsin(BETTA)
C|==================================================================|
C|      set      FZR and SFij to be zero                            |
C|==================================================================|
C|==================================================================|
C|  Calculatng the amplitude matrices SAij in the particle frame
C|  for all orientations defined by the Gauss quadratures           |
C|  for ZU and HI/PI.  (ÒHETA,PHI)  define the scattering difection |
C|  in the particle symmetry frame,  ZT and SZT are cos and sin THETA,
c   GAM is the angle between the scattering planes 
c   in the particle and the incident wave frames 
C|  For details: See scheme in Khebtsov, JQSRT, 2022                |
C|==================================================================|
C|          Loop over cos (ZU)                                      |
C|==================================================================|
	do i=1,4
	do j=1,4
	FSran(i,j)=0.q0
	end do
	end do


		IZU=1
		
60		DO 63 INZU=1,Norient
			ZU=0.5d0+0.5d0*IZU*ZUX1(INZU)

			theta_angle=180.q0*qacos(zu)/pi

			SZU=qsqrt(1.q0-ZU*ZU)
C|         Loop over azimuth          Hi                            |
C|==================================================================|
			IHI=1
61			DO 62 INHI=1,Norient
				HI=0.5d0+0.5d0*IHI*ZUX1(INHI)
				phi=pi*HI
				phi_angle=180.q0*HI

	call qAsca(Bet,theta_angle,phi_angle,Mmax, N,
     +			  	 A11, A12, A21, A22, ERR_Aij,M_max)
		

				CALL qMULLER(A11,A22,A21,A12,F)

C|==================================================================|
C|  Transformation of F(i,j),  into the incident wave frame
C|  from F(i,j) to FS(i,j) matrix                                    
C|==================================================================|
	

				CPHI=qcos(phi)
	
				ZT=CBET*ZU+SBET*SZU*CPHI
				if(abs(ZT).lt.1.q0) then
					SZT=qsqrt(1.q0-ZT*ZT)
				 	CGAM=(ZU*SBET-SZU*CBET*CPHI)/SZT
				else
					SZT=0.q0
					CGAM=1.q0
				end if

				CALL qSMG(CGAM,MG)
				CALL qSM(phi,MH)

c--------------------------------------------------------------------
				DO 843 ii=1,4
				DO 843 jj=1,4

				FS(ii,jj)=0.q0

					DO 841 K=1,4
					DO 841 L=1,4
841					FS(ii,jj)=FS(ii,jj)+MG(ii,K)*F(K,L)*MH(L,jj)
843				CONTINUE

c	Calcualte the random-averaged values

	do 842 i1=1,4
	do 842 j1=1,4
	FSran(i1,j1)=FSran(i1,j1)+FS(i1,j1)*ZUW(INHI)*ZUW(INZU)/4.q0
842     CONTINUE


c
62			continue
c---------------------------------------------------------------------
			IF (IHI.NE.1) GO TO 63
			IHI=-1
			GO TO 61
63		CONTINUE
		IF (IZU.NE.1) GO TO 64
		IZU=-1
		GO TO 60
64		CONTINUE

C==================================================================
C  Recalculate the Mueller matrix FS from (I,Q,U,V) representation    
C  into (Il,Ir,U,V) representation. For i,j=3,4 no recalculation
C  of  FLran(i,j) is needed, see R. Newton book
C==================================================================
c	Initial FSran(i,j) are in the  (I,Q,U,V) representation
        FLran(1,1) = (FSran(1,1)+FSran(2,2)+FSran(1,2)+FSran(2,1))/2.q0
        FLran(1,2) = (FSran(1,1)-FSran(1,2)+FSran(2,1)-FSran(2,2))/2.q0
        FLran(2,1) = (FSran(1,1)+FSran(1,2)-FSran(2,1)-FSran(2,2))/2.q0
        FLran(2,2) = (FSran(1,1)-FSran(1,2)-FSran(2,1)+FSran(2,2))/2.q0
c
c	For the incident VV polarization, the Stockes vector is (1,-1,0,0)
c	The total scattered intensity Isca=F11(S)*1+F12(S)*(-1)
	  Ftot=FSran(1,1)-FSran(1,2)
c
C===================================================================
c	Normalization to the scattering intensity at wavel=wavelbeg and betsca=90

c	if(wavel.eq.wavelbeg) then
c        if (BET.EQ.90.q0) then
c	      F22_norm=FZR(2,2)
c	      Ftot_norm=Ftot
c	   end if
c	end if

c----------------------------------------------------------------

	  Dvh=FLran(1,2)/FLran(2,2)

c----------------------------------------------------------------
c	Normalization
		x=2.q0*pi*refmed*a0/wavel
		anor=(8.q0/3.q0)/x**2.
		scaI_vv =anor*FLran(2,2)
		phaseF=(4./Qsca/abs(AL)**2)*FSran(1,1)	
c----------------------------------------------------------------
c	Print in output file
c92		format (1h ,f6.1,2x,e11.6)
93		format (1x,f6.1,3(2x,e11.6),2x,f6.1)
c		write (*,92) wavel,Dvh 
		write (16,93) wavel,Dvh,scaI_vv,phaseF,bet

c
c	end loop over  beta
99	Continue
c	end loop over wavelength
999   Continue

	CLOSE (UNIT=16)
	STOP
	END

c*******************************************************************
C                                                MUELLER           
C-------------------------------------------------------------------
C  Calculation of Mueller matrix FSij for   (I,Q,U,V) Stockes vextor
c  For details: See Bohren and Huffmen book
C-------------------------------------------------------------------
c
	SUBROUTINE qMULLER(A1,A2,A3,A4,F)
	IMPLICIT real*16(A-H,O-Z)
	complex*32 A1,A2,A3,A4,DA1,DA2,DA3,DA4
	real*16 F(4,4)
	DA1=conjg(A1)
	DA2=conjg(A2)
	DA3=conjg(A3)
	DA4=conjg(A4)
	F(1,1)=0.5d0*(A1*DA1+A2*DA2+A3*DA3+A4*DA4)
	F(1,2)=0.5d0*(A2*DA2-A1*DA1+A4*DA4-A3*DA3)
	F(1,3)=real(A2*DA3+A1*DA4)
	F(1,4)=imag(A2*DA3-A1*DA4)
	F(2,1)=0.5d0*(A2*DA2-A1*DA1+A3*DA3-A4*DA4)
	F(2,2)=0.5d0*(A1*DA1+A2*DA2-A3*DA3-A4*DA4)
	F(2,3)=real(A2*DA3-A1*DA4)
	F(2,4)=imag(A2*DA3+A1*DA4)
	F(3,1)=real(A2*DA4+A1*DA3)
	F(3,2)=real(A2*DA4-A1*DA3)
	F(3,3)=real(A1*DA2+A3*DA4)
	F(3,4)=imag(DA1*A2+DA3*A4)
	F(4,1)=imag(DA2*A4+A1*DA3)
	F(4,2)=imag(DA2*A4-A1*DA3)
	F(4,3)=imag(A1*DA2-A3*DA4)
	F(4,4)=real(A1*DA2-A3*DA4)
	RETURN
	END
C
C********************************************************************
C
C-------------------------------------------------------------------
c                                                qSMG               
C-------------------------------------------------------------------
C	Calculation of ÌG(I,J) matrix that transforms the Mueller FSij 
c	elements from particle coordinate system to the coordinate system
c    related to the scattering plane
C-------------------------------------------------------------------
	SUBROUTINE qSMG(CGAM,M)
	IMPLICIT real*16(A-H,O-Z)
	real*16 M(4,4)
	CG2=CGAM*CGAM
	SG2=1.q0-CG2
	SG=qsqrt(abs(SG2))
	S2G=2.q0*SG*CGAM
	C2G=CG2-SG2
	M(1,1)=1.q0
	M(1,2)=0.q0
	M(1,3)=0.q0
	M(1,4)=0.q0
	M(2,1)=0.q0
	M(2,2)= C2G
	M(2,3)=-S2G
	M(2,4)=0.q0
	M(3,1)=0.q0
	M(3,2)= S2G
	M(3,3)= C2G
	M(3,4)=0.q0
	M(4,1)=0.q0
	M(4,2)=0.q0
	M(4,3)=0.q0
	M(4,4)=1.q0
	RETURN
	END
C-------------------------------------------------------------------
C*******************************************************************
C                                                qSM                
C-------------------------------------------------------------------
C	Transformation matrix Ì(I,J) from ccordinate system of the incident wave
c	to the particle coordinate system
C-------------------------------------------------------------------

	SUBROUTINE qSM(ALPHA,M)
	IMPLICIT real*16(A-H,O-Z)
	real*16 M(4,4)
	CA=qcos(ALPHA)
	C2A=qcos(2.q0*ALPHA)
	SA=qsin(ALPHA)
	S2A=qsin(2.q0*ALPHA)
	M(1,1)=1.q0
	M(1,2)=0.q0
	M(1,3)=0.q0
	M(1,4)=0.q0
	M(2,1)=0.q0
	M(2,2)= C2A
	M(2,3)=-S2A
	M(2,4)=0.q0
	M(3,1)=0.q0
	M(3,2)= S2A
	M(3,3)= C2A
	M(3,4)=0.q0
	M(4,1)=0.q0
	M(4,2)=0.q0
	M(4,3)=0.q0
	M(4,4)=1.q0
	RETURN
	END
c******************************************************************

c********************************************************************
C                                                         qT
C-------------------------------------------------------------------
c	Subroutine qT to calculate T_m-matrixes for m=0,1,...Mmax
c	The T-matrix order equals N, 
c     AL=k*Rev, k - the complex wave number in medium,
c     epsiln = diameter/rotation length
c	refm - is the complex relative refractive index,
c	NGauss - the number of Gauss quadrature points,
c	TE_0mTM_0, TE_m, and TM_m - output T-matrix arrays
c	They are saved as one-dimensional arrays of length NLxNL
c	and sum_over_m(NmxNm), where Nm=N-m+1, m=1,2,...N
c	TM dimension (M>0): 
c	N=100 Dim_TM=338350 N=60 Dim_TM=73810 N=50 Dim_TM=42925
c	Here N=60 and Dim_TM is set 74000
c	NGauss (max) is set 192
c------------------------------------------------------------ 
        SUBROUTINE qT(AL,EPSILN,
     +	  Nshape,par_c,par_hi,REFM,Mmax,N,NGauss)
      
	 IMPLICIT real*16 (A-H,O-Z)
	 COMPLEX*32 TE_0,TM_0,TE_m,TM_m 
	 COMPLEX*32 AL,refm,aim
	 complex*32 SE11,SM11,SE22,SM22,SE,SM,
     +      TE(60,60),TM(60,60),
     +      TE011re(30,30), TM011re(30,30),
     +      TE022re(30,30), TM022re(30,30),
     +      TE011(30,30),   TM011(30,30), 
     +      TE022(30,30),   TM022(30,30), 
     +      TE011_1(30,30), TM011_1(30,30), 
     +      TE022_1(30,30), TM022_1(30,30), 
     +      TEre(60,60),  TMre(60,60),
     +      TEm(60,60),   TMm(60,60),
     +      TEm_1(60,60), TMm_1(60,60)

	  common /T/ TE_0(1000),TM_0(1000),TE_m(74000),TM_m(74000)

        COMMON /CGAU/ W(192),W1(192)
c------------------------------------------------------------ 
        CALL qGAU(NGauss)
C------------------------------------------------------------
c	Set all initial indexes to be equal zero for
c    arrays TE0(I0), TM0(I0) and TE_m(Im),TM_m(Im)
	I0=0
	Im=0
c------------------------------------------------------------

	do M=0,Mmax

       IF (M.GT.0) GO TO 10
C------------------------------------------------------------
c	Case M=0
c------------------------------------------------------------
	NL=N/2
	NL2=NL*NL

	aim=(0.q0,1.q0)
C------------------------------------------------------------

c	Calcualte TE0re TM0re
        KL=1

	  CALL qTQ(AL,EPSILN,Nshape,par_c,par_hi,
     +           REFM,0,N,NGauss,KL,TE,TM)

C------------------------------------------------------------
c	NOTE: TE011Re and TM011 (odd/odd)   are saved in the left upper submatrix (i,j)
c	NOTE: TE022Re and TM022 (even/even) are saved in the left bottom submatrix (i+NL,j)
C------------------------------------------------------------
	
	do i=1,NL
		iN=i+NL
		do j=1,NL
			TE011re(i,j)=TE(i,j)
			TE022re(i,j)=TE(iN,j)
			TM011re(i,j)=TM(i,j)
			TM022re(i,j)=TM(iN,j)
		end do
	end do

c	Calcualte TE0im and TMim 
	  KL=2			

	  CALL qTQ(AL,EPSILN,Nshape,par_c,par_hi,
     +           REFM,0,N,NGauss,KL,TE,TM)

C------------------------------------------------------------
c	NOTE: TE011Im and TM011Im are saved in the left upper corner (i,j)
c	NOTE  TE022Im and TM022Im are saved in the left bottom corner (i+NL,j)
C------------------------------------------------------------

c	Calculate TE011,TE022,TM011,TM022

	do i=1,NL
		iN=i+NL
		do j=1,NL
			TE011(i,j)=TE011re(i,j)   + aim*TE(i,j)
			TE022(i,j)=TE022re(i,j) + aim*TE(iN,j)
			TM011(i,j)=TM011re(i,j)   + aim*TM(i,j)
			TM022(i,j)=TM022re(i,j) + aim*TM(iN,j)
		end do
	end do
C------------------------------------------------------------
c Calculate inverse matrix TE011_1,TE022_1,TM011_1 and TM022_1 
        CALL qGGG0(NL,TE011,TE011_1)
        CALL qGGG0(NL,TE022,TE022_1)
        CALL qGGG0(NL,TM011,TM011_1)
        CALL qGGG0(NL,TM022,TM022_1)
C------------------------------------------------------------
c	Calcultae TE_011=-TE011re*TE011_1, TE_022=-TE022_re*TE022_1
c     and       TM_011=-TM011re*TM011_1, TM_022=-TM022_re*TM022_1
c	In usual notation: TEii=-RE(Qii)*Qii**(-1), i=1 or 2
c
c	Size I2max=2*NLxNL=N*N/2
C------------------------------------------------------------
	do 2 i=1,NL
		do 2 j=1,NL
			
			I0=I0+1
			I2=I0+NL2
			
			SE11=(0.q0,0.q0)
			SM11=(0.q0,0.q0)
			SE22=(0.q0,0.q0)
			SM22=(0.q0,0.q0)

			do k=1,NL
		       SE11=SE11-TE011re(i,k)*TE011_1(k,j)
		       SE22=SE22-TE022re(i,k)*TE022_1(k,j)
		       SM11=SM11-TM011re(i,k)*TM011_1(k,j)
		       SM22=SM22-TM022re(i,k)*TM022_1(k,j)
			end do
C------------------------------------------------------------
c	Write first TE022, and TM022 (even/even) and then TE011,TM011 (odd/odd) 
C------------------------------------------------------------
			TE_0(I0)  = SE22
			TM_0(I0)  = SM22
			TE_0(I2)  = SE11
			TM_0(I2)  = SM11
2	continue
	go to 30
C------------------------------------------------------------
c	Case M>0
10	continue
C------------------------------------------------------------
c	Calcualte TEre TMre
C------------------------------------------------------------
        KL=1

	  CALL qTQ(AL,EPSILN,Nshape,par_c,par_hi,
     +   	       REFM,M,N,NGauss,KL,TEre,TMre)
c	Calcualte TEim and TMim 

	  KL=2			
	  CALL qTQ(AL,EPSILN,Nshape,par_c,par_hi,
     +           REFM,M,N,NGauss,KL,TE,TM)
C------------------------------------------------------------
c	Calculate TEm=TEre+i*TEim and TMm=TMre+i*TMim

	NM=N-M+1
	do i=1,NM
		do j=1,NM
			TEm(i,j)=TEre(i,j)+aim*TE(i,j)
			TMm(i,j)=TMre(i,j)+aim*TM(i,j)

		end do
	end do


C------------------------------------------------------------
c Calculate inverse matrix TEm_1 and TMm_1
        CALL qGGG(NM,TEm,TEm_1)
        CALL qGGG(NM,TMm,TMm_1)
C------------------------------------------------------------
c	Calcultae TE_m=-TEre*TEm_1 and TM_m=- TMre*TMm_1 matrices
c	In usual notation: TEm=-RE(Qm)*Qm**(-1)
c
c	Size Immax=NLxNL=N*N/4
	Im=(N+1)*(M-1)*(N+1-M) + ((M-1)*M*(2*M-1))/6

	do 22 i=1,NM
		do 22 j=1,NM
			Im=Im+1
			SE=(0.q0,0.q0)
			SM=(0.q0,0.q0)
			do k=1,NM
		       SE=SE-TEre(i,k)*TEm_1(k,j)
		       SM=SM-TMre(i,k)*TMm_1(k,j)
			end do
			TE_m(Im)=SE
			TM_m(Im)=SM
22	continue
c-------------------------------------------------------------

c	END loop over M
30      end do

        RETURN
        END
C************************************************************
C************************************************************
C                                                qTQ        
C------------------------------------------------------------
C  Subroutine for T-matrix calculation                      
c      We consider an axilly symmetrical particle (a,b,b)
c  Nshape=1  spheroid: length aL=2*a, thickness d=2*b
c            semiaxis (a,b,b), a is rotation axis
c		   eqivolume radius = b*(a/b)**(1/3)
c            axis ratio e=a/b (<1 or >1), epsiln=b/a (>1 or <1)
c
c  Nshape=2  right circular cylinder
c            length aL= 2a, and thickness d=2*b (rotation axis), 
c            axis ratio e=a/b (<1 or >1), epsiln = (b/a >1 or <1
c		   eqivolume radius = b*[(3a/2b)]**(1/3)
c
c  Nshape=3  PROLATE s-cylinder with semi-spherical ends,
c            total length aL= 2a (rotation axis),total thickness d=2*b
c            Equivolume radius = b*(1+3(a-b)/2b)**1/3,
c            axis ratio e=a/b>1, epsiln = b/a<1
c
c            OBLATE s-cylinder with length L=2a (rotation axis)
c            radius of rotation is  b=d/2
c            equivolume radius = a*[1+(3/2)*((b-a)/a)**2+(3*pi/4)*((b-a)/a)]**1/3
c            axis ratio e=a/b<1, epsiln=b/a>1
c  Nshape=4  PROLATE Dog-bone rods with elliptical caps
c		   Four parameters have to be defeined:
c			Total length L=2(a+b_c)=2a*(1+par_c)
c			Cental thickness d=2*b
c 
c			OBLATE disks with eeliptical edges:
c			the disk thickness L=2a
c			the total disk diameter d=2b; 
c			The ellipce axes are a and b_c=a*par_c, a>b_c
c
c    Calculate   FW = r(theta)    FW1(M>0) = FW1(M=0)*sin(theta)*sin(theta)
c                                 
c
c   For the case M=0              FW1=[dr(theta)/d(theta)]/r(theta)/sin(theta)             
c
c               Y  = kr(theta) and  Y1 = Y*ref(m)
c                    
C------------------------------------------------------------
        SUBROUTINE qTQ(AL,EPS,Nshape,par_c,par_hi,
     +    	         REFM,M,N,NGauss,KLUCH,TE,TM)
c
        IMPLICIT real*16 (A-H,O-Z)
 	  real*16 qz,Q(60),S(60) 
	  complex*32 AL,refm	  
	  complex*32 X,Y, YC,Y1,I,J,K,L,JJ,JJ1,KMJ,
     +             JMK,IML,LMI
	  complex*32 HH,HH1,FK,FL,FI,FJ 
        complex*32 TE(60,60),TM(60,60)
        complex*32 qBES(101),qNEU(101)

        COMMON /qBES/ qBES
        COMMON /qNEU/ qNEU
        COMMON /CGAU/ W(192),W1(192)
        COMMON /QS/ Q,S
C------------------------------------------------------------
       NI=1  
C------------------------------------------------------------
C  Set all T-matrices =0                                     
C------------------------------------------------------------
        M2=M+1
        IF (M.GT.0) M2=M
        NM=N-M2+1
        DO 10 N1=1,NM
        DO 10 N2=1,NM
        TE(N1,N2)=(0.q0,0.q0)
10      TM(N1,N2)=(0.q0,0.q0)
C------------------------------------------------------------
C  IV is a switch for Gauss absc.  "+" (IV=1) and "-" for   |
C  (IV=-1); KI is a counter for number of subintervals   ,  |
C  NG is a counter for Gauss points                         |
C------------------------------------------------------------
        IV=1
c       KI=1
101     NG=1
11      Z=(1.0d0+IV*W1(NG))/2.0d0
C------------------------------------------------------------
c	Calculation of the angular functions Q[m,m+n-1] and S[m,m+n-1]
C------------------------------------------------------------
	  qz=Z
        CALL qQSM(M,NM,qz)
C------------------------------------------------------------
c	Calculateing the parameters Y=k1*r(theta), Y1=k2*r(theta),
c	k1=(2pi/lambda)*n1, k2=(2pi/lambda)*n2=k1*refm, n1=refmed, refm=n2/n1
c	V=RT**2*W(NG), RT(M=0)=[1/sin(theta)/r(theta)]*dr(theta)/dtheta
c	W(NG) is the Gauss quadrature weight, 
c	r(theta) is the radius-vector defining the particle shape
c	see, e.g. Msishchenko et al book, 2002 Eq.(5.270) for spheroids.
	call qRTFW1(Z,AL,refm,Nshape,
     +	eps,par_c, par_hi,NG,V,Y,Y1,FW1)

C------------------------------------------------------------

C------------------------------------------------------------
C  KLUCH is a switch for calc. Re T or Im T, equals 1 or 2,    |
C          respectively                                     |
C------------------------------------------------------------
        IF (KLUCH.EQ.2) GO TO 201
        YC=Y
c       CALL cBESSEL(YC,NM,M)
        CALL qBESSEL(YC,NM,M)
        DO 200 II=1,NM+1
200     qNEU(II)=qBES(II)
        GO TO 202
c201     CALL cNEYMAN(Y,NM,M)
c202     CALL cBESSEL(Y1,NM,M)
201     CALL qNEUMAN(Y,NM,M)
202     CALL qBESSEL(Y1,NM,M)
        IF (M.NE.0) GO TO 50
C------------------------------------------------------------
C  Calculation  T-matrix for the case M = 0                 |
C------------------------------------------------------------
        NL=N/2
C------------------------------------------------------------
C  KL is a switch for odd/odd (1) or even/even (2), cases      |
C------------------------------------------------------------
        KL=1
12      DO 14 NL1=1,NL
        N1=NL1*2-1
        IF (KL.EQ.2) N1=NL1*2
        N21=2*N1+1
        HH=N1*(N1+1)*N21*qNEU(N1)/Y
        HH1=HH/N1-N21*qNEU(N1+1)
        FK=-(FW1*HH*Q(N1)+HH1*S(N1))*V/N21
        FJ=qNEU(N1)*S(N1)*V
        DO 14 NL2=1,NL
        N2=NL2*2-1
        IF (KL.EQ.2) N2=NL2*2
        N12=N2+1
        JJ=N12*qBES(N2)/Y1
        K=FK*qBES(N2)*S(N2)
        J=FJ*((JJ-qBES(N2+1))*S(N2)+FW1*JJ*N2*Q(N2))
        IF (KL.EQ.2) GO TO 13
        TE(NL1,NL2)=TE(NL1,NL2)+(K+REFM*J)
        TM(NL1,NL2)=TM(NL1,NL2)+(K*REFM+J)
        GO TO 14
13      TE(NL1+NL,NL2)=TE(NL1+NL,NL2)+(K+REFM*J)
        TM(NL1+NL,NL2)=TM(NL1+NL,NL2)+(K*REFM+J)
14      CONTINUE
        IF (KL.EQ.2) GO TO 15
        KL=2
        GO TO 12
15      NG=NG+1
        IF (NG.LE.NGauss) GO TO 11
        IF (IV.NE.1) GO TO 16
        IV=-1
        GO TO 101
c16      WRITE (21) TM(N1,N2)
16	  RETURN
50      CONTINUE
C------------------------------------------------------------
C  Calculation of T-matrix for M>0                           |
C------------------------------------------------------------
        FW1=FW1*(1.0d0-Z*Z)
        DO 100 N1=1,NM
        N1M=N1+M-1
        N21=2*N1M+1
        LN1=((N1+1)/2-N1/2)*2
        FJ=qNEU(N1)*V
        FK=V/N21
        FI=FJ*M
        FL=FK*M
        HH=N1M*(N1M+1)*N21*qNEU(N1)/Y
        HH1=HH/N1M-N21*qNEU(N1+1)
        DO 100 N2=1,NM
        N2M=N2+M-1
        QQ=Q(N1)*Q(N2)
        QS=Q(N1)*S(N2)
        SQ=S(N1)*Q(N2)
        JJ=(N2M+1)*N2M*qBES(N2)/Y1
        JJ1=JJ/N2M-qBES(N2+1)
        LN=LN1+(N2+1)/2-N2/2+1
        GO TO (51,52,52,51),LN
51      QQSS=M*M*QQ+S(N1)*S(N2)
        K=FK*qBES(N2)*(-FW1*HH*QS-HH1*QQSS)
        J=FJ*(JJ1*QQSS+FW1*JJ*SQ)
        KMJ=(K+REFM*J)
        JMK=(J+REFM*K)
         GO TO 53
52      SQQS=SQ+QS
        I=FI*qBES(N2)*SQQS
        L=FL*(HH1*JJ1*SQQS+FW1*QQ*(HH*JJ1+JJ*HH1))
        LMI=(L+REFM*I)
        IML=(I+REFM*L)
53      GO TO (61,62,63,64),LN
61      TE(N1,N2)=TE(N1,N2)+JMK
        TM(N1,N2)=TM(N1,N2)+KMJ
        GO TO 100
62      TE(N1,N2)=TE(N1,N2)-IML
        TM(N1,N2)=TM(N1,N2)+LMI
        GO TO 100
63      TE(N1,N2)=TE(N1,N2)+LMI
        TM(N1,N2)=TM(N1,N2)-IML
        GO TO 100
64      TE(N1,N2)=TE(N1,N2)+KMJ
        TM(N1,N2)=TM(N1,N2)+JMK
100     CONTINUE
        NG=NG+1
        IF (NG.LE.NGauss) GO TO 11
        IF (IV.NE.1) GO TO 65
        IV=-1
        GO TO 101
65	  RETURN
        END
C
c*************************************************************
c	Subroutine qRTFW1 to calculate the parameters V, Y=k1*r(theta), Y1=k2*r(theta), FW1
c	k1=(2pi/lambda)*n1, k2=(2pi/lambda)*n2=k1*refm, n1=refmed, refm=n2/n1
c	V=RT**2*W(NG), RT(M=0)=[1/sin(theta)/r(theta)]*dr(theta)/dtheta
c	W(NG) is the Gauss quadrature weight, 
c	r(theta) is the radius-vector defining the particle shape
c	see, e.g. Msishchenko et al book, 2002 Eq.(5.270) for spheroids.

	subroutine qRTFW1(Z,AL,refm,Nshape,
     +	eps,par_c, par_hi,NG,V,Y,Y1,FW1)
c
	IMPLICIT REAL*16 (A-H,O-Z)
	complex*32 AL,refm,X,Y,Y1
      COMMON /CGAU/ W(768),W1(768)
c
C------------------------------------------------------------
c	Case 1: Spheroids
c
      if(Nshape.eq.1) then
        EPS2=(EPS*EPS-1.0D0)
        RT=1.q0/qsqrt(1.0D0+EPS2*Z*Z)
        RT2=RT*RT
        FW1=EPS2*Z*RT2
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
        X=AL*EPS**(1.q0/3.q0)
        V=W(NG)*RT2
        Y=X*RT
        Y1=REFM*Y
        go to 4321
      end if
C------------------------------------------------------------
c
c	Case 2: Right circular cylinders
c
      if(Nshape.eq.2) then
        X=AL*(2.q0*EPS/3.q0)**(1.q0/3.q0)
        SZ=qsqrt(1.q0-Z*Z)
        
        if(SZ/Z.lt.eps) then
          RT = 1.q0/EPS/Z
          FW1= 1.q0/Z
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
        else
          RT = 1.q0/SZ
          FW1= -z/(SZ*SZ)
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
        end if
      
        RT2=RT*RT
        V=W(NG)*RT2
        Y=X*RT
        Y1=REFM*Y
        go to 4321
      end if
C------------------------------------------------------------

c	Case 3: s-cylinders      
c
      if(Nshape.eq.3) then
        if(eps.eq.1.q0) then
          V=W(NG)
          Y=AL
          Y1=Refm*Y
	    FW1=0.q0
          go to 4321
        end if   
	  if(eps.lt.1.q0) then
c
c		PROLATE s-cylinder
c
	     SZ2=(1.q0-Z)*(1.q0+Z)
		 SZ=qsqrt(SZ2)
           ab=(1.q0-eps)/eps
           ba=1.q0/ab
           ba2=ba*ba
           X=AL/(1.q0+3.q0*ab/2.q0)**(1.q0/3.q0)
        
           if(SZ/Z.lt.ba) then
			RT1=qsqrt((ba-SZ)*(ba+SZ))
		    RT=ab*(Z+RT1)
			FW1=-1.q0/RT1

C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
		 else
			RT = 1.q0/SZ
			FW1= -z/SZ2

C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
		 end if
           RT2=RT*RT
           V=W(NG)*RT2
           Y=X*RT
           Y1=REFM*Y
	  else

c	     OBLATE s-cylinder e=b/(a+b)<1, b is rotation axis
c          eps=1/e=(a+b)/b
c
	     SZ2=(1.q0-Z)*(1.q0+Z)
		 SZ=qsqrt(SZ2)
	     ab=eps-1.q0
           ba=1.q0/ab
           ab2=ab*ab
		 ba2=ba*ba
           pi=qatan(1.q0)*4.q0
		 X=AL/(1.q0+3.q0*ab2/2.q0+3.q0*pi*ab/4.q0)**(1.q0/3.q0)
		        
          if(SZ/Z.gt.ab) then
            RT1=qsqrt((ba-Z)*(ba+Z))
		  RT=ab*(SZ+RT1)
            FW1=Z/SZ/RT1
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
          else
            RT = 1.q0/Z
            FW1= RT
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
          end if
          RT2=RT*RT
          V=W(NG)*RT2
          Y=X*RT
          Y1=REFM*Y
	  end if
 
	  go to 4321
	end if
C------------------------------------------------------------
c	Case 4: Dog-bone rods with an elliptic cap
c	Calcualtion of the equivolume radius:
c	Radius of a sphere for side rod curve
c		R=(a**2.q0+(b1-b)**2.q0)/2.q0/(b1-b)
c		p=a/b
c		zz=p*xx=2*p*hi/[p**2+hi**2] xx=b/R=2*hi/[p**2+hi**2]=zz/p
c       zp=zz*p
c
	if (Nshape.eq.4) then
c	if eps=1/e>1 then goto 4320 Case of disks with an elliptical edge
	if (eps.gt.1.) go to 4320
C	Case of dog-bone rods		
c	e=1/eps=aLi/2b aspect rato, b1b=b1/b, bc=par_c*b1, b1=b*(1+par_hi) [b=b_aver/(1+par_hi/2)]
		e =  1.q0/eps
	  b1b = 1.q0+par_hi
		p = e - par_c*b1b
		zz=2.q0*p*par_hi/(p*p+par_hi**2)
		xx=zz/p
		zp=zz*p
c		f(x)=1 + 2*(1+x)/[x**2]*fi(z) - p**2/3
c		fi(zz)=1 - (1/2)*[arcsin(zz)/zz+sqrt(1-zz**2)]     		
			if(zp.lt.0.001) then
			ff=1.q0+(zp/3.q0)*(1.q0+(3.q0/20.q0)*zz*zz)+zp*zp/20.q0
			else
			fi=1.q0-( qasin(zz)/zz+qsqrt((1.q0-zz)*(1.q0+zz)))/2.q0	
			ff=1.q0+(2.q0*(1.q0+xx)/xx/xx)*fi-p**2/3.q0
			end if

c		Recalulate "b" by using a0 and the scale multiplier
	scale_b = ( (b1b)**3*par_c + 1.5d0*p*ff )**(1.q0/3.q0)
         b = AL/scale_b
		b1 = b*b1b
		bc = b1*par_c
	  b2 = bc*bc/b1/b1
		a  = b*p
		a2 = 1.q0-bc*bc/a/a
	   X = b
       SZ2 = (1.q0-z)*(1.q0+z)
	 SZ = qsqrt(SZ2)

		if(sz/z.lt.b1/a) then
			g = z*z+sz2*b2
	       dg = 2.q0*z*(1.q0-b2)
		  RT1 = qsqrt(z*z-a2*g)
		   RT = p*(z+RT1)/g
             FW2 = dg*(1.q0/g+a2/2.q0/RT1/(z+RT1))
	      FW1 = -1.q0/RT1 +FW2
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
 		   RT2 = RT*RT
		 	V = W(NG)*RT2
			Y = X*RT
	       Y1 = REFM*Y
		else
			xx1=1.q0+xx
c*******************************************
			szb=qsqrt((1.q0/xx1-z)*(1.q0/xx1+z))
c			szb=qsqrt((xx1)**(-2)-z*z)   
c*******************************************
			RT=(2.q0+xx)/xx1/(sz+szb)
			FW1=-z/sz/szb
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
		    RT2 = RT*RT
		    V = W(NG)*RT2
		   	Y = X*RT
   	        Y1 = REFM*Y
		end if
c
		goto 4321
C------------------------------------------------------------
c	Case of disks with elliptic edge
c
4320	ec=par_c
c	NOTE!!!: here the symbol 'e' stands for the e=eps=b/a>1
c	Y=X*RT X=k1*a  FW1=d(RT)/d(theta)/RT/sin(theta)
c
c	if z=cos(theta)>=z1=a/sqrt(a**2+b**2) then
c	RT=1/z FW1=T/sz=RT, sz=sin(theta)
c
c	if z<z1 RT=p*[1+det] FW1=(1/sz/p)*dp/d(theta)-[1/2/sz/(1-q-det)]*dq/d(theta)
c	det=dsqrt(1.-q), q=q1*[1+ec**2*z*z/sz/sz], q1=(e*e-2*e*ec)/(e-ec)**2
c	dp/d(theta)=dpdt=(e-ec)*z/p1*[1-2*sz*sz(1-ec**2)/p1]
c	p1=1+(ec**1-1)*z*z
c
	e=eps
	
	eec2=(e-ec)*(e-ec)
	ec2=ec*ec
	
	pi=qatan(1.q0)*4.q0
	sc_a=(ec2+1.5q0*eec2+0.75q0*pi*ec*(e-ec))**(1.q0/3.q0)
c
	  X=AL/sc_a
c	  X==ka=AL/sc_a 'a' is the symmetry axis, L=2a is the disk thickness
	
		z1=1./qsqrt(1.q0+eec2)
		if (z.ge.z1) then
			RT=1.q0/z
			FW1=RT
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
		else
	
			sz2=(1.q0-Z)*(1.q0+Z)
			z2=z*z
			sz=qsqrt(sz2)
			ctg2=z*z/sz2
			p=(e-ec)*sz/(1.q0+(ec2-1.q0)*z2)
			q1=(e*e-2.q0*e*ec)/eec2
			qq=q1*(1.q0+ec2*ctg2)
			det=qsqrt(1.q0-qq)
			p1=1.q0+(ec2-1.q0)*z2
			
			dpdt=(e-ec)*z/p1*(1.q0-2.q0*sz2*(1.q0-ec2)/p1)
			dqdt=-2.q0*q1*ec2*z/sz/sz2

			RT=p*(1.q0+det)
			FW1=dpdt/sz/p-dqdt/2.q0/sz/(1.q0-qq+det)
C                                  IF M>0     FW1=FW1(M=0)*(1-Z*Z)
	    end if
c	 			
        RT2=RT*RT
        V=W(NG)*RT2
        Y=X*RT
        Y1=REFM*Y
c	end of Case Nshape=4 e=a/b<1 disck with elliptic side
C------------------------------------------------------------
	end if     

4321  continue        
	
	RETURN
	END
c******************************************************************
C
C                                                         qQSM     
C-------------------------------------------------------------------
c	Calculation of the angular function Q[m,m+n-1] and S[m,m+n-1]
c	These are saved in Q(n=1,...NM), S(n=1,...NM)
c	Case m=0: 
c	CMMON /qQSM/ Q,S
C-------------------------------------------------------------------
        SUBROUTINE qQSM(M,NM,Z)
        IMPLICIT REAL*16(A-H,O-Z)
        real*16 Z, Q(60),S(60)
	  COMMON /QS/ Q,S
C-------------------------------------------------------------------
C                                                   Case m>0
C-------------------------------------------------------------------
        IF (Z.NE.1.q0) GO TO 7
        Z12=0.q0
        GO TO 8
7       Z12=qSQRT(qABS((1.q0-Z)*(1.q0+Z)))
8       IF (M.EQ.0) GO TO 30
C-------------------------------------------------------------------
C                                                Calculate Q(1)
C-------------------------------------------------------------------
        U=1.0q0
        IF(M.EQ.1) GO TO 11
        DO 10 N=2,M
10      U=U*Z12*(2*N-1)
11      V=0.q0
        Q(1)=U
        S(1)=M*Z*U
C-------------------------------------------------------------------
C                                                Recursion
C-------------------------------------------------------------------
        DO 20 N=2,NM
        N1=N+M-1
        W=((2*N1-1)*Z*U-(N1+M-1)*V)/(N1-M)
        Q(N)=W
        S(N)=N1*Z*W-(N1+M)*U
        V=U
        U=W
20      CONTINUE
        RETURN
C-------------------------------------------------------------------
C                              End calculations for M>0
C-------------------------------------------------------------------
C                                                   Case m=0
C-------------------------------------------------------------------
30      V=Z12
        U=Z*Z12
        V1=0.q0
        U1=-Z12
        Q(1)=U
        S(1)=U1
        DO 40 N=1,NM-1
        ZN=Z*(2*N+1)
        W=(ZN*U-N*V)/(N+1)
        W1=(ZN*U1-(N+1)*V1)/N
        Q(N+1)=W
        S(N+1)=W1
        V=U
        U=W
        V1=U1
        U1=W1
40      CONTINUE
        RETURN
        END
c********************************************************************

c*************************************************************
c										qBESSEL(z,Nm,m)
c------------------------------------------------------------
C     Calculate the comples spherical Bessel functions j_n(z) of a complex argument z
c	m - Tm matrix order, Nm=[N for m=0, N-m-1 for m>0] Tm matrix dimension	     
C     Output: COMMON /qBES/  array     BES(NM+1).                    
C------------------------------------------------------------
        SUBROUTINE qBESSEL(z,NM,M)
	  implicit real*16 (A-H,O-Z)
	  COMPLEX*32 cqS,u,v,w
	  complex*32 qBES(61)
	  complex*32 z
c
	  COMMON /qBES/ qBES
c
	  IF (abs(z).NE.0.) GO TO 9
        DO 5 N1=1,NM+1
5       qBES(N1)=(0.q0,0.q0)
        RETURN

9       N=NM+M
        IF (M.EQ.0) N=N+1
        az=abs(z)
	  nmax=400.0d0/(3.q0-log10(az))

c        NMAX=80.0D0/(2.44D0-CDLOG(X)/2.3)
        U=(0.q0,0.q0)
c        V=(1.0D-36,0.0D0)
        V=(1.q-300,0.q0)
        N1=NMAX
10      N2=2*N1+1
        W=N2*V/z-U
        NL=N1-1
        IF (NL.EQ.0) GO TO 21
        IF (NL.GT.N) GO TO 20
        qBES(NL)=W
20      N1=NL
        U=V
        V=W
        IF (N1.GT.0) GO TO 10
21      zr=real(z)
	  zi=imag(z)
        CR=qsin(zr)*qcosh(zi)
	  CI=qcos(zr)*qsinh(zi)
        cqS=cmplx(CR,CI,kind=16)/z/W
c        C=cqSIN(z)/z/W
        DO 50 N1=1,NM+1
        N2=M+N1-1
        IF (M.EQ.0) N2=N1
50      qBES(N1)=qBES(N2)*cqS
        RETURN
        END
c************************************************************
c
c************************************************************
C                                         qNEUMAN(z, Nm,m)   
C------------------------------------------------------------
C     Calculate the complex Neuman functions of a complex argument  z                     
c	m - Tm matrix order, Nm=[N for m=0, N-m-1 for m>0] Tm matrix dimension	     
C     Output: COMMON /CNEU/  array  NEU(Nm+1).                        
C------------------------------------------------------------
        SUBROUTINE qNEUMAN(z,NM,M)
c
        IMPLICIT REAL*16 (A-H,O-Z)
        complex*32 DC,DS,U,V,W 
	  complex*32 z    
	  complex*32 qNEU(61)
        COMMON /qNEU/ qNEU
5       FORMAT (1H ,'ERROR !   Z = 0  ! ! !')
        IF (abs(z).NE.0.) GO TO 9
        WRITE (*,5)
	  STOP
        RETURN
9       N=NM+M
        IF (M.EQ.0) N=N+1
	  zr=real(z)
	  zi=imag(z)
	  CZr=qcos(zr)*qcosh(zi)
	  CZi=-qsin(zr)*qsinh(zi)
	  DC=-cmplx(CZr,CZi,kind=16)/z	
c       DC=-qcos(z)/z
	  SZr=qsin(zr)*qcosh(zi)
	  SZi=qcos(zr)*qsinh(zi)
        DS=cmplx(SZr,SZi,kind=16)/z		
c        DS=qsin(z)/z
        U=DC
        V=DC/z-DS
        qNEU(1)=V
        DO 40 N1=2,N
        W=(2*N1-1)/z*V-U
        U=V
        V=W
40      qNEU(N1)=W
        DO 50 N1=1,NM+1
        N2=M+N1-1
        IF (M.EQ.0) N2=N1
50      qNEU(N1)=qNEU(N2)
        RETURN
        END
c********************************************************************

c********************************************************************
C                                                   qGGG(N,A,B)      
C--------------------------------------------------------------------
C   Subroutine to calculate the inverse matrix by Gauss elimination
C   with a partial choise of the leading element
c   see Forsait book, 1980, p.65
c   Input: A(N,N) matrix, Output A**(-1)=B(N,N), IPVT(permutation vector)
C--------------------------------------------------------------------
        SUBROUTINE qGGG(N,A,B)
        IMPLICIT real*16(A-H,O-Z)
        complex*32 A(60,60),B(60,60),T
        INTEGER IPVT(60)
C--------------------------------------------------------------------
C     Define the unit matrix B=E
C********************************************************************
        DO 3 I=1,N
        DO 2 J=1,N
2       B(I,J)=(0.q0,0.q0)
3       B(I,I)=(1.q0,0.q0)
C--------------------------------------------------------------------
C    Triangular decomposition, the initial matrix A will be modified!!!
C--------------------------------------------------------------------
        IPVT(N)=1
        IF (N.NE.1) GO TO 14
C--------------------------------------------------------------------
C                                   Case N=1
C--------------------------------------------------------------------
        IF (A(1,1).NE.(0.q0,0.q0)) GO TO 141
        WRITE (*,1000)
        RETURN
141     B(1,1)=1.q0/A(1,1)
        RETURN
C--------------------------------------------------------------------
C                                    Case N>1
C--------------------------------------------------------------------
14      NM1=N-1
        DO 35 K=1,NM1
        KP1=K+1
        M=K
        DO 15 I=KP1,N
        IF (abs(A(I,K)).GT.abs(A(M,K))) M=I
15      CONTINUE
        IPVT(K)=M
        IF (M.NE.K) IPVT(N)=-IPVT(N)
        T=A(M,K)
        A(M,K)=A(K,K)
        A(K,K)=T
        IF (T.EQ.(0.q0,0.q0)) GO TO 35
        DO 20 I=KP1,N
        A(I,K)=-A(I,K)/T
20      CONTINUE
        DO 30 J=KP1,N
        T=A(M,J)
        A(M,J)=A(K,J)
        A(K,J)=T
        IF (T.EQ.(0.q0,0.q0)) GO TO 30
        DO 25 I=KP1,N
        A(I,J)=A(I,J)+A(I,K)*T
25      CONTINUE
30      CONTINUE
35      CONTINUE
C--------------------------------------------------------------------
C	Calculation of the inverse matrix: Direct calculations
C--------------------------------------------------------------------
42      CONTINUE
        DO 60 J=1,N
        DO 22 K=1,NM1
        KP1=K+1
        M=IPVT(K)
        T=B(M,J)
        B(M,J)=B(K,J)
        B(K,J)=T
        DO 12 I=KP1,N
        B(I,J)=B(I,J)+A(I,K)*T
12      CONTINUE
22      CONTINUE
60      CONTINUE
C--------------------------------------------------------------------
C       Inverse substitution
C--------------------------------------------------------------------
        DO 61 J=1,N
        DO 45 KB=1,NM1
        KM1=N-KB
        K=KM1+1
        B(K,J)=B(K,J)/A(K,K)
        T=-B(K,J)
        DO 32 I=1,KM1
        B(I,J)=B(I,J)+A(I,K)*T
32      CONTINUE
45      CONTINUE
C--------------------------------------------------------------------
C            Exact degeneracy
C--------------------------------------------------------------------
40      IF (A(1,1).NE.(0.q0,0.q0)) GO TO 41
1000    FORMAT (1H ,'Exact degeneracy!!!')
        WRITE (*,1000)
        RETURN
C--------------------------------------------------------------------
41      B(1,J)=B(1,J)/A(1,1)
61      CONTINUE
        RETURN
        END
c********************************************************************

C
c********************************************************************
C                                                  qGGG0(N,A,B)     
C	The same as qGGG but for Nmax=50
C--------------------------------------------------------------------
        SUBROUTINE qGGG0(N,A,B)
        IMPLICIT real*16(A-H,O-Z)
        complex*32 A(30,30), B(30,30), T
        INTEGER IPVT(50)

C--------------------------------------------------------------------
C     Define unit matrix B=E
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
        DO 3 I=1,N
        DO 2 J=1,N
2       B(I,J)=(0.q0,0.q0)
3       B(I,I)=(1.q0,0.q0)
C--------------------------------------------------------------------
C    Triangular decomposition, the initial matrix A will be modified!!!
C--------------------------------------------------------------------
        IPVT(N)=1
        IF (N.NE.1) GO TO 14
C--------------------------------------------------------------------
C                                   Case N=1
C--------------------------------------------------------------------
        IF (A(1,1).NE.(0.q0,0.q0)) GO TO 141
        WRITE (*,1000)
        RETURN
141     B(1,1)=1.q0/A(1,1)
        RETURN
C--------------------------------------------------------------------
C                                    Case N>1
C--------------------------------------------------------------------
14      NM1=N-1
        DO 35 K=1,NM1
        KP1=K+1
        M=K
        DO 15 I=KP1,N
        IF (abs(A(I,K)).GT.abs(A(M,K))) M=I
15      CONTINUE
        IPVT(K)=M
        IF (M.NE.K) IPVT(N)=-IPVT(N)
        T=A(M,K)
        A(M,K)=A(K,K)
        A(K,K)=T
        IF (T.EQ.(0.q0,0.q0)) GO TO 35
        DO 20 I=KP1,N
        A(I,K)=-A(I,K)/T
20      CONTINUE
        DO 30 J=KP1,N
        T=A(M,J)
        A(M,J)=A(K,J)
        A(K,J)=T
        IF (T.EQ.(0.q0,0.q0)) GO TO 30
        DO 25 I=KP1,N
        A(I,J)=A(I,J)+A(I,K)*T
25      CONTINUE
30      CONTINUE
35      CONTINUE
C--------------------------------------------------------------------
C	Calculation of the inverse matrix: Direct calculations
C--------------------------------------------------------------------
42      CONTINUE
        DO 60 J=1,N
        DO 22 K=1,NM1
        KP1=K+1
        M=IPVT(K)
        T=B(M,J)
        B(M,J)=B(K,J)
        B(K,J)=T
        DO 12 I=KP1,N
        B(I,J)=B(I,J)+A(I,K)*T
12      CONTINUE
22      CONTINUE
60      CONTINUE
C--------------------------------------------------------------------
C       Inverse substitution
C--------------------------------------------------------------------
        DO 61 J=1,N
        DO 45 KB=1,NM1
        KM1=N-KB
        K=KM1+1
        B(K,J)=B(K,J)/A(K,K)
        T=-B(K,J)
        DO 32 I=1,KM1
        B(I,J)=B(I,J)+A(I,K)*T
32      CONTINUE
45      CONTINUE
C--------------------------------------------------------------------
C            Exact degeneracy
C--------------------------------------------------------------------
40      IF (A(1,1).NE.(0.q0,0.q0)) GO TO 41
1000    FORMAT (1H ,'Exact degeneracy!!!')
        WRITE (*,1000)
        RETURN
C--------------------------------------------------------------------
41      B(1,J)=B(1,J)/A(1,1)
61      CONTINUE
        RETURN
        END
c********************************************************************
c												qQesa
c********************************************************************
c	Subroutine for calculation of the forfard-scattering amplitudes
c	A1 and A2 and the extinxtion, scattering, and absorption efficiencies
c     for random orientations and for TE and TM polarizations at longitudinal
c     perpendicular orientations of the symmetry axis
C--------------------------------------------------------------------

	subroutine qQesa(AL,Mmax,N,Z,A1,A2,Qext,Qsca,Qabs,
     +                 QexE,QexM,QscE,QscM,QabE,QabM,ERR_Q,M_max)
c	NOTE: M_max is actual M<=Mmax that ensures ERR_Q accuracy of calculations
c
	IMPLICIT real*16(A-H,O-Z)
        complex*32 PE,PM,QE,QM,AI1,AIM,
     +             DS,KOCM,KOCE,KOCX,A1,A2,A1m,A2m
	complex*32 AL

	A1=(0.q0,0.q0)
	A2=(0.q0,0.q0)
	Qext=0.
	Qsca=0.
	QexE=0.
	QexM=0.
	QscE=0.
	QscM=0.
c	Start summation over M
c	
	M=0
50	continue
	call qQesa_m(AL,M,N,Z,A1m,A2m,Qe_m,
     +  	Qs_m,Qe_Mm,Qs_Mm,Qe_Em,Qs_Em)

	A1=A1+A1m
	A2=A2+A2m
	Qext=Qext+Qe_m
	Qsca=Qsca+Qs_m
	QexE=QexE+Qe_Em
	QexM=QexM+Qe_Mm
	QscE=QscE+Qs_Em
	QscM=QscM+Qs_Mm
	Qabs=Qext-Qsca
	QabE=QexE-QscE
	QabM=QexM-QscM


c	Check accuracy for M>1	
	IF (M.LE.1) GO TO 80
	
	T1=A1*conjg(A1)
	T2=A2*conjg(A2)
	
	e1=A1m*conjg(A1m)
	e2=A2m*conjg(A2m)

	if(T1.ne.0.q0) err1= abs(e1/T1)
	if(T2.ne.0.q0) err2= abs(e2/T2)
	if(Qsca.ne.0.q0) err3= abs(Qs_m/Qsca)
	if(Qabs.ne.0.q0) err4= abs((Qe_m-Qs_m)/Qabs)

	err_max=max(err1,err2,err3,err4)
	 
	if(err_max.lt.ERR_Q) then
	  M_max=M
	  go to 81
	end if

80	continue
	M=M+1
c
	if (M.GT.Mmax) then
	   M_max=M-1
	   go to 81
	end if

	GO TO 50

81	continue

	RETURN
	END

c*******************************************************************
c
c												qQesa_m
c-------------------------------------------------------------------
c	Subroutine to calculate the partial 'm' forward scattering amplitudes A1 and A2
c     and the extinction and scattering efficiencies for random orientations 
c-------------------------------------------------------------------
	subroutine qQesa_m(AL,M,N,Z,A1,A2,OCX,SCX,OCM,SCM,OCE,SCE)
c
	IMPLICIT real*16(A-H,O-Z)
        complex*32 PE,PM,QE,QM,AI1,AIM,
     +             DS,KOCM,KOCE,KOCX,A1,A2
		complex*32 AL
        real*16 KSCM,KSCE,KSCX

	  real*16 qz,Q(60),S(60)		
        COMMON /QS/ Q,S

        COMMON /CPQM/ PE(60),PM(60),QE(60),QM(60)
        COMMON /CDMN/ D(60)
        COMMON /CKKK/ KSCM,KOCM,KSCE,KOCE,KOCX,KSCX

        AI1=(0.q0,1.q0)
c        PI=3.14159265358D0
c        THETA=NTHE*PI/180.q0
c        Z=qcos(THETA)

        A1=(0.q0,0.q0)
        A2=(0.q0,0.q0)
40      NM=N-M+1
        IF (M.EQ.0) NM=N
        M2=M*2
        AIM=1.q0
        
	  KOCX=(0.q0,0.q0)
        KOCE=(0.q0,0.q0)
        KOCM=(0.q0,0.q0)
c        
	  KSCX=0.q0
	  KSCE=0.q0
        KSCM=0.q0

        qz=Z
	  CALL qQSM(M,NM,qz)
        CALL qDMN(M,NM)
c****   CALL qPQM(M,N)
        CALL qPQM(M,NM,N)
        IF (M.NE.0) GO TO 42
C===================================================================
C                                                         M=0      
C====================================================================
        AIM=AI1
        DO 41 N1=1,N
        DS=D(N1)*S(N1)/AIM
        A1=A1-DS*PE(N1)*AI1
        A2=A2+DS*QE(N1)
41      AIM=AIM*AI1
        GO TO 45
C===================================================================
C³                                                         M>0      
C===================================================================
42      DO 43 N1=1,M
43      AIM=AIM*AI1
        DO 44 N1=1,NM
        DS=D(N1)/AIM*10.q0**(-f)
        A1=A1-DS*(S(N1)*PE(N1)*AI1+M*Q(N1)*QM(N1))
        A2=A2+DS*(S(N1)*QE(N1)-M*Q(N1)*PM(N1)*AI1)
44      AIM=AIM*AI1
45      SCE=KSCE/abs(AL)**2
        SCM=KSCM/abs(AL)**2
        SCX=2.q0*KSCX/abs(AL)**2
c        OCX=-2.q0*DREAL(KOCX)/AL/AL
c        OCE=-DREAL(KOCE)/AL/AL
c        OCM=-DREAL(KOCM)/AL/AL
        OCX=-2.q0*real(KOCX/AL)/real(AL)
        OCE=-real(KOCE/AL)/real(AL)
        OCM=-real(KOCM/AL)/real(AL)
 
c	Change sings of the scattering amplitudes
c	The optical theorem: Qext1,2= 4*Imag[A1,2(0)/alpha_0]/real(alpha_0)
c	alpha_0=[2*pi*refmed/lambda_vac]*a0=k_med*a_ev 
	A1=-A1
	A2=-A2

        RETURN
        END
C
c*******************************************************************

c*******************************************************************
c												qAsca
c-------------------------------------------------------------------
c	Subroutine to calculate complex amplitude scattering matrix Aij
c	by using calculated T-matrix ans three angles:
c	beta_sca, theta_angle, and phi_angle=hi-phi_sca
c	ERR_Aij - accuracy of SAij calculations, 
c     M_max - the actual number of the Tm matrix used
c-------------------------------------------------------------------
	subroutine qAsca(Beta_sca, theta_angle,phi_angle,
     +                 Mmax, N, A11, A12, A21, A22, ERR_Aij,M_max)
c
	IMPLICIT real*16(A-H,O-Z)
 
 	real*16 qzu,Q(60),S(60)

      COMPLEX*32 TE_0,TM_0,TE_m,TM_m
c      complex*32 SA11, SA12, SA21, SA22, A11, A12, A21, A22
      complex*32 A11, A12, A21, A22, A11m, A12m, A21m, A22m
	COMMON /QS/ Q,S
	COMMON /CPQM/ PE(60),PM(60),QE(60),QM(60)
	COMMON /CDMN/ D(60)
c-------------------------------------------------------------------
	pi=qacos(-1.q0)
	
	BETTA=beta_sca*PI/180.q0
	
	CBET=qcos(BETTA)
	SBET=qsin(BETTA)

c    Define the cos(polar orientational angle)
	ZU=qcos(PI*theta_angle/180.q0)
	
	SZU=qsin(PI*theta_angle/180.q0)

c    Define an azimuthal orientational angle

	phi=pi*phi_angle/180.q0

	CPHI=qcos(phi)
	SPHI=qsin(phi)
	ZT=CBET*ZU+SBET*SZU*CPHI
	if(abs(ZT).lt.1.q0) then
	   SZT=qsqrt(1.q0-ZT*ZT)
	   CF=(CBET*SZU-SBET*ZU*CPHI)/SZT
	   SF=SBET*SPHI/SZT
	else
	   SZT=0.q0
         CF=1.q0
	   SF=0.q0
	end if

c	Calculate amplitude matrix for a given M

c-------------------------------------------------------------------
c	Start: Set M=0
	M=0
C|==================================================================
C|      Set Mueller matrices FSij and amplitudes SFij to be zero     
C|==================================================================
c	Set all sums of amplitude matrices over M to be zero
	A11=(0.q0,0.q0)
	A12=(0.q0,0.q0)
	A21=(0.q0,0.q0)
      A22=(0.q0,0.q0)
c
50	 CONTINUE
c------------------------------------------------------------------
C  Calculatng the amplitude matrices SAij in the particle frame
C  for all orientations of the ncident wave vector
c  defined by the Gauss quadratures for ZU and HI/PI.        
C  The angles(ÒHETA,PHI)  define the scattering difection 
C  in the particle symmetry frame,  ZT and SZT are cos and sin THETA,
c   GAM is the angle between the scattering planes 
c   in the particle and the incident wave frames 
C  For details: see scheme in Khebtsov, JQSRT, 2022 https://doi.org/10.1016/j.jqsrt.2022.108069
c------------------------------------------------------------------
C	Calculate the size of T_m matrix for given N=n and M=m
	NM=N-M+1
	if(M.eq.0) NM=N
c------------------------------------------------------------------
c	Calculate angular function for a given orientation angle theta
c	Calculate normalization factors DMN
	qzu=ZU
	CALL qQSM(M,NM,qzu)

	CALL qDMN(M,NM)
C==================================================================
C  Caluclate Amplitudes Pnm and Qnm by using qPQM routine for given N and M
C==================================================================
c***	CALL qPQM(M,N)
	CALL qPQM(M,NM,N)
C==================================================================
c
	CALL qAm(M,N,ZT,CF,SF,A11m,A12m,A21m,A22m)

c	Add partial M-amplitudes to the sums over M
c
      A11=A11 + A11m
	A12=A12 + A12m
	A21=A21 + A21m
	A22=A22 + A22m

c	Check accuracy for M>1	
	IF (M.LE.1) GO TO 80
	
	T11=A11*conjg(A11)
	T12=A12*conjg(A12)
	T21=A21*conjg(A21)
	T22=A22*conjg(A22)

	e11=A11m*conjg(A11m)
	e12=A12m*conjg(A12m)
	e21=A21m*conjg(A21m)
	e22=A22m*conjg(A22m)

	if(T11.ne.0.q0) err11= abs(e11/T11)
	if(T12.ne.0.q0) err12= abs(e12/T12)
	if(T21.ne.0.q0) err21= abs(e21/T21)
	if(T22.ne.0.q0) err22= abs(e22/T22)

	err_max=max(err11,err12,err21,err22)
	 
	if(err_max.lt.ERR_Aij) then
	  M_max=M
	  go to 81
	end if

80	continue
	M=M+1

c	write(*,*) 'M=',M

	if (M.GT.Mmax) then
	   M_max=M-1
	   go to 81
	end if

	GO TO 50

81	continue

	Return
	end
c*******************************************************************

C------------------------------------------------------------
C                  Subroutine       qAm
C------------------------------------------------------------
C Calculation of the partial'm' amplitude matrix for a given M
c input: M,N, ZT=cos(THETA), CF=cos(FI),SF=sin(FI) in the particle frame 
C------------------------------------------------------------
	SUBROUTINE qAm(M,N,ZT,CF,SF,A11m,A12m,A21m,A22m)
	IMPLICIT real*16(A-H,O-Z)
	complex*32 PE,PM,QE,QM,AI1,AIM,CS,CSM,
     +             DA,A11m,A12m,A21m,A22m
	real*16 qzt,Q(60),S(60)
	COMMON /QS/ Q,S
	COMMON /CPQM/ PE(60),PM(60),QE(60),QM(60)
	COMMON /CDMN/ D(60)
c
	AI1=(0.q0,1.q0)
	PI=qacos(-1.q0)
c	CALL qDMN(M,NM)
C------------------------------------------------------------
	A11m=(0.q0,0.q0)
	A12m=(0.q0,0.q0)
	A21m=(0.q0,0.q0)
	A22m=(0.q0,0.q0)
40      NM=N-M+1
	IF (M.EQ.0) NM=N
	M2=M*2
	AIM=1.q0
	qzt=ZT
	CALL qQSM(M,NM,qzt)
	IF (M.NE.0) GO TO 42
C------------------------------------------------------------
C                                                      M=0
C------------------------------------------------------------
	AIM=AI1
	DO 41 N1=1,N
	DA=D(N1)/AIM
	DS=0.q0
	DC=1.q0
	A11m=A11m+DA*S(N1)*PE(N1)*DC/AI1
	A12m=A12m+DA*S(N1)*PM(N1)*DS/AI1
	A21m=A21m+DA*S(N1)*QM(N1)*DS
	A22m=A22m+DA*S(N1)*QE(N1)*DC
41      AIM=AIM*AI1
	GO TO 45
C------------------------------------------------------------
C                                                      M>0
C------------------------------------------------------------
42      CSM=(1.q0,0.q0)
	CS=CF+AI1*SF
	DO 43 N1=1,M
	AIM=AIM*AI1
43      CSM=CSM*CS
	DO 44 N1=1,NM
	DA=D(N1)/AIM*10.q0**(-M2)
	DS=imag(CSM)
	DC=real(CSM)
C------------------------------------------------------------
C       A11=A1, A12=A4, A21=A3, A22=A2
C------------------------------------------------------------
	A11m=A11m+DA*(S(N1)*PE(N1)/AI1-M*Q(N1)*QM(N1))*DC
	A12m=A12m+DA*(S(N1)*PM(N1)/AI1+M*Q(N1)*QE(N1))*DS
	A21m=A21m+DA*(-M*Q(N1)*PE(N1)/AI1+S(N1)*QM(N1))*DS
	A22m=A22m+DA*(M*Q(N1)*PM(N1)/AI1+S(N1)*QE(N1))*DC
44      AIM=AIM*AI1
45      CONTINUE
C------------------------------------------------------------
C Change signs of amplitude matrix according to Mischchenko book
C------------------------------------------------------------
	A11m=-A11m
	A12m=-A12m
	A21m=-A21m
	A22m=-A22m
	RETURN
	END
C
c*******************************************************************

c--------------------------------------------------------------------
c                                                         qPQM     
c--------------------------------------------------------------------
c	Calculation of the partial scattering amplitudes
c	PE=P_even, PM=P_odd, QE=Q_even, QM=Q_odd
c	m is the Tm matrix order, 
c     Nm=N-m+1(m>0) or Nm=N(m=0) is the Tm matriz dimension
c	INPUT: COMMON /T/ TE_0 , TM_0, TE_m, TM_m 
c	Output:COMMON /CPQM/ PE(N), PM(N), QE(N), QM(N)
c     COMMON /CKKK/ KSCM,KOCM,KSCE,KOCE,KOCX,KSCX
c	These are the extinction and scattering efficiencies for TE and TM polarizations
c	TE: the electric vector is perpendicular to the plane (k,a), a - symmetry axis
c	TM: the electric vector is parallel to the plane (k,a)	
c--------------------------------------------------------------------
        SUBROUTINE qPQM(M,NM,NN)
        IMPLICIT real*16(A-H,O-Z)
        complex*32 PE,PM,QE,QM,
     +             AE(60),AM(60),BE(60),BM(60),
     +             AI1,AIM,AIN,KOCM,KOCE,KOCX
        real*16 KSCM,KSCE,KSCX

        COMPLEX*32 TE_0,TM_0,TE_m,TM_m

        complex*32 TE(60,60),TM(60,60)
 
        common /T/ TE_0(1000),TM_0(1000),TE_m(74000),TM_m(74000)
	  
	  real*16 Q(60),S(60)
        COMMON /QS/ Q,S

        COMMON /CPQM/ PE(60),PM(60),QE(60),QM(60)
        COMMON /CDMN/ D(60)
        COMMON /CKKK/ KSCM,KOCM,KSCE,KOCE,KOCX,KSCX
c
        AI1=(0.q0,1.q0)

c	NM for M>0	  
c****	  NM=NN-M+1
c****	  if(M.eq.0) NM=NN
	   
	  IF (M.NE.0) GO TO 12
        AIM=AI1
        NL=NN/2
	  NL2=NL*NL	

c--------------------------------------------------------------------
c	Case m=0
c--------------------------------------------------------------------
		IM0=0
        DO 20 N1=1,NL
        DO 20 N2=1,NL
		IM0=IM0+1
          TE(N1,N2)=TE_0(IM0)
		TM(N1,N2)=TM_0(IM0)

20      continue

        DO 201 N1=1,NL
        DO 201 N2=1,NL
        KNL=N2+NL
		IM0=IM0+1
				
        TE(N1,KNL)=TE_0(IM0)
        TM(N1,KNL)=TM_0(IM0)

201		continue

        DO 10 L=1,NM
        AE(L)=-4.q0*AIM*S(L)
        BE(L)=-AI1*AE(L)
10      AIM=AIM*AI1
        DO 11 L=1,NL
        L2=2*L
        L1=L2-1
        PE(L1)=(0.q0,0.q0)
        PE(L2)=(0.q0,0.q0)
        QE(L1)=(0.q0,0.q0)
        QE(L2)=(0.q0,0.q0)
        DO 11 K=1,NL

        K2=2*K
        K1=K2-1
        KNL=K+NL
        IF (K.EQ.L) KOCX=KOCX+TE(L,K)+TM(L,K)+
     +                   TE(L,KNL)+TM(L,KNL)
        KSCX=KSCX+(TE(L,K)*conjg(TE(L,K))+TM(L,K)*
     +       conjg(TM(L,K)))*D(2*L)/D(2*K)
        KSCX=KSCX+(TE(L,KNL)*conjg(TE(L,KNL))+TM(L,KNL)*
     +       conjg(TM(L,KNL)))*D(2*L-1)/D(2*K-1)
        PE(L1)=PE(L1)+TE(L,KNL)*AE(K1)
        PE(L2)=PE(L2)+TE(L,K)*AE(K2)
        QE(L1)=QE(L1)+TM(L,KNL)*BE(K1)
11      QE(L2)=QE(L2)+TM(L,K)*BE(K2)
        DO 100 N=1,NM
        KSCE=KSCE+D(N)*PE(N)*conjg(PE(N))
        KSCM=KSCM+D(N)*QE(N)*conjg(QE(N))
        KOCE=KOCE+D(N)*conjg(AE(N))*PE(N)
        KOCM=KOCM+D(N)*conjg(BE(N))*QE(N)
100     CONTINUE
        RETURN

c--------------------------------------------------------------------
c	Case M>0
c--------------------------------------------------------------------
12      AIM=1.q0
        M2=M*2
        DO 13 L=1,M
13      AIM=AIM*AI1
        AIM=4.q0*AIM
        AIN=AIM
        DO 14 L=1,NM,2
        AE(L)=-AIM*S(L)
        AM(L)=-AIM*M*Q(L)
        BE(L)=-AE(L)*AI1
        BM(L)=AM(L)*AI1
14      AIM=-AIM
        DO 15 L=2,NM,2
        AE(L)=AIN*M*Q(L)
        AM(L)=-AIN*S(L)
        BE(L)=-AI1*AE(L)
        BM(L)=AI1*AM(L)
15      AIN=-AIN
c--------------------------------------------------------------------
c	Initial index for TE_m and TM_m arrays

		IM=(NN+1)*(M-1)*(NN+1-M)+((M-1)*M*(2*M-1))/6
 
      DO 202 N1=1,NM
        DO 202 N2=1,NM
		IM=IM+1
		TE(N1,N2)=TE_m(IM)	
		TM(N1,N2)=TM_m(IM)

202	continue

        DO 18 N=1,NM
        L1=2*((N+1)/2-N/2)+1
        PE(N)=(0.q0,0.q0)
        PM(N)=(0.q0,0.q0)
        QE(N)=(0.q0,0.q0)
        QM(N)=(0.q0,0.q0)
        DO 18 K=1,NM
        L=L1+(K+1)/2-K/2
        IF (N.EQ.K) KOCX=KOCX+2.q0*(TE(N,K)+TM(N,K))
        KSCX=KSCX+2.q0*(TE(N,K)*conjg(TE(N,K))+TM(N,K)*
     +       conjg(TM(N,K)))*D(N)/D(K)
        PE(N)=PE(N)+TE(N,K)*AE(K)
        QM(N)=QM(N)+TM(N,K)*BM(K)
        GO TO (16,17,17,16),L
16      PM(N)=PM(N)+TE(N,K)*AM(K)
        QE(N)=QE(N)+TM(N,K)*BE(K)
        GO TO 18
17      PM(N)=PM(N)-TE(N,K)*AM(K)
        QE(N)=QE(N)-TM(N,K)*BE(K)

18      NC=NC+1
c--------------------------------------------------------------------
c	Rewriting the even elements PE into QM, QM into PE and otherwise
        DO 19 N=2,NM,2
        AIM=PE(N)
        AIN=PM(N)
        PE(N)=QM(N)
        PM(N)=QE(N)
        QM(N)=AIM
        QE(N)=AIN
        AIM=AE(N)
        AIN=AM(N)
        AE(N)=BM(N)
        AM(N)=BE(N)
        BM(N)=AIM
19      BE(N)=AIN
c-------------------------------------------------------------------
c	Calculating the light scattering efficiencies
        DO 101 N=1,NM
        KSCM=KSCM+D(N)*(PM(N)*conjg(PM(N))+
     +       QE(N)*conjg(QE(N)))*10.q0**(-M2)
        KSCE=KSCE+D(N)*(PE(N)*conjg(PE(N))+
     +       QM(N)*conjg(QM(N)))*10.q0**(-M2)
        KOCM=KOCM+D(N)*(PM(N)*conjg(AM(N))+
     +       QE(N)*conjg(BE(N)))*10.q0**(-M2)
        KOCE=KOCE+D(N)*(PE(N)*conjg(AE(N))+
     +       QM(N)*conjg(BM(N)))*10.q0**(-M2)
101     CONTINUE
        RETURN
        END
C************************************************************

C************************************************************
C                                                qDMN
C------------------------------------------------------------
c	Calculation of the normalization constants
C	D=C*(2*N+1)*(N-M)[*10^(2*M)/(4*N*(N+1)*(N+M)]), 
c	where C=1 for M=0 and C=2 for M>0, N=NM+M-1 (M>0), N=NM(M=0)
c	OUTPUT: COMMON /CDMN/ D(N)	
C------------------------------------------------------------
	SUBROUTINE qDMN(M,NM)
	IMPLICIT real*16(A-H,O-Z)
	COMMON /CDMN/ D(60)
	IF (M.EQ.0) GO TO 30
C------------------------------------------------------------
C                                                      M>0
C------------------------------------------------------------
	M2=M*2
	FM=1.q0
	DO 10 I=1,M2
10      FM=FM*I
	D(1)=(M2+1.q0)/M2/(M+1.q0)/FM*10.q0**(M2)
	NM1=NM-1
	DO 20 N=1,NM1
	MPN=M+N
20      D(N+1)=D(N)*N*(2.q0*MPN+1.q0)*(MPN-1.q0)/(MPN+1.q0)/
     +         (MPN+M)/(2.q0*MPN-1.q0)
	RETURN
C------------------------------------------------------------
C                                                      M=0
C------------------------------------------------------------
30      DO 40 N=1,NM
40      D(N)=(2*N+1)/4.q0/N/(N+1)
	RETURN
	END

c

c***************************************************************
c           subroutine g34jCh1new(wavel,Rm,ref,eps)
c     new version: the data for 620-700 nm has been corrected slightly to ensure 
c	smooth transition from Otter to Johnson-Christy data
c***************************************************************
C     Calculations of real(n)=refre and imag(n)=imag(n)
c     for  GOLD particles according to G34J-Ch spline:
C      210.4-49.9 G03(Irani), 480-640 G04(Otter), 657.4-1605, Johnson and Christy 
c     with account for size-limiting effect.
c     wavel - is the input wavelength nanometers,
c     Rm    - is the radius of the gold particle in nanometers
c     refre and refim - are ouput parameters: real(n) and imag(n)
c     Correction to the size limiting effect:
c     eps1=eps1(bulk),    eps2=eps2(corrected) (ONLY eps2 is corrected!!!)
c     eps2=eps2(bulk)+(5./Rm)*)(wavel/520.)**3
c     The corrected refre and refim are calculated by:
c     refre=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)+eps1)/2.q0)
c     refim=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)-eps1)/2.q0)
c	For details, see Khlebtsov et al., J Phys Chem C, 2011.
C------------------------------------------------------------
       subroutine g34JCh1new(wavel,Rm,ref,eps)
      IMPLICIT real*16(A-H,O-Z)
      real*16 MN,MK
	complex*32 ref,eps,aim
      DIMENSION XN(100),YN(100),MN(100),XK(100),YK(100),MK(100)
c
	aim=(0.q0,1.q0)
c     
      xx=wavel

c     Spline G34JCh parameters
C
	XN(	1	)=	230.4
	XN(	2	)=	239.4
	XN(	3	)=	251.3
	XN(	4	)=	255.4
	XN(	5	)=	259.7
	XN(	6	)=	264.2
	XN(	7	)=	268.5
	XN(	8	)=	273.7
	XN(	9	)=	280.1
	XN(	10	)=	289.9
	XN(	11	)=	295.2
	XN(	12	)=	303.5
	XN(	13	)=	312.5
	XN(	14	)=	320.2
	XN(	15	)=	332.7
	XN(	16	)=	341.4
	XN(	17	)=	353.8
	XN(	18	)=	358.4
	XN(	19	)=	363.8
	XN(	20	)=	374.2
	XN(	21	)=	384
	XN(	22	)=	391
	XN(	23	)=	403.8
	XN(	24	)=	414.1
	XN(	25	)=	424.7
	XN(	26	)=	430
	XN(	27	)=	440
	XN(	28	)=	450.5
	XN(	29	)=	460.4
	XN(	30	)=	469.9
	XN(	31	)=	480
	XN(	32	)=	500
	XN(	33	)=	520
	XN(	34	)=	540
	XN(	35	)=	560
	XN(	36	)=	580
	XN(	37	)=	600
	XN(	38	)=	620
	XN(	39	)=	640
	XN(	40	)=	657.4
	XN(	41	)=	680
	XN(	42	)=	702.3
	XN(	43	)=	753.7
	XN(	44	)=	818.5
	XN(	45	)=	889.2
	XN(	46	)=	981
	XN(	47	)=	1084.2
	XN(	48	)=	1211.8
	XN(	49	)=	1388.8
	XN(	50	)=	1605.2

C
	YN(	1	)=	1.293
	YN(	2	)=	1.296
	YN(	3	)=	1.313
	YN(	4	)=	1.323
	YN(	5	)=	1.34
	YN(	6	)=	1.363
	YN(	7	)=	1.39
	YN(	8	)=	1.432
	YN(	9	)=	1.487
	YN(	10	)=	1.543
	YN(	11	)=	1.573
	YN(	12	)=	1.611
	YN(	13	)=	1.646
	YN(	14	)=	1.679
	YN(	15	)=	1.708
	YN(	16	)=	1.707
	YN(	17	)=	1.661
	YN(	18	)=	1.64
	YN(	19	)=	1.625
	YN(	20	)=	1.613
	YN(	21	)=	1.615
	YN(	22	)=	1.613
	YN(	23	)=	1.609
	YN(	24	)=	1.602
	YN(	25	)=	1.59
	YN(	26	)=	1.573
	YN(	27	)=	1.545
	YN(	28	)=	1.501
	YN(	29	)=	1.442
	YN(	30	)=	1.364
	YN(	31	)=	1.254
	YN(	32	)=	0.8031
	YN(	33	)=	0.5264
	YN(	34	)=	0.3772
	YN(	35	)=	0.3054
	YN(	36	)=	0.2524
	YN(	37	)=	0.2113
	YN(	38	)=	0.175
	YN(	39	)=	0.151
	YN(	40	)=	0.14
	YN(	41	)=	0.13
	YN(	42	)=	0.13
	YN(	43	)=	0.14
	YN(	44	)=	0.16
	YN(	45	)=	0.17
	YN(	46	)=	0.22
	YN(	47	)=	0.27
	YN(	48	)=	0.35
	YN(	49	)=	0.43
	YN(	50	)=	0.56
          
C

	MN(	1	)=	0.0001082668
	MN(	2	)=	0.0001082668
	MN(	3	)=	9.003636q-005
	MN(	4	)=	0.000461763
	MN(	5	)=	0.0002232546
	MN(	6	)=	0.0002290781
	MN(	7	)=	0.00045846
	MN(	8	)=	0.0002098571
	MN(	9	)=	-0.0006486994
	MN(	10	)=	0.0002446948
	MN(	11	)=	-0.0002558822
	MN(	12	)=	-9.989061q-005
	MN(	13	)=	0.0001603902
	MN(	14	)=	-0.0002697623
	MN(	15	)=	-0.0001704703
	MN(	16	)=	-0.0004608753
	MN(	17	)=	-5.132801q-005
	MN(	18	)=	0.0005058203
	MN(	19	)=	0.0001563588
	MN(	20	)=	0.0001991589
	MN(	21	)=	-0.0001555762
	MN(	22	)=	4.812338q-005
	MN(	23	)=	-7.636107q-005
	MN(	24	)=	6.886322q-005
	MN(	25	)=	-0.0004534667
	MN(	26	)=	0.0002334526
	MN(	27	)=	-0.0002294858
	MN(	28	)=	-0.0001208013
	MN(	29	)=	-0.0003309605
	MN(	30	)=	5.59582q-005
	MN(	31	)=	-0.00149829
	MN(	32	)=	0.0009854208
	MN(	33	)=	0.000169607
	MN(	34	)=	0.0002486505
	MN(	35	)=	-3.208495q-006
	MN(	36	)=	4.618334q-005
	MN(	37	)=	-3.024983q-006
	MN(	38	)=	3.791668q-005
	MN(	39	)=	3.585824q-005
	MN(	40	)=	-1.932519q-006
	MN(	41	)=	2.9597q-005
	MN(	42	)=	1.826824q-006
	MN(	43	)=	4.630904q-006
	MN(	44	)=	-7.493589q-006
	MN(	45	)=	1.028974q-005
	MN(	46	)=	-4.303401q-006
	MN(	47	)=	3.611746q-006
	MN(	48	)=	-2.886351q-006
	MN(	49	)=	1.398981q-006
	MN(	50	)=	1.398981q-006
      
C
	XK(	1	)=	230.4
	XK(	2	)=	239.4
	XK(	3	)=	251.3
	XK(	4	)=	255.4
	XK(	5	)=	259.7
	XK(	6	)=	264.2
	XK(	7	)=	268.5
	XK(	8	)=	273.7
	XK(	9	)=	280.1
	XK(	10	)=	289.9
	XK(	11	)=	295.2
	XK(	12	)=	303.5
	XK(	13	)=	312.5
	XK(	14	)=	320.2
	XK(	15	)=	332.7
	XK(	16	)=	341.4
	XK(	17	)=	353.8
	XK(	18	)=	358.4
	XK(	19	)=	363.8
	XK(	20	)=	374.2
	XK(	21	)=	384
	XK(	22	)=	391
	XK(	23	)=	403.8
	XK(	24	)=	414.1
	XK(	25	)=	424.7
	XK(	26	)=	430
	XK(	27	)=	440
	XK(	28	)=	450.5
	XK(	29	)=	460.4
	XK(	30	)=	469.9
	XK(	31	)=	480
	XK(	32	)=	500
	XK(	33	)=	520
	XK(	34	)=	540
	XK(	35	)=	560
	XK(	36	)=	580
	XK(	37	)=	600
	XK(	38	)=	620
	XK(	39	)=	640
	XK(	40	)=	657.4
	XK(	41	)=	680
	XK(	42	)=	702.3
	XK(	43	)=	753.7
	XK(	44	)=	818.5
	XK(	45	)=	889.2
	XK(	46	)=	981
	XK(	47	)=	1084.2
	XK(	48	)=	1211.8
	XK(	49	)=	1388.8
	XK(	50	)=	1605.2
     
C
	YK(	1	)=	1.477
	YK(	2	)=	1.543
	YK(	3	)=	1.641
	YK(	4	)=	1.676
	YK(	5	)=	1.715
	YK(	6	)=	1.755
	YK(	7	)=	1.79
	YK(	8	)=	1.824
	YK(	9	)=	1.843
	YK(	10	)=	1.864
	YK(	11	)=	1.87
	YK(	12	)=	1.874
	YK(	13	)=	1.877
	YK(	14	)=	1.87
	YK(	15	)=	1.84
	YK(	16	)=	1.811
	YK(	17	)=	1.792
	YK(	18	)=	1.804
	YK(	19	)=	1.823
	YK(	20	)=	1.854
	YK(	21	)=	1.869
	YK(	22	)=	1.877
	YK(	23	)=	1.887
	YK(	24	)=	1.882
	YK(	25	)=	1.866
	YK(	26	)=	1.853
	YK(	27	)=	1.838
	YK(	28	)=	1.806
	YK(	29	)=	1.767
	YK(	30	)=	1.728
	YK(	31	)=	1.73
	YK(	32	)=	1.818
	YK(	33	)=	2.1277
	YK(	34	)=	2.452
	YK(	35	)=	2.7501
	YK(	36	)=	3.0106
	YK(	37	)=	3.2411
	YK(	38	)=	3.44
	YK(	39	)=	3.62
	YK(	40	)=	3.76
	YK(	41	)=	3.94
	YK(	42	)=	4.13
	YK(	43	)=	4.542
	YK(	44	)=	5.083
	YK(	45	)=	5.663
	YK(	46	)=	6.35
	YK(	47	)=	7.15
	YK(	48	)=	8.145
	YK(	49	)=	9.519
	YK(	50	)=	11.21
C
	MK(	1	)=	0.0001090988
	MK(	2	)=	0.0001090988
	MK(	3	)=	-1.096293q-005
	MK(	4	)=	0.0002098452
	MK(	5	)=	-6.549076q-005
	MK(	6	)=	-0.0001855193
	MK(	7	)=	-0.0002177071
	MK(	8	)=	-0.0008985389
	MK(	9	)=	8.750596q-005
	MK(	10	)=	-0.0002081646
	MK(	11	)=	-0.0001199258
	MK(	12	)=	5.594176q-005
	MK(	13	)=	-0.0002035245
	MK(	14	)=	-0.0001506855
	MK(	15	)=	-0.0001032507
	MK(	16	)=	7.601288q-005
	MK(	17	)=	0.0006852496
	MK(	18	)=	0.0001314306
	MK(	19	)=	-5.959464q-005
	MK(	20	)=	-0.0001974043
	MK(	21	)=	-1.081461q-005
	MK(	22	)=	-4.084141q-006
	MK(	23	)=	-0.000150955
	MK(	24	)=	-5.570012q-005
	MK(	25	)=	-0.0002132856
	MK(	26	)=	0.0003230994
	MK(	27	)=	-0.0003039385
	MK(	28	)=	-5.264124q-006
	MK(	29	)=	-0.000196416
	MK(	30	)=	0.0007029396
	MK(	31	)=	1.290514q-005
	MK(	32	)=	0.0008667648
	MK(	33	)=	-0.0001544616
	MK(	34	)=	-2.992337q-005
	MK(	35	)=	-0.0001188422
	MK(	36	)=	-5.870466q-005
	MK(	37	)=	-9.634229q-005
	MK(	38	)=	-2.992621q-005
	MK(	39	)=	-6.745523q-005
	MK(	40	)=	-4.594196q-006
	MK(	41	)=	4.659774q-005
	MK(	42	)=	-3.3507q-005
	MK(	43	)=	1.69652q-005
	MK(	44	)=	-3.412724q-006
	MK(	45	)=	-1.478193q-005
	MK(	46	)=	7.900981q-006
	MK(	47	)=	-1.111366q-006
	MK(	48	)=	-2.133222q-007
	MK(	49	)=	3.458122q-007
	MK(	50	)=	3.458122q-007

C     Basic calculations
        J=2
10      IF (XX.LE.XN(J)) GO TO 20
        J=J+1
        GO TO 10
20      I=J
        XN0=XN(I)-XX
        XN1=XX-XN(I-1)
        H=XN(I)-XN(I-1)
        H2=H*H/6.0
        refre=MN(I-1)*XN0*XN0*XN0/6./H+MN(I)*XN1*XN1*XN1/6./H+
     +     (YN(I-1)-MN(I-1)*H2)*XN0/H+(YN(I)-MN(I)*H2)*XN1/H

         J=2
30      IF (XX.LE.XK(J)) GO TO 40
        J=J+1
        GO TO 30
40      I=J
        XN0=XK(I)-XX
        XN1=XX-XK(I-1)
        H=XK(I)-XK(I-1)
        H2=H*H/6.0
        refim=MK(I-1)*XN0*XN0*XN0/6./H+MK(I)*XN1*XN1*XN1/6./H+
     +    (YK(I-1)-MK(I-1)*H2)*XN0/H+(YK(I)-MK(I)*H2)*XN1/H

c     Correction to the size limiting effect
c     Step1: calculation of eps=eps1+ieps2
       eps1=refre*refre-refim*refim
       eps2=2.q0*refre*refim
c      Step 2: Calculation of the coorection term
c      by using ONLY imaginary contribution!!!! 
       deps=(5.q0/Rm)*(wavel/520.q0)**3.q0 
       eps2=eps2+deps
c      Step 3: Recalculation of tefre and refim: Bohren-Huffman, p.281
       refre=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)+eps1)/2.q0)
       refim=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)-eps1)/2.q0)
	 ref=refre+aim*refim
	 eps=eps1+aim*eps2
c---------------------------------------------------------------------
       RETURN
       END
c***************************************************************
C
c***************************************************************

C                                                         qGAU_short
c---------------------------------------------------------------------
c	A short version of qGAU for orientation averaging of the Mueller matrix
c---------------------------------------------------------------------
        SUBROUTINE qGAU_short(NG,W,X1)
        real*16 W(48),X1(48)
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 6                                    ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
6       IF (NG.GT.6) GO TO 12
        NG=6
C
	 w( 1)=0.4717533638651182719461596148501712q-01
	 w( 2)=0.1069393259953184309602547181939963q+00
	 w( 3)=0.1600783285433462263346525295433591q+00
	 w( 4)=0.2031674267230659217490644558097983q+00
	 w( 5)=0.2334925365383548087608498989248780q+00
	 w( 6)=0.2491470458134027850005624360429512q+00
c 
	 X1( 1)=0.9815606342467192506905490901492808q+00
	 X1( 2)=0.9041172563704748566784658661190963q+00
	 X1( 3)=0.7699026741943046870368938332128181q+00
	 X1( 4)=0.5873179542866174472967024189405343q+00
	 X1( 5)=0.3678314989981801937526915366437175q+00
	 X1( 6)=0.1252334085114689154724413694638531q+00     

C
        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 12                                   ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
12      IF (NG.GT.12) GO TO 24
        NG=12
C
	 w( 1)=0.1234122979998719954680566707003725q-01
	 w( 2)=0.2853138862893366318130781595187828q-01
	 w( 3)=0.4427743881741980616860274821133821q-01
	 w( 4)=0.5929858491543678074636775850010859q-01
	 w( 5)=0.7334648141108030573403361525311648q-01
	 w( 6)=0.8619016153195327591718520298374266q-01
	 w( 7)=0.9761865210411388826988066446424718q-01
	 w( 8)=0.1074442701159656347825773424466063q+00
	 w( 9)=0.1155056680537256013533444839067836q+00
	 w(10)=0.1216704729278033912044631534762624q+00
	 w(11)=0.1258374563468282961213753825111837q+00
	 w(12)=0.1279381953467521569740561652246955q+00
c
	 X1( 1)=0.9951872199970213601799974097007367q+00
	 X1( 2)=0.9747285559713094981983919930081690q+00
	 X1( 3)=0.9382745520027327585236490017087215q+00
	 X1( 4)=0.8864155270044010342131543419821967q+00
	 X1( 5)=0.8200019859739029219539498726697452q+00
	 X1( 6)=0.7401241915785543642438281030999784q+00
	 X1( 7)=0.6480936519369755692524957869107476q+00
	 X1( 8)=0.5454214713888395356583756172183724q+00
	 X1( 9)=0.4337935076260451384870842319133498q+00
	 X1(10)=0.3150426796961633743867932913198102q+00
	 X1(11)=0.1911188674736163091586398207570697q+00
	 X1(12)=0.6405689286260562608504308262474504q-01     


C
        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 24                                   ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
24      IF (NG.GT.24) GO TO 48
        NG=24
C
       w( 1)=0.3153346052305838632677311543891515q-02
	 w( 2)=0.7327553901276262102383979621786534q-02
	 w( 3)=0.1147723457923453948959266760909167q-01
	 w( 4)=0.1557931572294384872817695583446033q-01
	 w( 5)=0.1961616045735552781446071965221279q-01
	 w( 6)=0.2357076083932437914051930137844923q-01
	 w( 7)=0.2742650970835694820007383626250591q-01
	 w( 8)=0.3116722783279808890206575684635440q-01
	 w( 9)=0.3477722256477043889254858596380243q-01
	 w(10)=0.3824135106583070631721725652371569q-01
	 w(11)=0.4154508294346474921405882236106477q-01
	 w(12)=0.4467456085669428041944858712585040q-01
	 w(13)=0.4761665849249047482590662347892982q-01
	 w(14)=0.5035903555385447495780761908786555q-01
	 w(15)=0.5289018948519366709550505626469886q-01
	 w(16)=0.5519950369998416286820349519163536q-01
	 w(17)=0.5727729210040321570515023468470064q-01
	 w(18)=0.5911483969839563574647481743351990q-01
	 w(19)=0.6070443916589388005296923202782048q-01
	 w(20)=0.6203942315989266390419778413759853q-01
	 w(21)=0.6311419228625402565712602275023334q-01
	 w(22)=0.6392423858464818662390620182551554q-01
	 w(23)=0.6446616443595008220650419365770497q-01
	 w(24)=0.6473769681268392250302493873659151q-01
c 
	 X1( 1)=0.9987710072524261186005414915631136q+00
	 X1( 2)=0.9935301722663507575479287508490740q+00
	 X1( 3)=0.9841245837228268577445836000265988q+00
	 X1( 4)=0.9705915925462472504614119838006601q+00
	 X1( 5)=0.9529877031604308607229606660257181q+00
	 X1( 6)=0.9313866907065543331141743801016012q+00
	 X1( 7)=0.9058791367155696728220748356710117q+00
	 X1( 8)=0.8765720202742478859056935548050967q+00
	 X1( 9)=0.8435882616243935307110898445196559q+00
	 X1(10)=0.8070662040294426270825530430245384q+00
	 X1(11)=0.7671590325157403392538554375229691q+00
	 X1(12)=0.7240341309238146546744822334936653q+00
	 X1(13)=0.6778723796326639052118512806759091q+00
	 X1(14)=0.6288673967765136239951649330699947q+00
	 X1(15)=0.5772247260839727038178092385404788q+00
	 X1(16)=0.5231609747222330336782258691375086q+00
	 X1(17)=0.4669029047509584045449288616507985q+00
	 X1(18)=0.4086864819907167299162254958146334q+00
	 X1(19)=0.3487558862921607381598179372704079q+00
	 X1(20)=0.2873624873554555767358864613167976q+00
	 X1(21)=0.2247637903946890612248654401746923q+00
	 X1(22)=0.1612223560688917180564373907834977q+00
	 X1(23)=0.9700469920946269893005395585362452q-01
	 X1(24)=0.3238017096286936203332224315213445q-01     

C
        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 48                                   ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
C
48       if (NG.GT.48) then
	      write(*,*) '!!!   NG>48   !!'
		  stop
		end if		
C
        NG=48
C
	 w( 1)=0.7967920655520124294381434969436714q-03
	 w( 2)=0.1853960788946921732335925350893890q-02
	 w( 3)=0.2910731817934946408410617989400720q-02
	 w( 4)=0.3964554338444686673733415767419639q-02
	 w( 5)=0.5014202742927517692470194969030949q-02
	 w( 6)=0.6058545504235961683316742031729087q-02
	 w( 7)=0.7096470791153865269144160812143404q-02
	 w( 8)=0.8126876925698759217382427707855918q-02
	 w( 9)=0.9148671230783386632584602665207933q-02
	 w(10)=0.1016077053500841575758763695382536q-01
	 w(11)=0.1116210209983849859121326382856207q-01
	 w(12)=0.1215160467108831963518135273666725q-01
	 w(13)=0.1312822956696157263706366690259911q-01
	 w(14)=0.1409094177231486091586162472463661q-01
	 w(15)=0.1503872102699493800587627522209756q-01
	 w(16)=0.1597056290256229138061645679149306q-01
	 w(17)=0.1688547986424517245047754060686194q-01
	 w(18)=0.1778250231604526083761422648607104q-01
	 w(19)=0.1866067962741146738515675862213225q-01
	 w(20)=0.1951908114014502241008522028120468q-01
	 w(21)=0.2035679715433332459524521541727163q-01
	 w(22)=0.2117293989219129898767386719100365q-01
	 w(23)=0.2196664443874434919475638680156228q-01
	 w(24)=0.2273706965832937400134784197749033q-01
	 w(25)=0.2348339908592621984223593266761256q-01
	 w(26)=0.2420484179236469128226733787267704q-01
	 w(27)=0.2490063322248361028838218086833276q-01
	 w(28)=0.2557003600534936149879716794360011q-01
	 w(29)=0.2621234073567241391345796396446336q-01
	 w(30)=0.2682686672559176219805672871415660q-01
	 w(31)=0.2741296272602924282342108748909122q-01
	 w(32)=0.2797000761684833443981857658902251q-01
	 w(33)=0.2849741106508538564559951294580563q-01
	 w(34)=0.2899461415055523654267878127968158q-01
	 w(35)=0.2946108995816790597043633218285850q-01
	 w(36)=0.2989634413632838598438807579440005q-01
	 w(37)=0.3029991542082759379408876420650097q-01
	 w(38)=0.3067137612366914901422883035620429q-01
	 w(39)=0.3101033258631383742324977997063272q-01
	 w(40)=0.3131642559686135581278426671506315q-01
	 w(41)=0.3158933077072716855802074616996422q-01
	 w(42)=0.3182875889441100653475373988553221q-01
	 w(43)=0.3203445623199266321813897747021117q-01
	 w(44)=0.3220620479403025066866711455723255q-01
	 w(45)=0.3234382256857592842877483882894328q-01
	 w(46)=0.3244716371406426936401278844884587q-01
	 w(47)=0.3251611871386883598720549144777837q-01
	 w(48)=0.3255061449236316624196141829728569q-01
c
	 X1( 1)=0.9996895038832307668276901057843654q+00
	 X1( 2)=0.9983643758631816777241494395267193q+00
	 X1( 3)=0.9959818429872092906503990848793565q+00
	 X1( 4)=0.9925439003237626245718922977016011q+00
	 X1( 5)=0.9880541263296237994807627724184511q+00
	 X1( 6)=0.9825172635630146774470458318825480q+00
	 X1( 7)=0.9759391745851364664526010342418798q+00
	 X1( 8)=0.9683268284632642121736593664436509q+00
	 X1( 9)=0.9596882914487425393000680258847348q+00
	 X1(10)=0.9500327177844376357560989484438634q+00
	 X1(11)=0.9393703397527552169318573794316151q+00
	 X1(12)=0.9277124567223086909646904726965465q+00
	 X1(13)=0.9150714231208980742058844661504215q+00
	 X1(14)=0.9014606353158523413192326973096868q+00
	 X1(15)=0.8868945174024204160568774333937861q+00
	 X1(16)=0.8713885059092965028737747644930427q+00
	 X1(17)=0.8549590334346014554627869698935691q+00
	 X1(18)=0.8376235112281871214943028167646929q+00
	 X1(19)=0.8194003107379316755389996242244780q+00
	 X1(20)=0.8003087441391408172287961439713277q+00
	 X1(21)=0.7803690438674332176036044558893177q+00
	 X1(22)=0.7596023411766474987029704411196217q+00
	 X1(23)=0.7380306437444001328511657310964972q+00
	 X1(24)=0.7156768123489676262251441480576248q+00
	 X1(25)=0.6925645366421715613442457698182425q+00
	 X1(26)=0.6687183100439161539525572075727743q+00
	 X1(27)=0.6441634037849671067984123500633016q+00
	 X1(28)=0.6189258401254685703863692869878628q+00
	 X1(29)=0.5930323647775720806835557550437583q+00
	 X1(30)=0.5665104185613971684042501950834627q+00
	 X1(31)=0.5393881083243574362268025973255607q+00
	 X1(32)=0.5116941771546676735855097454288543q+00
	 X1(33)=0.4834579739205963597684056093638952q+00
	 X1(34)=0.4547094221677430086356761480863078q+00
	 X1(35)=0.4254789884073005453648192035699980q+00
	 X1(36)=0.3957976498289086032850002431351642q+00
	 X1(37)=0.3656968614723136350308955939954732q+00
	 X1(38)=0.3352085228926254226163256248056287q+00
	 X1(39)=0.3043649443544963530239297792952033q+00
	 X1(40)=0.2731988125910491414872721557245648q+00
	 X1(41)=0.2417431561638400123279319037406277q+00
	 X1(42)=0.2100313104605672036028471857122957q+00
	 X1(43)=0.1780968823676186027594026251764915q+00
	 X1(44)=0.1459737146548969419891073333433265q+00
	 X1(45)=0.1136958501106659209112080954143281q+00
	 X1(46)=0.8129749546442555899447126297475285q-01
	 X1(47)=0.4881298513604973111195820426127506q-01
	 X1(48)=0.1627674484960296957913456369523845q-01

56	RETURN
	END
C
C******************************************************************
C                                                         qGAU     
C******************************************************************
        SUBROUTINE qGAU(NG)
        real*16 W,X1
        COMMON /CGAU/ W(192),X1(192)
C
C-------------------------------------------------------------------
C                       N --- 6                                    
C-------------------------------------------------------------------
C
6       IF (NG.GT.6) GO TO 12
        NG=6
C
	 w( 1)=0.4717533638651182719461596148501712q-01
	 w( 2)=0.1069393259953184309602547181939963q+00
	 w( 3)=0.1600783285433462263346525295433591q+00
	 w( 4)=0.2031674267230659217490644558097983q+00
	 w( 5)=0.2334925365383548087608498989248780q+00
	 w( 6)=0.2491470458134027850005624360429512q+00
c 
	 X1( 1)=0.9815606342467192506905490901492808q+00
	 X1( 2)=0.9041172563704748566784658661190963q+00
	 X1( 3)=0.7699026741943046870368938332128181q+00
	 X1( 4)=0.5873179542866174472967024189405343q+00
	 X1( 5)=0.3678314989981801937526915366437175q+00
	 X1( 6)=0.1252334085114689154724413694638531q+00     

C
        GO TO 56
C
C-------------------------------------------------------------------
C³                       N --- 12                                   ³
C-------------------------------------------------------------------
C
12      IF (NG.GT.12) GO TO 24
        NG=12
C
	 w( 1)=0.1234122979998719954680566707003725q-01
	 w( 2)=0.2853138862893366318130781595187828q-01
	 w( 3)=0.4427743881741980616860274821133821q-01
	 w( 4)=0.5929858491543678074636775850010859q-01
	 w( 5)=0.7334648141108030573403361525311648q-01
	 w( 6)=0.8619016153195327591718520298374266q-01
	 w( 7)=0.9761865210411388826988066446424718q-01
	 w( 8)=0.1074442701159656347825773424466063q+00
	 w( 9)=0.1155056680537256013533444839067836q+00
	 w(10)=0.1216704729278033912044631534762624q+00
	 w(11)=0.1258374563468282961213753825111837q+00
	 w(12)=0.1279381953467521569740561652246955q+00
c
	 X1( 1)=0.9951872199970213601799974097007367q+00
	 X1( 2)=0.9747285559713094981983919930081690q+00
	 X1( 3)=0.9382745520027327585236490017087215q+00
	 X1( 4)=0.8864155270044010342131543419821967q+00
	 X1( 5)=0.8200019859739029219539498726697452q+00
	 X1( 6)=0.7401241915785543642438281030999784q+00
	 X1( 7)=0.6480936519369755692524957869107476q+00
	 X1( 8)=0.5454214713888395356583756172183724q+00
	 X1( 9)=0.4337935076260451384870842319133498q+00
	 X1(10)=0.3150426796961633743867932913198102q+00
	 X1(11)=0.1911188674736163091586398207570697q+00
	 X1(12)=0.6405689286260562608504308262474504q-01     


C
        GO TO 56
C
C-------------------------------------------------------------------
C³                       N --- 24                                   ³
C-------------------------------------------------------------------
C
24      IF (NG.GT.24) GO TO 48
        NG=24
C
       w( 1)=0.3153346052305838632677311543891515q-02
	 w( 2)=0.7327553901276262102383979621786534q-02
	 w( 3)=0.1147723457923453948959266760909167q-01
	 w( 4)=0.1557931572294384872817695583446033q-01
	 w( 5)=0.1961616045735552781446071965221279q-01
	 w( 6)=0.2357076083932437914051930137844923q-01
	 w( 7)=0.2742650970835694820007383626250591q-01
	 w( 8)=0.3116722783279808890206575684635440q-01
	 w( 9)=0.3477722256477043889254858596380243q-01
	 w(10)=0.3824135106583070631721725652371569q-01
	 w(11)=0.4154508294346474921405882236106477q-01
	 w(12)=0.4467456085669428041944858712585040q-01
	 w(13)=0.4761665849249047482590662347892982q-01
	 w(14)=0.5035903555385447495780761908786555q-01
	 w(15)=0.5289018948519366709550505626469886q-01
	 w(16)=0.5519950369998416286820349519163536q-01
	 w(17)=0.5727729210040321570515023468470064q-01
	 w(18)=0.5911483969839563574647481743351990q-01
	 w(19)=0.6070443916589388005296923202782048q-01
	 w(20)=0.6203942315989266390419778413759853q-01
	 w(21)=0.6311419228625402565712602275023334q-01
	 w(22)=0.6392423858464818662390620182551554q-01
	 w(23)=0.6446616443595008220650419365770497q-01
	 w(24)=0.6473769681268392250302493873659151q-01
c 
	 X1( 1)=0.9987710072524261186005414915631136q+00
	 X1( 2)=0.9935301722663507575479287508490740q+00
	 X1( 3)=0.9841245837228268577445836000265988q+00
	 X1( 4)=0.9705915925462472504614119838006601q+00
	 X1( 5)=0.9529877031604308607229606660257181q+00
	 X1( 6)=0.9313866907065543331141743801016012q+00
	 X1( 7)=0.9058791367155696728220748356710117q+00
	 X1( 8)=0.8765720202742478859056935548050967q+00
	 X1( 9)=0.8435882616243935307110898445196559q+00
	 X1(10)=0.8070662040294426270825530430245384q+00
	 X1(11)=0.7671590325157403392538554375229691q+00
	 X1(12)=0.7240341309238146546744822334936653q+00
	 X1(13)=0.6778723796326639052118512806759091q+00
	 X1(14)=0.6288673967765136239951649330699947q+00
	 X1(15)=0.5772247260839727038178092385404788q+00
	 X1(16)=0.5231609747222330336782258691375086q+00
	 X1(17)=0.4669029047509584045449288616507985q+00
	 X1(18)=0.4086864819907167299162254958146334q+00
	 X1(19)=0.3487558862921607381598179372704079q+00
	 X1(20)=0.2873624873554555767358864613167976q+00
	 X1(21)=0.2247637903946890612248654401746923q+00
	 X1(22)=0.1612223560688917180564373907834977q+00
	 X1(23)=0.9700469920946269893005395585362452q-01
	 X1(24)=0.3238017096286936203332224315213445q-01     
C
        GO TO 56
C
C-------------------------------------------------------------------
C³                       N --- 48                                   ³
C-------------------------------------------------------------------
C
C
48      IF (NG.GT.48) GO TO 96
C
        NG=48
C
	 w( 1)=0.7967920655520124294381434969436714q-03
	 w( 2)=0.1853960788946921732335925350893890q-02
	 w( 3)=0.2910731817934946408410617989400720q-02
	 w( 4)=0.3964554338444686673733415767419639q-02
	 w( 5)=0.5014202742927517692470194969030949q-02
	 w( 6)=0.6058545504235961683316742031729087q-02
	 w( 7)=0.7096470791153865269144160812143404q-02
	 w( 8)=0.8126876925698759217382427707855918q-02
	 w( 9)=0.9148671230783386632584602665207933q-02
	 w(10)=0.1016077053500841575758763695382536q-01
	 w(11)=0.1116210209983849859121326382856207q-01
	 w(12)=0.1215160467108831963518135273666725q-01
	 w(13)=0.1312822956696157263706366690259911q-01
	 w(14)=0.1409094177231486091586162472463661q-01
	 w(15)=0.1503872102699493800587627522209756q-01
	 w(16)=0.1597056290256229138061645679149306q-01
	 w(17)=0.1688547986424517245047754060686194q-01
	 w(18)=0.1778250231604526083761422648607104q-01
	 w(19)=0.1866067962741146738515675862213225q-01
	 w(20)=0.1951908114014502241008522028120468q-01
	 w(21)=0.2035679715433332459524521541727163q-01
	 w(22)=0.2117293989219129898767386719100365q-01
	 w(23)=0.2196664443874434919475638680156228q-01
	 w(24)=0.2273706965832937400134784197749033q-01
	 w(25)=0.2348339908592621984223593266761256q-01
	 w(26)=0.2420484179236469128226733787267704q-01
	 w(27)=0.2490063322248361028838218086833276q-01
	 w(28)=0.2557003600534936149879716794360011q-01
	 w(29)=0.2621234073567241391345796396446336q-01
	 w(30)=0.2682686672559176219805672871415660q-01
	 w(31)=0.2741296272602924282342108748909122q-01
	 w(32)=0.2797000761684833443981857658902251q-01
	 w(33)=0.2849741106508538564559951294580563q-01
	 w(34)=0.2899461415055523654267878127968158q-01
	 w(35)=0.2946108995816790597043633218285850q-01
	 w(36)=0.2989634413632838598438807579440005q-01
	 w(37)=0.3029991542082759379408876420650097q-01
	 w(38)=0.3067137612366914901422883035620429q-01
	 w(39)=0.3101033258631383742324977997063272q-01
	 w(40)=0.3131642559686135581278426671506315q-01
	 w(41)=0.3158933077072716855802074616996422q-01
	 w(42)=0.3182875889441100653475373988553221q-01
	 w(43)=0.3203445623199266321813897747021117q-01
	 w(44)=0.3220620479403025066866711455723255q-01
	 w(45)=0.3234382256857592842877483882894328q-01
	 w(46)=0.3244716371406426936401278844884587q-01
	 w(47)=0.3251611871386883598720549144777837q-01
	 w(48)=0.3255061449236316624196141829728569q-01
c
	 X1( 1)=0.9996895038832307668276901057843654q+00
	 X1( 2)=0.9983643758631816777241494395267193q+00
	 X1( 3)=0.9959818429872092906503990848793565q+00
	 X1( 4)=0.9925439003237626245718922977016011q+00
	 X1( 5)=0.9880541263296237994807627724184511q+00
	 X1( 6)=0.9825172635630146774470458318825480q+00
	 X1( 7)=0.9759391745851364664526010342418798q+00
	 X1( 8)=0.9683268284632642121736593664436509q+00
	 X1( 9)=0.9596882914487425393000680258847348q+00
	 X1(10)=0.9500327177844376357560989484438634q+00
	 X1(11)=0.9393703397527552169318573794316151q+00
	 X1(12)=0.9277124567223086909646904726965465q+00
	 X1(13)=0.9150714231208980742058844661504215q+00
	 X1(14)=0.9014606353158523413192326973096868q+00
	 X1(15)=0.8868945174024204160568774333937861q+00
	 X1(16)=0.8713885059092965028737747644930427q+00
	 X1(17)=0.8549590334346014554627869698935691q+00
	 X1(18)=0.8376235112281871214943028167646929q+00
	 X1(19)=0.8194003107379316755389996242244780q+00
	 X1(20)=0.8003087441391408172287961439713277q+00
	 X1(21)=0.7803690438674332176036044558893177q+00
	 X1(22)=0.7596023411766474987029704411196217q+00
	 X1(23)=0.7380306437444001328511657310964972q+00
	 X1(24)=0.7156768123489676262251441480576248q+00
	 X1(25)=0.6925645366421715613442457698182425q+00
	 X1(26)=0.6687183100439161539525572075727743q+00
	 X1(27)=0.6441634037849671067984123500633016q+00
	 X1(28)=0.6189258401254685703863692869878628q+00
	 X1(29)=0.5930323647775720806835557550437583q+00
	 X1(30)=0.5665104185613971684042501950834627q+00
	 X1(31)=0.5393881083243574362268025973255607q+00
	 X1(32)=0.5116941771546676735855097454288543q+00
	 X1(33)=0.4834579739205963597684056093638952q+00
	 X1(34)=0.4547094221677430086356761480863078q+00
	 X1(35)=0.4254789884073005453648192035699980q+00
	 X1(36)=0.3957976498289086032850002431351642q+00
	 X1(37)=0.3656968614723136350308955939954732q+00
	 X1(38)=0.3352085228926254226163256248056287q+00
	 X1(39)=0.3043649443544963530239297792952033q+00
	 X1(40)=0.2731988125910491414872721557245648q+00
	 X1(41)=0.2417431561638400123279319037406277q+00
	 X1(42)=0.2100313104605672036028471857122957q+00
	 X1(43)=0.1780968823676186027594026251764915q+00
	 X1(44)=0.1459737146548969419891073333433265q+00
	 X1(45)=0.1136958501106659209112080954143281q+00
	 X1(46)=0.8129749546442555899447126297475285q-01
	 X1(47)=0.4881298513604973111195820426127506q-01
	 X1(48)=0.1627674484960296957913456369523845q-01

        GO TO 56
C
C-------------------------------------------------------------------
C³                       N --- 96                                   ³
C-------------------------------------------------------------------
C
96      IF (NG.GT.96) GO TO 192
C
        NG=96
C
	 w( 1)=0.2002510147471267303712312450255703q-03
	 w( 2)=0.4660944855912469603944432543932574q-03
	 w( 3)=0.7322075666250185328339640008748637q-03
	 w( 4)=0.9981598380233457554121719231494520q-03
	 w( 5)=0.1263853132096728285868613120015403q-02
	 w( 6)=0.1529211875019241493105816989986518q-02
	 w( 7)=0.1794164116955045626355485501631911q-02
	 w( 8)=0.2058638858115252727715815761111229q-02
	 w( 9)=0.2322565485048108214442097204508405q-02
	 w(10)=0.2585873625444375300279098771560957q-02
	 w(11)=0.2848493111476581340970402165390526q-02
	 w(12)=0.3110353977068882106680359244296400q-02
	 w(13)=0.3371386467342533589805998887040961q-02
	 w(14)=0.3631521052910285400855749274005996q-02
	 w(15)=0.3890688446257093249913484957572511q-02
	 w(16)=0.4148819619065216099139199471198278q-02
	 w(17)=0.4405845819973666428314679521316033q-02
	 w(18)=0.4661698592528197474902387832688645q-02
	 w(19)=0.4916309793197733507646532723064942q-02
	 w(20)=0.5169611609390087300715718380698600q-02
	 w(21)=0.5421536577428234926249655845129893q-02
	 w(22)=0.5672017600463228374204949713339389q-02
	 w(23)=0.5920987966307842199944703122938390q-02
	 w(24)=0.6168381365179540674609447146559622q-02
	 w(25)=0.6414131907343948572275024532122624q-02
	 w(26)=0.6658174140651556266610738434308557q-02
	 w(27)=0.6900443067961339719780742657597472q-02
	 w(28)=0.7140874164445576807569565961274595q-02
	 w(29)=0.7379403394770533802095458002934840q-02
	 w(30)=0.7615967230147961172624843647609506q-02
	 w(31)=0.7850502665252524421764359380945037q-02
	 w(32)=0.8082947235000432632340753774352802q-02
	 w(33)=0.8313239031184632745769256416255048q-02
	 w(34)=0.8541316718962022605943311383449181q-02
	 w(35)=0.8767119553188207673350963346520281q-02
	 w(36)=0.8990587394595389693608817807139734q-02
	 w(37)=0.9211660725809033543826169770456859q-02
	 w(38)=0.9430280667199013005803355133957693q-02
	 w(39)=0.9646388992560988688871555468912967q-02
	 w(40)=0.9859928144623822671848493861984853q-02
	 w(41)=0.1007084125037888528666258758886661q-01
	 w(42)=0.1027907213622716026310745432162761q-01
       w(43)=0.1048456534294010549872985896280295q-01
       w(44)=0.1068726614043027822182015353540765q-01
       w(45)=0.1088712054232778542665459065331425q-01
       w(46)=0.1108407532035867328136178184595929q-01
       w(47)=0.1127807801852142281099299668499557q-01
       w(48)=0.1146907696705777358968578140016646q-01
       w(49)=0.1165702129621415246785292850221955q-01
       w(50)=0.1184186094979004053271591814062855q-01
       w(51)=0.1202354669846966856400385119887345q-01
       w(52)=0.1220203015293348920546093199411790q-01
       w(53)=0.1237726377674593292623165356797658q-01
       w(54)=0.1254920089901601458978069148872322q-01
       w(55)=0.1271779572682741807446091689663609q-01
       w(56)=0.1288300335743474788966926602598000q-01
       w(57)=0.1304477979022269909350163852581688q-01
       w(58)=0.1320308193842496002934176412517347q-01
       w(59)=0.1335786764059972644620390859976834q-01
       w(60)=0.1350909567185877043631007883710086q-01
       w(61)=0.1365672575484707329799543649929899q-01
       w(62)=0.1380071857047009789665883850857899q-01
       w(63)=0.1394103576836584333464095632054565q-01
       w(64)=0.1407763997711889273569915932555554q-01
       w(65)=0.1421049481421373368383892660234473q-01
       w(66)=0.1433956489572470031198777092644869q-01
       w(67)=0.1446481584573995619537692547529471q-01
       w(68)=0.1458621430551700804926100230937411q-01
       w(69)=0.1470372794236731174222873844805851q-01
       w(70)=0.1481732545826760429607097933074648q-01
       w(71)=0.1492697659819566833198447984901299q-01
       w(72)=0.1503265215818830882160301087474154q-01
       w(73)=0.1513432399311939599056635198463324q-01
       w(74)=0.1523196502419590278248468769118197q-01
       w(75)=0.1532554924616994040247847083367136q-01
       w(76)=0.1541505173426487110205508422732436q-01
       w(77)=0.1550044865081365352085059342061589q-01
       w(78)=0.1558171725160765254549632202027095q-01
       w(79)=0.1565883589195422276120398725391178q-01
       w(80)=0.1573178403244145213710425099997847q-01
       w(81)=0.1580054224440853058129963261059210q-01
       w(82)=0.1586509221512028640526162020266563q-01
       w(83)=0.1592541675264451252875776745840061q-01
       w(84)=0.1598149979043078341497460420826010q-01
       w(85)=0.1603332639158954322984235013839654q-01
       w(86)=0.1608088275287032554860887339963451q-01
       w(87)=0.1612415620833804506520547322754080q-01
       w(88)=0.1616313523274638217456553266231246q-01
       w(89)=0.1619780944460736197339177074554857q-01
       w(90)=0.1622816960895631013944081316803146q-01
       w(91)=0.1625420763981144928166232902304553q-01
       w(92)=0.1627591660232748068189234778679601q-01
       w(93)=0.1629329071464257785160185083292900q-01
       w(94)=0.1630632534941829998274048170492025q-01
       w(95)=0.1631501703507201515824817210238839q-01
       w(96)=0.1631936345670150508355623929347410q-01
c
       X1( 1)=0.9999219686591948443225892454883059q+00
       X1( 2)=0.9995888803394138541093349296043197q+00
       X1( 3)=0.9989897220082228180635279012706708q+00
       X1( 4)=0.9981245201227147822147993855489200q+00
       X1( 5)=0.9969934888089658276499713406258741q+00
       X1( 6)=0.9955969254054171881427829799848505q+00
       X1( 7)=0.9939352005680611362729976453293206q+00
       X1( 8)=0.9920087563423258732188356447321776q+00
       X1( 9)=0.9898181055591399025125907309048817q+00
       X1(10)=0.9873638315366433759878785288998682q+00
       X1(11)=0.9846465878619336764167001102060532q+00
       X1(12)=0.9816670981895111523937086705610450q+00
       X1(13)=0.9784261560354051061458091997963069q+00
       X1(14)=0.9749246245590810488762592665419201q+00
       X1(15)=0.9711634363298709597893756810891852q+00
       X1(16)=0.9671435930764901577246260169273177q+00
       X1(17)=0.9628661654189837492787994912503302q+00
       X1(18)=0.9583322925828043949981517329007647q+00
       X1(19)=0.9535431820948995063904708574256973q+00
       X1(20)=0.9485001094617782425446726957406750q+00
       X1(21)=0.9432044178295800403958793265082119q+00
       X1(22)=0.9376575176261969636700289673768231q+00
       X1(23)=0.9318608861855216719820172701963348q+00
       X1(24)=0.9258160673539062609128154511394387q+00
       X1(25)=0.9195246710789272245880472369639542q+00
       X1(26)=0.9129883729805597709004952516054175q+00
       X1(27)=0.9062089139048714766238190927436892q+00
       X1(28)=0.8991880994603512666000259568610842q+00
       X1(29)=0.8919277995369952102259633387959993q+00
       X1(30)=0.8844299478082758137705466992072068q+00
       X1(31)=0.8766965412161264494836028717904297q+00
       X1(32)=0.8687296394390773650062052096541978q+00
       X1(33)=0.8605313643436843991581167574853416q+00
       X1(34)=0.8521038994193961183039219550285386q+00
       X1(35)=0.8434494891970095976588225636632542q+00
       X1(36)=0.8345704386508695147393385041118966q+00
       X1(37)=0.8254691125849696045538228486220810q+00
       X1(38)=0.8161479350031198524218644920022155q+00
       X1(39)=0.8066093884633470731986553853368341q+00
       X1(40)=0.7968560134167007467203242051792675q+00
       X1(41)=0.7868904075306401492697565606090533q+00
       X1(42)=0.7767152249971829400510635555749104q+00
       X1(43)=0.7663331758259994299497462391583613q+00
       X1(44)=0.7557470251226407768848703807442227q+00
       X1(45)=0.7449595923520933172914643282748366q+00
       X1(46)=0.7339737505878551560441256578153761q+00
       X1(47)=0.7227924257467349966987265710195044q+00
       X1(48)=0.7114185958095769994814961263379360q+00
       X1(49)=0.6998552900281192051475475113416287q+00
       X1(50)=0.6881055881181967577927191283377693q+00
       X1(51)=0.6761726194395047980470490663027341q+00
       X1(52)=0.6640595621621394789116155447266718q+00
       X1(53)=0.6517696424201390789275344473475309q+00
       X1(54)=0.6393061334522506504939562614901836q+00
       X1(55)=0.6266723547301510440942257233984526q+00
       X1(56)=0.6138716710743544910676773343071987q+00
       X1(57)=0.6009074917580422075109537861843085q+00
       X1(58)=0.5877832695990526990513018873113860q+00
       X1(59)=0.5745025000402745997622328877199733q+00
       X1(60)=0.5610687202186869675605952774257933q+00
       X1(61)=0.5474855080232949822198280472393228q+00
       X1(62)=0.5337564811422119498589125315046377q+00
       X1(63)=0.5198852960991414086385045279429992q+00
       X1(64)=0.5058756472795159536497884562071851q+00
       X1(65)=0.4917312659465521538697583777912980q+00
       X1(66)=0.4774559192474836198483994125790001q+00
       X1(67)=0.4630534092102368967759482341181050q+00
       X1(68)=0.4485275717308174030575108927536054q+00
       X1(69)=0.4338822755516751088216675844448257q+00
       X1(70)=0.4191214212313220512518038667756913q+00
       X1(71)=0.4042489401054761136151635308892819q+00
       X1(72)=0.3892687932400077517555031425089794q+00
       X1(73)=0.3741849703759685350105230464587519q+00
       X1(74)=0.3590014888669824774342224742254436q+00
       X1(75)=0.3437223926092831692859114825625485q+00
       X1(76)=0.3283517509646816774505792608815985q+00
       X1(77)=0.3128936576767520662584319958241016q+00
       X1(78)=0.2973522297805231965735868566629289q+00
       X1(79)=0.2817316065059671905423468209843095q+00
       X1(80)=0.2660359481755766015698183178924267q+00
       X1(81)=0.2502694350963239034900637063721692q+00
       X1(82)=0.2344362664462984090903798119956503q+00
       X1(83)=0.2185406591563171457463379882822288q+00
       X1(84)=0.2025868467868075545434961993931732q+00
       X1(85)=0.1865790784002611385477991734258360q+00
       X1(86)=0.1705216174295583655043303941593981q+00
       X1(87)=0.1544187405424662298791849492833018q+00
       X1(88)=0.1382747365026108985189998820297989q+00
       X1(89)=0.1220939050272288030156827800270713q+00
       X1(90)=0.1058805556420003998801369668227609q+00
       X1(91)=0.8963900653327159661982075744412597q-01
       X1(92)=0.7337358339796853757389944640916407q-01
       X1(93)=0.5708861829151205770128136652529276q-01
       X1(94)=0.4078844847403864527738161700471162q-01
       X1(95)=0.2447741525523520549315998396275543q-01
       X1(96)=0.8159862838095286144725462139639069q-02
C
        GO TO 56
C
C-------------------------------------------------------------------
C³                       N --- 192                                  ³
C-------------------------------------------------------------------
C
192    continue  
C
         NG=192
C
       w(  1)=0.5019410348692173752939580444469347q-04
       w(  2)=0.1168390665730188627954282573512990q-03
       w(  3)=0.1835749193551260444766282687846253q-03
       w(  4)=0.2503070890844147244367666054448446q-03
       w(  5)=0.3170242698112706309204038817114463q-03
       w(  6)=0.3837208020912924380769776054432363q-03
       w(  7)=0.4503919137716877613812631428790360q-03
       w(  8)=0.5170330453491546380708515921736194q-03
       w(  9)=0.5836397042629966772866241664826160q-03
       w( 10)=0.6502074240969915513807189141732146q-03
       w( 11)=0.7167317509947397544734551873122715q-03
       w( 12)=0.7832082385905052313623308973839090q-03
       w( 13)=0.8496324460039079108554845086356611q-03
       w( 14)=0.9159999370632673684792719266271580q-03
       w( 15)=0.9823062800663751210383405437425685q-03
       w( 16)=0.1048547047793687527743074799806691q-02
       w( 17)=0.1114717817647312658872806756262422q-02
       w( 18)=0.1180814171855886469925513966209686q-02
       w( 19)=0.1246831697715436938175984171901580q-02
       w( 20)=0.1312765987850601065545162646052231q-02
       w( 21)=0.1378612640487679042415226157425863q-02
       w( 22)=0.1444367259734716215028497908118924q-02
       w( 23)=0.1510025455865806445301950005132107q-02
       w( 24)=0.1575582845607929439091787931193488q-02
       w( 25)=0.1641035052429278657631432315109335q-02
       w( 26)=0.1706377706828417190445077796020022q-02
       w( 27)=0.1771606446623829624142789694803191q-02
       w( 28)=0.1836716917243580798316589604725254q-02
       w( 29)=0.1901704772014882631752333906403915q-02
       w( 30)=0.1966565672453428345731083190287523q-02
       w( 31)=0.2031295288552391503809476913047267q-02
       w( 32)=0.2095889299071012646298758209214139q-02
       w( 33)=0.2160343391822713433314413546749456q-02
       w( 34)=0.2224653263962689952606445977779841q-02
       w( 35)=0.2288814622274945002047155445794683q-02
       w( 36)=0.2352823183458724890143200545087301q-02
       w( 37)=0.2416674674414330378673324717230462q-02
       w( 38)=0.2480364832528274330294802477920565q-02
       w( 39)=0.2543889405957760764920170338479650q-02
       w( 40)=0.2607244153914461606545741921501410q-02
       w( 41)=0.2670424846947568578829181362191945q-02
       w( 42)=0.2733427267226098595794760803005978q-02
       w( 43)=0.2796247208820431672903378569127224q-02
       w( 44)=0.2858880477983060909682385985388926q-02
       w( 45)=0.2921322893428534508519486756552241q-02
       w( 46)=0.2983570286612570124238412542485457q-02
       w( 47)=0.3045618502010322106922345096791569q-02
       w( 48)=0.3107463397393782421712854698238788q-02
       w( 49)=0.3169100844108296215501620180239663q-02
       w( 50)=0.3230526727348173160118143165099009q-02
       w( 51)=0.3291736946431375841217790633522596q-02
       w( 52)=0.3352727415073266586429107133956190q-02
       w( 53)=0.3413494061659394239096272143057411q-02
       w( 54)=0.3474032829517302487949757945706370q-02
       w( 55)=0.3534339677187341460401017745200512q-02
       w( 56)=0.3594410578692464379538911645922593q-02
       w( 57)=0.3654241523806991173587274806303053q-02
       w( 58)=0.3713828518324321012561164471010274q-02
       w( 59)=0.3773167584323575830912870459951390q-02
       w( 60)=0.3832254760435156977699152627836126q-02
       w( 61)=0.3891086102105197217710156177243493q-02
       w( 62)=0.3949657681858890388459062191086982q-02
       w( 63)=0.4007965589562681099242087057301574q-02
       w( 64)=0.4066005932685296939881577617258525q-02
       w( 65)=0.4123774836557605748452690149827222q-02
       w( 66)=0.4181268444631280569420091916030623q-02
       w( 67)=0.4238482918736255016298372596434625q-02
       w( 68)=0.4295414439336951836297167614953607q-02
       w( 69)=0.4352059205787267558498685300423258q-02
       w( 70)=0.4408413436584296192005077629759063q-02
       w( 71)=0.4464473369620775026236987584419831q-02
       w( 72)=0.4520235262436235672203545086973304q-02
       w( 73)=0.4575695392466843571130698967835241q-02
       w( 74)=0.4630850057293909285354959836294116q-02
       w( 75)=0.4685695574891054975883778349808404q-02
       w( 76)=0.4740228283870019561507770496630028q-02
       w( 77)=0.4794444543725086145835985874347312q-02
       w( 78)=0.4848340735076115391122488071208516q-02
       w( 79)=0.4901913259910168611267227157949856q-02
       w( 80)=0.4955158541821704450910967740974816q-02
       w( 81)=0.5008073026251333113105562322508227q-02
       w( 82)=0.5060653180723112194628321228366439q-02
       w( 83)=0.5112895495080368285622576841672028q-02
       w( 84)=0.5164796481720028588884672211481845q-02
       w( 85)=0.5216352675825446913778534306791000q-02
       w( 86)=0.5267560635597708500439975819280702q-02
       w( 87)=0.5318416942485398231630515459523692q-02
       w( 88)=0.5368918201412816892310868929135109q-02
       w( 89)=0.5419061041006630240722915992774529q-02
       w( 90)=0.5468842113820935759491052843955629q-02
       w( 91)=0.5518258096560732060974190869931274q-02
       w( 92)=0.5567305690303776027812743496351508q-02
       w( 93)=0.5615981620720812877314953198997290q-02
       w( 94)=0.5664282638294164447007804661322833q-02
       w( 95)=0.5712205518534661108333280508549526q-02
       w( 96)=0.5759747062196902826094378399809767q-02
       w( 97)=0.5806904095492834992840479525720119q-02
       w( 98)=0.5853673470303624779921533980872358q-02
       w( 99)=0.5900052064389823860428155662763938q-02
       w(100)=0.5946036781599803473663010369179463q-02
       w(101)=0.5991624552076447916150622583933663q-02
       w(102)=0.6036812332462092660480589860761524q-02
       w(103)=0.6081597106101693420485741134637252q-02
       w(104)=0.6125975883244212599374467415858922q-02
       w(105)=0.6169945701242209676457655012639579q-02
       w(106)=0.6213503624749622208027636789813014q-02
       w(107)=0.6256646745917724238751533827339605q-02
       w(108)=0.6299372184589249041626393658557443q-02
       w(109)=0.6341677088490663227100668621788001q-02
       w(110)=0.6383558633422579385387769513699468q-02
       w(111)=0.6425014023448294550274553378309604q-02
       w(112)=0.6466040491080441897852466968525296q-02
       w(113)=0.6506635297465743219563407826493119q-02
       w(114)=0.6546795732567849835747855268282835q-02
       w(115)=0.6586519115348259743501071436292001q-02
       w(116)=0.6625802793945298921075722647056811q-02
       w(117)=0.6664644145851154840307606718594922q-02
       w(118)=0.6703040578086950368576715976102163q-02
       w(119)=0.6740989527375846372639982636434371q-02
       w(120)=0.6778488460314161468276050167038758q-02
       w(121)=0.6815534873540497492057541529323633q-02
       w(122)=0.6852126293902859404703748170145686q-02
       w(123)=0.6888260278623758469357583153393716q-02
       w(124)=0.6923934415463287682766115548614226q-02
       w(125)=0.6959146322880158572715066354226053q-02
       w(126)=0.6993893650190688611165282578892367q-02
       w(127)=0.7028174077725728629354348759431617q-02
       w(128)=0.7061985316985519758650027489941020q-02
       w(129)=0.7095325110792469559164976829021600q-02
       w(130)=0.7128191233441837137054958908735018q-02
       w(131)=0.7160581490850317191016269447963523q-02
       w(132)=0.7192493720702513068763074281317703q-02
       w(133)=0.7223925792595289055192382656499322q-02
       w(134)=0.7254875608179992255524118965100018q-02
       w(135)=0.7285341101302534578926731351462142q-02
       w(136)=0.7315320238141325470995510168094003q-02
       w(137)=0.7344811017343046186931751382647443q-02
       w(138)=0.7373811470156256541366517357587224q-02
       w(139)=0.7402319660562825215473405886890571q-02
       w(140)=0.7430333685407174847310783050492116q-02
       w(141)=0.7457851674523333277215671196279933q-02
       w(142)=0.7484871790859782466529175526518238q-02
       w(143)=0.7511392230602096754957207871825633q-02
       w(144)=0.7537411223293362269450512667035560q-02
       w(145)=0.7562927031952369445614768837178497q-02
       w(146)=0.7587937953189570771324946434885427q-02
       w(147)=0.7612442317320796011408216543945447q-02
       w(148)=0.7636438488478717321966589889492444q-02
       w(149)=0.7659924864722056813124101781926437q-02
       w(150)=0.7682899878142529269693742460534107q-02
       w(151)=0.7705361994969512890456393254497869q-02
       w(152)=0.7727309715672441058417678260561926q-02
       w(153)=0.7748741575060908306548754509131573q-02
       w(154)=0.7769656142382483796113485647677424q-02
       w(155)=0.7790052021418225777726989227644196q-02
       w(156)=0.7809927850575890658769000066558214q-02
       w(157)=0.7829282302980830454679607002281370q-02
       w(158)=0.7848114086564572555984424289025206q-02
       w(159)=0.7866421944151075897620850722052546q-02
       w(160)=0.7884204653540657772256421021438552q-02
       w(161)=0.7901461027591585684794010336408512q-02
       w(162)=0.7918189914299328801136433509273305q-02
       w(163)=0.7934390196873463700524380495810257q-02
       w(164)=0.7950060793812229297356218285761781q-02
       w(165)=0.7965200658974725955335514386937625q-02
       w(166)=0.7979808781650753974061721090304012q-02
       w(167)=0.7993884186628286785770804762619995q-02
       w(168)=0.8007425934258574357835190243047412q-02
       w(169)=0.8020433120518872454835676165684177q-02
       w(170)=0.8032904877072793572511401702041506q-02
       w(171)=0.8044840371328275514666928430531446q-02
       w(172)=0.8056238806493163743157443720547272q-02
       w(173)=0.8067099421628403790373377115032378q-02
       w(174)=0.8077421491698840183193714388325539q-02
       w(175)=0.8087204327621618487162344370528777q-02
       w(176)=0.8096447276312187239653214547491269q-02
       w(177)=0.8105149720727896701017221018990322q-02
       w(178)=0.8113311079909191513135920363982719q-02
       w(179)=0.8120930809018394515433615358296035q-02
       w(180)=0.8128008399376079129209410370631953q-02
       w(181)=0.8134543378495027882133720414389284q-02
       w(182)=0.8140535310111774805898703493707823q-02
       w(183)=0.8145983794215729601308411459789776q-02
       w(184)=0.8150888467075881626531352383873214q-02
       w(185)=0.8155249001265081925804850360173328q-02
       w(186)=0.8159065105681901677566290907904316q-02
       w(187)=0.8162336525570065602780258027416842q-02
       w(188)=0.8165063042535459036121901533080955q-02
       w(189)=0.8167244474560707524654813078338456q-02
       w(190)=0.8168880676017327980695423063927239q-02
       w(191)=0.8169971537675450577674640072765497q-02
       w(192)=0.8170516986711110739980316957776570q-02
c
       X1(  1)=0.9999804411726473547435916080674043q+00
       X1(  2)=0.9998969471378595991468054406187173q+00
       X1(  3)=0.9997467408113522753988972414026949q+00
       X1(  4)=0.9995297988558859110049282600075867q+00
       X1(  5)=0.9992461316671845033851407841964070q+00
       X1(  6)=0.9988957572063257575662544639878614q+00
       X1(  7)=0.9984786985384589563177774591129954q+00
       X1(  8)=0.9979949833727937972095028935222694q+00
       X1(  9)=0.9974446439389107953978233829595598q+00
       X1( 10)=0.9968277169440913106011587182604726q+00
       X1( 11)=0.9961442435551086373704026483342547q+00
       X1( 12)=0.9953942693885953365534400381277338q+00
       X1( 13)=0.9945778445047067583913757666163057q+00
       X1( 14)=0.9936950234020882231198241628454611q+00
       X1( 15)=0.9927458650133152938198183044691738q+00
       X1( 16)=0.9917304327004320079833462855178316q+00
       X1( 17)=0.9906487942504060643059327389568952q+00
       X1( 18)=0.9895010218704086785345938500939276q+00
       X1( 19)=0.9882871921828698242876652773997211q+00
       X1( 20)=0.9870073862202814985581032513524182q+00
       X1( 21)=0.9856616894197333315606559277141099q+00
       X1( 22)=0.9842501916171713291634871599288893q+00
       X1( 23)=0.9827729870413742449448393161947657q+00
       X1( 24)=0.9812301743076442733429848278256271q+00
       X1( 25)=0.9796218564112100920471705233487289q+00
       X1( 26)=0.9779481407203411185936260705153875q+00
       X1( 27)=0.9762091389691723840573414857694040q+00
       X1( 28)=0.9744049672502397817145518444257367q+00
       X1( 29)=0.9725357460067256891482271342125922q+00
       X1( 30)=0.9706016000244151299324478785733199q+00
       X1( 31)=0.9686026584233627613468938469142880q+00
       X1( 32)=0.9665390546492710636330171850153891q+00
       X1( 33)=0.9644109264645801744096617485385597q+00
       X1( 34)=0.9622184159392698657186713295280470q+00
       X1( 35)=0.9599616694413742051951783143330813q+00
       X1( 36)=0.9576408376272094800157457907014543q+00
       X1( 37)=0.9552560754313159945804566987849321q+00
       X1( 38)=0.9528075420561143817127760802366948q+00
       X1( 39)=0.9502954009612770934723309692376080q+00
       X1( 40)=0.9477198198528157621437301911075226q+00
       X1( 41)=0.9450809706718851450664253714904562q+00
       X1( 42)=0.9423790295833043890479544200210348q+00
       X1( 43)=0.9396141769637963714023014550010982q+00
       X1( 44)=0.9367865973899458953560325449694781q+00
       X1( 45)=0.9338964796258775377989448863308919q+00
       X1( 46)=0.9309440166106539672203921286354512q+00
       X1( 47)=0.9279294054453955692394478673227273q+00
       X1( 48)=0.9248528473801222364606223411731031q+00
       X1( 49)=0.9217145478003181985074972751087354q+00
       X1( 50)=0.9185147162132207870349702943350292q+00
       X1( 51)=0.9152535662338340493199450128510118q+00
       X1( 52)=0.9119313155706681426978188772843445q+00
       X1( 53)=0.9085481860112054606613366380471192q+00
       X1( 54)=0.9051044034070944598794091866618400q+00
       X1( 55)=0.9016001976590721757340055783620212q+00
       X1( 56)=0.8980358027016164322189825155509176q+00
       X1( 57)=0.8944114564873287702000339620129925q+00
       X1( 58)=0.8907274009710491361030051461123556q+00
       X1( 59)=0.8869838820937033910809182316524975q+00
       X1( 60)=0.8831811497658847186098080398417027q+00
       X1( 61)=0.8793194578511700262809281937597013q+00
       X1( 62)=0.8753990641491724552926944264179587q+00
       X1( 63)=0.8714202303783311288001759785935056q+00
       X1( 64)=0.8673832221584392878530512051251624q+00
       X1( 65)=0.8632883089929119811445175287540802q+00
       X1( 66)=0.8591357642507944922033244774248383q+00
       X1( 67)=0.8549258651485127049883796988593641q+00
       X1( 68)=0.8506588927313666260896490767762044q+00
       X1( 69)=0.8463351318547682988996322075658534q+00
       X1( 70)=0.8419548711652253621957711528218481q+00
       X1( 71)=0.8375184030810715225649140193318276q+00
       X1( 72)=0.8330260237729452270055301759391325q+00
       X1( 73)=0.8284780331440178388608488149881887q+00
       X1( 74)=0.8238747348099726369655259746972595q+00
       X1( 75)=0.8192164360787359745288728764424483q+00
       X1( 76)=0.8145034479299619508281182573903209q+00
       X1( 77)=0.8097360849942719652446330962187934q+00
       X1( 78)=0.8049146655322505395435109626909809q+00
       X1( 79)=0.8000395114131988105713566146269877q+00
       X1( 80)=0.7951109480936471117275694026534827q+00
       X1( 81)=0.7901293045956280776497929820694757q+00
       X1( 82)=0.7850949134847117225435133205577459q+00
       X1( 83)=0.7800081108478039584779970252953213q+00
       X1( 84)=0.7748692362707100357648461775913433q+00
       X1( 85)=0.7696786328154644032303802974686968q+00
       X1( 86)=0.7644366469974285017878192808063809q+00
       X1( 87)=0.7591436287621580202088147500339729q+00
       X1( 88)=0.7537999314620411573852465037346998q+00
       X1( 89)=0.7484059118327094506603551001816149q+00
       X1( 90)=0.7429619299692227449922151617960297q+00
       X1( 91)=0.7374683493020298927912658642704690q+00
       X1( 92)=0.7319255365727067892461097624838560q+00
       X1( 93)=0.7263338618094733628170787542357151q+00
       X1( 94)=0.7206936983024911553341626661454074q+00
       X1( 95)=0.7150054225789431407838238964062842q+00
       X1( 96)=0.7092694143778974464070093666021278q+00
       X1( 97)=0.7034860566249566541573538749467285q+00
       X1( 98)=0.6976557354066943748831886867027143q+00
       X1( 99)=0.6917788399448808017985746044446108q+00
       X1(100)=0.6858557625704989638962256126180645q+00
       X1(101)=0.6798868986975534139279404201613783q+00
       X1(102)=0.6738726467966730994350850362576510q+00
       X1(103)=0.6678134083685101790518474960354795q+00
       X1(104)=0.6617095879169365599265010826617543q+00
       X1(105)=0.6555615929220399456098575348715410q+00
       X1(106)=0.6493698338129211971445671305719373q+00
       X1(107)=0.6431347239402948233530363131371438q+00
       X1(108)=0.6368566795488944294646016464415135q+00
       X1(109)=0.6305361197496849662433452508668629q+00
       X1(110)=0.6241734664918836346756934003516616q+00
       X1(111)=0.6177691445347913140508466508936422q+00
       X1(112)=0.6113235814194363939162948824599147q+00
       X1(113)=0.6048372074400329029143303550023993q+00
       X1(114)=0.5983104556152548399027509859457489q+00
       X1(115)=0.5917437616593286250330175920622718q+00
       X1(116)=0.5851375639529456006011742680295714q+00
       X1(117)=0.5784923035139965235000503483323530q+00
       X1(118)=0.5718084239681300029848340165972527q+00
       X1(119)=0.5650863715191368492172486913231557q+00
       X1(120)=0.5583265949191623096754895981040859q+00
       X1(121)=0.5515295454387481820070139260283435q+00
       X1(122)=0.5446956768367068032584569273814345q+00
       X1(123)=0.5378254453298289266406103544919316q+00
       X1(124)=0.5309193095624275080757997207496409q+00
       X1(125)=0.5239777305757194357293929483507150q+00
       X1(126)=0.5170011717770472465458339376197121q+00
       X1(127)=0.5099900989089428844917984178933022q+00
       X1(128)=0.5029449800180355657541031658341907q+00
       X1(129)=0.4958662854238058265471594928099309q+00
       X1(130)=0.4887544876871878394533531875564727q+00
       X1(131)=0.4816100615790220943490704702881039q+00
       X1(132)=0.4744334840483605499584968708817867q+00
       X1(133)=0.4672252341906263719261265188289889q+00
       X1(134)=0.4599857932156303830064757348266593q+00
       X1(135)=0.4527156444154463605351490729288543q+00
       X1(136)=0.4454152731321473257685195731821442q+00
       X1(137)=0.4380851667254049789592289513884266q+00
       X1(138)=0.4307258145399544431708683301861818q+00
       X1(139)=0.4233377078729264888269560524147196q+00
       X1(140)=0.4159213399410494198360858966659831q+00
       X1(141)=0.4084772058477228108362860940704126q+00
       X1(142)=0.4010058025499652936566260976977899q+00
       X1(143)=0.3935076288252385995023630085956919q+00
       X1(144)=0.3859831852381500716308716351558133q+00
       X1(145)=0.3784329741070358713986999053061071q+00
       X1(146)=0.3708574994704271085247931815886368q+00
       X1(147)=0.3632572670534011342307052690318163q+00
       X1(148)=0.3556327842338202435849388160444399q+00
       X1(149)=0.3479845600084600408949216252276269q+00
       X1(150)=0.3403131049590297293560263589567875q+00
       X1(151)=0.3326189312180865933819875314974761q+00
       X1(152)=0.3249025524348469491045798886735075q+00
       X1(153)=0.3171644837408958454420247591082737q+00
       X1(154)=0.3094052417157978048948243288461299q+00
       X1(155)=0.3016253443526108998341368253182129q+00
       X1(156)=0.2938253110233064665009572795200358q+00
       X1(157)=0.2860056624440967652338280796185322q+00
       X1(158)=0.2781669206406729015881504068838585q+00
       X1(159)=0.2703096089133553290009916150615600q+00
       X1(160)=0.2624342518021592594911847510776239q+00
       X1(161)=0.2545413750517773145651051338829009q+00
       X1(162)=0.2466315055764817540234059925128776q+00
       X1(163)=0.2387051714249486257328317902941839q+00
       X1(164)=0.2307629017450061846396456234879560q+00
       X1(165)=0.2228052267483099343568580151017328q+00
       X1(166)=0.2148326776749466495559910286414503q+00
       X1(167)=0.2068457867579697421352268902362858q+00
       X1(168)=0.1988450871878683387191577282500035q+00
       X1(169)=0.1908311130769724414706666451301036q+00
       X1(170)=0.1828043994237965484624347057734921q+00
       X1(171)=0.1747654820773241139638987097839602q+00
       X1(172)=0.1667148977012352329489065511304805q+00
       X1(173)=0.1586531837380799379195669054759770q+00
       X1(174)=0.1505808783733994997726148961683054q+00
       X1(175)=0.1424985204997981279057728905596349q+00
       X1(176)=0.1344066496809674680728439437942368q+00
       X1(177)=0.1263058061156662996474135866673014q+00
       X1(178)=0.1181965306016578369458432859131389q+00
       X1(179)=0.1100793644996070420905163308066413q+00
       X1(180)=0.1019548496969403595638551772479879q+00
       X1(181)=0.9382352857167028511229019667736039q-01
       X1(182)=0.8568594395618718400695586041120788q-01
       X1(183)=0.7754263910102077585426492071875290q-01
       X1(184)=0.6939415763857370517451817128958850q-01
       X1(185)=0.6124104354682961883021169726308101q-01
       X1(186)=0.5308384111303817308758075577250051q-01
       X1(187)=0.4492309489737939463505594887932863q-01
       X1(188)=0.3675934969660982126060034610772658q-01
       X1(189)=0.2859315050769284910625367206999765q-01
       X1(190)=0.2042504249141571447255039484499058q-01
       X1(191)=0.1225557093599553903773214761107266q-01
       X1(192)=0.4085281220676868086781824858344715q-02
c


Cã==================================================================¬

56       CONTINUE
C
        RETURN
        END
Cã==================================================================¬

c***************************************************************
c           subroutine Olmon(wavel,Rm,ref,eps,KD,KS) Dec 2019
c	 2019 New Bulk Constants
c 	 Olmon et al. PHYSICAL REVIEW B 86, 235147 (2012)
c***************************************************************
C     Calculations of ref=real(n)+i*imag(n) and eps=real(eps)+i*imag(eps)
c     for  GOLD according to Olmon PRB 2012 spline for Bulk Single Cristalyne sample (SC data)

c     The size correction is considered 
c     wavel - is the input wavelength nanometers,
c     Rm=aLeff is the size-corrected length
c     HERE: aLeff=Rev/As, As is the elecctron scattring parameter
c	a small As=1.e-09 corresponds to the BULK optical constaants
c	Typically As -0.-1 
c     ref and eps are ouput complex ref=sqrt(eps)
c     Correction to the size limiting effect:
c	gamma_b=gam_b*10**(14) 1/s, gam_b=1.05-1.6  
c      1.3 - Scaffardi JApplPhys 2014,(1.32) Etchegoin-Le Ru(2006)gam_b=1.3
c	 gam_b=1.3
c	Plasma wavelength
c	(h/2*pi)*omega_p=9 eV  lambda_p=140 nm
c	parameters gamma_b/omega and gamma_a/omega = [gamma_b+A*v_F/Leff]/omega
c	a=1 v_F=14.1*10**(14) nm/s Fermi velocity, Leff=Rm effetive length
c	gamma_b/omega=ga_om=0.53*10**(-4)*gam_b*lambda
c	v_F/(Leff*omega)=7.48*10**(-4)*lambda/Leff
c	KD is the parameter for size correction
c	if KD=0 Full correction is used

c	 Deps=(wavel/wp)**2*[ 1.q0/(1.q0+aim*gamb_om)-1.q0/(1.q0+aim*gama_om)]
c
c	if KD=1 approximation a9-11 is used
c      Deps = Deps2*[B+ i]
c      Deps2=(5.q0/Rm)*(wavel/520.q0)**3.q0
c	 B=(gama+gamb)/2/omega
c
c  	 if KD=2 B=0 is used
c	KS is the parameter for sample choice
c	if KS=0 TS sample constants are used
c	if KS=1 CS crystalline sample constants are used
c	for details see Olmon et al PRB 2012	 
C--------------------------------------------------------
       subroutine Olmon(wavel,Rm,ref,eps,KD,KS)
      IMPLICIT real*16(A-H,O-Z)
      real*16 MN,MK
      DIMENSION XN(1000),YN(1000),MN(1000),XK(1000),YK(1000),MK(1000)

	complex*32 ref, eps,epsb,Deps,aim
C--------------------------------------------------------
	aim=(0.q0,1.q0)
	N=258
c     
      xx=wavel

	xn(	1)=	299.9859497
	xn(	2)=	310.037992
	xn(	3)=	320.0417992
	xn(	4)=	330.0084988
	xn(	5)=	339.9621415
	xn(	6)=	350.0400706
	xn(	7)=	360.0005604
	xn(	8)=	369.9916234
	xn(	9)=	379.9699448
	xn(	10)=	390.010044
	xn(	11)=	399.9490097
	xn(	12)=	410.0006382
	xn(	13)=	420.0006538
	xn(	14)=	430.0526986
	xn(	15)=	439.972296
	xn(	16)=	450.0333684
	xn(	17)=	460.0526642
	xn(	18)=	469.9931501
	xn(	19)=	480.0007472
	xn(	20)=	490.0560988
	xn(	21)=	499.9362621
	xn(	22)=	510.0131345
	xn(	23)=	520.0679237
	xn(	24)=	530.0735058
	xn(	25)=	540.0008406
	xn(	26)=	550.0629681
	xn(	27)=	560.0008717
	xn(	28)=	570.0422667
	xn(	29)=	579.9073573
	xn(	30)=	590.1199096
	xn(	31)=	600.1171007
	xn(	32)=	609.858303
	xn(	33)=	619.920965
	xn(	34)=	630.0009807
	xn(	35)=	640.0835983
	xn(	36)=	650.1530834
	xn(	37)=	659.8413677
	xn(	38)=	669.8227607
	xn(	39)=	680.110768
	xn(	40)=	689.9509905
	xn(	41)=	700.0801412
	xn(	42)=	710.1041982
	xn(	43)=	720.0011208
	xn(	44)=	730.1778151
	xn(	45)=	740.2041373
	xn(	46)=	750.055614
	xn(	47)=	760.1728571
	xn(	48)=	770.0881553
	xn(	49)=	779.7747987
	xn(	50)=	790.2115551
	xn(	51)=	799.8980194
	xn(	52)=	809.8249053
	xn(	53)=	820.0012765
	xn(	54)=	829.8808099
	xn(	55)=	840.0013076
	xn(	56)=	849.7888485
	xn(	57)=	859.8071637
	xn(	58)=	870.0645123
	xn(	59)=	879.9445919
	xn(	60)=	890.0516368
	xn(	61)=	899.7401524
	xn(	62)=	910.3097871
	xn(	63)=	919.764043
	xn(	64)=	930.113976
	xn(	65)=	939.9863002
	xn(	66)=	950.0704444
	xn(	67)=	960.3732998
	xn(	68)=	970.1423552
	xn(	69)=	980.1121976
	xn(	70)=	990.2890815
	xn(	71)=	999.8725242
	xn(	72)=	1009.643265
	xn(	73)=	1019.60685
	xn(	74)=	1029.769045
	xn(	75)=	1040.135847
	xn(	76)=	1049.823819
	xn(	77)=	1059.693957
	xn(	78)=	1069.75145
	xn(	79)=	1080.001681
	xn(	80)=	1090.450246
	xn(	81)=	1100.125936
	xn(	82)=	1109.97487
	xn(	83)=	1120.001743
	xn(	84)=	1130.211422
	xn(	85)=	1139.560597
	xn(	86)=	1150.13166
	xn(	87)=	1159.814715
	xn(	88)=	1169.662198
	xn(	89)=	1179.678335
	xn(	90)=	1189.867495
	xn(	91)=	1200.234201
	xn(	92)=	1209.601883
	xn(	93)=	1220.31686
	xn(	94)=	1230.001915
	xn(	95)=	1239.965927
	xn(	96)=	1249.96666
	xn(	97)=	1260.001961
	xn(	98)=	1269.939496
	xn(	99)=	1280.03503
	xn(	100)=	1290.023858
	xn(	101)=	1300.03348
	xn(	102)=	1310.061211
	xn(	103)=	1319.963728
	xn(	104)=	1330.017089
	xn(	105)=	1339.935081
	xn(	106)=	1350.002101
	xn(	107)=	1430.036828
	xn(	108)=	1440.002242
	xn(	109)=	1449.937937
	xn(	110)=	1460.011693
	xn(	111)=	1470.052087
	xn(	112)=	1480.054829
	xn(	113)=	1490.015539
	xn(	114)=	1499.929748
	xn(	115)=	1509.976775
	xn(	116)=	1519.972943
	xn(	117)=	1529.913537
	xn(	118)=	1539.985008
	xn(	119)=	1549.996162
	xn(	120)=	1559.942036
	xn(	121)=	1570.016373
	xn(	122)=	1580.020301
	xn(	123)=	1589.948615
	xn(	124)=	1600.002491
	xn(	125)=	1609.975237
	xn(	126)=	1620.073082
	xn(	127)=	1630.084052
	xn(	128)=	1640.002553
	xn(	129)=	1650.042494
	xn(	130)=	1659.98384
	xn(	131)=	1670.045703
	xn(	132)=	1680.002615
	xn(	133)=	1690.078967
	xn(	134)=	1700.043782
	xn(	135)=	1709.890953
	xn(	136)=	1720.091468
	xn(	137)=	1729.931533
	xn(	138)=	1739.88483
	xn(	139)=	1749.953324
	xn(	140)=	1759.889184
	xn(	141)=	1769.938515
	xn(	142)=	1780.103274
	xn(	143)=	1790.126956
	xn(	144)=	1800.002802
	xn(	145)=	1809.988219
	xn(	146)=	1820.085041
	xn(	147)=	1830.024989
	xn(	148)=	1840.074102
	xn(	149)=	1849.958117
	xn(	150)=	1859.94889
	xn(	151)=	1870.04816
	xn(	152)=	1879.9726
	xn(	153)=	1890.002942
	xn(	154)=	1900.140889
	xn(	155)=	1910.093868
	xn(	156)=	1920.151665
	xn(	157)=	1930.015458
	xn(	158)=	1939.981114
	xn(	159)=	1950.05022
	xn(	160)=	1959.914527
	xn(	161)=	1969.879139
	xn(	162)=	1979.945592
	xn(	163)=	1990.115457
	xn(	164)=	2000.06764
	xn(	165)=	2006.541398
	xn(	166)=	2018.95771
	xn(	167)=	2025.554534
	xn(	168)=	2038.208006
	xn(	169)=	2044.594212
	xn(	170)=	2057.487438
	xn(	171)=	2063.995222
	xn(	172)=	2077.135081
	xn(	173)=	2084.118222
	xn(	174)=	2097.516376
	xn(	175)=	2104.280261
	xn(	176)=	2117.939751
	xn(	177)=	2124.836213
	xn(	178)=	2139.133765
	xn(	179)=	2146.169171
	xn(	180)=	2153.251007
	xn(	181)=	2167.555822
	xn(	182)=	2174.77974
	xn(	183)=	2189.759679
	xn(	184)=	2197.132607
	xn(	185)=	2204.555352
	xn(	186)=	2219.552327
	xn(	187)=	2227.127591
	xn(	188)=	2234.75474
	xn(	189)=	2242.43431
	xn(	190)=	2258.364171
	xn(	191)=	2266.207147
	xn(	192)=	2274.104787
	xn(	193)=	2282.057666
	xn(	194)=	2298.131474
	xn(	195)=	2306.68266
	xn(	196)=	2314.865441
	xn(	197)=	2323.106483
	xn(	198)=	2339.765861
	xn(	199)=	2348.185473
	xn(	200)=	2356.6659
	xn(	201)=	2365.207802
	xn(	202)=	2373.811851
	xn(	203)=	2382.936633
	xn(	204)=	2391.67039
	xn(	205)=	2409.331384
	xn(	206)=	2418.260055
	xn(	207)=	2427.255149
	xn(	208)=	2436.31741
	xn(	209)=	2445.930026
	xn(	210)=	2455.132535
	xn(	211)=	2464.404552
	xn(	212)=	2473.746868
	xn(	213)=	2483.160284
	xn(	214)=	2492.645617
	xn(	215)=	2502.203693
	xn(	216)=	2512.344336
	xn(	217)=	2522.054373
	xn(	218)=	2531.839759
	xn(	219)=	2541.701374
	xn(	220)=	2551.640111
	xn(	221)=	2561.65688
	xn(	222)=	2571.752603
	xn(	223)=	2581.928217
	xn(	224)=	2592.726746
	xn(	225)=	2603.069347
	xn(	226)=	2613.494793
	xn(	227)=	2624.004085
	xn(	228)=	2634.598236
	xn(	229)=	2645.27828
	xn(	230)=	2656.045266
	xn(	231)=	2666.900258
	xn(	232)=	2678.422834
	xn(	233)=	2689.461887
	xn(	234)=	2700.592311
	xn(	235)=	2711.815245
	xn(	236)=	2723.131847
	xn(	237)=	2734.543295
	xn(	238)=	2746.050786
	xn(	239)=	2757.655538
	xn(	240)=	2769.977502
	xn(	241)=	2781.785798
	xn(	242)=	2793.695201
	xn(	243)=	2805.707015
	xn(	244)=	2817.822568
	xn(	245)=	2830.043209
	xn(	246)=	2842.370312
	xn(	247)=	2854.805273
	xn(	248)=	2868.012792
	xn(	249)=	2880.673629
	xn(	250)=	2893.446744
	xn(	251)=	2906.333638
	xn(	252)=	2919.335837
	xn(	253)=	2932.454896
	xn(	254)=	2945.692397
	xn(	255)=	2959.756338
	xn(	256)=	2973.242038
	xn(	257)=	2986.851192
	xn(	258)=	3000.585503
c
	do i=1,N
		xk(i)=xn(i)
	end do
c
	if(KS.eq.0) then

c     TS sample Spline  parameters Real n

	yn(	1)=	1.667
	yn(	2)=	1.719
	yn(	3)=	1.754
	yn(	4)=	1.767
	yn(	5)=	1.759
	yn(	6)=	1.732
	yn(	7)=	1.696
	yn(	8)=	1.664
	yn(	9)=	1.645
	yn(	10)=	1.642
	yn(	11)=	1.645
	yn(	12)=	1.641
	yn(	13)=	1.621
	yn(	14)=	1.59
	yn(	15)=	1.556
	yn(	16)=	1.52
	yn(	17)=	1.466
	yn(	18)=	1.372
	yn(	19)=	1.227
	yn(	20)=	1.041
	yn(	21)=	0.8463
	yn(	22)=	0.678
	yn(	23)=	0.5507
	yn(	24)=	0.46
	yn(	25)=	0.3953
	yn(	26)=	0.3476
	yn(	27)=	0.3111
	yn(	28)=	0.2822
	yn(	29)=	0.2587
	yn(	30)=	0.2392
	yn(	31)=	0.2228
	yn(	32)=	0.209
	yn(	33)=	0.1973
	yn(	34)=	0.1874
	yn(	35)=	0.1789
	yn(	36)=	0.1717
	yn(	37)=	0.1656
	yn(	38)=	0.1604
	yn(	39)=	0.1561
	yn(	40)=	0.1524
	yn(	41)=	0.1494
	yn(	42)=	0.147
	yn(	43)=	0.145
	yn(	44)=	0.1435
	yn(	45)=	0.1424
	yn(	46)=	0.1416
	yn(	47)=	0.1411
	yn(	48)=	0.141
	yn(	49)=	0.141
	yn(	50)=	0.1414
	yn(	51)=	0.1419
	yn(	52)=	0.1426
	yn(	53)=	0.1435
	yn(	54)=	0.1445
	yn(	55)=	0.1457
	yn(	56)=	0.1471
	yn(	57)=	0.1486
	yn(	58)=	0.1501
	yn(	59)=	0.1518
	yn(	60)=	0.1536
	yn(	61)=	0.1555
	yn(	62)=	0.1575
	yn(	63)=	0.1596
	yn(	64)=	0.1617
	yn(	65)=	0.1639
	yn(	66)=	0.1662
	yn(	67)=	0.1686
	yn(	68)=	0.171
	yn(	69)=	0.1734
	yn(	70)=	0.176
	yn(	71)=	0.1786
	yn(	72)=	0.1812
	yn(	73)=	0.1839
	yn(	74)=	0.1866
	yn(	75)=	0.1894
	yn(	76)=	0.1923
	yn(	77)=	0.1951
	yn(	78)=	0.198
	yn(	79)=	0.201
	yn(	80)=	0.204
	yn(	81)=	0.207
	yn(	82)=	0.2101
	yn(	83)=	0.2132
	yn(	84)=	0.2163
	yn(	85)=	0.2195
	yn(	86)=	0.2227
	yn(	87)=	0.226
	yn(	88)=	0.2292
	yn(	89)=	0.2325
	yn(	90)=	0.2359
	yn(	91)=	0.2392
	yn(	92)=	0.2426
	yn(	93)=	0.2461
	yn(	94)=	0.2495
	yn(	95)=	0.253
	yn(	96)=	0.2565
	yn(	97)=	0.26
	yn(	98)=	0.2636
	yn(	99)=	0.2672
	yn(	100)=	0.2708
	yn(	101)=	0.2744
	yn(	102)=	0.2781
	yn(	103)=	0.2813
	yn(	104)=	0.2845
	yn(	105)=	0.2878
	yn(	106)=	0.291
	yn(	107)=	0.3214
	yn(	108)=	0.3249
	yn(	109)=	0.3284
	yn(	110)=	0.3319
	yn(	111)=	0.3355
	yn(	112)=	0.339
	yn(	113)=	0.3426
	yn(	114)=	0.3462
	yn(	115)=	0.3498
	yn(	116)=	0.3535
	yn(	117)=	0.3571
	yn(	118)=	0.3608
	yn(	119)=	0.3646
	yn(	120)=	0.3683
	yn(	121)=	0.3721
	yn(	122)=	0.3759
	yn(	123)=	0.3797
	yn(	124)=	0.3835
	yn(	125)=	0.3874
	yn(	126)=	0.3913
	yn(	127)=	0.3952
	yn(	128)=	0.3991
	yn(	129)=	0.4031
	yn(	130)=	0.4071
	yn(	131)=	0.4112
	yn(	132)=	0.4152
	yn(	133)=	0.4193
	yn(	134)=	0.4234
	yn(	135)=	0.4276
	yn(	136)=	0.4317
	yn(	137)=	0.4359
	yn(	138)=	0.4402
	yn(	139)=	0.4444
	yn(	140)=	0.4487
	yn(	141)=	0.4531
	yn(	142)=	0.4574
	yn(	143)=	0.4618
	yn(	144)=	0.4663
	yn(	145)=	0.4708
	yn(	146)=	0.4752
	yn(	147)=	0.4798
	yn(	148)=	0.4844
	yn(	149)=	0.489
	yn(	150)=	0.4936
	yn(	151)=	0.4983
	yn(	152)=	0.503
	yn(	153)=	0.5078
	yn(	154)=	0.5126
	yn(	155)=	0.5175
	yn(	156)=	0.5223
	yn(	157)=	0.5273
	yn(	158)=	0.5322
	yn(	159)=	0.5373
	yn(	160)=	0.5424
	yn(	161)=	0.5474
	yn(	162)=	0.5526
	yn(	163)=	0.5578
	yn(	164)=	0.5622
	yn(	165)=	0.5667
	yn(	166)=	0.5736
	yn(	167)=	0.5771
	yn(	168)=	0.5842
	yn(	169)=	0.5878
	yn(	170)=	0.5952
	yn(	171)=	0.5989
	yn(	172)=	0.6064
	yn(	173)=	0.6102
	yn(	174)=	0.618
	yn(	175)=	0.6219
	yn(	176)=	0.6299
	yn(	177)=	0.6339
	yn(	178)=	0.6421
	yn(	179)=	0.6463
	  yn(	180)=	0.6505
	  yn(	181)=	0.6591
	  yn(	182)=	0.6634
	  yn(	183)=	0.6722
	  yn(	184)=	0.6767
	  yn(	185)=	0.6812
	  yn(	186)=	0.6903
	  yn(	187)=	0.695
	  yn(	188)=	0.6997
	  yn(	189)=	0.7044
	  yn(	190)=	0.7141
	  yn(	191)=	0.719
	  yn(	192)=	0.7239
	  yn(	193)=	0.7289
	  yn(	194)=	0.7391
	  yn(	195)=	0.7442
	  yn(	196)=	0.7494
	  yn(	197)=	0.7547
	  yn(	198)=	0.7654
	  yn(	199)=	0.7708
	  yn(	200)=	0.7763
	  yn(	201)=	0.7819
	  yn(	202)=	0.7875
	  yn(	203)=	0.7932
	  yn(	204)=	0.7989
	  yn(	205)=	0.8106
	  yn(	206)=	0.8165
	  yn(	207)=	0.8225
	  yn(	208)=	0.8286
	  yn(	209)=	0.8347
	  yn(	210)=	0.8409
	  yn(	211)=	0.8472
	  yn(	212)=	0.8535
	  yn(	213)=	0.8599
	  yn(	214)=	0.8664
	  yn(	215)=	0.873
	  yn(	216)=	0.8796
	  yn(	217)=	0.8863
	  yn(	218)=	0.8931
	  yn(	219)=	0.8999
	  yn(	220)=	0.9069
	  yn(	221)=	0.9139
	  yn(	222)=	0.921
	  yn(	223)=	0.9282
	  yn(	224)=	0.9355
	  yn(	225)=	0.9428
	  yn(	226)=	0.9503
	  yn(	227)=	0.9578
	  yn(	228)=	0.9655
	  yn(	229)=	0.9732
	  yn(	230)=	0.981
	  yn(	231)=	0.9889
	  yn(	232)=	0.9969
	  yn(	233)=	1.005
	  yn(	234)=	1.013
	  yn(	235)=	1.022
	  yn(	236)=	1.03
	  yn(	237)=	1.038
	  yn(	238)=	1.047
	  yn(	239)=	1.056
	  yn(	240)=	1.065
	  yn(	241)=	1.074
	  yn(	242)=	1.083
	  yn(	243)=	1.092
	  yn(	244)=	1.101
	  yn(	245)=	1.111
	  yn(	246)=	1.12
	  yn(	247)=	1.13
	  yn(	248)=	1.139
	  yn(	249)=	1.149
	  yn(	250)=	1.159
	  yn(	251)=	1.17
	  yn(	252)=	1.18
	  yn(	253)=	1.19
	  yn(	254)=	1.201
	  yn(	255)=	1.212
	  yn(	256)=	1.223
	  yn(	257)=	1.234
	  yn(	258)=	1.245
c
		  mn(	1)=-0.1524273019q-03
		  mn(	2)=-0.1524273019q-03
		  mn(	3)=-0.2399223450q-03
		  mn(	4)=-0.2065229207q-03
		  mn(	5)=-0.2038618210q-03
		  mn(	6)=-0.1021415376q-03
		  mn(	7)=	5.391957312q-005
		  mn(	8)=	0.000133552102
		  mn(	9)=	0.0001923957612
		  mn(	10)=	5.940299086q-005
		  mn(	11)=-0.7057478066q-04
		  mn(	12)=-0.1957361553q-03
		  mn(	13)=-0.1053251955q-03
		  mn(	14)=-0.3207085087q-04
		  mn(	15)=	2.803498725q-005
		  mn(	16)=-0.1695365121q-03
		  mn(	17)=-0.4333706503q-03
		  mn(	18)=-0.5433802081q-03
		  mn(	19)=-0.4206389112q-03
		  mn(	20)=-0.1725724564q-03
		  mn(	21)=	0.0003905884493
		  mn(	22)=	0.0004110716963
		  mn(	23)=	0.0003738358622
		  mn(	24)=	0.0002440987556
		  mn(	25)=	0.0001827119006
		  mn(	26)=	9.272259691q-005
		  mn(	27)=	8.644361086q-005
		  mn(	28)=	3.9107112q-005
		  mn(	29)=	5.582399758q-005
		  mn(	30)=	2.045619589q-005
		  mn(	31)=	2.168515224q-005
		  mn(	32)=	2.897205493q-005
		  mn(	33)=	1.639084117q-005
		  mn(	34)=	1.305458033q-005
		  mn(	35)=	1.418166918q-005
		  mn(	36)=	6.436813542q-006
		  mn(	37)=	1.189836085q-005
		  mn(	38)=	1.217334038q-005
		  mn(	39)=	5.625120807q-007
		  mn(	40)=	1.055297079q-005
		  mn(	41)=	5.13266774q-006
		  mn(	42)=	2.666896676q-006
		  mn(	43)=	6.703240244q-006
		  mn(	44)=	3.204760292q-006
		  mn(	45)=	2.832329986q-006
		  mn(	46)=	2.669423561q-006
		  mn(	47)=	5.554927862q-006
		  mn(	48)=-0.1367072964q-05
		  mn(	49)=	6.093788393q-006
		  mn(	50)=-0.1970564037q-06
		  mn(	51)=	2.486507343q-006
		  mn(	52)=	1.788494009q-006
		  mn(	53)=	1.076536017q-006
		  mn(	54)=	1.547939399q-006
		  mn(	55)=	3.118225136q-006
		  mn(	56)=	7.136710083q-007
		  mn(	57)=-0.1863466392q-05
		  mn(	58)=	4.629005982q-006
		  mn(	59)=-0.1250770117q-05
		  mn(	60)=	4.001654718q-006
		  mn(	61)=-0.3891119482q-05
		  mn(	62)=	7.338061616q-006
		  mn(	63)=-0.5853453691q-05
		  mn(	64)=	4.554228978q-006
		  mn(	65)=-0.3989695329q-06
		  mn(	66)=	2.357217756q-007
		  mn(	67)=	2.290409003q-006
		  mn(	68)=-0.1842857399q-05
		  mn(	69)=	2.075303662q-006
		  mn(	70)=	2.287708556q-006
		  mn(	71)=-0.1733206262q-05
		  mn(	72)=	1.428907365q-006
		  mn(	73)=-0.1018213069q-05
		  mn(	74)=-0.4949059562q-06
		  mn(	75)=	5.506115788q-006
		  mn(	76)=-0.4152897421q-05
		  mn(	77)=	1.536423846q-006
		  mn(	78)=	7.660805399q-007
		  mn(	79)=-0.2006094547q-05
		  mn(	80)=	4.006454642q-006
		  mn(	81)=-0.2775141494q-06
		  mn(	82)=	2.723930442q-008
		  mn(	83)=-0.3177846095q-05
		  mn(	84)=	9.31762294q-006
		  mn(	85)=-0.1071569811q-04
		  mn(	86)=	9.689506723q-006
		  mn(	87)=-0.5235758405q-05
		  mn(	88)=	1.585956179q-006
		  mn(	89)=	1.560154939q-006
		  mn(	90)=-0.5261896434q-05
		  mn(	91)=	1.044322901q-005
		  mn(	92)=-0.9596132540q-05
		  mn(	93)=	6.512041612q-006
		  mn(	94)=-0.1693813576q-05
		  mn(	95)=	4.757966665q-007
		  mn(	96)=-0.9859204715q-06
		  mn(	97)=	2.741953236q-006
		  mn(	98)=-0.1878797056q-05
		  mn(	99)=	1.387793512q-006
		  mn(	100)=-0.1393813877q-05
		  mn(	101)=	3.735767225q-006
		  mn(	102)=-0.7960032629q-05
		  mn(	103)=	4.917638956q-007
		  mn(	104)=	2.994570167q-006
		  mn(	105)=-0.3830659164q-05
		  mn(	106)=	3.402913423q-006
		  mn(	107)=-0.2534655876q-05
		  mn(	108)=	1.220741021q-006
		  mn(	109)=-0.1713654575q-05
		  mn(	110)=	2.728168942q-006
		  mn(	111)=-0.2569732393q-05
		  mn(	112)=	2.372683077q-006
		  mn(	113)=	6.643309589q-009
		  mn(	114)=-0.1384538121q-05
		  mn(	115)=	2.628337047q-006
		  mn(	116)=-0.2049675089q-05
		  mn(	117)=	7.557014498q-007
		  mn(	118)=	2.131368532q-006
		  mn(	119)=-0.1998186968q-05
		  mn(	120)=	1.311094605q-006
		  mn(	121)=-0.1517119164q-06
		  mn(	122)=	8.809150761q-007
		  mn(	123)=-0.1635932008q-05
		  mn(	124)=	2.780304358q-006
		  mn(	125)=-0.1634469609q-05
		  mn(	126)=	8.728103937q-007
		  mn(	127)=	1.510507322q-007
		  mn(	128)=	7.090971599q-007
		  mn(	129)=	1.415964489q-007
		  mn(	130)=	1.099441215q-006
		  mn(	131)=-0.1458696617q-05
		  mn(	132)=	1.290647555q-006
		  mn(	133)=-0.6166590866q-06
		  mn(	134)=	3.917620158q-006
		  mn(	135)=-0.5957273842q-05
		  mn(	136)=	5.177586804q-006
		  mn(	137)=	2.600721755q-007
		  mn(	138)=-0.3023698269q-05
		  mn(	139)=	2.904323436q-006
		  mn(	140)=	8.095759409q-007
		  mn(	141)=-0.3067922983q-05
		  mn(	142)=	2.659760748q-006
		  mn(	143)=	1.932717772q-006
		  mn(	144)=-0.3443067892q-06
		  mn(	145)=-0.3546215489q-05
		  mn(	146)=	5.606799317q-006
		  mn(	147)=-0.2705077369q-05
		  mn(	148)=	2.214035954q-006
		  mn(	149)=-0.1538359169q-05
		  mn(	150)=	9.435752797q-007
		  mn(	151)=	7.117802797q-007
		  mn(	152)=	1.123949355q-006
		  mn(	153)=-0.2203542743q-05
		  mn(	154)=	4.649228595q-006
		  mn(	155)=-0.5164061826q-05
		  mn(	156)=	6.955902358q-006
		  mn(	157)=-0.4788222946q-05
		  mn(	158)=	3.009281177q-006
		  mn(	159)=	1.589432214q-006
		  mn(	160)=-0.3099240585q-05
		  mn(	161)=	1.58476775q-006
		  mn(	162)=	5.577203083q-006
		  mn(	163)=-0.2686300413q-04
		  mn(	164)=	6.120818588q-005
		  mn(	165)=-0.3482661527q-04
		  mn(	166)=	6.696582969q-006
		  mn(	167)=	4.0624334q-006
		  mn(	168)=-0.1364679681q-05
		  mn(	169)=	2.53474697q-006
		  mn(	170)=-0.2143997068q-05
		  mn(	171)=	2.787550172q-006
		  mn(	172)=-0.6255001610q-05
		  mn(	173)=	7.936717321q-006
		  mn(	174)=-0.3868442197q-05
		  mn(	175)=	2.392981403q-006
		  mn(	176)=-0.1251098839q-05
		  mn(	177)=-0.2210865803q-05
		  mn(	178)=	4.437512573q-006
		  mn(	179)=-0.2414404935q-05
		  mn(	180)=	1.901494611q-006
		  mn(	181)=-0.1080520262q-05
		  mn(	182)=-0.2268088327q-05
		  mn(	183)=	4.123621892q-006
		  mn(	184)=-0.1768926988q-05
		  mn(	185)=-0.3550874595q-06
		  mn(	186)=	2.154980054q-006
		  mn(	187)=-0.1326924485q-05
		  mn(	188)=-0.1709332219q-06
		  mn(	189)=-0.1287072675q-05
		  mn(	190)=	2.732089752q-006
		  mn(	191)=-0.1827711411q-05
		  mn(	192)=	1.287024571q-006
		  mn(	193)=	2.920029567q-006
		  mn(	194)=-0.7175356279q-05
		  mn(	195)=	9.058908873q-006
		  mn(	196)=-0.9032067325q-06
		  mn(	197)=	1.688217796q-007
		  mn(	198)=-0.3610011570q-06
		  mn(	199)=	1.160119306q-006
		  mn(	200)=	8.23396094q-007
		  mn(	201)=	5.110984171q-007
		  mn(	202)=-0.6156611359q-05
		  mn(	203)=	6.224635636q-006
		  mn(	204)=	1.897806813q-007
		  mn(	205)=-0.3036211723q-06
		  mn(	206)=	3.011880417q-007
		  mn(	207)=	3.261762821q-006
		  mn(	208)=-0.9264942739q-05
		  mn(	209)=	8.868951069q-006
		  mn(	210)=-0.1065054013q-05
		  mn(	211)=-0.8473766156q-06
		  mn(	212)=	1.150155764q-006
		  mn(	213)=-0.2177408122q-06
		  mn(	214)=	3.134353779q-006
		  mn(	215)=-0.8979854359q-05
		  mn(	216)=	8.461929633q-006
		  mn(	217)=-0.1021665487q-05
		  mn(	218)=-0.1317566356q-05
		  mn(	219)=	2.995484704q-006
		  mn(	220)=-0.1709952767q-05
		  mn(	221)=	5.545309251q-007
		  mn(	222)=	2.125854162q-006
		  mn(	223)=-0.6481266053q-05
		  mn(	224)=	5.640608911q-006
		  mn(	225)=	9.952886822q-007
		  mn(	226)=-0.1748468704q-05
		  mn(	227)=	2.701846371q-006
		  mn(	228)=-0.1575401221q-05
		  mn(	229)=	3.122886446q-007
		  mn(	230)=	2.250004897q-006
		  mn(	231)=-0.7427713728q-05
		  mn(	232)=	9.293522365q-006
		  mn(	233)=-0.8782556495q-05
		  mn(	234)=	1.767860055q-005
		  mn(	235)=-0.1724386760q-04
		  mn(	236)=	7.874787961q-007
		  mn(	237)=	1.087440861q-005
		  mn(	238)=-0.1838195127q-05
		  mn(	239)=-0.6850361185q-05
		  mn(	240)=	6.354277015q-006
		  mn(	241)=-0.2677163186q-05
		  mn(	242)=	1.102921471q-006
		  mn(	243)=-0.4956912182q-05
		  mn(	244)=	1.547213116q-005
		  mn(	245)=-0.1966874760q-04
		  mn(	246)=	2.007226097q-005
		  mn(	247)=-0.2469556866q-04
		  mn(	248)=	2.122942164q-005
		  mn(	249)=-0.9614708173q-05
		  mn(	250)=	1.398569444q-005
		  mn(	251)=-0.1325532962q-04
		  mn(	252)=-0.5950502285q-07
		  mn(	253)=	1.024096979q-005
		  mn(	254)=-0.9572533393q-05
		  mn(	255)=	6.693742982q-006
		  mn(	256)=-0.2445001490q-05
		  mn(	257)=-0.1595892729q-06
		  mn(	258)=-0.1595892729q-06
C
C	Spline	coefficinets Imag n

	  yk(	1)=	1.934
	  yk(	2)=	1.937
	  yk(	3)=	1.921
	  yk(	4)=	1.895
	  yk(	5)=	1.871
	  yk(	6)=	1.858
	  yk(	7)=	1.862
	  yk(	8)=	1.885
	  yk(	9)=	1.918
	  yk(	10)=	1.945
	  yk(	11)=	1.956
	  yk(	12)=	1.95
	  yk(	13)=	1.936
	  yk(	14)=	1.925
	  yk(	15)=	1.915
	  yk(	16)=	1.893
	  yk(	17)=	1.85
	  yk(	18)=	1.795
	  yk(	19)=	1.756
	  yk(	20)=	1.764
	  yk(	21)=	1.836
	  yk(	22)=	1.962
	  yk(	23)=	2.115
	  yk(	24)=	2.274
	  yk(	25)=	2.427
	  yk(	26)=	2.573
	  yk(	27)=	2.71
	  yk(	28)=	2.841
	  yk(	29)=	2.966
	  yk(	30)=	3.087
	  yk(	31)=	3.203
	  yk(	32)=	3.316
	  yk(	33)=	3.426
	  yk(	34)=	3.534
	  yk(	35)=	3.639
	  yk(	36)=	3.742
	  yk(	37)=	3.844
	  yk(	38)=	3.943
	  yk(	39)=	4.042
	  yk(	40)=	4.138
	  yk(	41)=	4.234
	  yk(	42)=	4.328
	  yk(	43)=	4.422
	  yk(	44)=	4.514
	  yk(	45)=	4.605
	  yk(	46)=	4.696
	  yk(	47)=	4.786
	  yk(	48)=	4.875
	  yk(	49)=	4.963
	  yk(	50)=	5.05
	  yk(	51)=	5.137
	  yk(	52)=	5.224
	  yk(	53)=	5.31
	  yk(	54)=	5.395
	  yk(	55)=	5.48
	  yk(	56)=	5.565
	  yk(	57)=	5.649
	  yk(	58)=	5.732
	  yk(	59)=	5.815
	  yk(	60)=	5.898
	  yk(	61)=	5.981
	  yk(	62)=	6.063
	  yk(	63)=	6.145
	  yk(	64)=	6.227
	  yk(	65)=	6.308
	  yk(	66)=	6.389
	  yk(	67)=	6.47
	  yk(	68)=	6.55
	  yk(	69)=	6.63
	  yk(	70)=	6.711
	  yk(	71)=	6.79
	  yk(	72)=	6.87
	  yk(	73)=	6.949
	  yk(	74)=	7.029
	  yk(	75)=	7.108
	  yk(	76)=	7.186
	  yk(	77)=	7.265
	  yk(	78)=	7.344
	  yk(	79)=	7.422
	  yk(	80)=	7.5
	  yk(	81)=	7.578
	  yk(	82)=	7.656
	  yk(	83)=	7.734
	  yk(	84)=	7.811
	  yk(	85)=	7.889
	  yk(	86)=	7.966
	  yk(	87)=	8.043
	  yk(	88)=	8.121
	  yk(	89)=	8.198
	  yk(	90)=	8.274
	  yk(	91)=	8.351
	  yk(	92)=	8.428
	  yk(	93)=	8.504
	  yk(	94)=	8.581
	  yk(	95)=	8.657
	  yk(	96)=	8.733
	  yk(	97)=	8.809
	  yk(	98)=	8.885
	  yk(	99)=	8.961
	  yk(	100)=	9.037
	  yk(	101)=	9.113
	  yk(	102)=	9.189
	  yk(	103)=	9.262
	  yk(	104)=	9.336
	  yk(	105)=	9.41
	  yk(	106)=	9.483
	  yk(	107)=	10.08
	  yk(	108)=	10.15
	  yk(	109)=	10.23
	  yk(	110)=	10.3
	  yk(	111)=	10.37
	  yk(	112)=	10.44
	  yk(	113)=	10.51
	  yk(	114)=	10.59
	  yk(	115)=	10.66
	  yk(	116)=	10.73
	  yk(	117)=	10.8
	  yk(	118)=	10.87
	  yk(	119)=	10.94
	  yk(	120)=	11.02
	  yk(	121)=	11.09
	  yk(	122)=	11.16
	  yk(	123)=	11.23
	  yk(	124)=	11.3
	  yk(	125)=	11.37
	  yk(	126)=	11.44
	  yk(	127)=	11.51
	  yk(	128)=	11.58
	  yk(	129)=	11.66
	  yk(	130)=	11.73
	  yk(	131)=	11.8
	  yk(	132)=	11.87
	  yk(	133)=	11.94
	  yk(	134)=	12.01
	  yk(	135)=	12.08
	  yk(	136)=	12.15
	  yk(	137)=	12.22
	  yk(	138)=	12.29
	  yk(	139)=	12.36
	  yk(	140)=	12.43
	  yk(	141)=	12.5
	  yk(	142)=	12.57
	  yk(	143)=	12.64
	  yk(	144)=	12.71
	  yk(	145)=	12.78
	  yk(	146)=	12.85
	  yk(	147)=	12.92
	  yk(	148)=	12.99
	  yk(	149)=	13.06
	  yk(	150)=	13.13
	  yk(	151)=	13.2
	  yk(	152)=	13.27
	  yk(	153)=	13.34
	  yk(	154)=	13.41
	  yk(	155)=	13.48
	  yk(	156)=	13.55
	  yk(	157)=	13.62
	  yk(	158)=	13.69
	  yk(	159)=	13.75
	    yk(	160)=	13.82
	    yk(	161)=	13.89
	    yk(	162)=	13.96
	    yk(	163)=	14.03
	    yk(	164)=	14.09
	    yk(	165)=	14.15
	    yk(	166)=	14.24
	    yk(	167)=	14.28
	    yk(	168)=	14.38
	    yk(	169)=	14.42
	    yk(	170)=	14.52
	    yk(	171)=	14.56
	    yk(	172)=	14.66
	    yk(	173)=	14.71
	    yk(	174)=	14.8
	    yk(	175)=	14.85
	    yk(	176)=	14.95
	    yk(	177)=	15
	    yk(	178)=	15.1
	    yk(	179)=	15.15
	    yk(	180)=	15.2
	    yk(	181)=	15.31
	    yk(	182)=	15.36
	    yk(	183)=	15.46
	    yk(	184)=	15.52
	    yk(	185)=	15.57
	    yk(	186)=	15.68
	    yk(	187)=	15.73
	    yk(	188)=	15.79
	    yk(	189)=	15.84
	    yk(	190)=	15.95
	    yk(	191)=	16.01
	    yk(	192)=	16.07
	    yk(	193)=	16.13
	    yk(	194)=	16.24
	    yk(	195)=	16.3
	    yk(	196)=	16.36
	    yk(	197)=	16.42
	    yk(	198)=	16.54
	    yk(	199)=	16.6
	    yk(	200)=	16.66
	    yk(	201)=	16.72
	    yk(	202)=	16.78
	    yk(	203)=	16.84
	    yk(	204)=	16.91
	    yk(	205)=	17.03
	    yk(	206)=	17.1
	    yk(	207)=	17.16
	    yk(	208)=	17.23
	    yk(	209)=	17.29
	    yk(	210)=	17.36
	    yk(	211)=	17.43
	    yk(	212)=	17.49
	    yk(	213)=	17.56
	    yk(	214)=	17.63
	    yk(	215)=	17.7
	    yk(	216)=	17.77
	    yk(	217)=	17.84
	    yk(	218)=	17.91
	    yk(	219)=	17.98
	    yk(	220)=	18.05
	    yk(	221)=	18.12
	    yk(	222)=	18.19
	    yk(	223)=	18.26
	    yk(	224)=	18.34
	    yk(	225)=	18.41
	    yk(	226)=	18.49
	    yk(	227)=	18.56
	    yk(	228)=	18.64
	    yk(	229)=	18.71
	    yk(	230)=	18.79
	    yk(	231)=	18.87
	    yk(	232)=	18.95
	    yk(	233)=	19.03
	    yk(	234)=	19.1
	    yk(	235)=	19.18
	    yk(	236)=	19.27
	    yk(	237)=	19.35
	    yk(	238)=	19.43
	    yk(	239)=	19.51
	    yk(	240)=	19.6
	    yk(	241)=	19.68
	    yk(	242)=	19.76
	    yk(	243)=	19.85
	    yk(	244)=	19.94
	    yk(	245)=	20.02
	    yk(	246)=	20.11
	    yk(	247)=	20.2
	    yk(	248)=	20.29
	    yk(	249)=	20.38
	    yk(	250)=	20.47
	    yk(	251)=	20.56
	    yk(	252)=	20.65
	    yk(	253)=	20.75
	    yk(	254)=	20.84
	    yk(	255)=	20.94
	    yk(	256)=	21.03
	    yk(	257)=	21.13
	    yk(	258)=	21.23
C
	    mk(	1)=-0.2056974797q-03
	    mk(	2)=-0.2056974797q-03
	    mk(	3)=-0.1068065588q-03
	    mk(	4)=	2.688397435q-005
	    mk(	5)=	0.0001183978268
	    mk(	6)=	0.0001703125192
	    mk(	7)=	0.000213886397
	    mk(	8)=	0.0001172744509
	    mk(	9)=-0.7917934408q-04
	    mk(	10)=-0.1701008492q-03
	    mk(	11)=-0.1914573904q-03
	    mk(	12)=-0.8722156273q-04
	    mk(	13)=	6.038513654q-005
	    mk(	14)=	2.832088618q-005
	    mk(	15)=-0.1230924566q-03
	    mk(	16)=-0.2418463267q-03
	    mk(	17)=-0.1676036327q-03
	    mk(	18)=	0.0001676511814
	    mk(	19)=	0.0004789119822
	    mk(	20)=	0.000722128429
	    mk(	21)=	0.0005407566278
	    mk(	22)=	0.0002561071177
	    mk(	23)=	5.128294627q-005
	    mk(	24)=-0.5852862237q-04
	    mk(	25)=-0.1062370393q-03
	    mk(	26)=-0.5809402905q-04
	    mk(	27)=-0.9587203859q-04
	    mk(	28)=-0.2927795607q-05
	    mk(	29)=-0.1187081651q-03
	    mk(	30)=-0.1380908200q-04
	    mk(	31)=	3.011215267q-005
	    mk(	32)=-0.1097369471q-03
	    mk(	33)=	4.058991964q-006
	    mk(	34)=-0.3597852373q-04
	    mk(	35)=-0.3887003001q-04
	    mk(	36)=	8.135027812q-005
	    mk(	37)=-0.1060730336q-03
	    mk(	38)=-0.2741568969q-04
	    mk(	39)=	3.854510182q-005
	    mk(	40)=-0.4791596583q-04
	    mk(	41)=-0.1335478074q-04
	    mk(	42)=	4.216835303q-005
	    mk(	43)=-0.8320096898q-04
	    mk(	44)=	1.740474506q-005
	    mk(	45)=	3.575885323q-005
	    mk(	46)=-0.6391073138q-04
	    mk(	47)=	1.494621512q-005
	    mk(	48)=	5.342520799q-005
	    mk(	49)=-0.1642259652q-03
	    mk(	50)=	0.000153261479
	    mk(	51)=-0.5989278296q-04
	    mk(	52)=-0.4435751202q-04
	    mk(	53)=	4.905884956q-005
	    mk(	54)=-0.6075869421q-04
	    mk(	55)=	7.080490163q-005
	    mk(	56)=-0.5006182511q-04
	    mk(	57)=-0.5082416871q-04
	    mk(	58)=	7.85017941q-005
	    mk(	59)=-0.7959787112q-04
	    mk(	60)=	0.0001260874676
	    mk(	61)=-0.2125139603q-03
	    mk(	62)=	0.000239937889
	    mk(	63)=-0.1979182266q-03
	    mk(	64)=	0.0001031185805
	    mk(	65)=-0.4357142927q-04
	    mk(	66)=-0.3104021235q-04
	    mk(	67)=	6.618865489q-005
	    mk(	68)=-0.3827516211q-04
	    mk(	69)=-0.1255056031q-04
	    mk(	70)=	4.887492889q-005
	    mk(	71)=-0.1031202277q-04
	    mk(	72)=-0.4127333641q-04
	    mk(	73)=	1.773812258q-005
	    mk(	74)=-0.6318548379q-04
	    mk(	75)=	8.710424826q-005
	    mk(	76)=-0.2624160969q-04
	    mk(	77)=-0.1024047049q-04
	    mk(	78)=-0.2261527568q-04
	    mk(	79)=-0.4390308411q-04
	    mk(	80)=	0.0001131857357
	    mk(	81)=-0.5364249211q-04
	    mk(	82)=	1.510185882q-005
	    mk(	83)=-0.9128112172q-04
	    mk(	84)=	0.0002076084199
	    mk(	85)=-0.2548356662q-03
	    mk(	86)=	0.000175775122
	    mk(	87)=-0.4321478085q-04
	    mk(	88)=-0.2045299000q-04
	    mk(	89)=-0.1609070508q-04
	    mk(	90)=-0.5074319542q-04
	    mk(	91)=	0.0001989441456
	    mk(	92)=-0.2746986030q-03
	    mk(	93)=	0.0002247755234
	    mk(	94)=-0.1117551322q-03
	    mk(	95)=	2.78139445q-005
	    mk(	96)=-0.1650973092q-04
	    mk(	97)=	2.255585476q-005
	    mk(	98)=-0.2901001452q-04
	    mk(	99)=	2.179430865q-005
	    mk(	100)=-0.1001771327q-04
	    mk(	101)=	8.80581103q-006
	    mk(	102)=-0.3339589332q-04
	    mk(	103)=	1.557357936q-008
	    mk(	104)=	2.618454061q-005
	    mk(	105)=-0.4469140921q-04
	    mk(	106)=	2.661110254q-005
	    mk(	107)=-0.3871237640q-04
	    mk(	108)=	0.0002236356442
	    mk(	109)=-0.2365730980q-03
	    mk(	110)=	6.226602689q-005
	    mk(	111)=	1.679391785q-006
	    mk(	112)=-0.5348922556q-04
	    mk(	113)=	0.0002305097825
	    mk(	114)=-0.2400842347q-03
	    mk(	115)=	6.842699935q-005
	    mk(	116)=-0.1182175862q-04
	    mk(	117)=	2.240200245q-006
	    mk(	118)=-0.5174987641q-04
	    mk(	119)=	0.0002304669954
	    mk(	120)=-0.2385694852q-03
	    mk(	121)=	6.840085604q-005
	    mk(	122)=-0.4986672689q-05
	    mk(	123)=-0.1669400647q-04
	    mk(	124)=	1.873465661q-005
	    mk(	125)=-0.2433643417q-04
	    mk(	126)=	2.657102842q-005
	    mk(	127)=-0.4614250040q-04
	    mk(	128)=	0.0001980460975
	    mk(	129)=-0.1975890376q-03
	    mk(	130)=	3.485954556q-005
	    mk(	131)=	6.32749399q-006
	    mk(	132)=-0.1648188110q-04
	    mk(	133)=	9.663241093q-006
	    mk(	134)=	2.461637308q-005
	    mk(	135)=-0.5769647750q-04
	    mk(	136)=	5.818428873q-005
	    mk(	137)=-0.2391347114q-04
	    mk(	138)=-0.1119725147q-04
	    mk(	139)=	2.022202635q-005
	    mk(	140)=-0.1403739152q-04
	    mk(	141)=-0.1165710249q-04
	    mk(	142)=	1.355163327q-005
	    mk(	143)=	1.52492754q-005
	    mk(	144)=-0.1169590712q-04
	    mk(	145)=-0.1528974344q-04
	    mk(	146)=	2.642432002q-005
	    mk(	147)=-0.2495404711q-04
	    mk(	148)=	2.746008004q-005
	    mk(	149)=-0.1475548351q-04
	    mk(	150)=-0.1390843013q-04
	    mk(	151)=	2.521323699q-005
	    mk(	152)=-0.1376991226q-04
	    mk(	153)=-0.1470449367q-04
	    mk(	154)=	2.82899601q-005
	    mk(	155)=-0.2187864294q-04
	    mk(	156)=	1.533923347q-005
	    mk(	157)=	4.361476917q-005
	    mk(	158)=-0.2324228631q-03
	    mk(	159)=	0.000246954308
	    mk(	160)=-0.6895262338q-04
	    mk(	161)=-0.1305734226q-04
	    mk(	162)=	7.785949066q-005
	    mk(	163)=-0.3386467637q-03
	    mk(	164)=	0.0007748249087
	    mk(	165)=-0.4090410153q-03
	    mk(	166)=-0.1353338039q-03
	    mk(	167)=	0.0004721961759
	    mk(	168)=-0.4939666099q-03
	    mk(	169)=	0.000469479856
	    mk(	170)=-0.4648155946q-03
	    mk(	171)=	0.0003573443869
	    mk(	172)=-0.1699676120q-03
	    mk(	173)=-0.7975126030q-04
	    mk(	174)=	0.0001329450659
	    mk(	175)=-0.3595832188q-04
	    mk(	176)=	1.038513955q-005
	    mk(	177)=-0.5230534190q-04
	    mk(	178)=	4.267573689q-005
	    mk(	179)=-0.5638982736q-04
	    mk(	180)=	0.0001429472027
	    mk(	181)=-0.1355180711q-03
	    mk(	182)=-0.1134252412q-03
	    mk(	183)=	0.0003031247883
	    mk(	184)=-0.4175554897q-03
	    mk(	185)=	0.0002303931365
	    mk(	186)=-0.2426309411q-03
	    mk(	187)=	0.0004081644099
	    mk(	188)=-0.3900467590q-03
	    mk(	189)=	9.016395449q-005
	    mk(	190)=	6.935937211q-005
	    mk(	191)=-0.3375169639q-04
	    mk(	192)=	2.543160308q-005
	    mk(	193)=-0.1076663665q-03
	    mk(	194)=	4.762030851q-005
	    mk(	195)=	4.960001629q-005
	    mk(	196)=-0.2099677794q-04
	    mk(	197)=-0.3301101301q-05
	    mk(	198)=-0.7651052891q-05
	    mk(	199)=-0.2711406005q-05
	    mk(	200)=-0.1775324660q-04
	    mk(	201)=	3.768255664q-005
	    mk(	202)=-0.1679410543q-03
	    mk(	203)=	0.0003553848668
	    mk(	204)=-0.2890638089q-03
	    mk(	205)=	0.0002737233476
	    mk(	206)=-0.3561114657q-03
	    mk(	207)=	0.0003673192324
	    mk(	208)=-0.4124951753q-03
	    mk(	209)=	0.0003310844259
	    mk(	210)=-0.3310944402q-04
	    mk(	211)=-0.2335613916q-03
	    mk(	212)=	0.0002396543562
	    mk(	213)=-0.7701317672q-04
	    mk(	214)=	3.338422991q-005
	    mk(	215)=-0.9185894828q-04
	    mk(	216)=	7.647561007q-005
	    mk(	217)=-0.2759584847q-04
	    mk(	218)=	3.501512179q-008
	    mk(	219)=-0.6400008719q-05
	    mk(	220)=-0.7786476221q-05
	    mk(	221)=	4.509850835q-006
	    mk(	222)=-0.4272364763q-04
	    mk(	223)=	0.0001336504311
	    mk(	224)=-0.1848693440q-03
	    mk(	225)=	0.0002447843083
	    mk(	226)=-0.2707702816q-03
	    mk(	227)=	0.0002577210848
	    mk(	228)=-0.2537834158q-03
	    mk(	229)=	0.0001952617005
	    mk(	230)=-0.3808913961q-04
	    mk(	231)=-0.7523675594q-04
	    mk(	232)=	0.0001057710681
	    mk(	233)=-0.1885282923q-03
	    mk(	234)=	0.0001297305257
	    mk(	235)=	0.0001188394073
	    mk(	236)=-0.1648188147q-03
	    mk(	237)=	4.317438174q-005
	    mk(	238)=-0.3904107566q-04
	    mk(	239)=	8.257237207q-005
	    mk(	240)=-0.8411750069q-04
	    mk(	241)=-0.1123723625q-04
	    mk(	242)=	9.918423324q-005
	    mk(	243)=	3.335442452q-006
	    mk(	244)=-0.1433911623q-03
	    mk(	245)=	0.0001346720585
	    mk(	246)=-0.2687995515q-04
	    mk(	247)=-0.5700657294q-04
	    mk(	248)=	5.433887629q-005
	    mk(	249)=-0.2314161276q-04
	    mk(	250)=	8.946712047q-006
	    mk(	251)=-0.4165567041q-04
	    mk(	252)=	0.0001284365424
	    mk(	253)=-0.1497577798q-03
	    mk(	254)=	9.574634766q-005
	    mk(	255)=-0.9787508130q-04
	    mk(	256)=	0.0001057708793
	    mk(	257)=-0.2690996341q-04
	    mk(	258)=-0.2690996341q-04

	    else
c       	SC Single Crystalline Au sample Spline	coefficinets
c        	n=Real n
		    yn(	1)=	1.562
		    yn(	2)=	1.609
		    yn(	3)=	1.642
		    yn(	4)=	1.656
		    yn(	5)=	1.65
		    yn(	6)=	1.627
		    yn(	7)=	1.595
		    yn(	8)=	1.566
		    yn(	9)=	1.55
		    yn(	10)=	1.549
		    yn(	11)=	1.554
		    yn(	12)=	1.552
		    yn(	13)=	1.533
		    yn(	14)=	1.503
		    yn(	15)=	1.473
		    yn(	16)=	1.444
		    yn(	17)=	1.399
		    yn(	18)=	1.315
		    yn(	19)=	1.179
		    yn(	20)=	1.003
		    yn(	21)=	0.8181
		    yn(	22)=	0.6582
		    yn(	23)=	0.5376
		    yn(	24)=	0.4515
		    yn(	25)=	0.3897
		    yn(	26)=	0.3437
		    yn(	27)=	0.3078
		    yn(	28)=	0.279
		    yn(	29)=	0.2553
		    yn(	30)=	0.2353
		    yn(	31)=	0.2184
		    yn(	32)=	0.204
		    yn(	33)=	0.1916
		    yn(	34)=	0.1809
		    yn(	35)=	0.1717
		    yn(	36)=	0.1638
		    yn(	37)=	0.157
		    yn(	38)=	0.1511
		    yn(	39)=	0.146
		    yn(	40)=	0.1417
		    yn(	41)=	0.138
		    yn(	42)=	0.1348
		    yn(	43)=	0.1322
		    yn(	44)=	0.13
		    yn(	45)=	0.1282
		    yn(	46)=	0.1268
		    yn(	47)=	0.1257
		    yn(	48)=	0.1248
		    yn(	49)=	0.1243
		    yn(	50)=	0.1239
		    yn(	51)=	0.1238
		    yn(	52)=	0.1239
		    yn(	53)=	0.1241
		    yn(	54)=	0.1245
		    yn(	55)=	0.1251
		    yn(	56)=	0.1258
		    yn(	57)=	0.1266
		    yn(	58)=	0.1275
		    yn(	59)=	0.1286
		    yn(	60)=	0.1297
		    yn(	61)=	0.1309
		    yn(	62)=	0.1323
		    yn(	63)=	0.1337
		    yn(	64)=	0.1351
		    yn(	65)=	0.1367
		    yn(	66)=	0.1383
		    yn(	67)=	0.14
		    yn(	68)=	0.1417
		    yn(	69)=	0.1435
		    yn(	70)=	0.1453
		    yn(	71)=	0.1472
		    yn(	72)=	0.1491
		    yn(	73)=	0.1511
		    yn(	74)=	0.1531
		    yn(	75)=	0.1552
		    yn(	76)=	0.1573
		    yn(	77)=	0.1594
		    yn(	78)=	0.1616
		    yn(	79)=	0.1638
		    yn(	80)=	0.166
		    yn(	81)=	0.1683
		    yn(	82)=	0.1705
		    yn(	83)=	0.1729
		    yn(	84)=	0.1752
		    yn(	85)=	0.1776
		    yn(	86)=	0.18
		    yn(	87)=	0.1824
		    yn(	88)=	0.1849
		    yn(	89)=	0.1874
		    yn(	90)=	0.1899
		    yn(	91)=	0.1924
		    yn(	92)=	0.1949
		    yn(	93)=	0.1975
		    yn(	94)=	0.2001
		    yn(	95)=	0.2027
		    yn(	96)=	0.2053
		    yn(	97)=	0.2079
		    yn(	98)=	0.2106
		    yn(	99)=	0.2133
		    yn(	100)=	0.216
		    yn(	101)=	0.2187
		    yn(	102)=	0.2214
		    yn(	103)=	0.2241
		    yn(	104)=	0.2269
		    yn(	105)=	0.2297
		    yn(	106)=	0.2325
		    yn(	107)=	0.2553
		    yn(	108)=	0.2582
		    yn(	109)=	0.2611
		    yn(	110)=	0.2644
		    yn(	111)=	0.2678
		    yn(	112)=	0.2712
		    yn(	113)=	0.2747
		    yn(	114)=	0.2782
		    yn(	115)=	0.2818
		    yn(	116)=	0.2854
		    yn(	117)=	0.2891
		    yn(	118)=	0.2928
		    yn(	119)=	0.2965
		    yn(	120)=	0.3003
		    yn(	121)=	0.3042
		    yn(	122)=	0.3081
		    yn(	123)=	0.312
		    yn(	124)=	0.3161
		    yn(	125)=	0.3201
		    yn(	126)=	0.3242
		    yn(	127)=	0.3284
		    yn(	128)=	0.3326
		    yn(	129)=	0.3369
		    yn(	130)=	0.3412
		    yn(	131)=	0.3457
		    yn(	132)=	0.3501
		    yn(	133)=	0.3546
		    yn(	134)=	0.3592
		    yn(	135)=	0.3639
		    yn(	136)=	0.3686
		    yn(	137)=	0.3733
		    yn(	138)=	0.3782
			  yn(	139)=	0.3831
			  yn(	140)=	0.3881
			  yn(	141)=	0.3931
			  yn(	142)=	0.3982
			  yn(	143)=	0.4034
			  yn(	144)=	0.4087
			  yn(	145)=	0.414
			  yn(	146)=	0.4194
			  yn(	147)=	0.4249
			  yn(	148)=	0.4305
			  yn(	149)=	0.4361
			  yn(	150)=	0.4418
			  yn(	151)=	0.4476
			  yn(	152)=	0.4535
			  yn(	153)=	0.4595
			  yn(	154)=	0.4655
			  yn(	155)=	0.4717
			  yn(	156)=	0.4779
			  yn(	157)=	0.4842
			  yn(	158)=	0.4906
			  yn(	159)=	0.4971
			  yn(	160)=	0.5037
			  yn(	161)=	0.5104
			  yn(	162)=	0.5172
			  yn(	163)=	0.5241
			  yn(	164)=	0.5294
			  yn(	165)=	0.5345
			  yn(	166)=	0.541
			  yn(	167)=	0.5443
			  yn(	168)=	0.551
			  yn(	169)=	0.5544
			  yn(	170)=	0.5613
			  yn(	171)=	0.5647
			  yn(	172)=	0.5718
			  yn(	173)=	0.5754
			  yn(	174)=	0.5827
			  yn(	175)=	0.5864
			  yn(	176)=	0.5939
			  yn(	177)=	0.5977
			  yn(	178)=	0.6054
			  yn(	179)=	0.6093
			  yn(	180)=	0.6133
			  yn(	181)=	0.6213
			  yn(	182)=	0.6254
			  yn(	183)=	0.6336
			  yn(	184)=	0.6378
			  yn(	185)=	0.6421
			  yn(	186)=	0.6507
			  yn(	187)=	0.655
			  yn(	188)=	0.6594
			  yn(	189)=	0.6639
			  yn(	190)=	0.6729
			  yn(	191)=	0.6775
			  yn(	192)=	0.6822
			  yn(	193)=	0.6869
			  yn(	194)=	0.6964
			  yn(	195)=	0.7012
			  yn(	196)=	0.7061
			  yn(	197)=	0.711
			  yn(	198)=	0.7211
			  yn(	199)=	0.7262
			  yn(	200)=	0.7313
			  yn(	201)=	0.7365
			  yn(	202)=	0.7418
			  yn(	203)=	0.7471
			  yn(	204)=	0.7525
			  yn(	205)=	0.7634
			  yn(	206)=	0.7689
			  yn(	207)=	0.7745
			  yn(	208)=	0.7802
			  yn(	209)=	0.786
			  yn(	210)=	0.7918
			  yn(	211)=	0.7976
			  yn(	212)=	0.8036
			  yn(	213)=	0.8096
			  yn(	214)=	0.8156
			  yn(	215)=	0.8218
			  yn(	216)=	0.828
			  yn(	217)=	0.8342
			  yn(	218)=	0.8406
			  yn(	219)=	0.847
			  yn(	220)=	0.8535
			  yn(	221)=	0.86
			  yn(	222)=	0.8667
			  yn(	223)=	0.8734
			  yn(	224)=	0.8802
			  yn(	225)=	0.8871
			  yn(	226)=	0.894
			  yn(	227)=	0.9011
			  yn(	228)=	0.9082
			  yn(	229)=	0.9154
			  yn(	230)=	0.9227
			  yn(	231)=	0.9301
			  yn(	232)=	0.9376
			  yn(	233)=	0.9451
			  yn(	234)=	0.9528
			  yn(	235)=	0.9605
			  yn(	236)=	0.9683
			  yn(	237)=	0.9763
			  yn(	238)=	0.9843
			  yn(	239)=	0.9924
			  yn(	240)=	1.001
			  yn(	241)=	1.009
			  yn(	242)=	1.017
			  yn(	243)=	1.026
			  yn(	244)=	1.035
			  yn(	245)=	1.043
			  yn(	246)=	1.052
			  yn(	247)=	1.061
			  yn(	248)=	1.07
			  yn(	249)=	1.08
			  yn(	250)=	1.089
			  yn(	251)=	1.098
			  yn(	252)=	1.108
			  yn(	253)=	1.118
			  yn(	254)=	1.128
			  yn(	255)=	1.138
			  yn(	256)=	1.148
			  yn(	257)=	1.158
			  yn(	258)=	1.168
c

			  mn(	1)=-0.1245561126q-03
			  mn(	2)=-0.1245561126q-03
			  mn(	3)=-0.2012568965q-03
			  mn(	4)=-0.2086908594q-03
			  mn(	5)=-0.1732610296q-03
			  mn(	6)=-0.1049726901q-03
			  mn(	7)=	3.716833379q-005
			  mn(	8)=	0.0001424316554
			  mn(	9)=	0.0001738580517
			  mn(	10)=	6.387356337q-005
			  mn(	11)=-0.6859806460q-04
			  mn(	12)=-0.2093654918q-03
			  mn(	13)=-0.1120377185q-03
			  mn(	14)=	7.95959396q-006
			  mn(	15)=	5.737959546q-005
			  mn(	16)=-0.1511173547q-03
			  mn(	17)=-0.4153923482q-03
			  mn(	18)=-0.5691245305q-03
			  mn(	19)=-0.3998182133q-03
			  mn(	20)=-0.1732426697q-03
			  mn(	21)=	0.0003705188885
			  mn(	22)=	0.00039697138
			  mn(	23)=	0.0003506137853
			  mn(	24)=	0.0002274908354
			  mn(	25)=	0.0001715019514
			  mn(	26)=	8.020038773q-005
			  mn(	27)=	8.264325156q-005
			  mn(	28)=	3.649838617q-005
			  mn(	29)=	5.183298811q-005
			  mn(	30)=	2.181534156q-005
			  mn(	31)=	1.963400408q-005
			  mn(	32)=	2.875700484q-005
			  mn(	33)=	1.447066802q-005
			  mn(	34)=	1.51095061q-005
			  mn(	35)=	1.379696645q-005
			  mn(	36)=	5.865069427q-006
			  mn(	37)=	1.293621391q-005
			  mn(	38)=	9.913243222q-006
			  mn(	39)=	4.011450864q-006
			  mn(	40)=	9.04153167q-006
			  mn(	41)=	2.923914511q-006
			  mn(	42)=	6.670597012q-006
			  mn(	43)=	4.452439979q-006
			  mn(	44)=	3.379751789q-006
			  mn(	45)=	3.794329437q-006
			  mn(	46)=	4.036785639q-006
			  mn(	47)=	1.69392801q-007
			  mn(	48)=	6.062427409q-006
			  mn(	49)=-0.4585804443q-06
			  mn(	50)=	3.782796168q-006
			  mn(	51)=	2.122178375q-006
			  mn(	52)=	2.514431197q-007
			  mn(	53)=	2.58461856q-006
			  mn(	54)=	1.900263696q-006
			  mn(	55)=	1.110794106q-006
			  mn(	56)=	1.016002781q-006
			  mn(	57)=-0.1110017853q-06
			  mn(	58)=	4.06069729q-006
			  mn(	59)=-0.2109946649q-05
			  mn(	60)=	2.891307593q-006
			  mn(	61)=-0.3103421283q-06
			  mn(	62)=	3.419512263q-006
			  mn(	63)=-0.4220785397q-05
			  mn(	64)=	5.600058803q-006
			  mn(	65)=-0.2227528568q-05
			  mn(	66)=	1.308567299q-006
			  mn(	67)=	6.924733636q-007
			  mn(	68)=	1.311887598q-006
			  mn(	69)=-0.1946027768q-05
			  mn(	70)=	4.254207103q-006
			  mn(	71)=-0.2087087851q-05
			  mn(	72)=	1.761883482q-006
			  mn(	73)=-0.1155178119q-05
			  mn(	74)=	5.31829998q-007
			  mn(	75)=	2.360817671q-006
			  mn(	76)=-0.1552544322q-05
			  mn(	77)=	1.403638175q-006
			  mn(	78)=-0.4714864649q-06
			  mn(	79)=-0.1916647019q-05
			  mn(	80)=	5.716862821q-006
			  mn(	81)=-0.4872569539q-05
			  mn(	82)=	4.969801675q-006
			  mn(	83)=-0.5353010929q-05
			  mn(	84)=	8.064836021q-006
			  mn(	85)=-0.7727012389q-05
			  mn(	86)=	5.147519974q-006
			  mn(	87)=-0.1971880920q-06
			  mn(	88)=-0.6137218928q-06
			  mn(	89)=	6.737655991q-008
			  mn(	90)=-0.2159772240q-05
			  mn(	91)=	6.066784999q-006
			  mn(	92)=-0.6698422523q-05
			  mn(	93)=	6.240665985q-006
			  mn(	94)=-0.2893421642q-05
			  mn(	95)=	8.199605866q-007
			  mn(	96)=-0.9658626330q-06
			  mn(	97)=	2.504224748q-006
			  mn(	98)=-0.1476178525q-05
			  mn(	99)=	8.663161872q-007
			  mn(	100)=-0.2757184830q-06
			  mn(	101)=-0.9937820269q-07
			  mn(	102)=	3.809138601q-007
			  mn(	103)=	6.302161968q-007
			  mn(	104)=	6.177404837q-007
			  mn(	105)=-0.8269604926q-06
			  mn(	106)=	1.838722039q-007
			  mn(	107)=	1.953304203q-007
			  mn(	108)=-0.1314023050q-05
			  mn(	109)=	5.593663053q-006
			  mn(	110)=	3.420028262q-007
			  mn(	111)=-0.3802321191q-06
			  mn(	112)=	1.945082501q-006
			  mn(	113)=-0.5034963908q-06
			  mn(	114)=	1.061914818q-006
			  mn(	115)=-0.5658039066q-06
			  mn(	116)=	2.295890398q-006
			  mn(	117)=-0.1353052415q-05
			  mn(	118)=	2.294584904q-007
			  mn(	119)=	1.767192934q-006
			  mn(	120)=	2.059537394q-007
			  mn(	121)=	4.469433179q-007
			  mn(	122)=-0.3673122810q-06
			  mn(	123)=	2.818782999q-006
			  mn(	124)=-0.1897982404q-05
			  mn(	125)=	7.442075011q-007
			  mn(	126)=	1.847841678q-006
			  mn(	127)=-0.7544378103q-07
			  mn(	128)=	8.041819273q-007
			  mn(	129)=-0.2313260706q-06
			  mn(	130)=	2.681357973q-006
			  mn(	131)=-0.1669101473q-05
			  mn(	132)=	7.906170673q-007
			  mn(	133)=	1.295970625q-006
			  mn(	134)=	3.039897132q-006
			  mn(	135)=-0.3995633004q-05
			  mn(	136)=	3.0460976q-006
			  mn(	137)=	2.025860059q-006
			  mn(	138)=-0.2231455894q-05
			  mn(	139)=	3.515536593q-006
			  mn(	140)=-0.1893974097q-05
			  mn(	141)=	6.647143614q-007
			  mn(	142)=	1.700775702q-006
			  mn(	143)=	2.673521752q-006
			  mn(	144)=-0.1630506933q-05
			  mn(	145)=	3.035713364q-007
			  mn(	146)=	2.810269876q-006
			  mn(	147)=-0.4704328427q-06
			  mn(	148)=	1.444382724q-006
			  mn(	149)=	3.030136079q-007
			  mn(	150)=-0.2593088346q-06
			  mn(	151)=	2.973159615q-006
			  mn(	152)=	4.745381979q-007
			  mn(	153)=-0.2620766273q-05
			  mn(	154)=	6.20024768q-006
			  mn(	155)=-0.3617900137q-05
			  mn(	156)=	4.387847274q-006
			  mn(	157)=-0.4930832898q-06
			  mn(	158)=-0.2698844084q-06
			  mn(	159)=	3.548302401q-006
			  mn(	160)=	2.532251367q-007
			  mn(	161)=-0.2533038688q-05
			  mn(	162)=	1.169675679q-005
			  mn(	163)=-0.4229311773q-04
			  mn(	164)=	7.109217955q-005
			  mn(	165)=-0.5917844237q-04
			  mn(	166)=	1.528557243q-005
			  mn(	167)=	2.113164455q-006
			  mn(	168)=-0.5251909602q-06
			  mn(	169)=	1.667892348q-006
			  mn(	170)=-0.3440146113q-05
			  mn(	171)=	5.485588319q-006
			  mn(	172)=-0.6532299769q-05
			  mn(	173)=	6.005809968q-006
			  mn(	174)=-0.1735594453q-05
			  mn(	175)=	3.769913715q-007
			  mn(	176)=	6.307864669q-007
			  mn(	177)=-0.2820630331q-05
			  mn(	178)=	2.831958311q-006
			  mn(	179)=	2.020014881q-006
			  mn(	180)=-0.1982616624q-05
			  mn(	181)=	2.590724496q-006
			  mn(	182)=-0.4616440967q-05
			  mn(	183)=	4.361113175q-006
			  mn(	184)=	1.044879778q-006
			  mn(	185)=-0.6978464396q-06
			  mn(	186)=-0.7717524888q-06
			  mn(	187)=	1.377384647q-006
			  mn(	188)=	2.552012673q-006
			  mn(	189)=-0.4444099871q-05
			  mn(	190)=	4.035500645q-006
			  mn(	191)=	1.037337116q-006
			  mn(	192)=-0.1607164126q-05
			  mn(	193)=	2.257709291q-006
			  mn(	194)=-0.5938416448q-05
			  mn(	195)=	9.120221636q-006
			  mn(	196)=-0.3604777080q-05
			  mn(	197)=	2.230133555q-006
			  mn(	198)=-0.6767697440q-06
			  mn(	199)=-0.7632395082q-06
			  mn(	200)=	6.406535801q-007
			  mn(	201)=	3.387399765q-006
			  mn(	202)=-0.9098115769q-05
			  mn(	203)=	9.044961048q-006
			  mn(	204)=-0.1753262229q-05
			  mn(	205)=	3.900952182q-007
			  mn(	206)=	3.474521734q-007
			  mn(	207)=	2.609460355q-006
			  mn(	208)=-0.6493106182q-05
			  mn(	209)=	6.784664862q-006
			  mn(	210)=-0.3429290194q-05
			  mn(	211)=	3.874486809q-006
			  mn(	212)=-0.1310150847q-05
			  mn(	213)=-0.1716342439q-05
			  mn(	214)=	5.082704026q-006
			  mn(	215)=-0.8436910855q-05
			  mn(	216)=	5.938629742q-006
			  mn(	217)=	1.283718812q-006
			  mn(	218)=-0.1490539180q-05
			  mn(	219)=	1.589383751q-006
			  mn(	220)=-0.1819922668q-05
			  mn(	221)=	2.622593725q-006
			  mn(	222)=	1.138230663q-007
			  mn(	223)=-0.6127826509q-05
			  mn(	224)=	7.73851902q-006
			  mn(	225)=-0.3525366195q-05
			  mn(	226)=	3.317356911q-006
			  mn(	227)=-0.1868776662q-05
			  mn(	228)=	1.089536181q-006
			  mn(	229)=-0.2546114578q-06
			  mn(	230)=	2.075693316q-006
			  mn(	231)=-0.5962883145q-05
			  mn(	232)=	5.157775461q-006
			  mn(	233)=	6.370254464q-007
			  mn(	234)=-0.9733084086q-06
			  mn(	235)=	1.967705243q-007
			  mn(	236)=	1.855634162q-006
			  mn(	237)=-0.1383811813q-05
			  mn(	238)=	6.212574018q-007
			  mn(	239)=	3.40426093q-007
			  mn(	240)=-0.1931148351q-05
			  mn(	241)=-0.2854066021q-05
			  mn(	242)=	1.038482945q-005
			  mn(	243)=	1.914624707q-007
			  mn(	244)=-0.1423563482q-04
			  mn(	245)=	1.319591991q-005
			  mn(	246)=-0.1710125306q-05
			  mn(	247)=-0.9326241537q-05
			  mn(	248)=	1.859140527q-005
			  mn(	249)=-0.1486778801q-04
			  mn(	250)=	7.450749045q-007
			  mn(	251)=	8.872967398q-006
			  mn(	252)=-0.3440007895q-05
			  mn(	253)=	1.771570592q-006
			  mn(	254)=-0.6736671156q-05
			  mn(	255)=	5.549128839q-006
			  mn(	256)=-0.2082301635q-05
			  mn(	257)=-0.1730124334q-06
			  mn(	258)=-0.1730124334q-06
c
c	k=Imag n
		    yk(	1)=	1.827
		    yk(	2)=	1.836
		    yk(	3)=	1.825
		    yk(	4)=	1.804
		    yk(	5)=	1.784
		    yk(	6)=	1.773
		    yk(	7)=	1.778
		    yk(	8)=	1.8
		    yk(	9)=	1.832
		    yk(	10)=	1.858
		    yk(	11)=	1.869
		    yk(	12)=	1.863
		    yk(	13)=	1.85
		    yk(	14)=	1.841
		    yk(	15)=	1.834
		    yk(	16)=	1.815
		    yk(	17)=	1.774
		    yk(	18)=	1.719
		    yk(	19)=	1.677
		    yk(	20)=	1.68
		    yk(	21)=	1.745
		    yk(	22)=	1.863
		    yk(	23)=	2.007
		    yk(	24)=	2.155
		    yk(	25)=	2.298
		    yk(	26)=	2.434
		    yk(	27)=	2.563
		    yk(	28)=	2.685
		    yk(	29)=	2.803
		    yk(	30)=	2.916
		    yk(	31)=	3.026
		    yk(	32)=	3.133
		    yk(	33)=	3.237
		    yk(	34)=	3.338
		    yk(	35)=	3.438
		    yk(	36)=	3.536
		    yk(	37)=	3.632
		    yk(	38)=	3.726
		    yk(	39)=	3.819
		    yk(	40)=	3.911
		    yk(	41)=	4.002
		    yk(	42)=	4.091
		    yk(	43)=	4.18
		    yk(	44)=	4.268
		    yk(	45)=	4.354
		    yk(	46)=	4.44
		    yk(	47)=	4.525
		    yk(	48)=	4.61
		    yk(	49)=	4.694
		    yk(	50)=	4.777
		    yk(	51)=	4.859
		    yk(	52)=	4.941
		    yk(	53)=	5.023
		    yk(	54)=	5.104
		    yk(	55)=	5.185
		    yk(	56)=	5.265
		    yk(	57)=	5.345
		    yk(	58)=	5.424
		    yk(	59)=	5.503
		    yk(	60)=	5.582
		    yk(	61)=	5.66
		    yk(	62)=	5.738
		    yk(	63)=	5.816
		    yk(	64)=	5.893
		    yk(	65)=	5.97
		    yk(	66)=	6.047
		    yk(	67)=	6.124
		    yk(	68)=	6.2
		    yk(	69)=	6.276
		    yk(	70)=	6.352
		    yk(	71)=	6.428
		    yk(	72)=	6.504
		    yk(	73)=	6.579
		    yk(	74)=	6.654
		    yk(	75)=	6.729
		    yk(	76)=	6.804
		    yk(	77)=	6.878
		    yk(	78)=	6.953
		    yk(	79)=	7.027
		    yk(	80)=	7.101
		    yk(	81)=	7.175
		    yk(	82)=	7.249
		    yk(	83)=	7.323
		    yk(	84)=	7.396
		    yk(	85)=	7.47
		    yk(	86)=	7.543
		    yk(	87)=	7.616
		    yk(	88)=	7.689
		    yk(	89)=	7.762
		    yk(	90)=	7.835
		    yk(	91)=	7.908
		    yk(	92)=	7.981
		    yk(	93)=	8.053
		    yk(	94)=	8.126
		    yk(	95)=	8.198
		    yk(	96)=	8.27
		    yk(	97)=	8.342
		    yk(	98)=	8.414
		    yk(	99)=	8.486
		    yk(	100)=	8.558
		    yk(	101)=	8.63
		    yk(	102)=	8.702
		    yk(	103)=	8.773
		    yk(	104)=	8.845
		    yk(	105)=	8.916
		    yk(	106)=	8.988
		    yk(	107)=	9.556
		    yk(	108)=	9.627
		    yk(	109)=	9.698
		    yk(	110)=	9.762
		    yk(	111)=	9.825
		    yk(	112)=	9.889
		    yk(	113)=	9.952
		    yk(	114)=	10.01
		    yk(	115)=	10.08
		    yk(	116)=	10.14
		    yk(	117)=	10.2
		    yk(	118)=	10.27
		    yk(	119)=	10.33
		    yk(	120)=	10.39
			  yk(	121)=	10.45
			  yk(	122)=	10.51
			  yk(	123)=	10.57
			  yk(	124)=	10.64
			  yk(	125)=	10.7
			  yk(	126)=	10.76
			  yk(	127)=	10.82
			  yk(	128)=	10.88
			  yk(	129)=	10.94
			  yk(	130)=	11
			  yk(	131)=	11.06
			  yk(	132)=	11.12
			  yk(	133)=	11.18
			  yk(	134)=	11.24
			  yk(	135)=	11.3
			  yk(	136)=	11.36
			  yk(	137)=	11.42
			  yk(	138)=	11.48
			  yk(	139)=	11.53
			  yk(	140)=	11.59
			  yk(	141)=	11.65
			  yk(	142)=	11.71
			  yk(	143)=	11.77
			  yk(	144)=	11.83
			  yk(	145)=	11.88
			  yk(	146)=	11.94
			  yk(	147)=	12
			  yk(	148)=	12.06
			  yk(	149)=	12.11
			  yk(	150)=	12.17
			  yk(	151)=	12.23
			  yk(	152)=	12.29
			  yk(	153)=	12.34
			  yk(	154)=	12.4
			  yk(	155)=	12.46
			  yk(	156)=	12.51
			  yk(	157)=	12.57
			  yk(	158)=	12.62
			  yk(	159)=	12.68
			  yk(	160)=	12.74
			  yk(	161)=	12.79
			  yk(	162)=	12.85
			  yk(	163)=	12.9
			  yk(	164)=	12.95
			  yk(	165)=	13
			  yk(	166)=	13.08
			  yk(	167)=	13.13
			  yk(	168)=	13.21
			  yk(	169)=	13.25
			  yk(	170)=	13.34
			  yk(	171)=	13.38
			  yk(	172)=	13.47
			  yk(	173)=	13.51
			  yk(	174)=	13.6
			  yk(	175)=	13.64
			  yk(	176)=	13.73
			  yk(	177)=	13.78
			  yk(	178)=	13.87
			  yk(	179)=	13.92
			  yk(	180)=	13.96
			  yk(	181)=	14.06
			  yk(	182)=	14.1
			  yk(	183)=	14.2
			  yk(	184)=	14.25
			  yk(	185)=	14.3
			  yk(	186)=	14.4
			  yk(	187)=	14.45
			  yk(	188)=	14.5
			  yk(	189)=	14.55
			  yk(	190)=	14.65
			  yk(	191)=	14.7
			  yk(	192)=	14.75
			  yk(	193)=	14.8
			  yk(	194)=	14.91
			  yk(	195)=	14.96
			  yk(	196)=	15.02
			  yk(	197)=	15.07
			  yk(	198)=	15.18
			  yk(	199)=	15.24
			  yk(	200)=	15.29
			  yk(	201)=	15.35
			  yk(	202)=	15.4
			  yk(	203)=	15.46
			  yk(	204)=	15.52
			  yk(	205)=	15.63
			  yk(	206)=	15.69
			  yk(	207)=	15.75
			  yk(	208)=	15.81
			  yk(	209)=	15.87
			  yk(	210)=	15.93
			  yk(	211)=	15.99
			  yk(	212)=	16.05
			  yk(	213)=	16.11
			  yk(	214)=	16.17
			  yk(	215)=	16.24
			  yk(	216)=	16.3
			  yk(	217)=	16.36
			  yk(	218)=	16.43
			  yk(	219)=	16.49
			  yk(	220)=	16.56
			  yk(	221)=	16.62
			  yk(	222)=	16.69
			  yk(	223)=	16.75
			  yk(	224)=	16.82
			  yk(	225)=	16.89
			  yk(	226)=	16.96
			  yk(	227)=	17.03
			  yk(	228)=	17.09
			  yk(	229)=	17.16
			  yk(	230)=	17.23
			  yk(	231)=	17.3
			  yk(	232)=	17.38
			  yk(	233)=	17.45
			  yk(	234)=	17.52
			  yk(	235)=	17.59
			  yk(	236)=	17.67
			  yk(	237)=	17.74
			  yk(	238)=	17.82
			  yk(	239)=	17.89
			  yk(	240)=	17.97
			  yk(	241)=	18.04
			  yk(	242)=	18.12
			  yk(	243)=	18.2
			  yk(	244)=	18.28
			  yk(	245)=	18.36
			  yk(	246)=	18.44
			  yk(	247)=	18.52
			  yk(	248)=	18.6
			  yk(	249)=	18.68
			  yk(	250)=	18.76
			  yk(	251)=	18.85
			  yk(	252)=	18.93
			  yk(	253)=	19.02
			  yk(	254)=	19.1
			  yk(	255)=	19.19
			  yk(	256)=	19.27
			  yk(	257)=	19.36
			  yk(	258)=	19.45
c
			  mk(	1)=-0.2184725378q-03
			  mk(	2)=-0.2184725378q-03
			  mk(	3)=-0.1009746558q-03
			  mk(	4)=	1.745580891q-005
			  mk(	5)=	9.013217706q-005
			  mk(	6)=	0.0001708878976
			  mk(	7)=	0.0001811027744
			  mk(	8)=	0.0001272354121
			  mk(	9)=-0.8629999345q-04
			  mk(	10)=-0.1512370398q-03
			  mk(	11)=-0.1999760698q-03
			  mk(	12)=-0.7199018306q-04
			  mk(	13)=	6.786485065q-005
			  mk(	14)=	4.239827116q-005
			  mk(	15)=-0.1247740365q-03
			  mk(	16)=-0.2515828531q-03
			  mk(	17)=-0.1859123197q-03
			  mk(	18)=	0.0001305039993
			  mk(	19)=	0.0004654615829
			  mk(	20)=	0.0006949400228
			  mk(	21)=	0.000535880421
			  mk(	22)=	0.0002512260722
			  mk(	23)=	1.532711276q-005
			  mk(	24)=-0.3195247352q-04
			  mk(	25)=-0.1210770415q-03
			  mk(	26)=-0.1730528312q-04
			  mk(	27)=-0.1310171955q-03
			  mk(	28)=	4.201119387q-005
			  mk(	29)=-0.1507348555q-03
			  mk(	30)=	2.536350854q-005
			  mk(	31)=	1.439024925q-005
			  mk(	32)=-0.9593952551q-04
			  mk(	33)=-0.2329699757q-04
			  mk(	34)=	1.137385227q-006
			  mk(	35)=-0.4181745294q-04
			  mk(	36)=	5.559796599q-005
			  mk(	37)=-0.7399659220q-04
			  mk(	38)=-0.5768552671q-04
			  mk(	39)=	7.871886153q-005
			  mk(	40)=-0.7287336315q-04
			  mk(	41)=-0.5588204713q-05
			  mk(	42)=	3.306056447q-005
			  mk(	43)=-0.5828671159q-04
			  mk(	44)=-0.5901823112q-05
			  mk(	45)=	4.118325351q-005
			  mk(	46)=-0.6747119137q-04
			  mk(	47)=	3.162579566q-005
			  mk(	48)=	4.459903188q-005
			  mk(	49)=-0.1514763427q-03
			  mk(	50)=	0.0001293496349
			  mk(	51)=-0.5661381330q-04
			  mk(	52)=-0.2642564237q-04
			  mk(	53)=	4.023087914q-005
			  mk(	54)=-0.5055937273q-04
			  mk(	55)=	4.482617587q-005
			  mk(	56)=-0.2580115079q-04
			  mk(	57)=-0.5454037546q-04
			  mk(	58)=	7.494068938q-005
			  mk(	59)=-0.7026625383q-04
			  mk(	60)=	9.805835428q-005
			  mk(	61)=-0.1822182243q-03
			  mk(	62)=	0.0002276274668
			  mk(	63)=-0.2079764961q-03
			  mk(	64)=	0.0001180701164
			  mk(	65)=-0.4692161956q-04
			  mk(	66)=-0.2735382143q-04
			  mk(	67)=	5.978241831q-005
			  mk(	68)=-0.2886741570q-04
			  mk(	69)=-0.3856270031q-04
			  mk(	70)=	8.952860388q-005
			  mk(	71)=-0.3872764839q-04
			  mk(	72)=-0.2773837181q-04
			  mk(	73)=-0.3240625925q-05
			  mk(	74)=-0.4682839908q-04
			  mk(	75)=	0.0001043358955
			  mk(	76)=-0.6790318180q-04
			  mk(	77)=	1.825016957q-005
			  mk(	78)=-0.2968524077q-04
			  mk(	79)=-0.3946653633q-04
			  mk(	80)=	0.0001067973708
			  mk(	81)=-0.5082207956q-04
			  mk(	82)=	1.462548411q-005
			  mk(	83)=-0.8785005922q-04
			  mk(	84)=	0.0001986725731
			  mk(	85)=-0.2443364409q-03
			  mk(	86)=	0.0001721776797
			  mk(	87)=-0.6112958120q-04
			  mk(	88)=-0.3524929943q-05
			  mk(	89)=-0.6914745129q-06
			  mk(	90)=-0.6667146362q-04
			  mk(	91)=	0.000194064391
			  mk(	92)=-0.2628679965q-03
			  mk(	93)=	0.0002147593076
			  mk(	94)=-0.1072421937q-03
			  mk(	95)=	2.671339631q-005
			  mk(	96)=-0.1572761982q-04
			  mk(	97)=	2.135330274q-005
			  mk(	98)=-0.2733358525q-04
			  mk(	99)=	2.006873707q-005
			  mk(	100)=-0.7315279057q-05
			  mk(	101)=	2.279410522q-007
			  mk(	102)=-0.1381154859q-05
			  mk(	103)=-0.8485344029q-06
			  mk(	104)=-0.1107589334q-06
			  mk(	105)=-0.5554785554q-06
			  mk(	106)=-0.1643503994q-05
			  mk(	107)=-0.3638884690q-06
			  mk(	108)=	3.646326051q-005
			  mk(	109)=-0.1328374414q-03
			  mk(	110)=	1.954190362q-005
			  mk(	111)=	8.078567103q-006
			  mk(	112)=	2.214381249q-005
			  mk(	113)=-0.1410858014q-03
			  mk(	114)=	0.0002561579763
			  mk(	115)=-0.2115479561q-03
			  mk(	116)=	1.169992928q-005
			  mk(	117)=	0.0001860548754
			  mk(	118)=-0.2061448555q-03
			  mk(	119)=	6.632145529q-005
			  mk(	120)=-0.3492763313q-04
			  mk(	121)=	2.752903188q-005
			  mk(	122)=-0.5018932075q-04
			  mk(	123)=	0.000201387471
			  mk(	124)=-0.2024121914q-03
			  mk(	125)=	4.07096111q-005
			  mk(	126)=-0.6213046265q-05
			  mk(	127)=	1.480154733q-005
			  mk(	128)=-0.1941022025q-04
			  mk(	129)=	1.882109066q-005
			  mk(	130)=-0.2028348040q-04
			  mk(	131)=	1.894512592q-005
			  mk(	132)=-0.1780677743q-04
			  mk(	133)=	9.552070267q-006
			  mk(	134)=	1.971480014q-005
			  mk(	135)=-0.4516569525q-04
			  mk(	136)=	3.435260853q-005
			  mk(	137)=	3.827277087q-005
			  mk(	138)=-0.2279975229q-03
			  mk(	139)=	0.0002359730529
			  mk(	140)=-0.7134831561q-04
			  mk(	141)=	9.762602466q-006
			  mk(	142)=-0.8310592001q-05
			  mk(	143)=	7.330437885q-005
			  mk(	144)=-0.2325381260q-03
			  mk(	145)=	0.0002107370061
			  mk(	146)=-0.5260881132q-04
			  mk(	147)=	5.464398353q-005
			  mk(	148)=-0.2045028466q-03
			  mk(	149)=	0.0002156634871
			  mk(	150)=-0.8708167100q-04
			  mk(	151)=	9.477879859q-005
			  mk(	152)=-0.2305662506q-03
			  mk(	153)=	0.0001890601274
			  mk(	154)=	2.836059556q-005
			  mk(	155)=-0.2407659185q-03
			  mk(	156)=	0.000299378404
			  mk(	157)=-0.2876279655q-03
			  mk(	158)=	0.0002067370639
			  mk(	159)=	2.304922047q-005
			  mk(	160)=-0.2289332745q-03
			  mk(	161)=	0.0002471709519
			  mk(	162)=-0.1952186896q-03
			  mk(	163)=-0.8363428677q-04
			  mk(	164)=	0.0006025169935
			  mk(	165)=-0.4270598614q-03
			  mk(	166)=	0.0003665920567
			  mk(	167)=-0.2758913618q-03
			  mk(	168)=	5.227445494q-005
			  mk(	169)=	0.0001796300758
			  mk(	170)=-0.2294777649q-03
			  mk(	171)=	0.0002434994648
			  mk(	172)=-0.2935776163q-03
			  mk(	173)=	0.0002703563238
			  mk(	174)=-0.2265061466q-03
			  mk(	175)=	0.0001019952717
			  mk(	176)=	0.000103684503
			  mk(	177)=-0.2448010275q-03
			  mk(	178)=	0.0002748521575
			  mk(	179)=-0.4767323190q-03
			  mk(	180)=	0.0003917978795
			  mk(	181)=-0.3724579994q-03
			  mk(	182)=	0.0002369276749
			  mk(	183)=-0.6676965816q-04
			  mk(	184)=	9.718158396q-006
			  mk(	185)=-0.9210702721q-05
			  mk(	186)=-0.4492788871q-05
			  mk(	187)=-0.8518499809q-05
			  mk(	188)=	3.098439266q-006
			  mk(	189)=-0.3885215752q-04
			  mk(	190)=	2.581202969q-005
			  mk(	191)=-0.2889695366q-05
			  mk(	192)=-0.4763846950q-04
			  mk(	193)=	0.0001595863175
			  mk(	194)=-0.2458283724q-03
			  mk(	195)=	0.000416801228
			  mk(	196)=-0.3587257086q-03
			  mk(	197)=	9.477398324q-005
			  mk(	198)=	8.707481517q-005
			  mk(	199)=-0.3333149163q-03
			  mk(	200)=	0.000371584864
			  mk(	201)=-0.3575582205q-03
			  mk(	202)=	0.0002103005432
			  mk(	203)=	2.250781923q-005
			  mk(	204)=-0.1095159532q-03
			  mk(	205)=	9.828610219q-005
			  mk(	206)=-0.3847922307q-04
			  mk(	207)=	2.268860268q-005
			  mk(	208)=-0.8495666013q-04
			  mk(	209)=	7.210329775q-005
			  mk(	210)=-0.2473547495q-04
			  mk(	211)=-0.4620578240q-05
			  mk(	212)=	1.168938078q-005
			  mk(	213)=-0.7291458205q-04
			  mk(	214)=	0.0002483833527
			  mk(	215)=-0.2908480887q-03
			  mk(	216)=	6.344881669q-005
			  mk(	217)=	0.0002064579096
			  mk(	218)=-0.2881794984q-03
			  mk(	219)=	0.000292800329
			  mk(	220)=-0.3017963686q-03
			  mk(	221)=	0.0002811063902
			  mk(	222)=-0.2597568199q-03
			  mk(	223)=	0.0001444815459
			  mk(	224)=	9.067872248q-006
			  mk(	225)=-0.2214577019q-04
			  mk(	226)=	4.828257986q-005
			  mk(	227)=-0.2009741102q-03
			  mk(	228)=	0.0001879765492
			  mk(	229)=-0.4908984277q-04
			  mk(	230)=-0.2038473680q-04
			  mk(	231)=	0.0001007653021
			  mk(	232)=-0.1148190884q-03
			  mk(	233)=	3.707988179q-005
			  mk(	234)=-0.6189515432q-04
			  mk(	235)=	0.0001820710515
			  mk(	236)=-0.2227488097q-03
			  mk(	237)=	0.0002150869394
			  mk(	238)=-0.2094660488q-03
			  mk(	239)=	0.000145411185
			  mk(	240)=-0.1432291149q-03
			  mk(	241)=	0.0001468414227
			  mk(	242)=-0.4518427879q-04
			  mk(	243)=	5.769432629q-006
			  mk(	244)=-0.6422902007q-05
			  mk(	245)=-0.8016803387q-05
			  mk(	246)=	1.077815432q-005
			  mk(	247)=-0.6213939901q-04
			  mk(	248)=	6.018391543q-005
			  mk(	249)=-0.5716523743q-04
			  mk(	250)=	0.0001419106992
			  mk(	251)=-0.1729332579q-03
			  mk(	252)=	0.0001645241643
			  mk(	253)=-0.1602238964q-03
			  mk(	254)=	0.0001047501143
			  mk(	255)=-0.1040425212q-03
			  mk(	256)=	0.0001080151449
			  mk(	257)=-0.2676916600q-04
			  mk(	258)=-0.2676916600q-04
c
		  end if

C     Final calculations
		    J=2
10      IF (XX.LE.XN(J)) GO TO 20
		    J=J+1
		    GO TO 10
20      I=J
		    XN0=XN(I)-XX
		    XN1=XX-XN(I-1)
		    H=XN(I)-XN(I-1)
		    H2=H*H/6.0
		    refre=MN(I-1)*XN0*XN0*XN0/6./H+MN(I)*XN1*XN1*XN1/6./H+
     +     (YN(I-1)-MN(I-1)*H2)*XN0/H+(YN(I)-MN(I)*H2)*XN1/H

		     J=2
30      IF (XX.LE.XK(J)) GO TO 40
		    J=J+1
		    GO TO 30
40      I=J
		    XN0=XK(I)-XX
		    XN1=XX-XK(I-1)
		    H=XK(I)-XK(I-1)
		    H2=H*H/6.0
		    refim=MK(I-1)*XN0*XN0*XN0/6./H+MK(I)*XN1*XN1*XN1/6./H+
     +    (YK(I-1)-MK(I-1)*H2)*XN0/H+(YK(I)-MK(I)*H2)*XN1/H


c     Correction to the size limiting effect
c     Step1: calculation of the experimental bulk epsb=eps1b+i*eps2b

		   eps1b=refre*refre-refim*refim
		   eps2b=2.d0*refre*refim

		   epsb=eps1b+aim*eps2b
C----------------------------------------------------------
c	SIZE CORRECTION 
c	gamma_b=gam_b*10**(14) 1/s, gam_b=1.05-1.6  
c      1.3 - Scaffardi JApplPhys 2014,(1.32) Etchegoin-Le Ru(2006)gam_b=1.3
c	 gam_b=1.3
c	Olmon PRB 2012 gam_b=1.4
c
		   gam_b=1.3
c
c	Plasma wavelength
c	(h/2*pi)*omega_p=9 eV  lambda_p=140 nm
c	Olmon PRB 2012 (h/2*pi)*omega_p=8.45 eV  lambda_p=131 nm
		   wp=140.
c
c	parameters gamma_b/omega and gamma_a/omega = [gamma_b+A*v_F/Leff]/omega
c	a=1 v_F=14.1*10**(14) nm/s Fermi velocity, Leff=Rm effetive length
c	gamma_b/omega=ga_om=0.53*10**(-4)*gam_b*lambda
c	v_F/(Leff*omega)=7.48*10**(-4)*lambda/Leff

c	if KD=0 Full correction
c	 Deps=(wavel/wp)**2*[ 1.d0/(1.d0+aim*gamb_om)-1.d0/(1.d0+aim*gama_om)]
c	if KD=1 approximation a9-11
c      Using approximationm for Deps = Deps2*[B+ i]
c       Deps2=(5.d0/Rm)*(wavel/520.d0)**3.d0
c	B=(gama+gamb)/2/omega
c	if KD=2 B=0	 
c
		   gamb_om=0.53d-04*gam_b*wavel

		   gama_om=gamb_om+7.48d-04*wavel/Rm

c      Step 2: 
		  IF(KD.eq.0) then 
c      2.1 Calculation of the FULL coorection term


		   Deps=1.d0/(1.d0+aim*gamb_om)-1.d0/(1.d0+aim*gama_om)

		   Deps=(wavel/wp)**2*Deps

		   go to 60
		  end if

c      Using approximationm for Deps = Deps2*[(gamma_a+gamma_b)/2/omega + i]
		   Deps2=(5.d0/Rm)*(wavel/520.d0)**3.d0
		   
		   B=(1.06*gam_b+7.48/Rm)*1.d-04*wavel

		  if(KD.eq.2) then
			  B=0.
		  end if
		   Deps=Deps2*(B+aim)

60	 continue
c      Step 3: Final correction
		   eps=epsb+Deps

		   eps1=real(eps)
		   eps2=imag(eps)
       refre=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)+eps1)/2.q0)
       refim=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)-eps1)/2.q0)
	 ref=refre+aim*refim
c---------------------------------------------------------------------
		   RETURN
		   END
c***************************************************************


c***************************************************************
c           subroutine Ag2021(wavel,Rm,ref,eps) 2021 July
c***************************************************************
C     2015: Calculations of complex ref=n'+i*n" eps=eps'+i*eps"
c     for  bilk SILVER: 
c      xx=wavel
c	if(xx.gt.600) Re(n)=0.0755 - 0.000155*xx + 2.15E-7*xx*xx
c     for wavel less then 600 nm
c     Re(n) - Irani+Schulz data was used;for 1000-1600 nm Jiang OE 2015
c     data was used as gide for eye to plot Re(n) points
c	Im(n) - JCh data was smoothed to plot Im(n) points.
c	Data from different sorces are in very good agreement.
c	Thus, for Im(n) splines Ag-2015 and Ag-JCH1-A,B are identical
c	
c     To ccount for size-limiting effect for silver NPs:
c     Rm    - is the radius of the gold particle in nanometers
c     refre and refim - are ouput parameters: real(n) and imag(n)
c     Correction to the size limiting effect:
c     eps1=eps1(bulk)+deps1,    eps2=eps2(bulk)+deps2
c        deps1=(WL/WLP)**2*(1.q0/(1.q0+(GB*WL/WLP)**2)-
c     +        1.q0/(1.q0+(GR*WL/WLP)**2))
c        deps2=(WL/WLP)**2*((GR*WL/WLP)/(1.q0+(GR*WL/WLP)**2)-
c     +        (GB*WL/WLP)/(1.q0+(GB*WL/WLP)**2))
c      WLP=136.1D0 - lambda_p  plasma wavelength, nm
c      GB=0.0019D0 - gamma(bulk)/omega_p
c      GR=GB+0.0047*136.1/(2.0*PI*Rm)
c      Gr=gamma_R/omega_p - size-corrected parameter
c      WL/WLP= lambda/lambda_p
c
c	2015 modification
c      Calculation of the coorection term
c      by using deps2 for IMAGINARY part of epsilon only!!!!
c	PARAMETERS:
c	hr=h/2/pi Reduced Plank constant=0.6582*10**(-15) eV*s
c	[hr]*omega_p=9 eV Average value for silver (vary from 9.5 to 8.7 eV
c	Lambda_p=138 nm (1240/9eV)
c	gamma_bulk~ 0.5-0.9 eV Here it is NOT USED
c	A=1 (Note that A=2-3 according to Quinten BOOK
c	v_F/c=0.0047 for Silver (Granquist)
c	hr*c/hr*omega_p=22 nm
c	deps2=(A/Rm)*(wavel/400)**3*(400/lambda_p)**3*(v_F/c)*(hr*c))/(hr*omega_p)
c	=2.5*(A/Rm)*(lambda/400)**3
c
c     The corrected refre and refim are calculated by:
c     refre=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)+eps1)/2.q0)
c     refim=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)-eps1)/2.q0)
c------------------------------------------------------------------
      subroutine Ag2021(wavel,Rm,ref,eps)
      IMPLICIT real*16(A-H,O-Z)
      real*16 MN,MK
      DIMENSION XN(100),YN(100),MN(100),XK(100),YK(100),MK(100)
	complex*32 ref,eps,im
c     
      xx=wavel
	if(xx.gt.600) go to 1
c	refre=0.0755 - 0.000155*xx + 2.15E-7*xx*xx
c	go to 2

c     Spline G1-4 parameters
C
	 XN(1)=  299.5
	 XN(2)=  302.9
	 XN(3)=  305.2
	 XN(4)=  308.6
	 XN(5)=  312.8
	 XN(6)=  315.3
	 XN(7)=  318.7
	 XN(8)=  322.7
	 XN(9)=  325.7
	 XN(10)=  328
	 XN(11)=  330
	 XN(12)=  333
	 XN(13)=  336.8
	 XN(14)=  341.5
	 XN(15)=  345.9
	 XN(16)=  349.4
	 XN(17)=  355
	 XN(18)=  360
	 XN(19)=  370
	 XN(20)=  380
	 XN(21)=  400
	 XN(22)=  420
	 XN(23)=  440
	 XN(24)=  460
	 XN(25)=  480
	 XN(26)=  500
	 XN(27)=  520
	 XN(28)=  540
	 XN(29)=  560
	 XN(30)=  580
	 XN(31)=  600
	 XN(32)=  625
	 XN(33)=  650
	 XN(34)=  675
	 XN(35)=  700
	 XN(36)=  750
	 XN(37)=  800
	 XN(38)=  850
	 XN(39)=  900
	 XN(40)=  950
	 XN(41)=  1000
	 XN(42)=  1050
	 XN(43)=  1100
	 XN(44)=  1150
	 XN(45)=  1200
	 XN(46)=  1250
	 XN(47)=  1300
	 XN(48)=  1350
	 XN(49)=  1400
	 XN(50)=  1450
	 XN(51)=  1500
	 XN(52)=  1550
	 XN(53)=  1600
C
	 YN(1)=  1.564
	 YN(2)=  1.556
	 YN(3)=  1.543
	 YN(4)=  1.47
	 YN(5)=  1.343
	 YN(6)=  1.248
	 YN(7)=  1.083
	 YN(8)=  0.877
	 YN(9)=  0.622
	 YN(10)=  0.389
	 YN(11)=  0.282
	 YN(12)=  0.207
	 YN(13)=  0.174
	 YN(14)=  0.161
	 YN(15)=  0.136
	 YN(16)=  0.119
	 YN(17)=  0.107
	 YN(18)=  0.098
	 YN(19)=  0.088
	 YN(20)=  0.081
	 YN(21)=  0.071
	 YN(22)=  0.0625
	 YN(23)=  0.056
	 YN(24)=  0.0528
	 YN(25)=  0.0508
	 YN(26)=  0.0505
	 YN(27)=  0.0512
	 YN(28)=  0.0528
	 YN(29)=  0.055
	 YN(30)=  0.0575
	 YN(31)=  0.0598
	 YN(32)=  0.0628
	 YN(33)=  0.0659
	 YN(34)=  0.0691
	 YN(35)=  0.0725
	 YN(36)=  0.0804
	 YN(37)=  0.089
	 YN(38)=  0.099
	 YN(39)=  0.11
	 YN(40)=  0.1225
	 YN(41)=  0.1355
	 YN(42)=  0.15
	 YN(43)=  0.165
	 YN(44)=  0.182
	 YN(45)=  0.199
	 YN(46)=  0.218
	 YN(47)=  0.237
	 YN(48)=  0.258
	 YN(49)=  0.28
	 YN(50)=  0.303
	 YN(51)=  0.3265
	 YN(52)=  0.352
	 YN(53)=  0.378
c
	 MN(1)=  -9.877E-005
	 MN(2)=  -9.877E-005
	 MN(3)=  -0.007971
	 MN(4)=  -0.001121
	 MN(5)=  -0.002014
	 MN(6)=  -0.00595
	 MN(7)=  0.00355
	 MN(8)=  -0.01253
	 MN(9)=  -0.01325
	 MN(10)=  0.03488
	 MN(11)=  0.008653
	 MN(12)=  0.0049
	 MN(13)=  0.001394
	 MN(14)=  -0.001449
	 MN(15)=  0.0005294
	 MN(16)=  0.0008458
	 MN(17)=  -0.0001715
	 MN(18)=  0.0001913
	 MN(19)=  -8.268E-006
	 MN(20)=  2.173E-005
	 MN(21)=  -1.057E-006
	 MN(22)=  4.999E-006
	 MN(23)=  1.106E-005
	 MN(24)=  2.58E-007
	 MN(25)=  5.907E-006
	 MN(26)=  1.613E-006
	 MN(27)=  2.642E-006
	 MN(28)=  1.318E-006
	 MN(29)=  1.087E-006
	 MN(30)=  -1.166E-006
	 MN(31)=  5.777E-007
	 MN(32)=  5.315E-008
	 MN(33)=  1.697E-007
	 MN(34)=  2.281E-007
	 MN(35)=  8.379E-007
	 MN(36)=  1.233E-008
	 MN(37)=  7.928E-007
	 MN(38)=  1.765E-007
	 MN(39)=  9.014E-007
	 MN(40)=  -1.82E-007
	 MN(41)=  1.027E-006
	 MN(42)=  -3.245E-007
	 MN(43)=  1.471E-006
	 MN(44)=  -7.612E-007
	 MN(45)=  1.573E-006
	 MN(46)=  -7.324E-007
	 MN(47)=  1.356E-006
	 MN(48)=  1.074E-007
	 MN(49)=  6.143E-007
	 MN(50)=  -1.644E-007
	 MN(51)=  1.243E-006
	 MN(52)=  -8.651E-009
	 MN(53)=  -8.651E-009
1	continue
C
       XK(1)=  248.9
       XK(2)=  255.2
       XK(3)=  261.6
       XK(4)=  268.9
       XK(5)=  276.1
       XK(6)=  284.4
       XK(7)=  292.5
       XK(8)=  300.9
       XK(9)=  310.7
       XK(10)=  320.4
       XK(11)=  331.5
       XK(12)=  342.5
       XK(13)=  354.3
       XK(14)=  367.9  
       XK(15)=  381.5
       XK(16)=  397.4
       XK(17)=  413.3
       XK(18)=  430.5
       XK(19)=  450.9
       XK(20)=  471.4
       XK(21)=  496
       XK(22)=  521
       XK(23)=  548.6
       XK(24)=  582.1
       XK(25)=  616.9
       XK(26)=  659.5
       XK(27)=  704.5
       XK(28)=  756.1
       XK(29)=  821.2
       XK(30)=  892
	 XK(31)=  984.1
	 XK(32)=  1087.7
	 XK(33)=  1215.5
	 XK(34)=  1393.2
       XK(35)=1610.2
C
       YK(1)= 1.38900E+00
       YK(2)= 1.39300E+00
       YK(3)= 1.38700E+00
       YK(4)= 1.37200E+00
       YK(5)= 1.33100E+00
       YK(6)= 1.26400E+00
       YK(7)= 1.16100E+00
       YK(8)= 9.64000E-01
       YK(9)= 6.16000E-01
       YK(10)=3.92000E-01
       YK(11)=8.29000E-01
       YK(12)=1.14200E+00
       YK(13)=1.41900E+00
       YK(14)=1.65700E+00
       YK(15)=1.86400E+00
       YK(16)=2.07000E+00
       YK(17)=2.27500E+00
       YK(18)=2.46200E+00
       YK(19)=2.65700E+00
       YK(20)=2.86900E+00
       YK(21)=3.09300E+00
       YK(22)=3.32400E+00
       YK(23)=3.58600E+00
       YK(24)=3.85800E+00
       YK(25)=4.15200E+00
       YK(26)=4.48300E+00
       YK(27)=4.83800E+00
       YK(28)=5.24200E+00
       YK(29)=5.72700E+00
       YK(30)=6.31200E+00
       YK(31)=6.99200E+00
       YK(32)=7.79500E+00
       YK(33)=8.82800E+00
       YK(34)=1.01000E+01
       YK(35)=1.18500E+01

C
	 MK(1)=  -0.0003008081
	 MK(2)=  -0.0003008081
	 MK(3)=  1.580026E-005
	 MK(4)=  -0.0007139049
	 MK(5)=  -0.0001736024
	 MK(6)=  -0.0004512451
	 MK(7)=  -0.001434666
	 MK(8)=  -0.001597511
	 MK(9)=  -0.0002189577
	 MK(10)=  0.01017516
	 MK(11)=  -0.004179221
	 MK(12)=  0.0005716763
	 MK(13)=  -0.0008454766
	 MK(14)=  2.622849E-005
	 MK(15)=  -0.0002650559
	 MK(16)=  0.0001065387
	 MK(17)=  -0.0001848276
	 MK(18)=  -9.212307E-005
	 MK(19)=  0.0001091794
	 MK(20)=  -0.0001149184
	 MK(21)=  2.897708E-005
	 MK(22)=  3.033301E-005
	 MK(23)=  -8.691511E-005
	 MK(24)=  4.607977E-005
	 MK(25)=  -4.050878E-005
	 MK(26)=  1.402204E-005
	 MK(27)=  -3.861147E-007
	 MK(28)=  -1.769247E-005
	 MK(29)=  2.877102E-005
	 MK(30)=  -2.531625E-005
	 MK(31)=  1.014666E-005
	 MK(32)=  5.466382E-006
	 MK(33)=  -1.243508E-005
	 MK(34)=  7.599362E-006
	 MK(35)=  7.599362E-006 

c-------------------------------------------------
C     Basic calculations
	if(xx.gt.600.) then
		refre=0.0755 - 0.000155*xx + 2.15E-7*xx*xx
		go to 29
	end if

        J=2
10      IF (XX.LE.XN(J)) GO TO 20
        J=J+1
        GO TO 10
20      I=J
        XN0=XN(I)-XX
        XN1=XX-XN(I-1)
        H=XN(I)-XN(I-1)
        H2=H*H/6.0
        refre=MN(I-1)*XN0*XN0*XN0/6./H+MN(I)*XN1*XN1*XN1/6./H+
     +     (YN(I-1)-MN(I-1)*H2)*XN0/H+(YN(I)-MN(I)*H2)*XN1/H

29      J=2
30      IF (XX.LE.XK(J)) GO TO 40
        J=J+1
        GO TO 30
40      I=J
        XN0=XK(I)-XX
        XN1=XX-XK(I-1)
        H=XK(I)-XK(I-1)
        H2=H*H/6.0
        refim=MK(I-1)*XN0*XN0*XN0/6./H+MK(I)*XN1*XN1*XN1/6./H+
     +    (YK(I-1)-MK(I-1)*H2)*XN0/H+(YK(I)-MK(I)*H2)*XN1/H
c----------------------------------------------------------------
c     Correction to the size limiting effect
c     Step1: calculation of the bulk parameter eps=eps1+ieps2
       eps1=refre*refre-refim*refim
       eps2=2.q0*refre*refim
C==================================================================
c      Step 2: Calculation of coorection terms
C    Ïîïðàâêà è ðåàëüíîé è ìíèìîé ÷àñòè EPSILON=EPS1+iEPS2
c      WL=WAVEL
c      WLP=136.1D0
c      GB=0.0019D0
c      pi=4.q0*datan(1.q0)
c      GR=GB+0.0047*136.1/(2.0*pi*Rm)
c        deps1=(WL/WLP)**2*(1.q0/(1.q0+(GB*WL/WLP)**2)-
c     +        1.q0/(1.q0+(GR*WL/WLP)**2))
c        deps2=(WL/WLP)**2*((GR*WL/WLP)/(1.q0+(GR*WL/WLP)**2)-
c     +        (GB*WL/WLP)/(1.q0+(GB*WL/WLP)**2))
c*********************************************
C        DEPS1=(OMPB/OM)**2*(1.q0/(1.q0+(GAMB/OM)**2)-
C     +        1.q0/(1.q0+(GAMJ/OM)**2))
C        DEPS2=(OMPB/OM)**2*((GAMJ/OM)/(1.q0+(GAMJ/OM)**2)-
C     +        (GAMB/OM)/(1.q0+(GAMB/OM)**2))

c	2015 modification
c      Step 2: Calculation of the coorection term
c      by using deps2 for IMAGINARY part of epsilon only!!!!
c	PARAMETERS:
c	hr=h/2/pi Reduced Plank constant=0.6582*10**(-15) eV*s
c	[hr]*omega_p=9 eV Average value for silver (vary from 9.5 to 8.7 eV
c	Lambda_p=138 nm (1240/9eV)
c	gamma_bulk~ 0.5-0.9 eV Here it is NOT USED
c	A=1 (Note that A=2-3 according to Quinten BOOK
c	v_F/c=0.0047 for Silver (Granquist)
c	hr*c/hr*omega_p=22 nm
c	deps2=(A/Rm)*(wavel/400)**3*(400/lambda_p)**3*(v_F/c)*(hr*c))/(hr*omega_p)
c	=2.5*(A/Rm)*(lambda/400)**3
c	 
       deps1=0.q0
	 As=1.q0
	 deps2=2.5d0*(As/Rm)*(wavel/400.q0)**3.q0 
       
       eps1=eps1  + deps1
       eps2=eps2  + deps2
c---------------------------------------------------------------------
c      Step 3: Recalculation of tefre and refim: Bohren-Huffman, p.281
       refre=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)+eps1)/2.q0)
       refim=qsqrt((qsqrt(eps1**2.q0+eps2**2.q0)-eps1)/2.q0)

	im=(0.q0,1.q0)
	eps=eps1+im*eps2
	ref=refre+im*refim
c---------------------------------------------------------------------
       RETURN
       END
c***************************************

C--------------------------------------------------------
c     Calculate size-correxted ref and eps
       subroutine g34JCh3(wavel,Rm,ref,KD)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MN,MK
      DIMENSION XN(100),YN(100),MN(100),XK(100),YK(100),MK(100)
      real*16 wavel,Rm
      complex*32 ref
	complex*16 eps,epsb,Deps,aim
C--------------------------------------------------------
c	Rm=aLeff=Rev/As !!!!!
	aim=(0.d0,1.d0)

c     
      xx=wavel
c     Spline G34JCh parameters
C
	XN(	1	)=	230.4
	XN(	2	)=	239.4
	XN(	3	)=	251.3
	XN(	4	)=	255.4
	XN(	5	)=	259.7
	XN(	6	)=	264.2
	XN(	7	)=	268.5
	XN(	8	)=	273.7
	XN(	9	)=	280.1
	XN(	10	)=	289.9
	XN(	11	)=	295.2
	XN(	12	)=	303.5
	XN(	13	)=	312.5
	XN(	14	)=	320.2
	XN(	15	)=	332.7
	XN(	16	)=	341.4
	XN(	17	)=	353.8
	XN(	18	)=	358.4
	XN(	19	)=	363.8
	XN(	20	)=	374.2
	XN(	21	)=	384
	XN(	22	)=	391
	XN(	23	)=	403.8
	XN(	24	)=	414.1
	XN(	25	)=	424.7
	XN(	26	)=	430
	XN(	27	)=	440
	XN(	28	)=	450.5
	XN(	29	)=	460.4
	XN(	30	)=	469.9
	XN(	31	)=	480
	XN(	32	)=	500
	XN(	33	)=	520
	XN(	34	)=	540
	XN(	35	)=	560
	XN(	36	)=	580
	XN(	37	)=	600
	XN(	38	)=	620
	XN(	39	)=	640
	XN(	40	)=	657.4
	XN(	41	)=	680
	XN(	42	)=	704.4
	XN(	43	)=	756
	XN(	44	)=	821
	XN(	45	)=	891.9
	XN(	46	)=	984
	XN(	47	)=	1087.5
	XN(	48	)=	1215.4
	XN(	49	)=	1393
	XN(	50	)=	1610


C
	YN(	1	)=	1.293
	YN(	2	)=	1.296
	YN(	3	)=	1.313
	YN(	4	)=	1.323
	YN(	5	)=	1.34
	YN(	6	)=	1.363
	YN(	7	)=	1.39
	YN(	8	)=	1.432
	YN(	9	)=	1.487
	YN(	10	)=	1.543
	YN(	11	)=	1.573
	YN(	12	)=	1.611
	YN(	13	)=	1.646
	YN(	14	)=	1.679
	YN(	15	)=	1.708
	YN(	16	)=	1.707
	YN(	17	)=	1.661
	YN(	18	)=	1.64
	YN(	19	)=	1.625
	YN(	20	)=	1.613
	YN(	21	)=	1.615
	YN(	22	)=	1.613
	YN(	23	)=	1.609
	YN(	24	)=	1.602
	YN(	25	)=	1.59
	YN(	26	)=	1.573
	YN(	27	)=	1.545
	YN(	28	)=	1.501
	YN(	29	)=	1.442
	YN(	30	)=	1.364
	YN(	31	)=	1.254
	YN(	32	)=	0.8031
	YN(	33	)=	0.5264
	YN(	34	)=	0.3772
	YN(	35	)=	0.3054
	YN(	36	)=	0.2524
	YN(	37	)=	0.2113
	YN(	38	)=	0.175
	YN(	39	)=	0.151
	YN(	40	)=	0.14
	YN(	41	)=	0.13
	YN(	42	)=	0.13
	YN(	43	)=	0.14
	YN(	44	)=	0.16
	YN(	45	)=	0.17
	YN(	46	)=	0.22
	YN(	47	)=	0.27
	YN(	48	)=	0.35
	YN(	49	)=	0.43
	YN(	50	)=	0.56


          
C

	MN(	1	)=	0.0001082668
	MN(	2	)=	0.0001082668
	MN(	3	)=	9.003636E-005
	MN(	4	)=	0.000461763
	MN(	5	)=	0.0002232546
	MN(	6	)=	0.0002290781
	MN(	7	)=	0.00045846
	MN(	8	)=	0.0002098571
	MN(	9	)=	-0.0006486994
	MN(	10	)=	0.0002446948
	MN(	11	)=	-0.0002558822
	MN(	12	)=	-9.989061E-005
	MN(	13	)=	0.0001603902
	MN(	14	)=	-0.0002697623
	MN(	15	)=	-0.0001704703
	MN(	16	)=	-0.0004608753
	MN(	17	)=	-5.132801E-005
	MN(	18	)=	0.0005058203
	MN(	19	)=	0.0001563588
	MN(	20	)=	0.0001991589
	MN(	21	)=	-0.0001555762
	MN(	22	)=	4.812338E-005
	MN(	23	)=	-7.636107E-005
	MN(	24	)=	6.886322E-005
	MN(	25	)=	-0.0004534667
	MN(	26	)=	0.0002334526
	MN(	27	)=	-0.0002294858
	MN(	28	)=	-0.0001208013
	MN(	29	)=	-0.0003309605
	MN(	30	)=	5.59582E-005
	MN(	31	)=	-0.00149829
	MN(	32	)=	0.0009854208
	MN(	33	)=	0.0001696069
	MN(	34	)=	0.0002486506
	MN(	35	)=	-3.20903E-006
	MN(	36	)=	4.618534E-005
	MN(	37	)=	-3.032441E-006
	MN(	38	)=	3.794452E-005
	MN(	39	)=	3.575436E-005
	MN(	40	)=	-1.517985E-006
	MN(	41	)=	2.820959E-005
	MN(	42	)=	1.535927E-006
	MN(	43	)=	4.670863E-006
	MN(	44	)=	-7.463649E-006
	MN(	45	)=	1.02274E-005
	MN(	46	)=	-4.276803E-006
	MN(	47	)=	3.597714E-006
	MN(	48	)=	-2.8772E-006
	MN(	49	)=	1.394114E-006
	MN(	50	)=	1.394114E-006
      
C
	XK(	1	)=	230.4
	XK(	2	)=	239.4
	XK(	3	)=	251.3
	XK(	4	)=	255.4
	XK(	5	)=	259.7
	XK(	6	)=	264.2
	XK(	7	)=	268.5
	XK(	8	)=	273.7
	XK(	9	)=	280.1
	XK(	10	)=	289.9
	XK(	11	)=	295.2
	XK(	12	)=	303.5
	XK(	13	)=	312.5
	XK(	14	)=	320.2
	XK(	15	)=	332.7
	XK(	16	)=	341.4
	XK(	17	)=	353.8
	XK(	18	)=	358.4
	XK(	19	)=	363.8
	XK(	20	)=	374.2
	XK(	21	)=	384.
	XK(	22	)=	391.
	XK(	23	)=	403.8
	XK(	24	)=	414.1
	XK(	25	)=	424.7
	XK(	26	)=	430.
	XK(	27	)=	440.
	XK(	28	)=	450.5
	XK(	29	)=	460.4
	XK(	30	)=	469.9
	XK(	31	)=	480.
	XK(	32	)=	500.
	XK(	33	)=	520.
	XK(	34	)=	540.
	XK(	35	)=	560.
	XK(	36	)=	580.
	XK(	37	)=	600.
	XK(	38	)=	620.
	XK(	39	)=	640.
	XK(	40	)=	657.4
	XK(	41	)=	680.
	XK(	42	)=	704.4
	XK(	43	)=	756.
	XK(	44	)=	821.
	XK(	45	)=	891.9
	XK(	46	)=	984.
	XK(	47	)=	1087.5
	XK(	48	)=	1215.4
	XK(	49	)=	1393.
	XK(	50	)=	1610.

      
C
	YK(	1	)=	1.477
	YK(	2	)=	1.543
	YK(	3	)=	1.641
	YK(	4	)=	1.676
	YK(	5	)=	1.715
	YK(	6	)=	1.755
	YK(	7	)=	1.79
	YK(	8	)=	1.824
	YK(	9	)=	1.843
	YK(	10	)=	1.864
	YK(	11	)=	1.87
	YK(	12	)=	1.874
	YK(	13	)=	1.877
	YK(	14	)=	1.87
	YK(	15	)=	1.84
	YK(	16	)=	1.811
	YK(	17	)=	1.792
	YK(	18	)=	1.804
	YK(	19	)=	1.823
	YK(	20	)=	1.854
	YK(	21	)=	1.869
	YK(	22	)=	1.877
	YK(	23	)=	1.887
	YK(	24	)=	1.882
	YK(	25	)=	1.866
	YK(	26	)=	1.853
	YK(	27	)=	1.838
	YK(	28	)=	1.806
	YK(	29	)=	1.767
	YK(	30	)=	1.728
	YK(	31	)=	1.73
	YK(	32	)=	1.818
	YK(	33	)=	2.1277
	YK(	34	)=	2.452
	YK(	35	)=	2.7501
	YK(	36	)=	3.0106
	YK(	37	)=	3.2411
	YK(	38	)=	3.44
	YK(	39	)=	3.62
	YK(	40	)=	3.76
	YK(	41	)=	3.94
	YK(	42	)=	4.13
	YK(	43	)=	4.542
	YK(	44	)=	5.083
	YK(	45	)=	5.663
	YK(	46	)=	6.35
	YK(	47	)=	7.15
	YK(	48	)=	8.145
	YK(	49	)=	9.519
	YK(	50	)=	11.21


C
	MK(	1	)=	0.0001090988
	MK(	2	)=	0.0001090988
	MK(	3	)=	-1.096293E-005
	MK(	4	)=	0.0002098452
	MK(	5	)=	-6.549076E-005
	MK(	6	)=	-0.0001855193
	MK(	7	)=	-0.0002177071
	MK(	8	)=	-0.0008985389
	MK(	9	)=	8.750596E-005
	MK(	10	)=	-0.0002081646
	MK(	11	)=	-0.0001199258
	MK(	12	)=	5.594176E-005
	MK(	13	)=	-0.0002035245
	MK(	14	)=	-0.0001506855
	MK(	15	)=	-0.0001032507
	MK(	16	)=	7.601288E-005
	MK(	17	)=	0.0006852496
	MK(	18	)=	0.0001314306
	MK(	19	)=	-5.959464E-005
	MK(	20	)=	-0.0001974043
	MK(	21	)=	-1.081461E-005
	MK(	22	)=	-4.084141E-006
	MK(	23	)=	-0.000150955
	MK(	24	)=	-5.570012E-005
	MK(	25	)=	-0.0002132856
	MK(	26	)=	0.0003230994
	MK(	27	)=	-0.0003039385
	MK(	28	)=	-5.264124E-006
	MK(	29	)=	-0.000196416
	MK(	30	)=	0.0007029396
	MK(	31	)=	1.290496E-005
	MK(	32	)=	0.0008667653
	MK(	33	)=	-0.0001544633
	MK(	34	)=	-2.99168E-005
	MK(	35	)=	-0.0001188668
	MK(	36	)=	-5.86132E-005
	MK(	37	)=	-9.668359E-005
	MK(	38	)=	-2.865242E-005
	MK(	39	)=	-7.220907E-005
	MK(	40	)=	1.437771E-005
	MK(	41	)=	-1.689947E-005
	MK(	42	)=	8.082517E-006
	MK(	43	)=	7.160682E-006
	MK(	44	)=	-8.522745E-007
	MK(	45	)=	-1.536105E-005
	MK(	46	)=	8.041652E-006
	MK(	47	)=	-1.062974E-006
	MK(	48	)=	-3.134448E-007
	MK(	49	)=	3.900975E-007
	MK(	50	)=	3.900975E-007

C     Basic calculations
        J=2
10      IF (XX.LE.XN(J)) GO TO 20
        J=J+1
        GO TO 10
20      I=J
        XN0=XN(I)-XX
        XN1=XX-XN(I-1)
        H=XN(I)-XN(I-1)
        H2=H*H/6.0
        refre=MN(I-1)*XN0*XN0*XN0/6./H+MN(I)*XN1*XN1*XN1/6./H+
     +     (YN(I-1)-MN(I-1)*H2)*XN0/H+(YN(I)-MN(I)*H2)*XN1/H

         J=2
30      IF (XX.LE.XK(J)) GO TO 40
        J=J+1
        GO TO 30
40      I=J
        XN0=XK(I)-XX
        XN1=XX-XK(I-1)
        H=XK(I)-XK(I-1)
        H2=H*H/6.0
        refim=MK(I-1)*XN0*XN0*XN0/6./H+MK(I)*XN1*XN1*XN1/6./H+
     +    (YK(I-1)-MK(I-1)*H2)*XN0/H+(YK(I)-MK(I)*H2)*XN1/H

c     Correction to the size limiting effect
c     Step1: calculation of experimental bulk epsb=eps1b+i*eps2b
       eps1b=refre*refre-refim*refim
       eps2b=2.d0*refre*refim

	 epsb=eps1b+aim*eps2b
C----------------------------------------------------------
c	SIZE CORRECTION 
c	gamma_b=gam_b*10**(14) 1/s, gam_b=1.05-1.6  
c      1.3 - Scaffardi JApplPhys 2014,(1.32) Etchegoin-Le Ru(2006)gam_b=1.3
c	 gam_b=1.3
c
	 gam_b=1.3
c
c	Plasma wavelength
c	(h/2*pi)*omega_p=9 eV  lambda_p=140 nm
	 wp=140.
c
c	parameters gamma_b/omega and gamma_a/omega = [gamma_b+A*v_F/Leff]/omega
c	a=1 v_F=14.1*10**(14) nm/s Fermi velocity, Leff=Rm effetive length
c	gamma_b/omega=ga_om=0.53*10**(-4)*gam_b*lambda
c	v_F/(Leff*omega)=7.48*10**(-4)*lambda/Leff

c	if KD=0 Full correction
c	 Deps=(wavel/wp)**2*[ 1.d0/(1.d0+aim*gamb_om)-1.d0/(1.d0+aim*gama_om)]
c	if KD=1 approximation a9-11
c      Using approximationm for Deps = Deps2*[B+ i]
c       Deps2=(5.d0/Rm)*(wavel/520.d0)**3.d0
c	B=(gama+gamb)/2/omega
c	if KD=2 B=0	 
c
	 gamb_om=0.53d-04*gam_b*wavel

	 gama_om=gamb_om+7.48d-04*wavel/Rm

c      Step 2: 
	IF(KD.eq.0) then 
c      2.1 Calculation of the FULL coorection term

	 Deps=1.d0/(1.d0+aim*gamb_om)-1.d0/(1.d0+aim*gama_om)

	 Deps=(wavel/wp)**2*Deps

	 go to 60
	end if

c      Using approximationm for Deps = Deps2*[(gamma_a+gamma_b)/2/omega + i]
       Deps2=(5.d0/Rm)*(wavel/520.d0)**3.d0
	 
	 B=(1.06*gam_b+7.48/Rm)*1.d-04*wavel

	if(KD.eq.2) then
		B=0.
	end if
	 Deps=Deps2*(B+aim)

60	 continue
c      Step 3: Final correction
       eps=epsb+Deps

	eps1=real(eps)
	eps2=imag(eps)
c      Step 3: Recalculation of tefre and refim: Bohren-Huffman, p.281
       refre=dsqrt((dsqrt(eps1**2.d0+eps2**2.d0)+eps1)/2.d0)
       refim=dsqrt((dsqrt(eps1**2.d0+eps2**2.d0)-eps1)/2.d0)

	 ref=refre+aim*refim

c---------------------------------------------------------------------
       RETURN
       END
c************************************************************************


c***************************************************************
c           subroutine P3HT2020DF(wavel,ref,eps) 2021
c	Hrostea et al.Coatings 2020, MDPI New formula Khlebtsov JQSRT, 2022
c***************************************************************
       subroutine P3HT2020DF(wavel,ref,eps)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*16 n,k, wavel
      dimension B(4), C(4), f(4), G(4), Om(4)

	complex*32 ref, eps
	complex*16 im
C--------------------------------------------------------
	im=(0.d0,1.d0)
	Omm=1240./wavel
C--------------------------------------------------------
	Omg=1.694058

	f(1)=-0.0005915
	f(2)=0.32539
	f(3)=0.1377056
	f(4)=-0.11381

	Om(1)=4.587146
	Om(2)=2.07346
	Om(3)=6.89049
	Om(4)=2.962906

	G(1)=0.070219
	G(2)=0.2504227
	G(3)=1.2102
	G(4)=-0.8820766 
	
	
	n=0.784217
	k=0.

	do j=1,4

	B(j)=f(j)/G(j)*(G(j)**2-(Om(j)-Omg)**2)

	C(j)=2.*f(j)*G(j)*(Om(j)-Omg)
	
	n= n + ( B(j)*(Omm-Om(j))+C(j) )/( (Omm-Om(j))**2+G(j)**2 )
		
	if(Omm.gt.Omg) then
	
		k=k+f(j)*(Omm-Omg)**2/( (Omm-Om(j))**2+G(j)**2 )
	else
		k=0.
	end if

	end do

	eps1=n**2-k**2
	eps2=2.*n*k

	ref=n+im*k
	eps=eps1+im*eps2

c---------------------------------------------------------------------
       RETURN
       END
c***************************************************************

c***************************************************************
c           subroutine PMMA(wavel,ref,eps) 2021
c	Khlebtsov JQSRT 2022	
c***************************************************************
       subroutine PMMA(wavel,ref,eps)
      IMPLICIT REAL*8(A-H,O-Z)
	complex*32 ref, eps
	complex*16 im
	real*16 wavel
	real a(0:6)
C--------------------------------------------------------
	x=wavel/1000.
	im=(0.d0,1.d0)
	a(0)=2.1778
	a(1)=6.1209d-03
	a(2)=-1.5004d-03
	a(3)=2.3678d-02
	a(4)=-4.2137d-03
	a(5)=7.3417d-04
	a(6)=-4.5042d-05
	ep=0.
	do i=0,2
	ep=ep+a(i)*x**(2*i)
	end do
	do i=3,6
	ep=ep+a(i)*x**(-2*(i-2))
	end do

	eps=ep+im*0.
	ref=dsqrt(ep)+im*0.

	return
	end
c***************************************************************
