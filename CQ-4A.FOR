c
C                                                                  
C                 Program       CQ-4A  (NGauss max<=768)           
c****
C	A new subroutine qRTFW1 has been implemented in a previous qCQ-3A version                                                                 
c*****
C                                                                 
C  This code is the double precision version of qCQ-4A.for 
c  with double precision T-matrix calculations
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
C|==================================================================

	IMPLICIT REAL*8(A-H,O-Z)
	COMPLEX*16 AI1,AL,refm,m2,eps2,refmed, epsmed,
     +           A1, A2, A11, A12, A21, A22

	REAL*8 F(4,4),FS(4,4),FSran(4,4), FLran(4,4),ZUW(48),ZUX1(48) 

      COMPLEX*16 TE_0,TM_0,TE_m,TM_m

	REAL*8 MG(4,4),MH(4,4)

      Character f_line*120,outfile*120
c
 	common/T/ TE_0(5000),TM_0(5000),TE_m(350000),TM_m(350000)

c      COMMON /CRES/ A1,A2,OCX,SCX,OCM,SCM,OCE,SCE

c######################################
c	for q programs
c      character (9) ::today
c######################################

C----- -------------------------------------------------------- 
      write (*,*)'Enter Output File name:'
      read (*,'(a)') outfile
      OPEN (UNIT=16, FILE=outfile)
	call CPU_TIME(atime0)
c######################################
c	for q programs
C     Write program-name and date
c      call date_and_time(today)
c      write (16,*) today
c	for d programs
1     FORMAT (I2, '/', I2.2, '/', I4.4,'CQ-4A.for')
      call getdat(year,month,day)
      WRITE (16,1) day,month,year

c######################################
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

	if (NRefm.gt.3) then
		write (16,*) 'Check NRefm=', NRefm
		STOP
	end if

	if(Nshape.gt.4) then
		write(16,*)'Check Nshape, Nshap=', Nshape
          stop
	end if
C--------------------------------------------------------
c	Averadge size parameters of rods
C     aLav and dav(central diameter)
	a=aLav/2.
	arod=a
	b_av=dav/2.
c---------------------------------------------------------
c	For Nshape=1-3  b=b_av=dav/2 
c	For Nshape=4 b=b_av/(1+par_hi/2)
c      and the aspect ratio e =a/b DEPENDS on the par_hi parameter!!! 
	b=b_av/(1.d0+par_hi/2.d0)
	e=a/b
C--------------------------------------------------------------   end
c	Parameter epsiln=1/e 
	epsiln=1.d0/e
C--------------------------------------------------------------   end
	AI1=(0.d0,1.d0)
      pi=dacos(-1.d0)
	CALL qGAU_short(Norient,ZUW,ZUX1)

C--------------------------------------------------------------   beg

		if (Nshape.eq.1) then
			a0=(a*b*b)**(1./3.)
			go to 1123
		end if

		if (Nshape.eq.2) then
			a0=b*(3.d0*a/2.d0/b)**(1./3.) 
			go to 1123
		end if

		if (Nshape.eq.3) then
		    if (e.ge.1.d0) then
c           Equivolume radius a0= b*(1+3(a-b)/2b)**1/3,
			   a0=b*(1.d0+3.d0*(a-b)/2.d0/b)**(1./3.)
	        else
c           Equivolume radius a0= a[1+(3/2)*((b-a)/a)**2+(3*pi/4)*((b-a))/a)]**1/3
	           a0=a*(1.d0+1.5d0*((b-a)/a)**2.d0+
     +	          (0.75d0*pi)*((b-a)/a))**(1.d0/3.d0)
			end if
			go to 1123
		 end if
c
c	For Nshape=4 
c
		if (Nshape.eq.4) then
c	
	if (e.lt.1.) goto 1120
c
c		Case 1: dog-bone rods		
c			end-cap diameter d1=2*b1 b1=b*(1+par_hi)			
			b1=b*(1.d0+par_hi)
c			end-cap secobd size bc=par_c*b1
			bc=b1*par_c
c			rods length without end-caps: 2a=aL-2*bc		
			ac=arod-bc
c		Calcualtion of the equivolume radius:
c		Radius of a sphere for side rod curve
		p=ac/b
		zz=2.d0*p*par_hi/(p*p+par_hi*par_hi)
		xx=zz/p
		zp=zz*p
			if(zp.lt.0.001) then			
			ff=1.d0+(zp/3.d0)*(1.d0+(3./20.)*zz*zz)+zp*zp/20.d0
			else
				fi=1.d0-( dasin(zz)/zz+dsqrt((1.d0-zz)*(1.d0+zz)) )/2.d0
				ff=1.d0+(2.d0*(1.d0+xx)/xx/xx)*fi-p**2/3.d0
			end if

c           The equivolume radius a0= b[(b1/b)**3*c+(3/2)*p*f(x)]**(1/3)
c		f(x)=1 + 2*(1+x)/[x**2]*fi(z) - p**2/3
c		fi(z)=1 - (1/2)*[arcsin(z)/z+sqrt(1-z**2)]
c		x=b/R  p=a/b  z=p*x  R=[a**2+(b1-b)**2]/[2*(b1-b)]

		a0=b*((b1/b)**3.d0*par_c + 1.5d0*p*ff )**(1./3.)
		end if
c---------------------------------------------------------
c		The disk with elliptical edge
c		This shape incldudes the right OBLATE cylinders (e.g.,par_c=0.001)
c		and OBLATE s-cylinders if par_c=1.
c
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
c*******************
	aLeff=a0/As
c*******************
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
C        Calculate refmed - water refr. index
c      W2=WAVEL*WAVEL
c      refmed_r=1.32334d0+3.479d+03/W2-5.111d+07/W2/W2
c	REFMED=(1.33,0.1)

	if(NRefm.eq.0) then
		refmed=dcmplx(refm1,refm2)
	end if

c	Complex ref and eps for P3HT polymer
	if(NRefm.eq.1) then
		call P3HT2020DF(wavel,refmed,epsmed)
	end if 

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
		m2=dcmplx(crm2,cim2)
		refm=m2/refmed
		goto 1111
	end if
C--------------------------------------------------------------------
C     Calculate the absolute refractive index m2 of Au or Ag particle            
C     with account for size limiting effect                         
c	Nmetal=1 Au
c	Nmetal=0 Ag
	if (Nmetal.eq.1) then 
		 if (NAu.eq.1) then
		    call Olmon(wavel,aLeff,m2,eps2,KD,KS)
			refm=m2/refmed
		 else
  			call g34JCh3(wavel,aLeff,m2,eps2,KD)
			refm=m2/refmed
		 end if
	else        
		 call Ag2021(wavel,aLeff,m2,eps2)
			refm=m2/refmed
	end if

C--------------------------------------------------------------------
1111	continue
c--------------------------------------------------------------------
	  AL0=2.d0*pi*a0/wavel	  
	  AL=AL0*REFMED
C--------------------------------------------------------------   end
c	Calculation of ALL T0 and Tm matrices for given N and Mmax
c	The particle parameters are specified by Al, EPSILN, Nshape, and Refm
c	The otuput matrices are saved in TE_0,TM_0, TE_m and TM_m one-dimensional arrays
c
        call qT(AL,EPSILN,Nshape,par_c,par_hi,REFM,Mmax,N,NGauss)

c	T-matrix has been calculated
C--------------------------------------------------------------------
c	T-matrix has been calculated
c     PRINT first T-matrix elements  
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
c     and for TE,TM polarizations with respect\
c      to the plane (symmetry axis, k). Theta is the angle between
c	the symmetry axis and the plane wave vector k.
c	Define Theta angles for given Ntheta
       i_th=0
1122         if(Ntheta.eq.0) then
                 Theta=0.d0
              else
                 Theta=i_th*90.d0/Ntheta
              end if   

	Z=dcos(theta*pi/180.d0)
c	Calculate A1(0), A2(0) light scattering amplitudes for 
c	TE [E perp. (a,k)] and TM [E parallel (a,k)] for different
c	orientation angles theta : cos (theta)=(ak)/|a||k|
c	NOTE: M_max is actual M that ensures ERR_Q accuracy of calculations
	call qQesa(AL,Mmax,N,Z,A1,A2,Qext,Qsca,Qabs,
     +           QexE,QexM,QscE,QscM,QabE,QabM,ERR_Q,M_max)

c********************************
c     TEST for the scattering amplitudes

C     Optical Theorem for the scattering amplitudes A1(0) and A2(0)
c	reads as follows:

c      Qext1=4.*imag(A1)/(AL)**2
c      Qext2=4.*imag(A2)/(AL)**2
c	write(16,876)wavel, Qext1, Qext2,A1,A2
c876	format (1x,f6.1,2(1x,e11.5),2(1x,E11.4),2(1x,E11.4))
c************************************

      if (Ntheta.eq.0) then
	    WRITE (16,991) WAVEL,Qext, Qsca, Qabs
c	    WRITE (16,991) WAVEL,Qext, Qsca, Qabs,M_max
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
c	The normalization coefficient anor=(8.d0/3.d0)/x**2
c	Itot_nor=Itot*anor  Ivv_nor=Ivv*anor

	if (NeedFij.eq.0) goto 999
C--------------------------------------------------------------------

c	write(16,192)
c192	format(2x,'Wavel',4x,'Ivh/Ivv[90]',3x,'Ivv_norm',
c     +       5x,'FS11_norm',6x,'sca_ang')

c-----------------------------------------------------------------------
c	WRITE (16,91) wavel,A0,E,REFM,Norient,N,MMAX,NGauss,NI
C|==================================================================|
C|  BETTA - the scattering angle, ZU is the cosine of the polar orienation angle
C|  of the symmetry axis,  HI is the azimuth in the incident wave refernce frame                             |
C|==================================================================|
c	start loop ober beta (scattering angle)
	DO 99 Bet=Bet_beg,Bet_end,Bet_step
		BETTA=Bet*PI/180.D0
		CBET=DCOS(BETTA)
		SBET=DSIN(BETTA)
C|==================================================================
C|  Calculatng the orientation averaged Mueller matrices FSran (i,j) and FL(i,j)
c   in (I,Q,U.V) and (Il,Ir,U,V) representations, respectively (see Newton book)	
C|  Calculations are prformed numerically by integration over all
c   random orientations of the inindent wave vector in the particle frame
c   defined by the Gauss quadratures for ZU and HI/PI.
c   The angles (ÒHETA,PHI)  define the scattering difection 
C|  in the particle symmetry frame,  ZT and SZT are cos and sin THETA,
c   GAM is the angle between the scattering planes 
c   in the particle and the incident wave frames 
C|  For details: see scheme in:
c   khlebtsov, JQSRT, 2022, https://doi.org/10.1016/j.jqsrt.2022.108069              
C|==================================================================
C|      set      FZR and SFij to be zero                            
C|==================================================================
	do i=1,4
	do j=1,4
	FSran(i,j)=0.d0
	end do
	end do
C|==================================================================
C|          Loop over cos (ZU)                                      
C|==================================================================
		IZU=1
		
60		DO 63 INZU=1,Norient
			ZU=0.5d0+0.5d0*IZU*ZUX1(INZU)

			theta_angle=180.d0*dacos(zu)/pi

			SZU=DSQRT(1.d0-ZU*ZU)
C|         loop over azimuth 
C|==================================================================
			IHI=1
61			DO 62 INHI=1,Norient
				HI=0.5d0+0.5d0*IHI*ZUX1(INHI)
				phi=pi*HI
				phi_angle=180.d0*HI

	call qAsca(Bet,theta_angle,phi_angle,Mmax, N,
     +			  	 A11, A12, A21, A22, ERR_Aij,M_max)
		

				CALL qMULLER(A11,A22,A21,A12,F)

C|==================================================================
C|  Transformation of F(i,j),  into the incident wave frame
C|  from F(i,j) to FS(i,j) matrix                                    
C|==================================================================
				CPHI=DCOS(phi)
	
				ZT=CBET*ZU+SBET*SZU*CPHI
				if(dabs(ZT).lt.1.d0) then
					SZT=DSQRT(1.d0-ZT*ZT)
				 	CGAM=(ZU*SBET-SZU*CBET*CPHI)/SZT
				else
					SZT=0.d0
					CGAM=1.d0
				end if

				CALL qSMG(CGAM,MG)
				CALL qSM(phi,MH)

c--------------------------------------------------------------------
				DO 843 ii=1,4
				DO 843 jj=1,4

				FS(ii,jj)=0.d0

					DO 841 K=1,4
					DO 841 L=1,4
841					FS(ii,jj)=FS(ii,jj)+MG(ii,K)*F(K,L)*MH(L,jj)
843				CONTINUE

c	Calcualte the random-averaged values

	do 842 i1=1,4
	do 842 j1=1,4
	FSran(i1,j1)=FSran(i1,j1)+FS(i1,j1)*ZUW(INHI)*ZUW(INZU)/4.d0
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

C  Recalculate the Mueller matrix from (I,Q,U,V) representation    
C  into (Il,Ir,U,V) representation. For i,j=3,4 no recalculation
C  of  FLran(i,j) is needed; see R. Newton book.
c
C	Initial FSran(i,j) are in the  (I,Q,U,V) representation
C
        FLran(1,1) = (FSran(1,1)+FSran(2,2)+FSran(1,2)+FSran(2,1))/2.d0
        FLran(1,2) = (FSran(1,1)-FSran(1,2)+FSran(2,1)-FSran(2,2))/2.d0
        FLran(2,1) = (FSran(1,1)+FSran(1,2)-FSran(2,1)-FSran(2,2))/2.d0
        FLran(2,2) = (FSran(1,1)-FSran(1,2)-FSran(2,1)+FSran(2,2))/2.d0
c
c	For incident VV polarization, the Stockes vector = (1,-1,0,0)
c	The total scattered intensity Isca=F11(S)*1+F12(S)*(-1)
	  Fvv_tot=FSran(1,1)-FSran(1,2)

c	Phase function:
c	phaseF=(4*pi/Csca)*FSran(1,1)(theta)
C===================================================================
c	Normalization to scattering intensity at wavel=wavelbeg and betsca=90

c	if(wavel.eq.wavelbeg) then
c        if (BET.EQ.90.d0) then
c	      F22_norm=FZR(2,2)
c	      Ftot_norm=Ftot
c	   end if
c	end if

c----------------------------------------------------------------

	  Dvh=FLran(1,2)/FLran(2,2)

c----------------------------------------------------------------
c	Normalization
		x=2.d0*pi*refmed*a0/wavel
		anor=(8.d0/3.d0)/x**2.
c		scaI_tot=anor*Fvv_tot
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
C  Calculating the Mueller matrix FSij for (I,Q,U,V) Stockes vector
c  For details: see Bohren and Huffman book
C-------------------------------------------------------------------
c
	SUBROUTINE qMULLER(A1,A2,A3,A4,F)
	IMPLICIT REAL*8(A-H,O-Z)
	COMPLEX*16 A1,A2,A3,A4,DA1,DA2,DA3,DA4
	REAL*8 F(4,4)
	DA1=dCONJG(A1)
	DA2=dCONJG(A2)
	DA3=dCONJG(A3)
	DA4=dCONJG(A4)
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
C	Calculation of ÌG(I,J) matrix that trahcforms Mueller FSij 
c	from particle coordinate system to the coordinate system related 	
c	to the scattering plane
C-------------------------------------------------------------------
	SUBROUTINE qSMG(CGAM,M)
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 M(4,4)
	CG2=CGAM*CGAM
	SG2=1.d0-CG2
	SG=DSQRT(dabs(SG2))
	S2G=2.d0*SG*CGAM
	C2G=CG2-SG2
	M(1,1)=1.d0
	M(1,2)=0.d0
	M(1,3)=0.d0
	M(1,4)=0.d0
	M(2,1)=0.d0
	M(2,2)= C2G
	M(2,3)=-S2G
	M(2,4)=0.d0
	M(3,1)=0.d0
	M(3,2)= S2G
	M(3,3)= C2G
	M(3,4)=0.d0
	M(4,1)=0.d0
	M(4,2)=0.d0
	M(4,3)=0.d0
	M(4,4)=1.d0
	RETURN
	END
C*******************************************************************
C-------------------------------------------------------------------
C                                                qSM                
C-------------------------------------------------------------------
C	Transformation matrix Ì(I,J) from ccordinate system of incident wave
c	to the particle coordinate system
C-------------------------------------------------------------------

	SUBROUTINE qSM(ALPHA,M)
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 M(4,4)
	CA=DCOS(ALPHA)
	C2A=DCOS(2.d0*ALPHA)
	SA=DSIN(ALPHA)
	S2A=DSIN(2.d0*ALPHA)
	M(1,1)=1.d0
	M(1,2)=0.d0
	M(1,3)=0.d0
	M(1,4)=0.d0
	M(2,1)=0.d0
	M(2,2)= C2A
	M(2,3)=-S2A
	M(2,4)=0.d0
	M(3,1)=0.d0
	M(3,2)= S2A
	M(3,3)= C2A
	M(3,4)=0.d0
	M(4,1)=0.d0
	M(4,2)=0.d0
	M(4,3)=0.d0
	M(4,4)=1.d0
	RETURN
	END
c******************************************************************

c******************************************************************
C
C                                                         qQSM     
C-------------------------------------------------------------------
c	Calculation of angular function Q[m,m+n-1] and S[m,m+n-1]
c	These are saved in Q(n=1,...NM), S(n=1,...NM)
c	Case m=0: 
c	CMMON /CQSM/ Q,S
C-------------------------------------------------------------------
        SUBROUTINE qQSM(M,NM,Z)
        IMPLICIT REAL*8(A-H,O-Z)
        COMMON /CQSM/ Q(100),S(100)
C-------------------------------------------------------------------
C                                                   Case m>0
C-------------------------------------------------------------------
        IF (Z.NE.1.d0) GO TO 7
        Z12=0.d0
        GO TO 8
7       Z12=DSQRT(DABS((1.d0-Z)*(1.d0+Z)))
8       IF (M.EQ.0) GO TO 30
C-------------------------------------------------------------------
C                                                Calculate Q(1)
C-------------------------------------------------------------------
        U=1.0d0
        IF(M.EQ.1) GO TO 11
        DO 10 N=2,M
10      U=U*Z12*(2*N-1)
11      V=0.d0
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
        V1=0.d0
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

c********************************************************************
C                                                         qT
C-------------------------------------------------------------------
c	Subroutine qT to calculate T_m-matrces for m=0,1,...Mmax
c	The T-matrix order euqals N, 
c     AL=k*Rev, k - the complex wave number in medium,
c     epsiln = diameter/rotation length
c	refm - is the complex relative refractive index,
c	NGauss - the number of Gauss quadrature points,
c	TE_0mTM_0, TE_m, and TM_m - output T-matrix arrays
c	They are saved as one-dimensional arrays of length NLxNL
c	and sum_over_m(NmxNm), where Nm=N-m+1, m=1,2,...N
c------------------------------------------------------------ 
        SUBROUTINE qT(AL,EPSILN,
     +	  Nshape,par_c,par_hi,REFM,Mmax,N,NGauss)
      
	 IMPLICIT REAL*8 (A-H,O-Z)
	 COMPLEX*16 TE_0,TM_0,TE_m,TM_m 
	 COMPLEX*16 AL,refm,aim,SE11,SM11,SE22,SM22,SE,SM,
     +      TE(100,100),TM(100,100),
     +      TE011re(50,50), TM011re(50,50),
     +      TE022re(50,50), TM022re(50,50),
     +      TE011(50,50),   TM011(50,50), 
     +      TE022(50,50),   TM022(50,50), 
     +      TE011_1(50,50), TM011_1(50,50), 
     +      TE022_1(50,50), TM022_1(50,50), 
     +      TEre(100,100),  TMre(100,100),
     +      TEm(100,100),   TMm(100,100),
     +      TEm_1(100,100), TMm_1(100,100)

	  common /T/ TE_0(5000),TM_0(5000),TE_m(350000),TM_m(350000)

        COMMON /CGAU/ W(768),W1(768)
c------------------------------------------------------------ 
        CALL qGAU(NGauss)
C------------------------------------------------------------
c	Set all initialindexes to be equal zero for
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

	aim=(0.d0,1.d0)
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
			
			SE11=(0.d0,0.d0)
			SM11=(0.d0,0.d0)
			SE22=(0.d0,0.d0)
			SM22=(0.d0,0.d0)

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
			SE=(0.d0,0.d0)
			SM=(0.d0,0.d0)
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
C
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
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 REFM,YC,Y1,I,J,K,L,JJ,JJ1,KMJ,
     +             JMK,IML,LMI
	complex*16 HH,HH1,FK,FL,FI,FJ 

	  complex*16 AL,X, Y
        COMPLEX*16 TE(100,100),TM(100,100)
        complex*16 BES,NEY
        COMMON /CBES/ BES(1001)
        COMMON /CNEY/ NEY(1001)
        COMMON /CGAU/ W(768),W1(768)
        COMMON /CQSM/ Q(100),S(100)
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
        TE(N1,N2)=(0.d0,0.d0)
10      TM(N1,N2)=(0.d0,0.d0)
C------------------------------------------------------------
C  IV - Switch for Gauss absc.  "+" (IV=1) and "-" for      |
C  (IV=-1); KI - counter for number of subintervals   , NG -|
C  counter for Gauss absc.                                  |
C------------------------------------------------------------
        IV=1
c       KI=1
101     NG=1
11      Z=(1.0d0+IV*W1(NG))/2.0d0
C------------------------------------------------------------
c	Calculation of the angular functions Q[m,m+n-1] and S[m,m+n-1]
C------------------------------------------------------------
        CALL qQSM(M,NM,Z)
C------------------------------------------------------------
c	Calculateing the parameters Y=k1*r(theta), Y1=k2*r(theta),
c	k1=(2pi/lambda)*n1, k2=(2pi/lambda)*n2=k1*refm, n1=refmed, refm=n2/n1
c	V=RT**2*W(NG), RT(M=0)=[1/sin(theta)/r(theta)]*dr(theta)/dtheta
c	W(NG) is the Gauss quadrature weight, 
c	r(theta) is the radius-vector defining the particle shape
c	see, e.g. Msishchenko et al book, 2002 Eq.(5.270) for spheroids.
	call RTFW1(Z,AL,refm,Nshape,
     +	eps,par_c, par_hi,NG,V,Y,Y1,FW1)

C------------------------------------------------------------

C------------------------------------------------------------
C  KLUCH - switch for calc. Re T or Im T, equals 1 or 2,    
C          respectively                                     
C------------------------------------------------------------
        IF (KLUCH.EQ.2) GO TO 201
        YC=Y
        CALL cBESSEL(YC,NM,M)
        DO 200 II=1,NM+1
200     NEY(II)=BES(II)
        GO TO 202
201     CALL cNEYMAN(Y,NM,M)
202     CALL cBESSEL(Y1,NM,M)
        IF (M.NE.0) GO TO 50
C------------------------------------------------------------
C  Calculation  T-matrix for the case M = 0                 
C------------------------------------------------------------
        NL=N/2
C------------------------------------------------------------
C  KL - switch for odd/odd (1) or even/even (2), cases      
C------------------------------------------------------------
        KL=1
12      DO 14 NL1=1,NL
        N1=NL1*2-1
        IF (KL.EQ.2) N1=NL1*2
        N21=2*N1+1
        HH=N1*(N1+1)*N21*NEY(N1)/Y
        HH1=HH/N1-N21*NEY(N1+1)
        FK=-(FW1*HH*Q(N1)+HH1*S(N1))*V/N21
        FJ=NEY(N1)*S(N1)*V
        DO 14 NL2=1,NL
        N2=NL2*2-1
        IF (KL.EQ.2) N2=NL2*2
        N12=N2+1
        JJ=N12*BES(N2)/Y1
        K=FK*BES(N2)*S(N2)
        J=FJ*((JJ-BES(N2+1))*S(N2)+FW1*JJ*N2*Q(N2))
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
        FJ=NEY(N1)*V
        FK=V/N21
        FI=FJ*M
        FL=FK*M
        HH=N1M*(N1M+1)*N21*NEY(N1)/Y
        HH1=HH/N1M-N21*NEY(N1+1)
        DO 100 N2=1,NM
        N2M=N2+M-1
        QQ=Q(N1)*Q(N2)
        QS=Q(N1)*S(N2)
        SQ=S(N1)*Q(N2)
        JJ=(N2M+1)*N2M*BES(N2)/Y1
        JJ1=JJ/N2M-BES(N2+1)
        LN=LN1+(N2+1)/2-N2/2+1
        GO TO (51,52,52,51),LN
51      QQSS=M*M*QQ+S(N1)*S(N2)
        K=FK*BES(N2)*(-FW1*HH*QS-HH1*QQSS)
        J=FJ*(JJ1*QQSS+FW1*JJ*SQ)
        KMJ=(K+REFM*J)
        JMK=(J+REFM*K)
 
        GO TO 53
52      SQQS=SQ+QS
        I=FI*BES(N2)*SQQS
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

	subroutine RTFW1(Z,AL,refm,Nshape,
     +	eps,par_c, par_hi,NG,V,Y,Y1,FW1)
c
	IMPLICIT REAL*8 (A-H,O-Z)
	complex*16 AL,refm,X,Y,Y1
      COMMON /CGAU/ W(768),W1(768)
C------------------------------------------------------------
c
c	Case 1: Spheroids
c
      if(Nshape.eq.1) then
        EPS2=(EPS*EPS-1.0D0)
        RT=1.d0/dsqrt(1.0D0+EPS2*Z*Z)
        RT2=RT*RT
        FW1=EPS2*Z*RT2
C                                  IF M>0     FW1=FW1*(1-Z*Z)
        X=AL*EPS**(1.d0/3.d0)
        V=W(NG)*RT2
        Y=X*RT
        Y1=REFM*Y
        go to 4321
      end if
c
c	Case 2: Right circular cylinders
c
      if(Nshape.eq.2) then
        X=AL*(2.d0*EPS/3.d0)**(1.d0/3.d0)
        SZ=dsqrt(1.d0-Z*Z)
        
        if(SZ/Z.lt.eps) then
          RT = 1.d0/EPS/Z
          FW1= 1.d0/Z
C                                  IF M>0     FW1=FW1*(1-Z*Z)
        else
          RT = 1.d0/SZ
          FW1= -z/(SZ*SZ)
C                                  IF M>0     FW1=FW1*(1-Z*Z)
        end if
      
        RT2=RT*RT
        V=W(NG)*RT2
        Y=X*RT
        Y1=REFM*Y
        go to 4321
      end if

c	Case 3: s-cylinders      
c
      if(Nshape.eq.3) then
        if(eps.eq.1.d0) then
          V=W(NG)
          Y=AL
          Y1=Refm*Y
	    FW1=0.d0
          go to 4321
        end if   
	  if(eps.lt.1.d0) then
c
c		PROLATE s-cylinder
c
	     SZ2=(1.d0-Z)*(1.d0+Z)
		 SZ=dsqrt(SZ2)
           ab=(1.d0-eps)/eps
           ba=1.d0/ab
           ba2=ba*ba
           X=AL/(1.d0+3.d0*ab/2.d0)**(1.d0/3.d0)
        
           if(SZ/Z.lt.ba) then
			RT1=dsqrt((ba-SZ)*(ba+SZ))
		    RT=ab*(Z+RT1)
			FW1=-1.d0/RT1

C                                  IF M>0     FW1=FW1*(1-Z*Z)
		 else
			RT = 1.d0/SZ
			FW1= -z/SZ2

C                                  IF M>0     FW1=FW1*(1-Z*Z)
		 end if
           RT2=RT*RT
           V=W(NG)*RT2
           Y=X*RT
           Y1=REFM*Y
	  else

c	     OBLATE s-cylinder, i.e. disk with a circular edge
c          here a=L/2 is the rotation axis, b=d/2 is the radius of rotation
c          e=a/b<1 the aspect ratio <1, eps=epsiln=1/e=b/a >1
c		 the size parameter is X=ka  
c
	     SZ2=(1.d0-Z)*(1.d0+Z)
		 SZ=dsqrt(SZ2)
	     ab=eps-1.d0
           ba=1.d0/ab
           ab2=ab*ab
		 ba2=ba*ba
           pi=datan(1.d0)*4.d0
		 X=AL/(1.d0+3.d0*ab2/2.d0+3.d0*pi*ab/4.d0)**(1.d0/3.d0)
		        
          if(SZ/Z.gt.ab) then
            RT1=dsqrt((ba-Z)*(ba+Z))
		  RT=ab*(SZ+RT1)
            FW1=Z/SZ/RT1
C                                  IF M>0     FW1=FW1*(1-Z*Z)
          else
            RT = 1.d0/Z
            FW1= RT
C                                  IF M>0     FW1=FW1*(1-Z*Z)
          end if
          RT2=RT*RT
          V=W(NG)*RT2
          Y=X*RT
          Y1=REFM*Y
	  end if
 
	  go to 4321
	end if
c	Case 4: Dog-bone RODS with elliptic cap
c
c	Note: Nsape=4 should not be used for OBLATE particles in this version!!!
c
c	Calcualtion of the equivolume radius:
c	Radius of a sphere for side rod curve
c		R=(a**2.d0+(b1-b)**2.d0)/2.d0/(b1-b)
c		p=a/b
c		zz=p*xx=2*p*hi/[p**2+hi**2] xx=b/R=2*hi/[p**2+hi**2]=zz/p
c       zp=zz*p
c
	if (Nshape.eq.4) then
c
c	if eps=1/e>1 then goto 4320 Case of disks with elliptic side
	if (eps.gt.1.) go to 4320
c	
C	Case of dog-bone rods		
c	e=1/eps=aLi/2b aspect rato, b1b=b1/b, bc=par_c*b1, b1=b*(1+par_hi) [b=b_aver/(1+par_hi/2)]
		e =  1.d0/eps
	  b1b = 1.d0+par_hi
		p = e - par_c*b1b
		zz=2.d0*p*par_hi/(p*p+par_hi**2)
		xx=zz/p
		zp=zz*p
c		f(x)=1 + 2*(1+x)/[x**2]*fi(z) - p**2/3
c		fi(zz)=1 - (1/2)*[arcsin(zz)/zz+sqrt(1-zz**2)]     		
			if(zp.lt.0.001) then
			ff=1.d0+(zp/3.d0)*(1.d0+(3.d0/20.d0)*zz*zz)+zp*zp/20.d0
			else
			fi=1.d0-( dasin(zz)/zz+dsqrt((1.d0-zz)*(1.d0+zz)))/2.d0	
			ff=1.d0+(2.d0*(1.d0+xx)/xx/xx)*fi-p**2/3.d0
			end if

c		Recalulate b by using a0 and the scale multiplier
	scale_b = ( (b1b)**3*par_c + 1.5d0*p*ff )**(1.d0/3.d0)
         b = AL/scale_b
		b1 = b*b1b
		bc = b1*par_c
	  b2 = bc*bc/b1/b1
		a  = b*p
		a2 = 1.d0-bc*bc/a/a
	   X = b
       SZ2 = (1.d0-z)*(1.d0+z)
	 SZ = dsqrt(SZ2)

		if(sz/z.lt.b1/a) then
			g = z*z+sz2*b2
	       dg = 2.d0*z*(1.d0-b2)
		  RT1 = dsqrt(z*z-a2*g)
		   RT = p*(z+RT1)/g
             FW2 = dg*(1.d0/g+a2/2.d0/RT1/(z+RT1))
	      FW1 = -1.d0/RT1 +FW2
 		   RT2 = RT*RT
		 	V = W(NG)*RT2
			Y = X*RT
	       Y1 = REFM*Y
		else
			xx1=1.d0+xx
c*******************************************
			szb=dsqrt((1.d0/xx1-z)*(1.d0/xx1+z))
c			szb=dsqrt((xx1)**(-2)-z*z)   
c*******************************************
			RT=(2.d0+xx)/xx1/(sz+szb)
			FW1=-z/sz/szb
 		    RT2 = RT*RT
		    V = W(NG)*RT2
		   	Y = X*RT
   	        Y1 = REFM*Y
		end if
c
		goto 4321
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
	
	pi=datan(1.d0)*4.d0
	sc_a=(ec2+1.5*eec2+0.75*pi*ec*(e-ec))**(1./3.)
c
	  X=AL/sc_a
c	  X==ka=AL/sc_a 'a' is the symmetry axis, L=2a is the disk thickness
	
		z1=1./dsqrt(1.+eec2)
		if (z.ge.z1) then
			RT=1./z
			FW1=RT
		else
	
			sz2=(1.d0-Z)*(1.d0+Z)
			z2=z*z
			sz=dsqrt(sz2)
			ctg2=z*z/sz2
			p=(e-ec)*sz/(1.d0+(ec2-1.d0)*z2)
			q1=(e*e-2.*e*ec)/eec2
			qq=q1*(1.+ec2*ctg2)
			det=dsqrt(1.d0-qq)
			p1=1.d0+(ec2-1.d0)*z2
			
			dpdt=(e-ec)*z/p1*(1.d0-2.d0*sz2*(1.-ec2)/p1)
			dqdt=-2.d0*q1*ec2*z/sz/sz2

			RT=p*(1.d0+det)
			FW1=dpdt/sz/p-dqdt/2.d0/sz/(1.d0-qq+det)
	    end if
c	 			
        RT2=RT*RT
        V=W(NG)*RT2
        Y=X*RT
        Y1=REFM*Y
c	end of Case Nshape=4 e=a/b<1 disck with elliptic side
	end if     
4321  continue        

	RETURN
	END


c*************************************************************
c										cBESSEL(x,Nm,m)
c------------------------------------------------------------
C     Calculate spherical Bessels j_n(z) of COMPLEX argument z=X
c	m - Tm matrix order, Nm=[N for m=0, N-m-1 for m>0] Tm matrix dimension	     
C     Output: COMMON /CBES/  array     BES(NM+1).                    
C------------------------------------------------------------
        SUBROUTINE cBESSEL(X,NM,M)
	  COMPLEX*16 x,c,u,v,w
	  COMPLEX*16 bes(1001)
        REAL*8 xrr,xii,cr,ci
        REAL*8 ax
	  COMMON /CBES/ BES
	  IF (cdabs(X).NE.0.d0) GO TO 9
        DO 5 N1=1,NM+1
5       BES(N1)=0.0d0
        RETURN
9       N=NM+M
        IF (M.EQ.0) N=N+1
        ax=cdabs(x)
	  nmax=400.0d0/(3.d0-log10(ax))

c        NMAX=80.0D0/(2.44D0-CDLOG(X)/2.3)
        U=(0.d0,0.d0)
c        V=(1.0D-36,0.0D0)
        V=(1.d-300,0.d0)
        N1=NMAX
10      N2=2*N1+1
        W=N2*V/X-U
        NL=N1-1
        IF (NL.EQ.0) GO TO 21
        IF (NL.GT.N) GO TO 20
        BES(NL)=W
20      N1=NL
        U=V
        V=W
        IF (N1.GT.0) GO TO 10

c21      xrr=dreal(X)
c	  xii=dimag(X)
c        CR=DSIN(xrr)*dcosh(xii)
c	  CI=DCOS(xrr)*dsinh(xii)
c	!!kind=16 for q accuracy
c        C=cmplx(CR,CI,kind=8)/X/W
21	  C=cdsin(x)/x/W

        DO 50 N1=1,NM+1
        N2=M+N1-1
        IF (M.EQ.0) N2=N1
50      BES(N1)=BES(N2)*C
        RETURN
        END

c********************************************************************
C                                  cNEUMAN(z, Nm,m)   
C------------------------------------------------------------
C     Calculate COMPLEX Neuman functions of a COMPLEX argument  z                     
c	m - Tm matrix order, Nm=[N for m=0, N-m-1 for m>0] Tm matrix dimension	     
C     Output: COMMON /CNEY/  array  NEY(Nm+1).                        
C------------------------------------------------------------
        SUBROUTINE cNEYMAN(z,NM,M)
c
        IMPLICIT REAL*8 (A-H,O-Z)
        complex*16 z,DC,DS,U,V,W     
	  complex*16 NEY(1001)
        COMMON /CNEY/ NEY
5       FORMAT (1H ,'ERROR !   Z = 0  ! ! !')
        IF (abs(z).NE.0.) GO TO 9
        WRITE (*,5)
	  STOP
        RETURN
9       N=NM+M
        IF (M.EQ.0) N=N+1
c	  zr=real(z)
c	  zi=imag(z)
c	  CZr=dcos(zr)*dcosh(zi)
c	  CZi=-dsin(zr)*dsinh(zi)
c	  DC=-cmplx(CZr,CZi,kind=8)/z	
       DC=-cDCOS(z)/z
c	  SZr=dsin(zr)*dcosh(zi)
c	  SZi=dcos(zr)*dsinh(zi)
c        DS=cmplx(SZr,SZi,kind=8)/z		
        DS=cDSIN(z)/z
        U=DC
        V=DC/z-DS
        NEY(1)=V
        DO 40 N1=2,N
        W=(2*N1-1)/z*V-U
        U=V
        V=W
40      NEY(N1)=W
        DO 50 N1=1,NM+1
        N2=M+N1-1
        IF (M.EQ.0) N2=N1
50      NEY(N1)=NEY(N2)
        RETURN
        END

c********************************************************************


c********************************************************************
C                                                   qGGG(N,A,B)      
C--------------------------------------------------------------------
C   Subroutine to calculate inverse matrix by Gauss elimination
C   with a partial choise of the leading element
c   Forsait, 1980, p.65
c   Input: A(N,N) matrix, Output A**(-1)=B(N,N), IPVT(permutation vector)
C--------------------------------------------------------------------
        SUBROUTINE qGGG(N,A,B)
        IMPLICIT REAL*8(A-H,O-Z)
        COMPLEX*16 A(100,100),B(100,100),T
        INTEGER IPVT(100)
C--------------------------------------------------------------------
C     Define unit matrix B=E
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
        DO 3 I=1,N
        DO 2 J=1,N
2       B(I,J)=(0.d0,0.d0)
3       B(I,I)=(1.d0,0.d0)
C--------------------------------------------------------------------
C    Triangular decomposition, the initial matrix A will be modified!!!
C--------------------------------------------------------------------
        IPVT(N)=1
        IF (N.NE.1) GO TO 14
C--------------------------------------------------------------------
C                                   Case N=1
C--------------------------------------------------------------------
        IF (A(1,1).NE.(0.d0,0.d0)) GO TO 141
        WRITE (*,1000)
        RETURN
141     B(1,1)=1.d0/A(1,1)
        RETURN
C--------------------------------------------------------------------
C                                    Case N>1
C--------------------------------------------------------------------
14      NM1=N-1
        DO 35 K=1,NM1
        KP1=K+1
        M=K
        DO 15 I=KP1,N
        IF (CDABS(A(I,K)).GT.CDABS(A(M,K))) M=I
15      CONTINUE
        IPVT(K)=M
        IF (M.NE.K) IPVT(N)=-IPVT(N)
        T=A(M,K)
        A(M,K)=A(K,K)
        A(K,K)=T
        IF (T.EQ.(0.d0,0.d0)) GO TO 35
        DO 20 I=KP1,N
        A(I,K)=-A(I,K)/T
20      CONTINUE
        DO 30 J=KP1,N
        T=A(M,J)
        A(M,J)=A(K,J)
        A(K,J)=T
        IF (T.EQ.(0.d0,0.d0)) GO TO 30
        DO 25 I=KP1,N
        A(I,J)=A(I,J)+A(I,K)*T
25      CONTINUE
30      CONTINUE
35      CONTINUE
C--------------------------------------------------------------------
C	Calculation of inverse matrix: Direct calculations
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
40      IF (A(1,1).NE.(0.d0,0.d0)) GO TO 41
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
C	The same as qGGG for Nmax=500
C--------------------------------------------------------------------
        SUBROUTINE qGGG0(N,A,B)
        IMPLICIT REAL*8(A-H,O-Z)
        COMPLEX*16 A(50,50), B(50,50), T
        INTEGER IPVT(50)

C--------------------------------------------------------------------
C     Define unit matrix B=E
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
        DO 3 I=1,N
        DO 2 J=1,N
2       B(I,J)=(0.d0,0.d0)
3       B(I,I)=(1.d0,0.d0)
C--------------------------------------------------------------------
C    Triangular decomposition, the initial matrix A will be modified!!!
C--------------------------------------------------------------------
        IPVT(N)=1
        IF (N.NE.1) GO TO 14
C--------------------------------------------------------------------
C                                   Case N=1
C--------------------------------------------------------------------
        IF (A(1,1).NE.(0.d0,0.d0)) GO TO 141
        WRITE (*,1000)
        RETURN
141     B(1,1)=1.d0/A(1,1)
        RETURN
C--------------------------------------------------------------------
C                                    Case N>1
C--------------------------------------------------------------------
14      NM1=N-1
        DO 35 K=1,NM1
        KP1=K+1
        M=K
        DO 15 I=KP1,N
        IF (CdABS(A(I,K)).GT.CdABS(A(M,K))) M=I
15      CONTINUE
        IPVT(K)=M
        IF (M.NE.K) IPVT(N)=-IPVT(N)
        T=A(M,K)
        A(M,K)=A(K,K)
        A(K,K)=T
        IF (T.EQ.(0.d0,0.d0)) GO TO 35
        DO 20 I=KP1,N
        A(I,K)=-A(I,K)/T
20      CONTINUE
        DO 30 J=KP1,N
        T=A(M,J)
        A(M,J)=A(K,J)
        A(K,J)=T
        IF (T.EQ.(0.d0,0.d0)) GO TO 30
        DO 25 I=KP1,N
        A(I,J)=A(I,J)+A(I,K)*T
25      CONTINUE
30      CONTINUE
35      CONTINUE
C--------------------------------------------------------------------
C	Calculation of inverse matrix: Direct calculations
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
40      IF (A(1,1).NE.(0.d0,0.d0)) GO TO 41
1000    FORMAT (1H ,'Exact degeneracy!!!')
        WRITE (*,1000)
        RETURN
C--------------------------------------------------------------------
41      B(1,J)=B(1,J)/A(1,1)
61      CONTINUE
        RETURN
        END
c********************************************************************

c#################################################################
	subroutine qQesa(AL,Mmax,N,Z,A1,A2,Qext,Qsca,Qabs,
     +                 QexE,QexM,QscE,QscM,QabE,QabM,ERR_Q,M_max)
c	NOTE: M_max is actual M<=Mmax that ensures ERR_Q accuracy of calculations
c
	IMPLICIT REAL*8(A-H,O-Z)
        COMPLEX*16 PE,PM,QE,QM,AI1,AIM,
     +             DS,KOCM,KOCE,KOCX,A1,A2,A1m,A2m
	complex*16 AL

	A1=(0.d0,0.d0)
	A2=(0.d0,0.d0)
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
	
	T1=A1*dconjg(A1)
	T2=A2*dconjg(A2)
	
	e1=A1m*dconjg(A1m)
	e2=A2m*dconjg(A2m)

	if(T1.ne.0.d0) err1= dabs(e1/T1)
	if(T2.ne.0.d0) err2= dabs(e2/T2)
	if(Qsca.ne.0.d0) err3= dabs(Qs_m/Qsca)
	if(Qabs.ne.0.d0) err4= dabs((Qe_m-Qs_m)/Qabs)

	err_max=dmax1(err1,err2,err3,err4)
	 
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
c	Subroutine to calculate complex amplitude scattering matrix SAij
c	by using calculated T-matrix ans three angles:
c	beta_sca, theta_angle, and phi_angle=hi-phi_sca
c	ERR_Aij - accuracy of SAij calculations, M_max - the actual number of Tm matrix used
c-------------------------------------------------------------------
	subroutine qQesa_m(AL,M,N,Z,A1,A2,OCX,SCX,OCM,SCM,OCE,SCE)
c
	IMPLICIT REAL*8(A-H,O-Z)
        COMPLEX*16 PE,PM,QE,QM,AI1,AIM,
     +             DS,KOCM,KOCE,KOCX,A1,A2
		complex*16 AL
        REAL*8 KSCM,KSCE,KSCX
        COMMON /CQSM/ Q(100),S(100)
        COMMON /CPQM/ PE(100),PM(100),QE(100),QM(100)
        COMMON /CDMN/ D(100)
        COMMON /CKKK/ KSCM,KOCM,KSCE,KOCE,KOCX,KSCX
c        COMMON /CRES/ A1,A2,OCX,SCX,OCM,SCM,OCE,SCE
        AI1=(0.D0,1.D0)
c        PI=3.14159265358D0
c        THETA=NTHE*PI/180.D0
c        Z=DCOS(THETA)

        A1=(0.D0,0.D0)
        A2=(0.D0,0.D0)
40      NM=N-M+1
        IF (M.EQ.0) NM=N
        M2=M*2
        AIM=1.D0
        
	  KOCX=(0.D0,0.D0)
        KOCE=(0.D0,0.D0)
        KOCM=(0.D0,0.D0)
c        
	  KSCX=0.D0
	  KSCE=0.D0
        KSCM=0.D0
        
	  CALL qQSM(M,NM,Z)
        CALL qDMN(M,NM)
c****   CALL qPQM(M,N)
        CALL qPQM(M,NM,N)
        IF (M.NE.0) GO TO 42
C===================================================================
C³                                                         M=0      ³
C===================================================================
        AIM=AI1
        DO 41 N1=1,N
        DS=D(N1)*S(N1)/AIM
        A1=A1-DS*PE(N1)*AI1
        A2=A2+DS*QE(N1)
41      AIM=AIM*AI1
        GO TO 45
C===================================================================
C³                                                         M>0      ³
C===================================================================
42      DO 43 N1=1,M
43      AIM=AIM*AI1
        DO 44 N1=1,NM
        DS=D(N1)/AIM*10.D0**(-M2)
        A1=A1-DS*(S(N1)*PE(N1)*AI1+M*Q(N1)*QM(N1))
        A2=A2+DS*(S(N1)*QE(N1)-M*Q(N1)*PM(N1)*AI1)
44      AIM=AIM*AI1
45      SCE=KSCE/abs(AL)**2
        SCM=KSCM/abs(AL)**2
        SCX=2.D0*KSCX/abs(AL)**2
c        OCX=-2.D0*DREAL(KOCX)/AL/AL
c        OCE=-DREAL(KOCE)/AL/AL
c        OCM=-DREAL(KOCM)/AL/AL
        OCX=-2.D0*DREAL(KOCX/AL)/dreal(AL)
        OCE=-DREAL(KOCE/AL)/dreal(AL)
        OCM=-DREAL(KOCM/AL)/dreal(AL)
 
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
c
c												qAsca
c-------------------------------------------------------------------
c	Subroutine to calculate complex amplitude scattering matrix SAij
c	by using calculated T-matrix ans three angles:
c	beta_sca, theta_angle, and phi_angle=hi-phi_sca
c	ERR_Aij - accuracy of SAij calculations, M_max - the actual number of Tm matrix used
c-------------------------------------------------------------------
	subroutine qAsca(Beta_sca, theta_angle,phi_angle,
     +                 Mmax, N, A11, A12, A21, A22, ERR_Aij,M_max)
c
	IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 TE_0,TM_0,TE_m,TM_m
c      COMPLEX*16 SA11, SA12, SA21, SA22, A11, A12, A21, A22
      COMPLEX*16 A11, A12, A21, A22, A11m, A12m, A21m, A22m
	
c      common /T/ TE_0(5000),TM_0(5000),TE_m(350000),TM_m(350000)
	COMMON /CQSM/ Q(100),S(100)
	COMMON /CPQM/ PE(100),PM(100),QE(100),QM(100)
	COMMON /CDMN/ D(100)
c-------------------------------------------------------------------
	pi=Dacos(-1.d0)
	
	BETTA=beta_sca*PI/180.d0
	
	CBET=DCOS(BETTA)
	SBET=DSIN(BETTA)

c    Define the cos(polar orientational angle)
	ZU=dcos(PI*theta_angle/180.d0)
	
	SZU=dsin(PI*theta_angle/180.d0)

c    Define an azimuthal orientational angle

	phi=pi*phi_angle/180.d0

	CPHI=DCOS(phi)
	SPHI=DSIN(phi)
	ZT=CBET*ZU+SBET*SZU*CPHI
	if(dabs(ZT).lt.1.d0) then
	   SZT=DSQRT(1.D0-ZT*ZT)
	   CF=(CBET*SZU-SBET*ZU*CPHI)/SZT
	   SF=SBET*SPHI/SZT
	else
	   SZT=0.d0
         CF=1.d0
	   SF=0.d0
	end if

c	Calculate amplitude matrix for a given M

c-------------------------------------------------------------------
c	Start: Set M=0
	M=0
C|==================================================================|
C|      Set Mueller matrices FSij and amplitudes SFij to be zero     |
C|==================================================================|
c	Set all sums of amplitude matrices over M to be zero
	A11=(0.d0,0.d0)
	A12=(0.d0,0.d0)
	A21=(0.d0,0.d0)
      A22=(0.d0,0.d0)
c
50	 CONTINUE
C==================================================================
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
	CALL qQSM(M,NM,ZU)
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
	
	T11=A11*dconjg(A11)
	T12=A12*dconjg(A12)
	T21=A21*dconjg(A21)
	T22=A22*dconjg(A22)

	e11=A11m*dconjg(A11m)
	e12=A12m*dconjg(A12m)
	e21=A21m*dconjg(A21m)
	e22=A22m*dconjg(A22m)

	if(T11.ne.0.d0) err11= dabs(e11/T11)
	if(T12.ne.0.d0) err12= dabs(e12/T12)
	if(T21.ne.0.d0) err21= dabs(e21/T21)
	if(T22.ne.0.d0) err22= dabs(e22/T22)

	err_max=dmax1(err11,err12,err21,err22)
	 
c*********************************************
c	write(16,2003) e11,T11,e12,T12,e21,T21,e22,T22
c2003	format(1h , 'eijTij=',8(2x, q11.4))
c	write(16,2000) M,err_max
c2000	format(1h , 'M=', I4,'  ERR_max=',e11.4)
c*********************************************

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
C Calculation of the amplitude matrix for a particular M
c input: M,N, ZT=cos(THETA), CF=cos(FI),SF=sin(FI) in the particle system 
C------------------------------------------------------------
	SUBROUTINE qAm(M,N,ZT,CF,SF,A11m,A12m,A21m,A22m)
	IMPLICIT REAL*8(A-H,O-Z)
	COMPLEX*16 PE,PM,QE,QM,AI1,AIM,CS,CSM,
     +             DA,A11m,A12m,A21m,A22m
	COMMON /CQSM/ Q(100),S(100)
	COMMON /CPQM/ PE(100),PM(100),QE(100),QM(100)
	COMMON /CDMN/ D(100)
c
	AI1=(0.d0,1.d0)
	PI=dacos(-1.d0)
c	CALL qDMN(M,NM)
C------------------------------------------------------------
	A11m=(0.d0,0.d0)
	A12m=(0.d0,0.d0)
	A21m=(0.d0,0.d0)
	A22m=(0.d0,0.d0)
40      NM=N-M+1
	IF (M.EQ.0) NM=N
	M2=M*2
	AIM=1.d0
	CALL qQSM(M,NM,ZT)
	IF (M.NE.0) GO TO 42
C------------------------------------------------------------
C                                                      M=0
C------------------------------------------------------------
	AIM=AI1
	DO 41 N1=1,N
	DA=D(N1)/AIM
	DS=0.d0
	DC=1.d0
	A11m=A11m+DA*S(N1)*PE(N1)*DC/AI1
	A12m=A12m+DA*S(N1)*PM(N1)*DS/AI1
	A21m=A21m+DA*S(N1)*QM(N1)*DS
	A22m=A22m+DA*S(N1)*QE(N1)*DC
41      AIM=AIM*AI1
	GO TO 45
C------------------------------------------------------------
C                                                      M>0
C------------------------------------------------------------
42      CSM=(1.d0,0.d0)
	CS=CF+AI1*SF
	DO 43 N1=1,M
	AIM=AIM*AI1
43      CSM=CSM*CS
	DO 44 N1=1,NM
	DA=D(N1)/AIM*10.d0**(-M2)
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
        IMPLICIT REAL*8(A-H,O-Z)
        COMPLEX*16 PE,PM,QE,QM,
     +             AE(100),AM(100),BE(100),BM(100),
     +             AI1,AIM,AIN,KOCM,KOCE,KOCX
        REAL*8 KSCM,KSCE,KSCX

        COMPLEX*16 TE_0,TM_0,TE_m,TM_m

        COMPLEX*16 TE(100,100),TM(100,100)
 
        common /T/ TE_0(5000),TM_0(5000),TE_m(350000),TM_m(350000)

        COMMON /CQSM/ Q(100),S(100)
        COMMON /CPQM/ PE(100),PM(100),QE(100),QM(100)
        COMMON /CDMN/ D(100)
        COMMON /CKKK/ KSCM,KOCM,KSCE,KOCE,KOCX,KSCX
c
        AI1=(0.d0,1.d0)

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
        AE(L)=-4.d0*AIM*S(L)
        BE(L)=-AI1*AE(L)
10      AIM=AIM*AI1
        DO 11 L=1,NL
        L2=2*L
        L1=L2-1
        PE(L1)=(0.d0,0.d0)
        PE(L2)=(0.d0,0.d0)
        QE(L1)=(0.d0,0.d0)
        QE(L2)=(0.d0,0.d0)
        DO 11 K=1,NL

        K2=2*K
        K1=K2-1
        KNL=K+NL
        IF (K.EQ.L) KOCX=KOCX+TE(L,K)+TM(L,K)+
     +                   TE(L,KNL)+TM(L,KNL)
        KSCX=KSCX+(TE(L,K)*DCONJG(TE(L,K))+TM(L,K)*
     +       DCONJG(TM(L,K)))*D(2*L)/D(2*K)
        KSCX=KSCX+(TE(L,KNL)*DCONJG(TE(L,KNL))+TM(L,KNL)*
     +       DCONJG(TM(L,KNL)))*D(2*L-1)/D(2*K-1)
        PE(L1)=PE(L1)+TE(L,KNL)*AE(K1)
        PE(L2)=PE(L2)+TE(L,K)*AE(K2)
        QE(L1)=QE(L1)+TM(L,KNL)*BE(K1)
11      QE(L2)=QE(L2)+TM(L,K)*BE(K2)
        DO 100 N=1,NM
        KSCE=KSCE+D(N)*PE(N)*DCONJG(PE(N))
        KSCM=KSCM+D(N)*QE(N)*DCONJG(QE(N))
        KOCE=KOCE+D(N)*DCONJG(AE(N))*PE(N)
        KOCM=KOCM+D(N)*DCONJG(BE(N))*QE(N)
100     CONTINUE
        RETURN

c--------------------------------------------------------------------
c	Case M>0
c--------------------------------------------------------------------
12      AIM=1.d0
        M2=M*2
        DO 13 L=1,M
13      AIM=AIM*AI1
        AIM=4.d0*AIM
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
        PE(N)=(0.d0,0.d0)
        PM(N)=(0.d0,0.d0)
        QE(N)=(0.d0,0.d0)
        QM(N)=(0.d0,0.d0)
        DO 18 K=1,NM
        L=L1+(K+1)/2-K/2
        IF (N.EQ.K) KOCX=KOCX+2.d0*(TE(N,K)+TM(N,K))
        KSCX=KSCX+2.d0*(TE(N,K)*DCONJG(TE(N,K))+TM(N,K)*
     +       DCONJG(TM(N,K)))*D(N)/D(K)
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
        KSCM=KSCM+D(N)*(PM(N)*DCONJG(PM(N))+
     +       QE(N)*DCONJG(QE(N)))*10.d0**(-M2)
        KSCE=KSCE+D(N)*(PE(N)*DCONJG(PE(N))+
     +       QM(N)*DCONJG(QM(N)))*10.d0**(-M2)
        KOCM=KOCM+D(N)*(PM(N)*DCONJG(AM(N))+
     +       QE(N)*DCONJG(BE(N)))*10.d0**(-M2)
        KOCE=KOCE+D(N)*(PE(N)*DCONJG(AE(N))+
     +       QM(N)*DCONJG(BM(N)))*10.d0**(-M2)
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
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /CDMN/ D(100)
	IF (M.EQ.0) GO TO 30
C------------------------------------------------------------
C                                                      M>0
C------------------------------------------------------------
	M2=M*2
	FM=1.d0
	DO 10 I=1,M2
10      FM=FM*I
	D(1)=(M2+1.d0)/M2/(M+1.d0)/FM*10.d0**(M2)
	NM1=NM-1
	DO 20 N=1,NM1
	MPN=M+N
20      D(N+1)=D(N)*N*(2.d0*MPN+1.d0)*(MPN-1.d0)/(MPN+1.d0)/
     +         (MPN+M)/(2.d0*MPN-1.d0)
	RETURN
C------------------------------------------------------------
C                                                      M=0
C------------------------------------------------------------
30      DO 40 N=1,NM
40      D(N)=(2*N+1)/4.d0/N/(N+1)
	RETURN
	END

c

c***************************************************************
c           subroutine P3HT2020DF(wavel,ref,eps) 2021
c	Hrostea et al. Coatings 2020, MDPI New Formula see Khlebtsov JQSRT, 2022.	
c***************************************************************
       subroutine P3HT2020DF(wavel,ref,eps)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 n,k
      dimension B(4), C(4), f(4), G(4), Om(4)

	complex*16 ref, eps,im
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
c	see  Khlebtsov JQSRT 2022
c***************************************************************
       subroutine PMMA(wavel,ref,eps)
      IMPLICIT REAL*8(A-H,O-Z)
	complex*16 ref, eps,im
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

C
c***************************************************************

C                                                         qGAU_short
c---------------------------------------------------------------------
c	A short version of qGAU with arrays in the parameter list
c---------------------------------------------------------------------
        SUBROUTINE qGAU_short(NG,W,X1)
        REAL*8 W(48),X1(48)
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 6                                    ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
6       IF (NG.GT.6) GO TO 12
        NG=6
C
	 w( 1)=0.4717533638651182719461596148501712d-01
	 w( 2)=0.1069393259953184309602547181939963d+00
	 w( 3)=0.1600783285433462263346525295433591d+00
	 w( 4)=0.2031674267230659217490644558097983d+00
	 w( 5)=0.2334925365383548087608498989248780d+00
	 w( 6)=0.2491470458134027850005624360429512d+00
c 
	 X1( 1)=0.9815606342467192506905490901492808d+00
	 X1( 2)=0.9041172563704748566784658661190963d+00
	 X1( 3)=0.7699026741943046870368938332128181d+00
	 X1( 4)=0.5873179542866174472967024189405343d+00
	 X1( 5)=0.3678314989981801937526915366437175d+00
	 X1( 6)=0.1252334085114689154724413694638531d+00     

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
	 w( 1)=0.1234122979998719954680566707003725d-01
	 w( 2)=0.2853138862893366318130781595187828d-01
	 w( 3)=0.4427743881741980616860274821133821d-01
	 w( 4)=0.5929858491543678074636775850010859d-01
	 w( 5)=0.7334648141108030573403361525311648d-01
	 w( 6)=0.8619016153195327591718520298374266d-01
	 w( 7)=0.9761865210411388826988066446424718d-01
	 w( 8)=0.1074442701159656347825773424466063d+00
	 w( 9)=0.1155056680537256013533444839067836d+00
	 w(10)=0.1216704729278033912044631534762624d+00
	 w(11)=0.1258374563468282961213753825111837d+00
	 w(12)=0.1279381953467521569740561652246955d+00
c
	 X1( 1)=0.9951872199970213601799974097007367d+00
	 X1( 2)=0.9747285559713094981983919930081690d+00
	 X1( 3)=0.9382745520027327585236490017087215d+00
	 X1( 4)=0.8864155270044010342131543419821967d+00
	 X1( 5)=0.8200019859739029219539498726697452d+00
	 X1( 6)=0.7401241915785543642438281030999784d+00
	 X1( 7)=0.6480936519369755692524957869107476d+00
	 X1( 8)=0.5454214713888395356583756172183724d+00
	 X1( 9)=0.4337935076260451384870842319133498d+00
	 X1(10)=0.3150426796961633743867932913198102d+00
	 X1(11)=0.1911188674736163091586398207570697d+00
	 X1(12)=0.6405689286260562608504308262474504d-01     


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
       w( 1)=0.3153346052305838632677311543891515d-02
	 w( 2)=0.7327553901276262102383979621786534d-02
	 w( 3)=0.1147723457923453948959266760909167d-01
	 w( 4)=0.1557931572294384872817695583446033d-01
	 w( 5)=0.1961616045735552781446071965221279d-01
	 w( 6)=0.2357076083932437914051930137844923d-01
	 w( 7)=0.2742650970835694820007383626250591d-01
	 w( 8)=0.3116722783279808890206575684635440d-01
	 w( 9)=0.3477722256477043889254858596380243d-01
	 w(10)=0.3824135106583070631721725652371569d-01
	 w(11)=0.4154508294346474921405882236106477d-01
	 w(12)=0.4467456085669428041944858712585040d-01
	 w(13)=0.4761665849249047482590662347892982d-01
	 w(14)=0.5035903555385447495780761908786555d-01
	 w(15)=0.5289018948519366709550505626469886d-01
	 w(16)=0.5519950369998416286820349519163536d-01
	 w(17)=0.5727729210040321570515023468470064d-01
	 w(18)=0.5911483969839563574647481743351990d-01
	 w(19)=0.6070443916589388005296923202782048d-01
	 w(20)=0.6203942315989266390419778413759853d-01
	 w(21)=0.6311419228625402565712602275023334d-01
	 w(22)=0.6392423858464818662390620182551554d-01
	 w(23)=0.6446616443595008220650419365770497d-01
	 w(24)=0.6473769681268392250302493873659151d-01
c 
	 X1( 1)=0.9987710072524261186005414915631136d+00
	 X1( 2)=0.9935301722663507575479287508490740d+00
	 X1( 3)=0.9841245837228268577445836000265988d+00
	 X1( 4)=0.9705915925462472504614119838006601d+00
	 X1( 5)=0.9529877031604308607229606660257181d+00
	 X1( 6)=0.9313866907065543331141743801016012d+00
	 X1( 7)=0.9058791367155696728220748356710117d+00
	 X1( 8)=0.8765720202742478859056935548050967d+00
	 X1( 9)=0.8435882616243935307110898445196559d+00
	 X1(10)=0.8070662040294426270825530430245384d+00
	 X1(11)=0.7671590325157403392538554375229691d+00
	 X1(12)=0.7240341309238146546744822334936653d+00
	 X1(13)=0.6778723796326639052118512806759091d+00
	 X1(14)=0.6288673967765136239951649330699947d+00
	 X1(15)=0.5772247260839727038178092385404788d+00
	 X1(16)=0.5231609747222330336782258691375086d+00
	 X1(17)=0.4669029047509584045449288616507985d+00
	 X1(18)=0.4086864819907167299162254958146334d+00
	 X1(19)=0.3487558862921607381598179372704079d+00
	 X1(20)=0.2873624873554555767358864613167976d+00
	 X1(21)=0.2247637903946890612248654401746923d+00
	 X1(22)=0.1612223560688917180564373907834977d+00
	 X1(23)=0.9700469920946269893005395585362452d-01
	 X1(24)=0.3238017096286936203332224315213445d-01     

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
	 w( 1)=0.7967920655520124294381434969436714d-03
	 w( 2)=0.1853960788946921732335925350893890d-02
	 w( 3)=0.2910731817934946408410617989400720d-02
	 w( 4)=0.3964554338444686673733415767419639d-02
	 w( 5)=0.5014202742927517692470194969030949d-02
	 w( 6)=0.6058545504235961683316742031729087d-02
	 w( 7)=0.7096470791153865269144160812143404d-02
	 w( 8)=0.8126876925698759217382427707855918d-02
	 w( 9)=0.9148671230783386632584602665207933d-02
	 w(10)=0.1016077053500841575758763695382536d-01
	 w(11)=0.1116210209983849859121326382856207d-01
	 w(12)=0.1215160467108831963518135273666725d-01
	 w(13)=0.1312822956696157263706366690259911d-01
	 w(14)=0.1409094177231486091586162472463661d-01
	 w(15)=0.1503872102699493800587627522209756d-01
	 w(16)=0.1597056290256229138061645679149306d-01
	 w(17)=0.1688547986424517245047754060686194d-01
	 w(18)=0.1778250231604526083761422648607104d-01
	 w(19)=0.1866067962741146738515675862213225d-01
	 w(20)=0.1951908114014502241008522028120468d-01
	 w(21)=0.2035679715433332459524521541727163d-01
	 w(22)=0.2117293989219129898767386719100365d-01
	 w(23)=0.2196664443874434919475638680156228d-01
	 w(24)=0.2273706965832937400134784197749033d-01
	 w(25)=0.2348339908592621984223593266761256d-01
	 w(26)=0.2420484179236469128226733787267704d-01
	 w(27)=0.2490063322248361028838218086833276d-01
	 w(28)=0.2557003600534936149879716794360011d-01
	 w(29)=0.2621234073567241391345796396446336d-01
	 w(30)=0.2682686672559176219805672871415660d-01
	 w(31)=0.2741296272602924282342108748909122d-01
	 w(32)=0.2797000761684833443981857658902251d-01
	 w(33)=0.2849741106508538564559951294580563d-01
	 w(34)=0.2899461415055523654267878127968158d-01
	 w(35)=0.2946108995816790597043633218285850d-01
	 w(36)=0.2989634413632838598438807579440005d-01
	 w(37)=0.3029991542082759379408876420650097d-01
	 w(38)=0.3067137612366914901422883035620429d-01
	 w(39)=0.3101033258631383742324977997063272d-01
	 w(40)=0.3131642559686135581278426671506315d-01
	 w(41)=0.3158933077072716855802074616996422d-01
	 w(42)=0.3182875889441100653475373988553221d-01
	 w(43)=0.3203445623199266321813897747021117d-01
	 w(44)=0.3220620479403025066866711455723255d-01
	 w(45)=0.3234382256857592842877483882894328d-01
	 w(46)=0.3244716371406426936401278844884587d-01
	 w(47)=0.3251611871386883598720549144777837d-01
	 w(48)=0.3255061449236316624196141829728569d-01
c
	 X1( 1)=0.9996895038832307668276901057843654d+00
	 X1( 2)=0.9983643758631816777241494395267193d+00
	 X1( 3)=0.9959818429872092906503990848793565d+00
	 X1( 4)=0.9925439003237626245718922977016011d+00
	 X1( 5)=0.9880541263296237994807627724184511d+00
	 X1( 6)=0.9825172635630146774470458318825480d+00
	 X1( 7)=0.9759391745851364664526010342418798d+00
	 X1( 8)=0.9683268284632642121736593664436509d+00
	 X1( 9)=0.9596882914487425393000680258847348d+00
	 X1(10)=0.9500327177844376357560989484438634d+00
	 X1(11)=0.9393703397527552169318573794316151d+00
	 X1(12)=0.9277124567223086909646904726965465d+00
	 X1(13)=0.9150714231208980742058844661504215d+00
	 X1(14)=0.9014606353158523413192326973096868d+00
	 X1(15)=0.8868945174024204160568774333937861d+00
	 X1(16)=0.8713885059092965028737747644930427d+00
	 X1(17)=0.8549590334346014554627869698935691d+00
	 X1(18)=0.8376235112281871214943028167646929d+00
	 X1(19)=0.8194003107379316755389996242244780d+00
	 X1(20)=0.8003087441391408172287961439713277d+00
	 X1(21)=0.7803690438674332176036044558893177d+00
	 X1(22)=0.7596023411766474987029704411196217d+00
	 X1(23)=0.7380306437444001328511657310964972d+00
	 X1(24)=0.7156768123489676262251441480576248d+00
	 X1(25)=0.6925645366421715613442457698182425d+00
	 X1(26)=0.6687183100439161539525572075727743d+00
	 X1(27)=0.6441634037849671067984123500633016d+00
	 X1(28)=0.6189258401254685703863692869878628d+00
	 X1(29)=0.5930323647775720806835557550437583d+00
	 X1(30)=0.5665104185613971684042501950834627d+00
	 X1(31)=0.5393881083243574362268025973255607d+00
	 X1(32)=0.5116941771546676735855097454288543d+00
	 X1(33)=0.4834579739205963597684056093638952d+00
	 X1(34)=0.4547094221677430086356761480863078d+00
	 X1(35)=0.4254789884073005453648192035699980d+00
	 X1(36)=0.3957976498289086032850002431351642d+00
	 X1(37)=0.3656968614723136350308955939954732d+00
	 X1(38)=0.3352085228926254226163256248056287d+00
	 X1(39)=0.3043649443544963530239297792952033d+00
	 X1(40)=0.2731988125910491414872721557245648d+00
	 X1(41)=0.2417431561638400123279319037406277d+00
	 X1(42)=0.2100313104605672036028471857122957d+00
	 X1(43)=0.1780968823676186027594026251764915d+00
	 X1(44)=0.1459737146548969419891073333433265d+00
	 X1(45)=0.1136958501106659209112080954143281d+00
	 X1(46)=0.8129749546442555899447126297475285d-01
	 X1(47)=0.4881298513604973111195820426127506d-01
	 X1(48)=0.1627674484960296957913456369523845d-01

56	RETURN
	END


C
CÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»
Cº                                                         qGAU      º
CÇÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¶
Cº      ‡ ¤ ­¨¥ ã§«®¢ ¨ ¢¥á®¢ ¤«ï ª¢ ¤à âãà­®© ä®à¬ã«ë ƒ ãáá        º
CÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼
        SUBROUTINE qGAU(NG)
        REAL*8 W,X1
        COMMON /CGAU/ W(768),X1(768)
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 6                                    ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
6       IF (NG.GT.6) GO TO 12
        NG=6
C
	 w( 1)=0.4717533638651182719461596148501712d-01
	 w( 2)=0.1069393259953184309602547181939963d+00
	 w( 3)=0.1600783285433462263346525295433591d+00
	 w( 4)=0.2031674267230659217490644558097983d+00
	 w( 5)=0.2334925365383548087608498989248780d+00
	 w( 6)=0.2491470458134027850005624360429512d+00
c 
	 X1( 1)=0.9815606342467192506905490901492808d+00
	 X1( 2)=0.9041172563704748566784658661190963d+00
	 X1( 3)=0.7699026741943046870368938332128181d+00
	 X1( 4)=0.5873179542866174472967024189405343d+00
	 X1( 5)=0.3678314989981801937526915366437175d+00
	 X1( 6)=0.1252334085114689154724413694638531d+00     

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
	 w( 1)=0.1234122979998719954680566707003725d-01
	 w( 2)=0.2853138862893366318130781595187828d-01
	 w( 3)=0.4427743881741980616860274821133821d-01
	 w( 4)=0.5929858491543678074636775850010859d-01
	 w( 5)=0.7334648141108030573403361525311648d-01
	 w( 6)=0.8619016153195327591718520298374266d-01
	 w( 7)=0.9761865210411388826988066446424718d-01
	 w( 8)=0.1074442701159656347825773424466063d+00
	 w( 9)=0.1155056680537256013533444839067836d+00
	 w(10)=0.1216704729278033912044631534762624d+00
	 w(11)=0.1258374563468282961213753825111837d+00
	 w(12)=0.1279381953467521569740561652246955d+00
c
	 X1( 1)=0.9951872199970213601799974097007367d+00
	 X1( 2)=0.9747285559713094981983919930081690d+00
	 X1( 3)=0.9382745520027327585236490017087215d+00
	 X1( 4)=0.8864155270044010342131543419821967d+00
	 X1( 5)=0.8200019859739029219539498726697452d+00
	 X1( 6)=0.7401241915785543642438281030999784d+00
	 X1( 7)=0.6480936519369755692524957869107476d+00
	 X1( 8)=0.5454214713888395356583756172183724d+00
	 X1( 9)=0.4337935076260451384870842319133498d+00
	 X1(10)=0.3150426796961633743867932913198102d+00
	 X1(11)=0.1911188674736163091586398207570697d+00
	 X1(12)=0.6405689286260562608504308262474504d-01     


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
       w( 1)=0.3153346052305838632677311543891515d-02
	 w( 2)=0.7327553901276262102383979621786534d-02
	 w( 3)=0.1147723457923453948959266760909167d-01
	 w( 4)=0.1557931572294384872817695583446033d-01
	 w( 5)=0.1961616045735552781446071965221279d-01
	 w( 6)=0.2357076083932437914051930137844923d-01
	 w( 7)=0.2742650970835694820007383626250591d-01
	 w( 8)=0.3116722783279808890206575684635440d-01
	 w( 9)=0.3477722256477043889254858596380243d-01
	 w(10)=0.3824135106583070631721725652371569d-01
	 w(11)=0.4154508294346474921405882236106477d-01
	 w(12)=0.4467456085669428041944858712585040d-01
	 w(13)=0.4761665849249047482590662347892982d-01
	 w(14)=0.5035903555385447495780761908786555d-01
	 w(15)=0.5289018948519366709550505626469886d-01
	 w(16)=0.5519950369998416286820349519163536d-01
	 w(17)=0.5727729210040321570515023468470064d-01
	 w(18)=0.5911483969839563574647481743351990d-01
	 w(19)=0.6070443916589388005296923202782048d-01
	 w(20)=0.6203942315989266390419778413759853d-01
	 w(21)=0.6311419228625402565712602275023334d-01
	 w(22)=0.6392423858464818662390620182551554d-01
	 w(23)=0.6446616443595008220650419365770497d-01
	 w(24)=0.6473769681268392250302493873659151d-01
c 
	 X1( 1)=0.9987710072524261186005414915631136d+00
	 X1( 2)=0.9935301722663507575479287508490740d+00
	 X1( 3)=0.9841245837228268577445836000265988d+00
	 X1( 4)=0.9705915925462472504614119838006601d+00
	 X1( 5)=0.9529877031604308607229606660257181d+00
	 X1( 6)=0.9313866907065543331141743801016012d+00
	 X1( 7)=0.9058791367155696728220748356710117d+00
	 X1( 8)=0.8765720202742478859056935548050967d+00
	 X1( 9)=0.8435882616243935307110898445196559d+00
	 X1(10)=0.8070662040294426270825530430245384d+00
	 X1(11)=0.7671590325157403392538554375229691d+00
	 X1(12)=0.7240341309238146546744822334936653d+00
	 X1(13)=0.6778723796326639052118512806759091d+00
	 X1(14)=0.6288673967765136239951649330699947d+00
	 X1(15)=0.5772247260839727038178092385404788d+00
	 X1(16)=0.5231609747222330336782258691375086d+00
	 X1(17)=0.4669029047509584045449288616507985d+00
	 X1(18)=0.4086864819907167299162254958146334d+00
	 X1(19)=0.3487558862921607381598179372704079d+00
	 X1(20)=0.2873624873554555767358864613167976d+00
	 X1(21)=0.2247637903946890612248654401746923d+00
	 X1(22)=0.1612223560688917180564373907834977d+00
	 X1(23)=0.9700469920946269893005395585362452d-01
	 X1(24)=0.3238017096286936203332224315213445d-01     
C
        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 48                                   ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
C
48      IF (NG.GT.48) GO TO 96
C
        NG=48
C
	 w( 1)=0.7967920655520124294381434969436714d-03
	 w( 2)=0.1853960788946921732335925350893890d-02
	 w( 3)=0.2910731817934946408410617989400720d-02
	 w( 4)=0.3964554338444686673733415767419639d-02
	 w( 5)=0.5014202742927517692470194969030949d-02
	 w( 6)=0.6058545504235961683316742031729087d-02
	 w( 7)=0.7096470791153865269144160812143404d-02
	 w( 8)=0.8126876925698759217382427707855918d-02
	 w( 9)=0.9148671230783386632584602665207933d-02
	 w(10)=0.1016077053500841575758763695382536d-01
	 w(11)=0.1116210209983849859121326382856207d-01
	 w(12)=0.1215160467108831963518135273666725d-01
	 w(13)=0.1312822956696157263706366690259911d-01
	 w(14)=0.1409094177231486091586162472463661d-01
	 w(15)=0.1503872102699493800587627522209756d-01
	 w(16)=0.1597056290256229138061645679149306d-01
	 w(17)=0.1688547986424517245047754060686194d-01
	 w(18)=0.1778250231604526083761422648607104d-01
	 w(19)=0.1866067962741146738515675862213225d-01
	 w(20)=0.1951908114014502241008522028120468d-01
	 w(21)=0.2035679715433332459524521541727163d-01
	 w(22)=0.2117293989219129898767386719100365d-01
	 w(23)=0.2196664443874434919475638680156228d-01
	 w(24)=0.2273706965832937400134784197749033d-01
	 w(25)=0.2348339908592621984223593266761256d-01
	 w(26)=0.2420484179236469128226733787267704d-01
	 w(27)=0.2490063322248361028838218086833276d-01
	 w(28)=0.2557003600534936149879716794360011d-01
	 w(29)=0.2621234073567241391345796396446336d-01
	 w(30)=0.2682686672559176219805672871415660d-01
	 w(31)=0.2741296272602924282342108748909122d-01
	 w(32)=0.2797000761684833443981857658902251d-01
	 w(33)=0.2849741106508538564559951294580563d-01
	 w(34)=0.2899461415055523654267878127968158d-01
	 w(35)=0.2946108995816790597043633218285850d-01
	 w(36)=0.2989634413632838598438807579440005d-01
	 w(37)=0.3029991542082759379408876420650097d-01
	 w(38)=0.3067137612366914901422883035620429d-01
	 w(39)=0.3101033258631383742324977997063272d-01
	 w(40)=0.3131642559686135581278426671506315d-01
	 w(41)=0.3158933077072716855802074616996422d-01
	 w(42)=0.3182875889441100653475373988553221d-01
	 w(43)=0.3203445623199266321813897747021117d-01
	 w(44)=0.3220620479403025066866711455723255d-01
	 w(45)=0.3234382256857592842877483882894328d-01
	 w(46)=0.3244716371406426936401278844884587d-01
	 w(47)=0.3251611871386883598720549144777837d-01
	 w(48)=0.3255061449236316624196141829728569d-01
c
	 X1( 1)=0.9996895038832307668276901057843654d+00
	 X1( 2)=0.9983643758631816777241494395267193d+00
	 X1( 3)=0.9959818429872092906503990848793565d+00
	 X1( 4)=0.9925439003237626245718922977016011d+00
	 X1( 5)=0.9880541263296237994807627724184511d+00
	 X1( 6)=0.9825172635630146774470458318825480d+00
	 X1( 7)=0.9759391745851364664526010342418798d+00
	 X1( 8)=0.9683268284632642121736593664436509d+00
	 X1( 9)=0.9596882914487425393000680258847348d+00
	 X1(10)=0.9500327177844376357560989484438634d+00
	 X1(11)=0.9393703397527552169318573794316151d+00
	 X1(12)=0.9277124567223086909646904726965465d+00
	 X1(13)=0.9150714231208980742058844661504215d+00
	 X1(14)=0.9014606353158523413192326973096868d+00
	 X1(15)=0.8868945174024204160568774333937861d+00
	 X1(16)=0.8713885059092965028737747644930427d+00
	 X1(17)=0.8549590334346014554627869698935691d+00
	 X1(18)=0.8376235112281871214943028167646929d+00
	 X1(19)=0.8194003107379316755389996242244780d+00
	 X1(20)=0.8003087441391408172287961439713277d+00
	 X1(21)=0.7803690438674332176036044558893177d+00
	 X1(22)=0.7596023411766474987029704411196217d+00
	 X1(23)=0.7380306437444001328511657310964972d+00
	 X1(24)=0.7156768123489676262251441480576248d+00
	 X1(25)=0.6925645366421715613442457698182425d+00
	 X1(26)=0.6687183100439161539525572075727743d+00
	 X1(27)=0.6441634037849671067984123500633016d+00
	 X1(28)=0.6189258401254685703863692869878628d+00
	 X1(29)=0.5930323647775720806835557550437583d+00
	 X1(30)=0.5665104185613971684042501950834627d+00
	 X1(31)=0.5393881083243574362268025973255607d+00
	 X1(32)=0.5116941771546676735855097454288543d+00
	 X1(33)=0.4834579739205963597684056093638952d+00
	 X1(34)=0.4547094221677430086356761480863078d+00
	 X1(35)=0.4254789884073005453648192035699980d+00
	 X1(36)=0.3957976498289086032850002431351642d+00
	 X1(37)=0.3656968614723136350308955939954732d+00
	 X1(38)=0.3352085228926254226163256248056287d+00
	 X1(39)=0.3043649443544963530239297792952033d+00
	 X1(40)=0.2731988125910491414872721557245648d+00
	 X1(41)=0.2417431561638400123279319037406277d+00
	 X1(42)=0.2100313104605672036028471857122957d+00
	 X1(43)=0.1780968823676186027594026251764915d+00
	 X1(44)=0.1459737146548969419891073333433265d+00
	 X1(45)=0.1136958501106659209112080954143281d+00
	 X1(46)=0.8129749546442555899447126297475285d-01
	 X1(47)=0.4881298513604973111195820426127506d-01
	 X1(48)=0.1627674484960296957913456369523845d-01

        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 96                                   ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
96      IF (NG.GT.96) GO TO 192
C
        NG=96
C
	 w( 1)=0.2002510147471267303712312450255703d-03
	 w( 2)=0.4660944855912469603944432543932574d-03
	 w( 3)=0.7322075666250185328339640008748637d-03
	 w( 4)=0.9981598380233457554121719231494520d-03
	 w( 5)=0.1263853132096728285868613120015403d-02
	 w( 6)=0.1529211875019241493105816989986518d-02
	 w( 7)=0.1794164116955045626355485501631911d-02
	 w( 8)=0.2058638858115252727715815761111229d-02
	 w( 9)=0.2322565485048108214442097204508405d-02
	 w(10)=0.2585873625444375300279098771560957d-02
	 w(11)=0.2848493111476581340970402165390526d-02
	 w(12)=0.3110353977068882106680359244296400d-02
	 w(13)=0.3371386467342533589805998887040961d-02
	 w(14)=0.3631521052910285400855749274005996d-02
	 w(15)=0.3890688446257093249913484957572511d-02
	 w(16)=0.4148819619065216099139199471198278d-02
	 w(17)=0.4405845819973666428314679521316033d-02
	 w(18)=0.4661698592528197474902387832688645d-02
	 w(19)=0.4916309793197733507646532723064942d-02
	 w(20)=0.5169611609390087300715718380698600d-02
	 w(21)=0.5421536577428234926249655845129893d-02
	 w(22)=0.5672017600463228374204949713339389d-02
	 w(23)=0.5920987966307842199944703122938390d-02
	 w(24)=0.6168381365179540674609447146559622d-02
	 w(25)=0.6414131907343948572275024532122624d-02
	 w(26)=0.6658174140651556266610738434308557d-02
	 w(27)=0.6900443067961339719780742657597472d-02
	 w(28)=0.7140874164445576807569565961274595d-02
	 w(29)=0.7379403394770533802095458002934840d-02
	 w(30)=0.7615967230147961172624843647609506d-02
	 w(31)=0.7850502665252524421764359380945037d-02
	 w(32)=0.8082947235000432632340753774352802d-02
	 w(33)=0.8313239031184632745769256416255048d-02
	 w(34)=0.8541316718962022605943311383449181d-02
	 w(35)=0.8767119553188207673350963346520281d-02
	 w(36)=0.8990587394595389693608817807139734d-02
	 w(37)=0.9211660725809033543826169770456859d-02
	 w(38)=0.9430280667199013005803355133957693d-02
	 w(39)=0.9646388992560988688871555468912967d-02
	 w(40)=0.9859928144623822671848493861984853d-02
	 w(41)=0.1007084125037888528666258758886661d-01
	 w(42)=0.1027907213622716026310745432162761d-01
       w(43)=0.1048456534294010549872985896280295d-01
       w(44)=0.1068726614043027822182015353540765d-01
       w(45)=0.1088712054232778542665459065331425d-01
       w(46)=0.1108407532035867328136178184595929d-01
       w(47)=0.1127807801852142281099299668499557d-01
       w(48)=0.1146907696705777358968578140016646d-01
       w(49)=0.1165702129621415246785292850221955d-01
       w(50)=0.1184186094979004053271591814062855d-01
       w(51)=0.1202354669846966856400385119887345d-01
       w(52)=0.1220203015293348920546093199411790d-01
       w(53)=0.1237726377674593292623165356797658d-01
       w(54)=0.1254920089901601458978069148872322d-01
       w(55)=0.1271779572682741807446091689663609d-01
       w(56)=0.1288300335743474788966926602598000d-01
       w(57)=0.1304477979022269909350163852581688d-01
       w(58)=0.1320308193842496002934176412517347d-01
       w(59)=0.1335786764059972644620390859976834d-01
       w(60)=0.1350909567185877043631007883710086d-01
       w(61)=0.1365672575484707329799543649929899d-01
       w(62)=0.1380071857047009789665883850857899d-01
       w(63)=0.1394103576836584333464095632054565d-01
       w(64)=0.1407763997711889273569915932555554d-01
       w(65)=0.1421049481421373368383892660234473d-01
       w(66)=0.1433956489572470031198777092644869d-01
       w(67)=0.1446481584573995619537692547529471d-01
       w(68)=0.1458621430551700804926100230937411d-01
       w(69)=0.1470372794236731174222873844805851d-01
       w(70)=0.1481732545826760429607097933074648d-01
       w(71)=0.1492697659819566833198447984901299d-01
       w(72)=0.1503265215818830882160301087474154d-01
       w(73)=0.1513432399311939599056635198463324d-01
       w(74)=0.1523196502419590278248468769118197d-01
       w(75)=0.1532554924616994040247847083367136d-01
       w(76)=0.1541505173426487110205508422732436d-01
       w(77)=0.1550044865081365352085059342061589d-01
       w(78)=0.1558171725160765254549632202027095d-01
       w(79)=0.1565883589195422276120398725391178d-01
       w(80)=0.1573178403244145213710425099997847d-01
       w(81)=0.1580054224440853058129963261059210d-01
       w(82)=0.1586509221512028640526162020266563d-01
       w(83)=0.1592541675264451252875776745840061d-01
       w(84)=0.1598149979043078341497460420826010d-01
       w(85)=0.1603332639158954322984235013839654d-01
       w(86)=0.1608088275287032554860887339963451d-01
       w(87)=0.1612415620833804506520547322754080d-01
       w(88)=0.1616313523274638217456553266231246d-01
       w(89)=0.1619780944460736197339177074554857d-01
       w(90)=0.1622816960895631013944081316803146d-01
       w(91)=0.1625420763981144928166232902304553d-01
       w(92)=0.1627591660232748068189234778679601d-01
       w(93)=0.1629329071464257785160185083292900d-01
       w(94)=0.1630632534941829998274048170492025d-01
       w(95)=0.1631501703507201515824817210238839d-01
       w(96)=0.1631936345670150508355623929347410d-01
c
       X1( 1)=0.9999219686591948443225892454883059d+00
       X1( 2)=0.9995888803394138541093349296043197d+00
       X1( 3)=0.9989897220082228180635279012706708d+00
       X1( 4)=0.9981245201227147822147993855489200d+00
       X1( 5)=0.9969934888089658276499713406258741d+00
       X1( 6)=0.9955969254054171881427829799848505d+00
       X1( 7)=0.9939352005680611362729976453293206d+00
       X1( 8)=0.9920087563423258732188356447321776d+00
       X1( 9)=0.9898181055591399025125907309048817d+00
       X1(10)=0.9873638315366433759878785288998682d+00
       X1(11)=0.9846465878619336764167001102060532d+00
       X1(12)=0.9816670981895111523937086705610450d+00
       X1(13)=0.9784261560354051061458091997963069d+00
       X1(14)=0.9749246245590810488762592665419201d+00
       X1(15)=0.9711634363298709597893756810891852d+00
       X1(16)=0.9671435930764901577246260169273177d+00
       X1(17)=0.9628661654189837492787994912503302d+00
       X1(18)=0.9583322925828043949981517329007647d+00
       X1(19)=0.9535431820948995063904708574256973d+00
       X1(20)=0.9485001094617782425446726957406750d+00
       X1(21)=0.9432044178295800403958793265082119d+00
       X1(22)=0.9376575176261969636700289673768231d+00
       X1(23)=0.9318608861855216719820172701963348d+00
       X1(24)=0.9258160673539062609128154511394387d+00
       X1(25)=0.9195246710789272245880472369639542d+00
       X1(26)=0.9129883729805597709004952516054175d+00
       X1(27)=0.9062089139048714766238190927436892d+00
       X1(28)=0.8991880994603512666000259568610842d+00
       X1(29)=0.8919277995369952102259633387959993d+00
       X1(30)=0.8844299478082758137705466992072068d+00
       X1(31)=0.8766965412161264494836028717904297d+00
       X1(32)=0.8687296394390773650062052096541978d+00
       X1(33)=0.8605313643436843991581167574853416d+00
       X1(34)=0.8521038994193961183039219550285386d+00
       X1(35)=0.8434494891970095976588225636632542d+00
       X1(36)=0.8345704386508695147393385041118966d+00
       X1(37)=0.8254691125849696045538228486220810d+00
       X1(38)=0.8161479350031198524218644920022155d+00
       X1(39)=0.8066093884633470731986553853368341d+00
       X1(40)=0.7968560134167007467203242051792675d+00
       X1(41)=0.7868904075306401492697565606090533d+00
       X1(42)=0.7767152249971829400510635555749104d+00
       X1(43)=0.7663331758259994299497462391583613d+00
       X1(44)=0.7557470251226407768848703807442227d+00
       X1(45)=0.7449595923520933172914643282748366d+00
       X1(46)=0.7339737505878551560441256578153761d+00
       X1(47)=0.7227924257467349966987265710195044d+00
       X1(48)=0.7114185958095769994814961263379360d+00
       X1(49)=0.6998552900281192051475475113416287d+00
       X1(50)=0.6881055881181967577927191283377693d+00
       X1(51)=0.6761726194395047980470490663027341d+00
       X1(52)=0.6640595621621394789116155447266718d+00
       X1(53)=0.6517696424201390789275344473475309d+00
       X1(54)=0.6393061334522506504939562614901836d+00
       X1(55)=0.6266723547301510440942257233984526d+00
       X1(56)=0.6138716710743544910676773343071987d+00
       X1(57)=0.6009074917580422075109537861843085d+00
       X1(58)=0.5877832695990526990513018873113860d+00
       X1(59)=0.5745025000402745997622328877199733d+00
       X1(60)=0.5610687202186869675605952774257933d+00
       X1(61)=0.5474855080232949822198280472393228d+00
       X1(62)=0.5337564811422119498589125315046377d+00
       X1(63)=0.5198852960991414086385045279429992d+00
       X1(64)=0.5058756472795159536497884562071851d+00
       X1(65)=0.4917312659465521538697583777912980d+00
       X1(66)=0.4774559192474836198483994125790001d+00
       X1(67)=0.4630534092102368967759482341181050d+00
       X1(68)=0.4485275717308174030575108927536054d+00
       X1(69)=0.4338822755516751088216675844448257d+00
       X1(70)=0.4191214212313220512518038667756913d+00
       X1(71)=0.4042489401054761136151635308892819d+00
       X1(72)=0.3892687932400077517555031425089794d+00
       X1(73)=0.3741849703759685350105230464587519d+00
       X1(74)=0.3590014888669824774342224742254436d+00
       X1(75)=0.3437223926092831692859114825625485d+00
       X1(76)=0.3283517509646816774505792608815985d+00
       X1(77)=0.3128936576767520662584319958241016d+00
       X1(78)=0.2973522297805231965735868566629289d+00
       X1(79)=0.2817316065059671905423468209843095d+00
       X1(80)=0.2660359481755766015698183178924267d+00
       X1(81)=0.2502694350963239034900637063721692d+00
       X1(82)=0.2344362664462984090903798119956503d+00
       X1(83)=0.2185406591563171457463379882822288d+00
       X1(84)=0.2025868467868075545434961993931732d+00
       X1(85)=0.1865790784002611385477991734258360d+00
       X1(86)=0.1705216174295583655043303941593981d+00
       X1(87)=0.1544187405424662298791849492833018d+00
       X1(88)=0.1382747365026108985189998820297989d+00
       X1(89)=0.1220939050272288030156827800270713d+00
       X1(90)=0.1058805556420003998801369668227609d+00
       X1(91)=0.8963900653327159661982075744412597d-01
       X1(92)=0.7337358339796853757389944640916407d-01
       X1(93)=0.5708861829151205770128136652529276d-01
       X1(94)=0.4078844847403864527738161700471162d-01
       X1(95)=0.2447741525523520549315998396275543d-01
       X1(96)=0.8159862838095286144725462139639069d-02
C
        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 192                                  ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
192      IF (NG.GT.192) GO TO 384
C
         NG=192
C
       w(  1)=0.5019410348692173752939580444469347d-04
       w(  2)=0.1168390665730188627954282573512990d-03
       w(  3)=0.1835749193551260444766282687846253d-03
       w(  4)=0.2503070890844147244367666054448446d-03
       w(  5)=0.3170242698112706309204038817114463d-03
       w(  6)=0.3837208020912924380769776054432363d-03
       w(  7)=0.4503919137716877613812631428790360d-03
       w(  8)=0.5170330453491546380708515921736194d-03
       w(  9)=0.5836397042629966772866241664826160d-03
       w( 10)=0.6502074240969915513807189141732146d-03
       w( 11)=0.7167317509947397544734551873122715d-03
       w( 12)=0.7832082385905052313623308973839090d-03
       w( 13)=0.8496324460039079108554845086356611d-03
       w( 14)=0.9159999370632673684792719266271580d-03
       w( 15)=0.9823062800663751210383405437425685d-03
       w( 16)=0.1048547047793687527743074799806691d-02
       w( 17)=0.1114717817647312658872806756262422d-02
       w( 18)=0.1180814171855886469925513966209686d-02
       w( 19)=0.1246831697715436938175984171901580d-02
       w( 20)=0.1312765987850601065545162646052231d-02
       w( 21)=0.1378612640487679042415226157425863d-02
       w( 22)=0.1444367259734716215028497908118924d-02
       w( 23)=0.1510025455865806445301950005132107d-02
       w( 24)=0.1575582845607929439091787931193488d-02
       w( 25)=0.1641035052429278657631432315109335d-02
       w( 26)=0.1706377706828417190445077796020022d-02
       w( 27)=0.1771606446623829624142789694803191d-02
       w( 28)=0.1836716917243580798316589604725254d-02
       w( 29)=0.1901704772014882631752333906403915d-02
       w( 30)=0.1966565672453428345731083190287523d-02
       w( 31)=0.2031295288552391503809476913047267d-02
       w( 32)=0.2095889299071012646298758209214139d-02
       w( 33)=0.2160343391822713433314413546749456d-02
       w( 34)=0.2224653263962689952606445977779841d-02
       w( 35)=0.2288814622274945002047155445794683d-02
       w( 36)=0.2352823183458724890143200545087301d-02
       w( 37)=0.2416674674414330378673324717230462d-02
       w( 38)=0.2480364832528274330294802477920565d-02
       w( 39)=0.2543889405957760764920170338479650d-02
       w( 40)=0.2607244153914461606545741921501410d-02
       w( 41)=0.2670424846947568578829181362191945d-02
       w( 42)=0.2733427267226098595794760803005978d-02
       w( 43)=0.2796247208820431672903378569127224d-02
       w( 44)=0.2858880477983060909682385985388926d-02
       w( 45)=0.2921322893428534508519486756552241d-02
       w( 46)=0.2983570286612570124238412542485457d-02
       w( 47)=0.3045618502010322106922345096791569d-02
       w( 48)=0.3107463397393782421712854698238788d-02
       w( 49)=0.3169100844108296215501620180239663d-02
       w( 50)=0.3230526727348173160118143165099009d-02
       w( 51)=0.3291736946431375841217790633522596d-02
       w( 52)=0.3352727415073266586429107133956190d-02
       w( 53)=0.3413494061659394239096272143057411d-02
       w( 54)=0.3474032829517302487949757945706370d-02
       w( 55)=0.3534339677187341460401017745200512d-02
       w( 56)=0.3594410578692464379538911645922593d-02
       w( 57)=0.3654241523806991173587274806303053d-02
       w( 58)=0.3713828518324321012561164471010274d-02
       w( 59)=0.3773167584323575830912870459951390d-02
       w( 60)=0.3832254760435156977699152627836126d-02
       w( 61)=0.3891086102105197217710156177243493d-02
       w( 62)=0.3949657681858890388459062191086982d-02
       w( 63)=0.4007965589562681099242087057301574d-02
       w( 64)=0.4066005932685296939881577617258525d-02
       w( 65)=0.4123774836557605748452690149827222d-02
       w( 66)=0.4181268444631280569420091916030623d-02
       w( 67)=0.4238482918736255016298372596434625d-02
       w( 68)=0.4295414439336951836297167614953607d-02
       w( 69)=0.4352059205787267558498685300423258d-02
       w( 70)=0.4408413436584296192005077629759063d-02
       w( 71)=0.4464473369620775026236987584419831d-02
       w( 72)=0.4520235262436235672203545086973304d-02
       w( 73)=0.4575695392466843571130698967835241d-02
       w( 74)=0.4630850057293909285354959836294116d-02
       w( 75)=0.4685695574891054975883778349808404d-02
       w( 76)=0.4740228283870019561507770496630028d-02
       w( 77)=0.4794444543725086145835985874347312d-02
       w( 78)=0.4848340735076115391122488071208516d-02
       w( 79)=0.4901913259910168611267227157949856d-02
       w( 80)=0.4955158541821704450910967740974816d-02
       w( 81)=0.5008073026251333113105562322508227d-02
       w( 82)=0.5060653180723112194628321228366439d-02
       w( 83)=0.5112895495080368285622576841672028d-02
       w( 84)=0.5164796481720028588884672211481845d-02
       w( 85)=0.5216352675825446913778534306791000d-02
       w( 86)=0.5267560635597708500439975819280702d-02
       w( 87)=0.5318416942485398231630515459523692d-02
       w( 88)=0.5368918201412816892310868929135109d-02
       w( 89)=0.5419061041006630240722915992774529d-02
       w( 90)=0.5468842113820935759491052843955629d-02
       w( 91)=0.5518258096560732060974190869931274d-02
       w( 92)=0.5567305690303776027812743496351508d-02
       w( 93)=0.5615981620720812877314953198997290d-02
       w( 94)=0.5664282638294164447007804661322833d-02
       w( 95)=0.5712205518534661108333280508549526d-02
       w( 96)=0.5759747062196902826094378399809767d-02
       w( 97)=0.5806904095492834992840479525720119d-02
       w( 98)=0.5853673470303624779921533980872358d-02
       w( 99)=0.5900052064389823860428155662763938d-02
       w(100)=0.5946036781599803473663010369179463d-02
       w(101)=0.5991624552076447916150622583933663d-02
       w(102)=0.6036812332462092660480589860761524d-02
       w(103)=0.6081597106101693420485741134637252d-02
       w(104)=0.6125975883244212599374467415858922d-02
       w(105)=0.6169945701242209676457655012639579d-02
       w(106)=0.6213503624749622208027636789813014d-02
       w(107)=0.6256646745917724238751533827339605d-02
       w(108)=0.6299372184589249041626393658557443d-02
       w(109)=0.6341677088490663227100668621788001d-02
       w(110)=0.6383558633422579385387769513699468d-02
       w(111)=0.6425014023448294550274553378309604d-02
       w(112)=0.6466040491080441897852466968525296d-02
       w(113)=0.6506635297465743219563407826493119d-02
       w(114)=0.6546795732567849835747855268282835d-02
       w(115)=0.6586519115348259743501071436292001d-02
       w(116)=0.6625802793945298921075722647056811d-02
       w(117)=0.6664644145851154840307606718594922d-02
       w(118)=0.6703040578086950368576715976102163d-02
       w(119)=0.6740989527375846372639982636434371d-02
       w(120)=0.6778488460314161468276050167038758d-02
       w(121)=0.6815534873540497492057541529323633d-02
       w(122)=0.6852126293902859404703748170145686d-02
       w(123)=0.6888260278623758469357583153393716d-02
       w(124)=0.6923934415463287682766115548614226d-02
       w(125)=0.6959146322880158572715066354226053d-02
       w(126)=0.6993893650190688611165282578892367d-02
       w(127)=0.7028174077725728629354348759431617d-02
       w(128)=0.7061985316985519758650027489941020d-02
       w(129)=0.7095325110792469559164976829021600d-02
       w(130)=0.7128191233441837137054958908735018d-02
       w(131)=0.7160581490850317191016269447963523d-02
       w(132)=0.7192493720702513068763074281317703d-02
       w(133)=0.7223925792595289055192382656499322d-02
       w(134)=0.7254875608179992255524118965100018d-02
       w(135)=0.7285341101302534578926731351462142d-02
       w(136)=0.7315320238141325470995510168094003d-02
       w(137)=0.7344811017343046186931751382647443d-02
       w(138)=0.7373811470156256541366517357587224d-02
       w(139)=0.7402319660562825215473405886890571d-02
       w(140)=0.7430333685407174847310783050492116d-02
       w(141)=0.7457851674523333277215671196279933d-02
       w(142)=0.7484871790859782466529175526518238d-02
       w(143)=0.7511392230602096754957207871825633d-02
       w(144)=0.7537411223293362269450512667035560d-02
       w(145)=0.7562927031952369445614768837178497d-02
       w(146)=0.7587937953189570771324946434885427d-02
       w(147)=0.7612442317320796011408216543945447d-02
       w(148)=0.7636438488478717321966589889492444d-02
       w(149)=0.7659924864722056813124101781926437d-02
       w(150)=0.7682899878142529269693742460534107d-02
       w(151)=0.7705361994969512890456393254497869d-02
       w(152)=0.7727309715672441058417678260561926d-02
       w(153)=0.7748741575060908306548754509131573d-02
       w(154)=0.7769656142382483796113485647677424d-02
       w(155)=0.7790052021418225777726989227644196d-02
       w(156)=0.7809927850575890658769000066558214d-02
       w(157)=0.7829282302980830454679607002281370d-02
       w(158)=0.7848114086564572555984424289025206d-02
       w(159)=0.7866421944151075897620850722052546d-02
       w(160)=0.7884204653540657772256421021438552d-02
       w(161)=0.7901461027591585684794010336408512d-02
       w(162)=0.7918189914299328801136433509273305d-02
       w(163)=0.7934390196873463700524380495810257d-02
       w(164)=0.7950060793812229297356218285761781d-02
       w(165)=0.7965200658974725955335514386937625d-02
       w(166)=0.7979808781650753974061721090304012d-02
       w(167)=0.7993884186628286785770804762619995d-02
       w(168)=0.8007425934258574357835190243047412d-02
       w(169)=0.8020433120518872454835676165684177d-02
       w(170)=0.8032904877072793572511401702041506d-02
       w(171)=0.8044840371328275514666928430531446d-02
       w(172)=0.8056238806493163743157443720547272d-02
       w(173)=0.8067099421628403790373377115032378d-02
       w(174)=0.8077421491698840183193714388325539d-02
       w(175)=0.8087204327621618487162344370528777d-02
       w(176)=0.8096447276312187239653214547491269d-02
       w(177)=0.8105149720727896701017221018990322d-02
       w(178)=0.8113311079909191513135920363982719d-02
       w(179)=0.8120930809018394515433615358296035d-02
       w(180)=0.8128008399376079129209410370631953d-02
       w(181)=0.8134543378495027882133720414389284d-02
       w(182)=0.8140535310111774805898703493707823d-02
       w(183)=0.8145983794215729601308411459789776d-02
       w(184)=0.8150888467075881626531352383873214d-02
       w(185)=0.8155249001265081925804850360173328d-02
       w(186)=0.8159065105681901677566290907904316d-02
       w(187)=0.8162336525570065602780258027416842d-02
       w(188)=0.8165063042535459036121901533080955d-02
       w(189)=0.8167244474560707524654813078338456d-02
       w(190)=0.8168880676017327980695423063927239d-02
       w(191)=0.8169971537675450577674640072765497d-02
       w(192)=0.8170516986711110739980316957776570d-02
c
       X1(  1)=0.9999804411726473547435916080674043d+00
       X1(  2)=0.9998969471378595991468054406187173d+00
       X1(  3)=0.9997467408113522753988972414026949d+00
       X1(  4)=0.9995297988558859110049282600075867d+00
       X1(  5)=0.9992461316671845033851407841964070d+00
       X1(  6)=0.9988957572063257575662544639878614d+00
       X1(  7)=0.9984786985384589563177774591129954d+00
       X1(  8)=0.9979949833727937972095028935222694d+00
       X1(  9)=0.9974446439389107953978233829595598d+00
       X1( 10)=0.9968277169440913106011587182604726d+00
       X1( 11)=0.9961442435551086373704026483342547d+00
       X1( 12)=0.9953942693885953365534400381277338d+00
       X1( 13)=0.9945778445047067583913757666163057d+00
       X1( 14)=0.9936950234020882231198241628454611d+00
       X1( 15)=0.9927458650133152938198183044691738d+00
       X1( 16)=0.9917304327004320079833462855178316d+00
       X1( 17)=0.9906487942504060643059327389568952d+00
       X1( 18)=0.9895010218704086785345938500939276d+00
       X1( 19)=0.9882871921828698242876652773997211d+00
       X1( 20)=0.9870073862202814985581032513524182d+00
       X1( 21)=0.9856616894197333315606559277141099d+00
       X1( 22)=0.9842501916171713291634871599288893d+00
       X1( 23)=0.9827729870413742449448393161947657d+00
       X1( 24)=0.9812301743076442733429848278256271d+00
       X1( 25)=0.9796218564112100920471705233487289d+00
       X1( 26)=0.9779481407203411185936260705153875d+00
       X1( 27)=0.9762091389691723840573414857694040d+00
       X1( 28)=0.9744049672502397817145518444257367d+00
       X1( 29)=0.9725357460067256891482271342125922d+00
       X1( 30)=0.9706016000244151299324478785733199d+00
       X1( 31)=0.9686026584233627613468938469142880d+00
       X1( 32)=0.9665390546492710636330171850153891d+00
       X1( 33)=0.9644109264645801744096617485385597d+00
       X1( 34)=0.9622184159392698657186713295280470d+00
       X1( 35)=0.9599616694413742051951783143330813d+00
       X1( 36)=0.9576408376272094800157457907014543d+00
       X1( 37)=0.9552560754313159945804566987849321d+00
       X1( 38)=0.9528075420561143817127760802366948d+00
       X1( 39)=0.9502954009612770934723309692376080d+00
       X1( 40)=0.9477198198528157621437301911075226d+00
       X1( 41)=0.9450809706718851450664253714904562d+00
       X1( 42)=0.9423790295833043890479544200210348d+00
       X1( 43)=0.9396141769637963714023014550010982d+00
       X1( 44)=0.9367865973899458953560325449694781d+00
       X1( 45)=0.9338964796258775377989448863308919d+00
       X1( 46)=0.9309440166106539672203921286354512d+00
       X1( 47)=0.9279294054453955692394478673227273d+00
       X1( 48)=0.9248528473801222364606223411731031d+00
       X1( 49)=0.9217145478003181985074972751087354d+00
       X1( 50)=0.9185147162132207870349702943350292d+00
       X1( 51)=0.9152535662338340493199450128510118d+00
       X1( 52)=0.9119313155706681426978188772843445d+00
       X1( 53)=0.9085481860112054606613366380471192d+00
       X1( 54)=0.9051044034070944598794091866618400d+00
       X1( 55)=0.9016001976590721757340055783620212d+00
       X1( 56)=0.8980358027016164322189825155509176d+00
       X1( 57)=0.8944114564873287702000339620129925d+00
       X1( 58)=0.8907274009710491361030051461123556d+00
       X1( 59)=0.8869838820937033910809182316524975d+00
       X1( 60)=0.8831811497658847186098080398417027d+00
       X1( 61)=0.8793194578511700262809281937597013d+00
       X1( 62)=0.8753990641491724552926944264179587d+00
       X1( 63)=0.8714202303783311288001759785935056d+00
       X1( 64)=0.8673832221584392878530512051251624d+00
       X1( 65)=0.8632883089929119811445175287540802d+00
       X1( 66)=0.8591357642507944922033244774248383d+00
       X1( 67)=0.8549258651485127049883796988593641d+00
       X1( 68)=0.8506588927313666260896490767762044d+00
       X1( 69)=0.8463351318547682988996322075658534d+00
       X1( 70)=0.8419548711652253621957711528218481d+00
       X1( 71)=0.8375184030810715225649140193318276d+00
       X1( 72)=0.8330260237729452270055301759391325d+00
       X1( 73)=0.8284780331440178388608488149881887d+00
       X1( 74)=0.8238747348099726369655259746972595d+00
       X1( 75)=0.8192164360787359745288728764424483d+00
       X1( 76)=0.8145034479299619508281182573903209d+00
       X1( 77)=0.8097360849942719652446330962187934d+00
       X1( 78)=0.8049146655322505395435109626909809d+00
       X1( 79)=0.8000395114131988105713566146269877d+00
       X1( 80)=0.7951109480936471117275694026534827d+00
       X1( 81)=0.7901293045956280776497929820694757d+00
       X1( 82)=0.7850949134847117225435133205577459d+00
       X1( 83)=0.7800081108478039584779970252953213d+00
       X1( 84)=0.7748692362707100357648461775913433d+00
       X1( 85)=0.7696786328154644032303802974686968d+00
       X1( 86)=0.7644366469974285017878192808063809d+00
       X1( 87)=0.7591436287621580202088147500339729d+00
       X1( 88)=0.7537999314620411573852465037346998d+00
       X1( 89)=0.7484059118327094506603551001816149d+00
       X1( 90)=0.7429619299692227449922151617960297d+00
       X1( 91)=0.7374683493020298927912658642704690d+00
       X1( 92)=0.7319255365727067892461097624838560d+00
       X1( 93)=0.7263338618094733628170787542357151d+00
       X1( 94)=0.7206936983024911553341626661454074d+00
       X1( 95)=0.7150054225789431407838238964062842d+00
       X1( 96)=0.7092694143778974464070093666021278d+00
       X1( 97)=0.7034860566249566541573538749467285d+00
       X1( 98)=0.6976557354066943748831886867027143d+00
       X1( 99)=0.6917788399448808017985746044446108d+00
       X1(100)=0.6858557625704989638962256126180645d+00
       X1(101)=0.6798868986975534139279404201613783d+00
       X1(102)=0.6738726467966730994350850362576510d+00
       X1(103)=0.6678134083685101790518474960354795d+00
       X1(104)=0.6617095879169365599265010826617543d+00
       X1(105)=0.6555615929220399456098575348715410d+00
       X1(106)=0.6493698338129211971445671305719373d+00
       X1(107)=0.6431347239402948233530363131371438d+00
       X1(108)=0.6368566795488944294646016464415135d+00
       X1(109)=0.6305361197496849662433452508668629d+00
       X1(110)=0.6241734664918836346756934003516616d+00
       X1(111)=0.6177691445347913140508466508936422d+00
       X1(112)=0.6113235814194363939162948824599147d+00
       X1(113)=0.6048372074400329029143303550023993d+00
       X1(114)=0.5983104556152548399027509859457489d+00
       X1(115)=0.5917437616593286250330175920622718d+00
       X1(116)=0.5851375639529456006011742680295714d+00
       X1(117)=0.5784923035139965235000503483323530d+00
       X1(118)=0.5718084239681300029848340165972527d+00
       X1(119)=0.5650863715191368492172486913231557d+00
       X1(120)=0.5583265949191623096754895981040859d+00
       X1(121)=0.5515295454387481820070139260283435d+00
       X1(122)=0.5446956768367068032584569273814345d+00
       X1(123)=0.5378254453298289266406103544919316d+00
       X1(124)=0.5309193095624275080757997207496409d+00
       X1(125)=0.5239777305757194357293929483507150d+00
       X1(126)=0.5170011717770472465458339376197121d+00
       X1(127)=0.5099900989089428844917984178933022d+00
       X1(128)=0.5029449800180355657541031658341907d+00
       X1(129)=0.4958662854238058265471594928099309d+00
       X1(130)=0.4887544876871878394533531875564727d+00
       X1(131)=0.4816100615790220943490704702881039d+00
       X1(132)=0.4744334840483605499584968708817867d+00
       X1(133)=0.4672252341906263719261265188289889d+00
       X1(134)=0.4599857932156303830064757348266593d+00
       X1(135)=0.4527156444154463605351490729288543d+00
       X1(136)=0.4454152731321473257685195731821442d+00
       X1(137)=0.4380851667254049789592289513884266d+00
       X1(138)=0.4307258145399544431708683301861818d+00
       X1(139)=0.4233377078729264888269560524147196d+00
       X1(140)=0.4159213399410494198360858966659831d+00
       X1(141)=0.4084772058477228108362860940704126d+00
       X1(142)=0.4010058025499652936566260976977899d+00
       X1(143)=0.3935076288252385995023630085956919d+00
       X1(144)=0.3859831852381500716308716351558133d+00
       X1(145)=0.3784329741070358713986999053061071d+00
       X1(146)=0.3708574994704271085247931815886368d+00
       X1(147)=0.3632572670534011342307052690318163d+00
       X1(148)=0.3556327842338202435849388160444399d+00
       X1(149)=0.3479845600084600408949216252276269d+00
       X1(150)=0.3403131049590297293560263589567875d+00
       X1(151)=0.3326189312180865933819875314974761d+00
       X1(152)=0.3249025524348469491045798886735075d+00
       X1(153)=0.3171644837408958454420247591082737d+00
       X1(154)=0.3094052417157978048948243288461299d+00
       X1(155)=0.3016253443526108998341368253182129d+00
       X1(156)=0.2938253110233064665009572795200358d+00
       X1(157)=0.2860056624440967652338280796185322d+00
       X1(158)=0.2781669206406729015881504068838585d+00
       X1(159)=0.2703096089133553290009916150615600d+00
       X1(160)=0.2624342518021592594911847510776239d+00
       X1(161)=0.2545413750517773145651051338829009d+00
       X1(162)=0.2466315055764817540234059925128776d+00
       X1(163)=0.2387051714249486257328317902941839d+00
       X1(164)=0.2307629017450061846396456234879560d+00
       X1(165)=0.2228052267483099343568580151017328d+00
       X1(166)=0.2148326776749466495559910286414503d+00
       X1(167)=0.2068457867579697421352268902362858d+00
       X1(168)=0.1988450871878683387191577282500035d+00
       X1(169)=0.1908311130769724414706666451301036d+00
       X1(170)=0.1828043994237965484624347057734921d+00
       X1(171)=0.1747654820773241139638987097839602d+00
       X1(172)=0.1667148977012352329489065511304805d+00
       X1(173)=0.1586531837380799379195669054759770d+00
       X1(174)=0.1505808783733994997726148961683054d+00
       X1(175)=0.1424985204997981279057728905596349d+00
       X1(176)=0.1344066496809674680728439437942368d+00
       X1(177)=0.1263058061156662996474135866673014d+00
       X1(178)=0.1181965306016578369458432859131389d+00
       X1(179)=0.1100793644996070420905163308066413d+00
       X1(180)=0.1019548496969403595638551772479879d+00
       X1(181)=0.9382352857167028511229019667736039d-01
       X1(182)=0.8568594395618718400695586041120788d-01
       X1(183)=0.7754263910102077585426492071875290d-01
       X1(184)=0.6939415763857370517451817128958850d-01
       X1(185)=0.6124104354682961883021169726308101d-01
       X1(186)=0.5308384111303817308758075577250051d-01
       X1(187)=0.4492309489737939463505594887932863d-01
       X1(188)=0.3675934969660982126060034610772658d-01
       X1(189)=0.2859315050769284910625367206999765d-01
       X1(190)=0.2042504249141571447255039484499058d-01
       X1(191)=0.1225557093599553903773214761107266d-01
       X1(192)=0.4085281220676868086781824858344715d-02
c
        GO TO 56
C
CÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
C³                       N --- 384                                  ³
CÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C
384      IF (NG.GT.384) GO TO 768
         NG=384
C
       w(  1)=0.1256492650122374769407672465620970d-04
       w(  2)=0.2924855339195397923090402019676296d-04
       w(  3)=0.4595639581651652537714570719159999d-04
       w(  4)=0.6266561599843194923302312584609243d-04
       w(  5)=0.7937421975252693487486411044483867d-04
       w(  6)=0.9608162609380794717612541320368149d-04
       w(  7)=0.1127874758170570632739999421182688d-03
       w(  8)=0.1294914625472838870266679925310064d-03
       w(  9)=0.1461932961959379392871611535521947d-03
       w( 10)=0.1628926926784371500832355258072985d-03
       w( 11)=0.1795893704439478210866983885762628d-03
       w( 12)=0.1962830491368474503771049130422675d-03
       w( 13)=0.2129734490251937594959759796672039d-03
       w( 14)=0.2296602907368082685391787520191135d-03
       w( 15)=0.2463432951299490661507339453159447d-03
       w( 16)=0.2630221832272330304788562310072830d-03
       w( 17)=0.2796966761811145625865776565522124d-03
       w( 18)=0.2963664952559300631695746581176249d-03
       w( 19)=0.3130313618190279569993428094186245d-03
       w( 20)=0.3296909973370752613991422344349431d-03
       w( 21)=0.3463451233754138504059705522189756d-03
       w( 22)=0.3629934615992673516559701188629916d-03
       w( 23)=0.3796357337761009941778966343798979d-03
       w( 24)=0.3962716617787168291446323348797247d-03
       w( 25)=0.4129009675888279577287316413204026d-03
       w( 26)=0.4295233733009507076610478131508901d-03
       w( 27)=0.4461386011265114315799542146715876d-03
       w( 28)=0.4627463733981003530736066528383651d-03
       w( 29)=0.4793464125738274805255273572646418d-03
       w( 30)=0.4959384412417501552038573091179351d-03
       w( 31)=0.5125221821243513267913444235122212d-03
       w( 32)=0.5290973580830539884038044940544978d-03
       w( 33)=0.5456636921227614832074885762343488d-03
       w( 34)=0.5622209073964163244064546124050348d-03
       w( 35)=0.5787687272095722013774187404686595d-03
       w( 36)=0.5953068750249752695855301356789180d-03
       w( 37)=0.6118350744671518327255165990608134d-03
       w( 38)=0.6283530493270002501559373089866779d-03
       w( 39)=0.6448605235663854272948677150293173d-03
       w( 40)=0.6613572213227346299763904636349361d-03
       w( 41)=0.6778428669136336463357320490805802d-03
       w( 42)=0.6943171848414225298166286736852358d-03
       w( 43)=0.7107798997977903142319974667635517d-03
       w( 44)=0.7272307366683682105579785289614903d-03
       w( 45)=0.7436694205373208854038048236150482d-03
       w( 46)=0.7600956766919354901743526340126280d-03
       w( 47)=0.7765092306272081631440800269956227d-03
       w( 48)=0.7929098080504277678910395271770850d-03
       w( 49)=0.8092971348857566636809837484172100d-03
       w( 50)=0.8256709372788083285902105236740425d-03
       w( 51)=0.8420309416012216760194724473897056d-03
       w( 52)=0.8583768744552319209921643628217264d-03
       w( 53)=0.8747084626782378651675338978273906d-03
       w( 54)=0.8910254333473654795352741729082663d-03
       w( 55)=0.9073275137840276718296596767097826d-03
       w( 56)=0.9236144315584801322243819958156021d-03
       w( 57)=0.9398859144943731561651381785110409d-03
       w( 58)=0.9561416906732993475164309080946744d-03
       w( 59)=0.9723814884393371087381062934307581d-03
       w( 60)=0.9886050364035898277201064180877550d-03
       w( 61)=0.1004812063448720673312542527548449d-02
       w( 62)=0.1021002298733482913589089007632623d-02
       w( 63)=0.1037175471697245672551770448982949d-02
       w( 64)=0.1053331312064515042385950141487114d-02
       w( 65)=0.1069469549849450469555069201082068d-02
       w( 66)=0.1085589915360376334025171607244356d-02
       w( 67)=0.1101692139204288641761603499449172d-02
       w( 68)=0.1117775952291356751370548202109099d-02
       w( 69)=0.1133841085839420056387452893928299d-02
       w( 70)=0.1149887271378479545260227283754957d-02
       w( 71)=0.1165914240755184161551524025938063d-02
       w( 72)=0.1181921726137311887303088611253613d-02
       w( 73)=0.1197909460018245472875685688211880d-02
       w( 74)=0.1213877175221442736908394502449039d-02
       w( 75)=0.1229824604904901360337684298367892d-02
       w( 76)=0.1245751482565618098685085398611283d-02
       w( 77)=0.1261657542044042337067005944910552d-02
       w( 78)=0.1277542517528523912605100961549805d-02
       w( 79)=0.1293406143559755129123754340570558d-02
       w( 80)=0.1309248155035206889215349329172888d-02
       w( 81)=0.1325068287213558868936319299974435d-02
       w( 82)=0.1340866275719123660569376188390366d-02
       w( 83)=0.1356641856546264809051403980930257d-02
       w( 84)=0.1372394766063808667823629233264517d-02
       w( 85)=0.1388124741019450000011985132484457d-02
       w( 86)=0.1403831518544151250992043621454577d-02
       w( 87)=0.1419514836156535418535330677541148d-02
       w( 88)=0.1435174431767272446872970428123720d-02
       w( 89)=0.1450810043683459071149029141376016d-02
       w( 90)=0.1466421410612992038870168139360581d-02
       w( 91)=0.1482008271668934635090709608604356d-02
       w( 92)=0.1497570366373876438203352882171661d-02
       w( 93)=0.1513107434664286233335879600009519d-02
       w( 94)=0.1528619216894858010483537448216614d-02
       w( 95)=0.1544105453842849974635638446703080d-02
       w( 96)=0.1559565886712416495283460395245769d-02
       w( 97)=0.1575000257138932922824982117966890d-02
       w( 98)=0.1590408307193313199510473223188977d-02
       w( 99)=0.1605789779386320192701634931168762d-02
       w(100)=0.1621144416672868678345969788943525d-02
       w(101)=0.1636471962456320902697448901050989d-02
       w(102)=0.1651772160592774650444436324516531d-02
       w(103)=0.1667044755395343747536300539979140d-02
       w(104)=0.1682289491638430927131261369887137d-02
       w(105)=0.1697506114561992987219847570665483d-02
       w(106)=0.1712694369875798168610928174103534d-02
       w(107)=0.1727854003763675682100675468854042d-02
       w(108)=0.1742984762887757313779059499673966d-02
       w(109)=0.1758086394392711037563598285982956d-02
       w(110)=0.1773158645909966564186125384990869d-02
       w(111)=0.1788201265561932755995313856023168d-02
       w(112)=0.1803214001966206837075636651213019d-02
       w(113)=0.1818196604239775328322368527647823d-02
       w(114)=0.1833148822003206637252161705625006d-02
       w(115)=0.1848070405384835232469672351943834d-02
       w(116)=0.1862961105024937332852691190573173d-02
       w(117)=0.1877820672079898041661250997274438d-02
       w(118)=0.1892648858226369855920256706609953d-02
       w(119)=0.1907445415665422481570319242930935d-02
       w(120)=0.1922210097126683885027679641287117d-02
       w(121)=0.1936942655872472511941392115415059d-02
       w(122)=0.1951642845701920604084299052117853d-02
       w(123)=0.1966310420955088545463782194658880d-02
       w(124)=0.1980945136517070168888816483391552d-02
       w(125)=0.1995546747822088954381489425296149d-02
       w(126)=0.2010115010857585050973882128989251d-02
       w(127)=0.2024649682168293053585040384740670d-02
       w(128)=0.2039150518860310466827697017144431d-02
       w(129)=0.2053617278605156787750441379656299d-02
       w(130)=0.2068049719643823139678169089741598d-02
       w(131)=0.2082447600790812389471885361452602d-02
       w(132)=0.2096810681438169680688278700872301d-02
       w(133)=0.2111138721559503315279928123840040d-02
       w(134)=0.2125431481713995916638556010271550d-02
       w(135)=0.2139688723050405806946389566313642d-02
       w(136)=0.2153910207311058531964445753337759d-02
       w(137)=0.2168095696835828466551406397563510d-02
       w(138)=0.2182244954566110434372700777086511d-02
       w(139)=0.2196357744048781275426460896120142d-02
       w(140)=0.2210433829440151295181158355925855d-02
       w(141)=0.2224472975509905529288969542058027d-02
       w(142)=0.2238474947645034758009245970566147d-02
       w(143)=0.2252439511853756204647887162654119d-02
       w(144)=0.2266366434769423852490922336605903d-02
       w(145)=0.2280255483654428314884202411901696d-02
       w(146)=0.2294106426404086193285783120765071d-02
       w(147)=0.2307919031550518858293341143871203d-02
       w(148)=0.2321693068266520588825805782083304d-02
       w(149)=0.2335428306369416004816306328126257d-02
       w(150)=0.2349124516324906728952527529910951d-02
       w(151)=0.2362781469250907213180629799551571d-02
       w(152)=0.2376398936921369665870024521545631d-02
       w(153)=0.2389976691770098015718495300027889d-02
       w(154)=0.2403514506894550848660420558904114d-02
       w(155)=0.2417012156059633254225178823882160d-02
       w(156)=0.2430469413701477517978202481832408d-02
       w(157)=0.2443886054931212596863585998835254d-02
       w(158)=0.2457261855538722314454647612560767d-02
       w(159)=0.2470596591996392213307386488665498d-02
       w(160)=0.2483890041462845001801367300349758d-02
       w(161)=0.2497141981786664533043198176700662d-02
       w(162)=0.2510352191510108253599442957934101d-02
       w(163)=0.2523520449872808060018521651874764d-02
       w(164)=0.2536646536815459501294900833672047d-02
       w(165)=0.2549730232983499265623655368628656d-02
       w(166)=0.2562771319730770889989291137019303d-02
       w(167)=0.2575769579123178631329552243861045d-02
       w(168)=0.2588724793942329438212792323639090d-02
       w(169)=0.2601636747689162962166364792580729d-02
       w(170)=0.2614505224587569547993378027225108d-02
       w(171)=0.2627330009587996142616065201878961d-02
       w(172)=0.2640110888371040062185931619922435d-02
       w(173)=0.2652847647351030557403761522955011d-02
       w(174)=0.2665540073679598117196488233385838d-02
       w(175)=0.2678187955249231451102852734496533d-02
       w(176)=0.2690791080696822090925693050204575d-02
       w(177)=0.2703349239407196552415616666846630d-02
       w(178)=0.2715862221516635997958707332948458d-02
       w(179)=0.2728329817916383341449802451700784d-02
       w(180)=0.2740751820256137736742744496514825d-02
       w(181)=0.2753128020947536391279855964918790d-02
       w(182)=0.2765458213167623646714708853339822d-02
       w(183)=0.2777742190862307268555052979277365d-02
       w(184)=0.2789979748749801887066529173636682d-02
       w(185)=0.2802170682324059531892519871802510d-02
       w(186)=0.2814314787858187203061177386443895d-02
       w(187)=0.2826411862407851421267315569158806d-02
       w(188)=0.2838461703814669700534450065225257d-02
       w(189)=0.2850464110709588886580822321694497d-02
       w(190)=0.2862418882516250304432739292214943d-02
       w(191)=0.2874325819454341659049000743466717d-02
       w(192)=0.2886184722542935632941565542015290d-02
       w(193)=0.2897995393603815124999923603945699d-02
       w(194)=0.2909757635264785074949887623258803d-02
       w(195)=0.2921471250962970818101694542515389d-02
       w(196)=0.2933136044948102915267407257950862d-02
       w(197)=0.2944751822285788402953628512248948d-02
       w(198)=0.2956318388860768409162477556130036d-02
       w(199)=0.2967835551380162080361632173556793d-02
       w(200)=0.2979303117376696765413000267372858d-02
       w(201)=0.2990720895211924402479252579392644d-02
       w(202)=0.3002088694079424055158017442411283d-02
       w(203)=0.3013406324007990544325005886707692d-02
       w(204)=0.3024673595864809122399697090230560d-02
       w(205)=0.3035890321358616136980466194192126d-02
       w(206)=0.3047056313042845631030175013310114d-02
       w(207)=0.3058171384318761827028267246236166d-02
       w(208)=0.3069235349438577442741309515443362d-02
       w(209)=0.3080248023508557786500694000592415d-02
       w(210)=0.3091209222492110580113863623772193d-02
       w(211)=0.3102118763212861457773932732731248d-02
       w(212)=0.3112976463357715089571951028220216d-02
       w(213)=0.3123782141479901878456292098060331d-02
       w(214)=0.3134535617002010179724736343141721d-02
       w(215)=0.3145236710219003992376757284290258d-02
       w(216)=0.3155885242301226071896306184422531d-02
       w(217)=0.3166481035297386414279018553865192d-02
       w(218)=0.3177023912137536061362233360082810d-02
       w(219)=0.3187513696636026177761517554064120d-02
       w(220)=0.3197950213494452349963520757863564d-02
       w(221)=0.3208333288304584058371943520824363d-02
       w(222)=0.3218662747551279273351183321410847d-02
       w(223)=0.3228938418615384126560821329035567d-02
       w(224)=0.3239160129776617609123525693657901d-02
       w(225)=0.3249327710216441248419169634404999d-02
       w(226)=0.3259440990020913715548990672503500d-02
       w(227)=0.3269499800183530315765446805326954d-02
       w(228)=0.3279503972608047314416052040764090d-02
       w(229)=0.3289453340111291051202893284378150d-02
       w(230)=0.3299347736425951795813738862600705d-02
       w(231)=0.3309186996203362298235641726837766d-02
       w(232)=0.3318970955016260987317713356042740d-02
       w(233)=0.3328699449361539771406293286072797d-02
       w(234)=0.3338372316662976395133059757004342d-02
       w(235)=0.3347989395273951306694714885584509d-02
       w(236)=0.3357550524480148990221728727286378d-02
       w(237)=0.3367055544502243718093236266359034d-02
       w(238)=0.3376504296498569678315545425417360d-02
       w(239)=0.3385896622567775432342828268438405d-02
       w(240)=0.3395232365751462658980427319874729d-02
       w(241)=0.3404511370036809140273809962796693d-02
       w(242)=0.3413733480359175945549541822949861d-02
       w(243)=0.3422898542604698770038720493397373d-02
       w(244)=0.3432006403612863384778109493787745d-02
       w(245)=0.3441056911179065154749734564795144d-02
       w(246)=0.3450049914057152582485945835357768d-02
       w(247)=0.3458985261961954834633905619361508d-02
       w(248)=0.3467862805571793209241128138779127d-02
       w(249)=0.3476682396530976501792069859476097d-02
       w(250)=0.3485443887452280228294842879603472d-02
       w(251)=0.3494147131919409663986894432637192d-02
       w(252)=0.3502791984489446656498958549946604d-02
       w(253)=0.3511378300695280172587736751934951d-02
       w(254)=0.3519905937048020537819598771430602d-02
       w(255)=0.3528374751039397328860107215820433d-02
       w(256)=0.3536784601144140878297357191731248d-02
       w(257)=0.3545135346822347352200978682938862d-02
       w(258)=0.3553426848521827360893171312443124d-02
       w(259)=0.3561658967680438063683323446006602d-02
       w(260)=0.3569831566728398728593605808556765d-02
       w(261)=0.3577944509090589708379419277379175d-02
       w(262)=0.3585997659188834794425712666751387d-02
       w(263)=0.3593990882444166910377964496543157d-02
       w(264)=0.3601924045279077107645038300550199d-02
       w(265)=0.3609797015119746825190169326484103d-02
       w(266)=0.3617609660398263376306016845571757d-02
       w(267)=0.3625361850554818625350016052082300d-02
       w(268)=0.3633053456039890817697182007821238d-02
       w(269)=0.3640684348316409526449050579388924d-02
       w(270)=0.3648254399861903679719583122200977d-02
       w(271)=0.3655763484170632632601608070124337d-02
       w(272)=0.3663211475755700248200718868196870d-02
       w(273)=0.3670598250151151952407489103342676d-02
       w(274)=0.3677923683914054727363397499380814d-02
       w(275)=0.3685187654626560008860972893197955d-02
       w(276)=0.3692390040897949453204367634162397d-02
       w(277)=0.3699530722366663539342842274430044d-02
       w(278)=0.3706609579702312972376490159837028d-02
       w(279)=0.3713626494607672854820942796132125d-02
       w(280)=0.3720581349820659592305770850448786d-02
       w(281)=0.3727474029116290500669826540834169d-02
       w(282)=0.3734304417308626081705856145711898d-02
       w(283)=0.3741072400252694935096341599487118d-02
       w(284)=0.3747777864846401274372702790105688d-02
       w(285)=0.3754420699032415015020702390164044d-02
       w(286)=0.3761000791800044403146137977068052d-02
       w(287)=0.3767518033187091153406676962783524d-02
       w(288)=0.3773972314281688065207983584105842d-02
       w(289)=0.3780363527224119086455099015560195d-02
       w(290)=0.3786691565208621794443360665823825d-02
       w(291)=0.3792956322485172263766980003150705d-02
       w(292)=0.3799157694361252291417734915401084d-02
       w(293)=0.3805295577203598949541067727232303d-02
       w(294)=0.3811369868439936436612208643768696d-02
       w(295)=0.3817380466560690198090761631356739d-02
       w(296)=0.3823327271120683287908490638377190d-02
       w(297)=0.3829210182740814942441823651005611d-02
       w(298)=0.3835029103109721338917845410739642d-02
       w(299)=0.3840783934985418510500271725007911d-02
       w(300)=0.3846474582196927390600084203939314d-02
       w(301)=0.3852100949645880959254148972370310d-02
       w(302)=0.3857662943308113464714241445686560d-02
       w(303)=0.3863160470235231693688446622837407d-02
       w(304)=0.3868593438556168263976895534030418d-02
       w(305)=0.3873961757478716913544228470988324d-02
       w(306)=0.3879265337291049760372039403682023d-02
       w(307)=0.3884504089363216507735848521335864d-02
       w(308)=0.3889677926148625569852866092270037d-02
       w(309)=0.3894786761185507093148945774740489d-02
       w(310)=0.3899830509098357848695674080317694d-02
       w(311)=0.3904809085599367971671499836541973d-02
       w(312)=0.3909722407489829524004168153814834d-02
       w(313)=0.3914570392661526856655482503197028d-02
       w(314)=0.3919352960098108748313570980783901d-02
       w(315)=0.3924070029876442297562373587852114d-02
       w(316)=0.3928721523167948545902991304749588d-02
       w(317)=0.3933307362239919809306839784875856d-02
       w(318)=0.3937827470456818696286225541077080d-02
       w(319)=0.3942281772281558790774005432150614d-02
       w(320)=0.3946670193276766978410395967538528d-02
       w(321)=0.3950992660106027395141762313265623d-02
       w(322)=0.3955249100535106977343332775328879d-02
       w(323)=0.3959439443433162592985247825776378d-02
       w(324)=0.3963563618773929733669158283511038d-02
       w(325)=0.3967621557636892747670729922597088d-02
       w(326)=0.3971613192208436594431886406065917d-02
       w(327)=0.3975538455782980101255423877909035d-02
       w(328)=0.3979397282764090703263753629632188d-02
       w(329)=0.3983189608665580647992968824580796d-02
       w(330)=0.3986915370112584646303182142094886d-02
       w(331)=0.3990574504842618951596138218071169d-02
       w(332)=0.3994166951706621849641462727242397d-02
       w(333)=0.3997692650669975541623563688947344d-02
       w(334)=0.4001151542813509403332144890922481d-02
       w(335)=0.4004543570334484603730521018419190d-02
       w(336)=0.4007868676547560066447433947677230d-02
       w(337)=0.4011126805885739758049854507730859d-02
       w(338)=0.4014317903901301287266308622353622d-02
       w(339)=0.4017441917266705799642585899820217d-02
       w(340)=0.4020498793775489152424267222812503d-02
       w(341)=0.4023488482343134354773340480628548d-02
       w(342)=0.4026410933007925258739255053020374d-02
       w(343)=0.4029266096931781486718090767484862d-02
       w(344)=0.4032053926401074581447080573589550d-02
       w(345)=0.4034774374827425364895522868867169d-02
       w(346)=0.4037427396748482492727144027005445d-02
       w(347)=0.4040012947828682191323218972725337d-02
       w(348)=0.4042530984859989164670222367390990d-02
       w(349)=0.4044981465762618658730459859912641d-02
       w(350)=0.4047364349585739671229012660215573d-02
       w(351)=0.4049679596508159295105414145333274d-02
       w(352)=0.4051927167838988184193759045545077d-02
       w(353)=0.4054107026018287130010418711138043d-02
       w(354)=0.4056219134617694738844194757361071d-02
       w(355)=0.4058263458341036198659582750878851d-02
       w(356)=0.4060239963024913125639832257425484d-02
       w(357)=0.4062148615639274480512674236384285d-02
       w(358)=0.4063989384287968545117936159808615d-02
       w(359)=0.4065762238209275949992774064337211d-02
       w(360)=0.4067467147776423744066913725255839d-02
       w(361)=0.4069104084498080497877104980761935d-02
       w(362)=0.4070673021018832432026948636948904d-02
       w(363)=0.4072173931119640562935349053397411d-02
       w(364)=0.4073606789718278858234072146485073d-02
       w(365)=0.4074971572869753394492242851323050d-02
       w(366)=0.4076268257766702510263092750185028d-02
       w(367)=0.4077496822739777947765862299956569d-02
       w(368)=0.4078657247258006976833467566080494d-02
       w(369)=0.4079749511929135495074353286429461d-02
       w(370)=0.4080773598499952098514867134440952d-02
       w(371)=0.4081729489856593117306498913944260d-02
       w(372)=0.4082617170024828611400427784051899d-02
       w(373)=0.4083436624170329321410005165405628d-02
       w(374)=0.4084187838598914570200065401677652d-02
       w(375)=0.4084870800756781111060295223904513d-02
       w(376)=0.4085485499230712918638301270107063d-02
       w(377)=0.4086031923748271919126487027698552d-02
       w(378)=0.4086510065177969656515381269310362d-02
       w(379)=0.4086919915529419892044644020797167d-02
       w(380)=0.4087261467953472134301608009360245d-02
       w(381)=0.4087534716742326097735888065178992d-02
       w(382)=0.4087739657329627087677302766197494d-02
       w(383)=0.4087876286290542310263096396667234d-02
       w(384)=0.4087944601341818105999219709072160d-02
c

c
       X1(  1)=0.9999951039143946033831476342009955d+00
       X1(  2)=0.9999742029064620072299638133436188d+00
       X1(  3)=0.9999366007554937931553363074495621d+00
       X1(  4)=0.9998822897379444526114654385140236d+00
       X1(  5)=0.9998112697383989329639633946593274d+00
       X1(  6)=0.9997235416993873692257339295441196d+00
       X1(  7)=0.9996191070055611184133547274119944d+00
       X1(  8)=0.9994979673689068484926857183791214d+00
       X1(  9)=0.9993601247982200545826798467407850d+00
       X1( 10)=0.9992055815889295082589377489533662d+00
       X1( 11)=0.9990343403191145837052802064305259d+00
       X1( 12)=0.9988464038477390529126187726024920d+00
       X1( 13)=0.9986417753137795765855865248793286d+00
       X1( 14)=0.9984204581357500596801402080512267d+00
       X1( 15)=0.9981824560114139054078948650523083d+00
       X1( 16)=0.9979277729175902025857745514091375d+00
       X1( 17)=0.9976564131100084719319204281829502d+00
       X1( 18)=0.9973683811231888002010426185798071d+00
       X1( 19)=0.9970636817703349506884851789039380d+00
       X1( 20)=0.9967423201432335218513903152626904d+00
       X1( 21)=0.9964043016121551447298674936116798d+00
       X1( 22)=0.9960496318257553242558883556461144d+00
       X1( 23)=0.9956783167109734531067313345983190d+00
       X1( 24)=0.9952903624729290712501249551070601d+00
       X1( 25)=0.9948857755948147741094889607216418d+00
       X1( 26)=0.9944645628377853769377334947637377d+00
       X1( 27)=0.9940267312408430728359097488985648d+00
       X1( 28)=0.9935722881207184059134407798096178d+00
       X1( 29)=0.9931012410717469365211786448144234d+00
       X1( 30)=0.9926135979657415126741997293112952d+00
       X1( 31)=0.9921093669518600871230854826653932d+00
       X1( 32)=0.9915885564564690370593004899706116d+00
       X1( 33)=0.9910511751830019557310550255672409d+00
       X1( 34)=0.9904972321118138939784095334026751d+00
       X1( 35)=0.9899267365000310359782966952000499d+00
       X1( 36)=0.9893396978813957980636269125887754d+00
       X1( 37)=0.9887361260661073428492043544904975d+00
       X1( 38)=0.9881160311406575034066736914645285d+00
       X1( 39)=0.9874794234676621141214653771888950d+00
       X1( 40)=0.9868263136856877463055141584790664d+00
       X1( 41)=0.9861567127090738477501735398887924d+00
       X1( 42)=0.9854706317277502862700548453940144d+00
       X1( 43)=0.9847680822070502979725263423237299d+00
       X1( 44)=0.9840490758875188415345580553788553d+00
       X1( 45)=0.9833136247847163602117524721776670d+00
       X1( 46)=0.9825617411890179536684028093330464d+00
       X1( 47)=0.9817934376654079620206606830660696d+00
       X1( 48)=0.9810087270532699647411969105366585d+00
       X1( 49)=0.9802076224661721972935613750509449d+00
       X1( 50)=0.9793901372916483885557539721284314d+00
       X1( 51)=0.9785562851909740222614092489888565d+00
       X1( 52)=0.9777060800989380258381731518626460d+00
       X1( 53)=0.9768395362236098901599534562329871d+00
       X1( 54)=0.9759566680461022238555929806515035d+00
       X1( 55)=0.9750574903203287459333676908402415d+00
       X1( 56)=0.9741420180727577205902977146878754d+00
       X1( 57)=0.9732102666021608381789593119654931d+00
       X1( 58)=0.9722622514793575464033966106116951d+00
       X1( 59)=0.9712979885469548359107278864424280d+00
       X1( 60)=0.9703174939190824845368224942235018d+00
       X1( 61)=0.9693207839811237645535541023765002d+00
       X1( 62)=0.9683078753894416173520678220358400d+00
       X1( 63)=0.9672787850711003000816009863272522d+00
       X1( 64)=0.9662335302235825088469686707099779d+00
       X1( 65)=0.9651721283145019831501093193590496d+00
       X1( 66)=0.9640945970813115963422824708661830d+00
       X1( 67)=0.9630009545310069369337831883531485d+00
       X1( 68)=0.9618912189398253856875210503245164d+00
       X1( 69)=0.9607654088529406935016165676507646d+00
       X1( 70)=0.9596235430841530651643865932288343d+00
       X1( 71)=0.9584656407155747541427989327144732d+00
       X1( 72)=0.9572917210973111736427383311992758d+00
       X1( 73)=0.9561018038471375292562941032568975d+00
       X1( 74)=0.9548959088501709785877979634581957d+00
       X1( 75)=0.9536740562585383233265459496411786d+00
       X1( 76)=0.9524362664910392393100615686964721d+00
       X1( 77)=0.9511825602328050501974243076344623d+00
       X1( 78)=0.9499129584349530504476201587547198d+00
       X1( 79)=0.9486274823142363833730870502885741d+00
       X1( 80)=0.9473261533526894801136433851034275d+00
       X1( 81)=0.9460089932972690654508151567116250d+00
       X1( 82)=0.9446760241594907364572271618993273d+00
       X1( 83)=0.9433272682150611200502057536088887d+00
       X1( 84)=0.9419627480035056155930619981569922d+00
       X1( 85)=0.9405824863277917287616913945825030d+00
       X1( 86)=0.9391865062539480029681448010895071d+00
       X1( 87)=0.9377748311106785547066993151399710d+00
       X1( 88)=0.9363474844889732192616912241680764d+00
       X1( 89)=0.9349044902417133132899687828887823d+00
       X1( 90)=0.9334458724832730208642829203963566d+00
       X1( 91)=0.9319716555891164096372609952658708d+00
       X1( 92)=0.9304818641953900838588039470904505d+00
       X1( 93)=0.9289765231985114810528118312724878d+00
       X1( 94)=0.9274556577547528192320776581055625d+00
       X1( 95)=0.9259192932798207016029953115388187d+00
       X1( 96)=0.9243674554484313857844044942833241d+00
       X1( 97)=0.9228001701938817246374443340324525d+00
       X1( 98)=0.9212174637076157858757075199280613d+00
       X1( 99)=0.9196193624387871576972785004911871d+00
       X1(100)=0.9180058930938169477524021170612074d+00
       X1(101)=0.9163770826359474828325627116604821d+00
       X1(102)=0.9147329582847917167386577792052799d+00
       X1(103)=0.9130735475158783538577240904901973d+00
       X1(104)=0.9113988780601926960493172798809608d+00
       X1(105)=0.9097089779037132205141574911843565d+00
       X1(106)=0.9080038752869438963890330714388604d+00
       X1(107)=0.9062835987044422478832007108911500d+00
       X1(108)=0.9045481769043431718426330206740836d+00
       X1(109)=0.9027976388878785176994424525348403d+00
       X1(110)=0.9010320139088924378346527993163504d+00
       X1(111)=0.8992513314733525164531953439168511d+00
       X1(112)=0.8974556213388566851405750961055088d+00
       X1(113)=0.8956449135141359333410824968059141d+00
       X1(114)=0.8938192382585528220677164849899114d+00
       X1(115)=0.8919786260815958092241349030802829d+00
       X1(116)=0.8901231077423693949889568380374121d+00
       X1(117)=0.8882527142490800957826076196626968d+00
       X1(118)=0.8863674768585182554066197765926865d+00
       X1(119)=0.8844674270755357020148812259090712d+00
       X1(120)=0.8825525966525192596457542779127774d+00
       X1(121)=0.8806230175888601231132745999656794d+00
       X1(122)=0.8786787221304191051247770228806416d+00
       X1(123)=0.8767197427689877645612839054758616d+00
       X1(124)=0.8747461122417454249258306086006116d+00
       X1(125)=0.8727578635307120920335903765725019d+00
       X1(126)=0.8707550298621972800861964859634762d+00
       X1(127)=0.8687376447062447553410418010700677d+00
       X1(128)=0.8667057417760732066545637723728156d+00
       X1(129)=0.8646593550275128522465953275962200d+00
       X1(130)=0.8625985186584379921007779324024604d+00
       X1(131)=0.8605232671081955154837912364530575d+00
       X1(132)=0.8584336350570293731337530674380261d+00
       X1(133)=0.8563296574255010237355829876111641d+00
       X1(134)=0.8542113693739058643684010815892853d+00
       X1(135)=0.8520788063016856546771499984223863d+00
       X1(136)=0.8499320038468369445875814237376556d+00
       X1(137)=0.8477709978853155154505370085886901d+00
       X1(138)=0.8455958245304368445680772311407471d+00
       X1(139)=0.8434065201322726031204686174515521d+00
       X1(140)=0.8412031212770431975793291017951851d+00
       X1(141)=0.8389856647865063647583519702647759d+00
       X1(142)=0.8367541877173418307189797105567530d+00
       X1(143)=0.8345087273605320438141790945180099d+00
       X1(144)=0.8322493212407389922190768588211989d+00
       X1(145)=0.8299760071156771163626503356753354d+00
       X1(146)=0.8276888229754823267399282345814490d+00
       X1(147)=0.8253878070420771376492424048420895d+00
       X1(148)=0.8230729977685319274639807361756589d+00
       X1(149)=0.8207444338384223361130233031009507d+00
       X1(150)=0.8184021541651828105085973519224657d+00
       X1(151)=0.8160461978914563087246606938529941d+00
       X1(152)=0.8136766043884401737931164333163257d+00
       X1(153)=0.8112934132552281880491736586634247d+00
       X1(154)=0.8088966643181488190209976879687700d+00
       X1(155)=0.8064863976300996679224386325316571d+00
       X1(156)=0.8040626534698781318710873552402688d+00
       X1(157)=0.8016254723415082910171823029255019d+00
       X1(158)=0.7991748949735640318319781269497708d+00
       X1(159)=0.7967109623184884178670864231461757d+00
       X1(160)=0.7942337155519093193590092723750299d+00
       X1(161)=0.7917431960719513131156065008888122d+00
       X1(162)=0.7892394454985438641835666628789414d+00
       X1(163)=0.7867225056727258008580886365309640d+00
       X1(164)=0.7841924186559460946579243832087413d+00
       X1(165)=0.7816492267293609569506828136798493d+00
       X1(166)=0.7790929723931272639748488053556867d+00
       X1(167)=0.7765236983656923220663291932654264d+00
       X1(168)=0.7739414475830799849584979909967869d+00
       X1(169)=0.7713462631981731350856751653972953d+00
       X1(170)=0.7687381885799925408807359729357104d+00
       X1(171)=0.7661172673129721021181101520169632d+00
       X1(172)=0.7634835431962304954137911432353399d+00
       X1(173)=0.7608370602428392320541339708402052d+00
       X1(174)=0.7581778626790871403851754591899054d+00
       X1(175)=0.7555059949437412850539610766311056d+00
       X1(176)=0.7528215016873043354529078983676506d+00
       X1(177)=0.7501244277712683957775719651575565d+00
       X1(178)=0.7474148182673653091673196951795101d+00
       X1(179)=0.7446927184568134484573259946402949d+00
       X1(180)=0.7419581738295610061290353245892309d+00
       X1(181)=0.7392112300835257961047252363603423d+00
       X1(182)=0.7364519331238315800901038089501447d+00
       X1(183)=0.7336803290620409312269520348275741d+00
       X1(184)=0.7308964642153846478756885360466352d+00
       X1(185)=0.7281003851059877304053860835093661d+00
       X1(186)=0.7252921384600919339262062757761409d+00
       X1(187)=0.7224717712072749099564394504814907d+00
       X1(188)=0.7196393304796659500733404953137031d+00
       X1(189)=0.7167948636111583446537367443880545d+00
       X1(190)=0.7139384181366183698669506410422096d+00
       X1(191)=0.7110700417910909161389263746038943d+00
       X1(192)=0.7081897825090017713625753151786194d+00
       X1(193)=0.7052976884233565721852588393065489d+00
       X1(194)=0.7023938078649364367600081264613590d+00
       X1(195)=0.6994781893614902924025377815367432d+00
       X1(196)=0.6965508816369239116513427751391772d+00
       X1(197)=0.6936119336104856702831752688705647d+00
       X1(198)=0.6906613943959490408909784877904865d+00
       X1(199)=0.6876993133007918356859080016172937d+00
       X1(200)=0.6847257398253722122393956684718645d+00
       X1(201)=0.6817407236621014559353071724043921d+00
       X1(202)=0.6787443146946135529561096447110666d+00
       X1(203)=0.6757365629969315676806003991186895d+00
       X1(204)=0.6727175188326308384241504361182610d+00
       X1(205)=0.6696872326539990055055861897583901d+00
       X1(206)=0.6666457551011928856777691126442873d+00
       X1(207)=0.6635931370013922070116342372049473d+00
       X1(208)=0.6605294293679502183759149328805577d+00
       X1(209)=0.6574546833995411877070108230636944d+00
       X1(210)=0.6543689504793048033154483596649224d+00
       X1(211)=0.6512722821739874925271380082882936d+00
       X1(212)=0.6481647302330806720091475083908697d+00
       X1(213)=0.6450463465879559441809863796471949d+00
       X1(214)=0.6419171833509972541634318912320901d+00
       X1(215)=0.6387772928147300217677202420784099d+00
       X1(216)=0.6356267274509472630784778685978484d+00
       X1(217)=0.6324655399098327162340757571564457d+00
       X1(218)=0.6292937830190809860581535511169995d+00
       X1(219)=0.6261115097830147222458792699181757d+00
       X1(220)=0.6229187733816988458580837679641201d+00
       X1(221)=0.6197156271700518389257358256616224d+00
       X1(222)=0.6165021246769541120163031594946355d+00
       X1(223)=0.6132783196043534646623758424298413d+00
       X1(224)=0.6100442658263676536015108242088319d+00
       X1(225)=0.6068000173883840838245886213681964d+00
       X1(226)=0.6035456285061566374780550014769096d+00
       X1(227)=0.6002811535648996557132508116157587d+00
       X1(228)=0.5970066471183790886236111982708179d+00
       X1(229)=0.5937221638880008284578405395101957d+00
       X1(230)=0.5904277587618962413442406697251350d+00
       X1(231)=0.5871234867940049128081866357477758d+00
       X1(232)=0.5838094032031546224113054984664107d+00
       X1(233)=0.5804855633721385628872188080692583d+00
       X1(234)=0.5771520228467898191947575599634441d+00
       X1(235)=0.5738088373350531229553489127333365d+00
       X1(236)=0.5704560627060538977868059540233386d+00
       X1(237)=0.5670937549891646110910245740746575d+00
       X1(238)=0.5637219703730684478981042932644886d+00
       X1(239)=0.5603407652048203224141619371841827d+00
       X1(240)=0.5569501959889052429645976127014543d+00
       X1(241)=0.5535503193862940460688007676005228d+00
       X1(242)=0.5501411922134965154262494756846857d+00
       X1(243)=0.5467228714416119016376577439533786d+00
       X1(244)=0.5432954141953768585282628583315466d+00
       X1(245)=0.5398588777522108119835168435527848d+00
       X1(246)=0.5364133195412587772503522897293257d+00
       X1(247)=0.5329587971424316406998323758742957d+00
       X1(248)=0.5294953682854439220893671866221110d+00
       X1(249)=0.5260230908488490334047826645353418d+00
       X1(250)=0.5225420228590720504043640630831864d+00
       X1(251)=0.5190522224894400130285618655104417d+00
       X1(252)=0.5155537480592097708803441177709686d+00
       X1(253)=0.5120466580325933900222042993741981d+00
       X1(254)=0.5085310110177811373765875388103057d+00
       X1(255)=0.5050068657659620590569794891390588d+00
       X1(256)=0.5014742811703421689971108378566889d+00
       X1(257)=0.4979333162651602642856655613615334d+00
       X1(258)=0.4943840302247013836535419808601945d+00
       X1(259)=0.4908264823623079256001017706109198d+00
       X1(260)=0.4872607321293884426839526528019921d+00
       X1(261)=0.4836868391144241285426449325261027d+00
       X1(262)=0.4801048630419730142442196322677866d+00
       X1(263)=0.4765148637716718906118261337264269d+00
       X1(264)=0.4729169012972359732006292859449468d+00
       X1(265)=0.4693110357454563266439492575440737d+00
       X1(266)=0.4656973273751950651230213669717067d+00
       X1(267)=0.4620758365763783457519270923324711d+00
       X1(268)=0.4584466238689871717061308205087884d+00
       X1(269)=0.4548097499020460219596590275556735d+00
       X1(270)=0.4511652754526093245322788770803973d+00
       X1(271)=0.4475132614247457901840710735482361d+00
       X1(272)=0.4438537688485206235305466109497063d+00
       X1(273)=0.4401868588789756285869282165158498d+00
       X1(274)=0.4365125927951072257854042114090537d+00
       X1(275)=0.4328310319988423975440646075302030d+00
       X1(276)=0.4291422380140125795008459485109248d+00
       X1(277)=0.4254462724853255145601421051140025d+00
       X1(278)=0.4217431971773350869337823769364888d+00
       X1(279)=0.4180330739734091533918352645780969d+00
       X1(280)=0.4143159648746953889721655951797591d+00
       X1(281)=0.4105919319990851644308537501391043d+00
       X1(282)=0.4068610375801754727484780023732914d+00
       X1(283)=0.4031233439662289220398638720536747d+00
       X1(284)=0.3993789136191318122472174094451253d+00
       X1(285)=0.3956278091133503130285818713540369d+00
       X1(286)=0.3918700931348847602852888385851751d+00
       X1(287)=0.3881058284802220888035148954532547d+00
       X1(288)=0.3843350780552864185162030333864440d+00
       X1(289)=0.3805579048743878119224634284721283d+00
       X1(290)=0.3767743720591692202321306617952813d+00
       X1(291)=0.3729845428375516358334232908692962d+00
       X1(292)=0.3691884805426774687116264345267272d+00
       X1(293)=0.3653862486118521644763982014062182d+00
       X1(294)=0.3615779105854840816846858776479296d+00
       X1(295)=0.3577635301060226461753273015301693d+00
       X1(296)=0.3539431709168948001602063054212130d+00
       X1(297)=0.3501168968614397638453280174064692d+00
       X1(298)=0.3462847718818421273833797100617672d+00
       X1(299)=0.3424468600180632909872452908216747d+00
       X1(300)=0.3386032254067712710615459809518221d+00
       X1(301)=0.3347539322802688902365857669698158d+00
       X1(302)=0.3308990449654203692160873731509722d+00
       X1(303)=0.3270386278825763383768123451777719d+00
       X1(304)=0.3231727455444972870845669067226968d+00
       X1(305)=0.3193014625552754687172031114795175d+00
       X1(306)=0.3154248436092552794110320265816311d+00
       X1(307)=0.3115429534899521285725718181953260d+00
       X1(308)=0.3076558570689698192227582401099378d+00
       X1(309)=0.3037636193049164562656477301632951d+00
       X1(310)=0.2998663052423189007982436811827550d+00
       X1(311)=0.2959639800105357886023740616851320d+00
       X1(312)=0.2920567088226691309835430108179287d+00
       X1(313)=0.2881445569744745161453699209606820d+00
       X1(314)=0.2842275898432699293116164541329877d+00
       X1(315)=0.2803058728868432098308845240484319d+00
       X1(316)=0.2763794716423581635218461285707802d+00
       X1(317)=0.2724484517252593485393386568224387d+00
       X1(318)=0.2685128788281755530638265457162433d+00
       X1(319)=0.2645728187198219831385915517860734d+00
       X1(320)=0.2606283372439011790005690705468867d+00
       X1(321)=0.2566795003180026782719965170052486d+00
       X1(322)=0.2527263739325014444009814222753839d+00
       X1(323)=0.2487690241494550787597312525707188d+00
       X1(324)=0.2448075171014998348295136733094818d+00
       X1(325)=0.2408419189907454529214347230088801d+00
       X1(326)=0.2368722960876688339018327945059790d+00
       X1(327)=0.2328987147300065704104881154472553d+00
       X1(328)=0.2289212413216463540789402517031861d+00
       X1(329)=0.2249399423315172772749897073226533d+00
       X1(330)=0.2209548842924790479179336489463760d+00
       X1(331)=0.2169661338002101359272498325247549d+00
       X1(332)=0.2129737575120948698852966521880590d+00
       X1(333)=0.2089778221461095025121405668344182d+00
       X1(334)=0.2049783944797072635678546962611889d+00
       X1(335)=0.2009755413487024188145538274893673d+00
       X1(336)=0.1969693296461533536870411505389741d+00
       X1(337)=0.1929598263212447003372404737292341d+00
       X1(338)=0.1889470983781685267335741792325083d+00
       X1(339)=0.1849312128750046065121215029494607d+00
       X1(340)=0.1809122369225997882917535968560880d+00
       X1(341)=0.1768902376834464831804910000836637d+00
       X1(342)=0.1728652823705602892150653556387500d+00
       X1(343)=0.1688374382463567714900902165905225d+00
       X1(344)=0.1648067726215274167473553477191134d+00
       X1(345)=0.1607733528539147812095548102419475d+00
       X1(346)=0.1567372463473868504562410877675871d+00
       X1(347)=0.1526985205507106301529653457667393d+00
       X1(348)=0.1486572429564249864574173945457888d+00
       X1(349)=0.1446134810997127549389178321607606d+00
       X1(350)=0.1405673025572721368598389693771835d+00
       X1(351)=0.1365187749461874016794402793824842d+00
       X1(352)=0.1324679659227989146521980714747037d+00
       X1(353)=0.1284149431815725084039876666340251d+00
       X1(354)=0.1243597744539682173804393652517676d+00
       X1(355)=0.1203025275073083940724367601385810d+00
       X1(356)=0.1162432701436452259340572833296310d+00
       X1(357)=0.1121820701986276719182701104992017d+00
       X1(358)=0.1081189955403678375654055146224325d+00
       X1(359)=0.1040541140683068075887922987951735d+00
       X1(360)=0.9998749371207995491102589004945587d-01
       X1(361)=0.9591920243038174511307889007769378d-01
       X1(362)=0.9184930820983005526689820879242158d-01
       X1(363)=0.8777787906383002613024821168465461d-01
       X1(364)=0.8370498303143746669035745638526070d-01
       X1(365)=0.7963068817622183005040744736265877d-01
       X1(366)=0.7555506258512877966006527529167245d-01
       X1(367)=0.7147817436734236489810790959186666d-01
       X1(368)=0.6740009165314682502171416455502455d-01
       X1(369)=0.6332088259278804050321085226478460d-01
       X1(370)=0.5924061535533465078095226537123128d-01
       X1(371)=0.5515935812753885745658680124570850d-01
       X1(372)=0.5107717911269693197622115294615573d-01
       X1(373)=0.4699414652950944683793096392581819d-01
       X1(374)=0.4291032861094124937268708998980420d-01
       X1(375)=0.3882579360308119715006855643417078d-01
       X1(376)=0.3474060976400167406411686869977669d-01
       X1(377)=0.3065484536261790615835146466727391d-01
       X1(378)=0.2656856867754709625231272466110254d-01
       X1(379)=0.2248184799596739643502702539258966d-01
       X1(380)=0.1839475161247673749349778585213595d-01
       X1(381)=0.1430734782795153434671726127928119d-01
       X1(382)=0.1021970494840528655776595584357610d-01
       X1(383)=0.6131891283847092998319910930853611d-02
       X1(384)=0.2043975147140099741320754647262302d-02

	go to 56

C
C-------------------------------------------------------------------¬
C¦                       N --- 768                                  ¦
CL-------------------------------------------------------------------

C
768      NG=768
	 w(  1)=0.3143280544300424052208816662575444d-05
	 w(  2)=0.7316946032956578863012167907742587d-05
	 w(  3)=0.1149676101862034392206079095617816d-04
	 w(  4)=0.1567706472452444466259585704073367d-04
	 w(  5)=0.1985741065987670454180922553021230d-04
	 w(  6)=0.2403770585968651029362125448649260d-04
	 w(  7)=0.2821791283477775856531783853620180d-04
	 w(  8)=0.3239800730583026899006640713578990d-04
	 w(  9)=0.3657796906039584244398190130921543d-04
	 w( 10)=0.4075777937951431805179904877408765d-04
	 w( 11)=0.4493742016850046872377096968141648d-04
	 w( 12)=0.4911687362099842947804053242992730d-04
	 w( 13)=0.5329612207491258991138609235175838d-04
	 w( 14)=0.5747514794529235198792302614432616d-04
	 w( 15)=0.6165393369088629776479758590684672d-04
	 w( 16)=0.6583246179651924956589407987910883d-04
	 w( 17)=0.7001071476336385355183968864024948d-04
	 w( 18)=0.7418867510335657593333101233895482d-04
	 w( 19)=0.7836632533588686575709115336574925d-04
	 w( 20)=0.8254364798578161686958377524091682d-04
	 w( 21)=0.8672062558205289028801845070916012d-04
	 w( 22)=0.9089724065710895541350738427175353d-04
	 w( 23)=0.9507347574625413533010643929408480d-04
	 w( 24)=0.9924931338737301294015717009533744d-04
	 w( 25)=0.1034247361207348839177097825291373d-03
	 w( 26)=0.1075997264888781848155705320864356d-03
	 w( 27)=0.1117742670365490669137175039591069d-03
	 w( 28)=0.1159483403106772304064582601700427d-03
	 w( 29)=0.1201219288603777857435698102502636d-03
	 w( 30)=0.1242950152369715478801345224204701d-03
	 w( 31)=0.1284675819940185524193907234491157d-03
	 w( 32)=0.1326396116873611683673167887989347d-03
	 w( 33)=0.1368110868751742529199829765391187d-03
	 w( 34)=0.1409819901180205265866084569458080d-03
	 w( 35)=0.1451523039789098550231450551246585d-03
	 w( 36)=0.1493220110233614803866061661218102d-03
	 w( 37)=0.1534910938194684978768675748562775d-03
	 w( 38)=0.1576595349379640543953895781465469d-03
	 w( 39)=0.1618273169522888774745112335573217d-03
	 w( 40)=0.1659944224386598385144916988820898d-03
	 w( 41)=0.1701608339761393250416464060236082d-03
	 w( 42)=0.1743265341467052492293013699334018d-03
	 w( 43)=0.1784915055353215592707136093830105d-03
	 w( 44)=0.1826557307300091498868092794420082d-03
	 w( 45)=0.1868191923219170908186929057941004d-03
	 w( 46)=0.1909818729053941094212185451722315d-03
	 w( 47)=0.1951437550780602767692150333197093d-03
	 w( 48)=0.1993048214408788569880234701668598d-03
	 w( 49)=0.2034650545982282875461122256448156d-03
	 w( 50)=0.2076244371579742645365175769397026d-03
	 w( 51)=0.2117829517315419119279171147577503d-03
	 w( 52)=0.2159405809339880176885155531861651d-03
	 w( 53)=0.2200973073840733228066164830116653d-03
	 w( 54)=0.2242531137043348517262431159771565d-03
	 w( 55)=0.2284079825211582747190702945270584d-03
	 w( 56)=0.2325618964648502943290582359867539d-03
	 w( 57)=0.2367148381697110493340148771895750d-03
	 w( 58)=0.2408667902741065307315517956519949d-03
	 w( 59)=0.2450177354205410051246081746131944d-03
	 w( 60)=0.2491676562557294415924796408467345d-03
	 w( 61)=0.2533165354306699387175723165697932d-03
	 w( 62)=0.2574643556007161489200701344458866d-03
	 w( 63)=0.2616110994256496976515962683504457d-03
	 w( 64)=0.2657567495697525953301492566049252d-03
	 w( 65)=0.2699012887018796401744446093330107d-03
	 w( 66)=0.2740446994955308103262352885860330d-03
	 w( 67)=0.2781869646289236438422554287939199d-03
	 w( 68)=0.2823280667850656052996488643410248d-03
	 w( 69)=0.2864679886518264378954115596133892d-03
	 w( 70)=0.2906067129220105000358252261416418d-03
	 w( 71)=0.2947442222934290855096346512179442d-03
	 w( 72)=0.2988804994689727264217367315875138d-03
	 w( 73)=0.3030155271566834781348049742568666d-03
	 w( 74)=0.3071492880698271855265522767357458d-03
	 w( 75)=0.3112817649269657299218800506694907d-03
	 w( 76)=0.3154129404520292561033380453246622d-03
	 w( 77)=0.3195427973743883788412633547862087d-03
	 w( 78)=0.3236713184289263684176283152536017d-03
	 w( 79)=0.3277984863561113146457999407244359d-03
	 w( 80)=0.3319242839020682689127649861633284d-03
	 w( 81)=0.3360486938186513637914654808561384d-03
	 w( 82)=0.3401716988635159097891924051513972d-03
	 w( 83)=0.3442932818001904688138995147996761d-03
	 w( 84)=0.3484134253981489039541633393541364d-03
	 w( 85)=0.3525321124328824051806160435271366d-03
	 w( 86)=0.3566493256859714905872591147589828d-03
	 w( 87)=0.3607650479451579828003355738606860d-03
	 w( 88)=0.3648792620044169601905740165729783d-03
	 w( 89)=0.3689919506640286825317710151441833d-03
	 w( 90)=0.3731030967306504907549792485555268d-03
	 w( 91)=0.3772126830173886804531287840814310d-03
	 w( 92)=0.3813206923438703487958241636969629d-03
	 w( 93)=0.3854271075363152145184130181216723d-03
	 w( 94)=0.3895319114276074106532841644364708d-03
	 w( 95)=0.3936350868573672496747861860845774d-03
	 w( 96)=0.3977366166720229607322147129619028d-03
	 w( 97)=0.4018364837248823986480442862057897d-03
	 w( 98)=0.4059346708762047243610189902076279d-03
	 w( 99)=0.4100311609932720564958999115306117d-03
	 w(100)=0.4141259369504610937436273813544597d-03
	 w(101)=0.4182189816293147077374184276097423d-03
	 w(102)=0.4223102779186135061119080996854570d-03
	 w(103)=0.4263998087144473654338776218091322d-03
	 w(104)=0.4304875569202869336944104530937508d-03
	 w(105)=0.4345735054470551020534948709799979d-03
	 w(106)=0.4386576372131984455291623408268307d-03
	 w(107)=0.4427399351447586323242267269376850d-03
	 w(108)=0.4468203821754438014845809398314591d-03
	 w(109)=0.4508989612466999085838242462771552d-03
	 w(110)=0.4549756553077820391297434410508958d-03
	 w(111)=0.4590504473158256893888616817720210d-03
	 w(112)=0.4631233202359180143259064695159752d-03
	 w(113)=0.4671942570411690423556387303071122d-03
	 w(114)=0.4712632407127828566050332811942827d-03
	 w(115)=0.4753302542401287423843116457674584d-03
	 w(116)=0.4793952806208123005658052134638387d-03
	 w(117)=0.4834583028607465265700736695258147d-03
	 w(118)=0.4875193039742228546591236236004932d-03
	 w(119)=0.4915782669839821672369682566720984d-03
	 w(120)=0.4956351749212857688581431064568330d-03
	 w(121)=0.4996900108259863246451480698545146d-03
	 w(122)=0.5037427577465987628161233286065570d-03
	 w(123)=0.5077933987403711410243890007508569d-03
	 w(124)=0.5118419168733554762117864980031147d-03
	 w(125)=0.5158882952204785376780552733526929d-03
	 w(126)=0.5199325168656126030687631711662162d-03
	 w(127)=0.5239745649016461769845831091260814d-03
	 w(128)=0.5280144224305546719149743748889073d-03
	 w(129)=0.5320520725634710511995843529238052d-03
	 w(130)=0.5360874984207564337209368573258327d-03
	 w(131)=0.5401206831320706600322171994348870d-03
	 w(132)=0.5441516098364428196242023548080918d-03
	 w(133)=0.5481802616823417390356177353532917d-03
	 w(134)=0.5522066218277464305114306819490888d-03
	 w(135)=0.5562306734402165009138153794190396d-03
	 w(136)=0.5602523996969625205907449196393386d-03
	 w(137)=0.5642717837849163519073841166752214d-03
	 w(138)=0.5682888089008014371456717878850689d-03
	 w(139)=0.5723034582512030454776938994209966d-03
	 w(140)=0.5763157150526384787186595443123046d-03
	 w(141)=0.5803255625316272355655004585960466d-03
	 w(142)=0.5843329839247611340273219422350673d-03
	 w(143)=0.5883379624787743917541388700633125d-03
	 w(144)=0.5923404814506136639705351659265141d-03
	 w(145)=0.5963405241075080387210888637608324d-03
	 w(146)=0.6003380737270389891346078686027968d-03
	 w(147)=0.6043331135972102824144239189266906d-03
	 w(148)=0.6083256270165178452621941856443121d-03
	 w(149)=0.6123155972940195854428615561831824d-03
	 w(150)=0.6163030077494051691986260663011532d-03
	 w(151)=0.6202878417130657542199812691929434d-03
	 w(152)=0.6242700825261636778820706729370670d-03
	 w(153)=0.6282497135407021004548208266687754d-03
	 w(154)=0.6322267181195946029955092783948391d-03
	 w(155)=0.6362010796367347396326275412431960d-03
	 w(156)=0.6401727814770655439501014617686721d-03
	 w(157)=0.6441418070366489891811340493493364d-03
	 w(158)=0.6481081397227354019211389598941042d-03
	 w(159)=0.6520717629538328290694364852656686d-03
	 w(160)=0.6560326601597763577095881326863682d-03
	 w(161)=0.6599908147817973876384507324845720d-03
	 w(162)=0.6639462102725928562542365306076755d-03
	 w(163)=0.6678988300963944155140719436688452d-03
	 w(164)=0.6718486577290375606717546150268588d-03
	 w(165)=0.6757956766580307105066161437410759d-03
	 w(166)=0.6797398703826242387546063947548363d-03
	 w(167)=0.6836812224138794564529246664592597d-03
	 w(168)=0.6876197162747375449097332167618995d-03
	 w(169)=0.6915553355000884390106997547590573d-03
	 w(170)=0.6954880636368396605743275140367538d-03
	 w(171)=0.6994178842439851014682444557169871d-03
	 w(172)=0.7033447808926737561988370232481585d-03
	 w(173)=0.7072687371662784036868287037716979d-03
	 w(174)=0.7111897366604642379416194584945617d-03
	 w(175)=0.7151077629832574473474188814335654d-03
	 w(176)=0.7190227997551137422744237456047588d-03
	 w(177)=0.7229348306089868307285094105863587d-03
	 w(178)=0.7268438391903968417531244067993136d-03
	 w(179)=0.7307498091574986962972983903350803d-03
	 w(180)=0.7346527241811504252638955874054420d-03
	 w(181)=0.7385525679449814344524688284166638d-03
	 w(182)=0.7424493241454607161112933164989873d-03
	 w(183)=0.7463429764919650068133843916208306d-03
	 w(184)=0.7502335087068468913715297461468198d-03
	 w(185)=0.7541209045255028525075938272706467d-03
	 w(186)=0.7580051476964412659915805320543581d-03
	 w(187)=0.7618862219813503409661697672509124d-03
	 w(188)=0.7657641111551660051726740136436304d-03
	 w(189)=0.7696387990061397347945927078653937d-03
	 w(190)=0.7735102693359063286351750377382327d-03
	 w(191)=0.7773785059595516263456356439269530d-03
	 w(192)=0.7812434927056801704209027346220753d-03
	 w(193)=0.7851052134164828116800142542521490d-03
	 w(194)=0.7889636519478042579485150047802296d-03
	 w(195)=0.7928187921692105656604460016079540d-03
	 w(196)=0.7966706179640565740977568578966397d-03
	 w(197)=0.8005191132295532819852126333853020d-03
	 w(198)=0.8043642618768351661591083585059897d-03
	 w(199)=0.8082060478310274420283473535238277d-03
	 w(200)=0.8120444550313132655466835071454839d-03
	 w(201)=0.8158794674310008764151728609479130d-03
	 w(202)=0.8197110689975906822341261663269124d-03
	 w(203)=0.8235392437128422833241015405480123d-03
	 w(204)=0.8273639755728414379357249488567647d-03
	 w(205)=0.8311852485880669675683759812974926d-03
	 w(206)=0.8350030467834576021180272766035256d-03
	 w(207)=0.8388173541984787645747779718471973d-03
	 w(208)=0.8426281548871892949908747259572466d-03
	 w(209)=0.8464354329183081134402681781086798d-03
	 w(210)=0.8502391723752808216910081586536467d-03
	 w(211)=0.8540393573563462433120375708946619d-03
	 w(212)=0.8578359719746029019362026067204177d-03
	 w(213)=0.8616290003580754374015558479752984d-03
	 w(214)=0.8654184266497809594932888383805659d-03
	 w(215)=0.8692042350077953390088918877751430d-03
	 w(216)=0.8729864096053194358694011912300056d-03
	 w(217)=0.8767649346307452639998568099907002d-03
	 w(218)=0.8805397942877220927023596689391653d-03
	 w(219)=0.8843109727952224842453814760033837d-03
	 w(220)=0.8880784543876082673932483623002194d-03
	 w(221)=0.8918422233146964465999870773405596d-03
	 w(222)=0.8956022638418250465919917508829358d-03
	 w(223)=0.8993585602499188920642395514726311d-03
	 w(224)=0.9031110968355553222150550307914928d-03
	 w(225)=0.9068598579110298398446955420731718d-03
	 w(226)=0.9106048278044216947433038593787087d-03
	 w(227)=0.9143459908596594010940490018135081d-03
	 w(228)=0.9180833314365861886175522821005043d-03
	 w(229)=0.9218168339110253871839727515819446d-03
	 w(230)=0.9255464826748457446194045029431023d-03
	 w(231)=0.9292722621360266774335177169560066d-03
	 w(232)=0.9329941567187234541956558995188610d-03
	 w(233)=0.9367121508633323112868834493834321d-03
	 w(234)=0.9404262290265555007557605243604895d-03
	 w(235)=0.9441363756814662700059061335943353d-03
	 w(236)=0.9478425753175737730436954747969215d-03
	 w(237)=0.9515448124408879130147237572130483d-03
	 w(238)=0.9552430715739841157579561026045878d-03
	 w(239)=0.9589373372560680341067715967388044d-03
	 w(240)=0.9626275940430401826663991717625925d-03
	 w(241)=0.9663138265075605027975337344482792d-03
	 w(242)=0.9699960192391128575362128155976644d-03
	 w(243)=0.9736741568440694561803270008598454d-03
	 w(244)=0.9773482239457552082734315118176963d-03
	 w(245)=0.9810182051845120067168215373689662d-03
	 w(246)=0.9846840852177629397411302681032761d-03
	 w(247)=0.9883458487200764314690060594709290d-03
	 w(248)=0.9920034803832303108007237419608761d-03
	 w(249)=0.9956569649162758083548848071410413d-03
	 w(250)=0.9993062870456014811966620261490724d-03
	 w(251)=0.1002951431514997065086346000920292d-02
	 w(252)=0.1006592383085717253981254206961438d-02
	 w(253)=0.1010229126536545406524367258669418d-02
	 w(254)=0.1013861646663857179253362412891089d-02
	 w(255)=0.1017489928281684086264020722445945d-02
	 w(256)=0.1021113956221776985062291757509600d-02
	 w(257)=0.1024733715333669488339608427882921d-02
	 w(258)=0.1028349190484741301406354162052435d-02
	 w(259)=0.1031960366560281485018695528365822d-02
	 w(260)=0.1035567228463551643334305318385575d-02
	 w(261)=0.1039169761115849036732814151311535d-02
	 w(262)=0.1042767949456569619237142800044923d-02
	 w(263)=0.1046361778443271000272182682754918d-02
	 w(264)=0.1049951233051735330497608307451694d-02
	 w(265)=0.1053536298276032111451922899805609d-02
	 w(266)=0.1057116959128580928745156985011717d-02
	 w(267)=0.1060693200640214108537959331607377d-02
	 w(268)=0.1064265007860239297045140397522391d-02
	 w(269)=0.1067832365856501962802050244989076d-02
	 w(270)=0.1071395259715447821432495809961483d-02
	 w(271)=0.1074953674542185182657226422089871d-02
	 w(272)=0.1078507595460547219282341571752565d-02
	 w(273)=0.1082057007613154157907301109847314d-02
	 w(274)=0.1085601896161475391092545342663913d-02
	 w(275)=0.1089142246285891510727060846867143d-02
	 w(276)=0.1092678043185756262336557277078829d-02
	 w(277)=0.1096209272079458420073250969414884d-02
	 w(278)=0.1099735918204483582128582757260817d-02
	 w(279)=0.1103257966817475886310530109199084d-02
	 w(280)=0.1106775403194299645527507471975379d-02
	 w(281)=0.1110288212630100902921183552339488d-02
	 w(282)=0.1113796380439368906390880199147034d-02
	 w(283)=0.1117299891955997502252554549883887d-02
	 w(284)=0.1120798732533346447775704182389205d-02
	 w(285)=0.1124292887544302642341874161616277d-02
	 w(286)=0.1127782342381341276968785091389522d-02
	 w(287)=0.1131267082456586901944442570885277d-02
	 w(288)=0.1134747093201874412315930813582720d-02
	 w(289)=0.1138222360068809950977936611283228d-02
	 w(290)=0.1141692868528831729106394316066102d-02
	 w(291)=0.1145158604073270763682988067314589d-02
	 w(292)=0.1148619552213411531856594107778532d-02
	 w(293)=0.1152075698480552541888093710608508d-02
	 w(294)=0.1155527028426066820425335976959912d-02
	 w(295)=0.1158973527621462315855379559685293d-02
	 w(296)=0.1162415181658442217481493221356198d-02
	 w(297)=0.1165851976148965190272747042932456d-02
	 w(298)=0.1169283896725305524934379061367904d-02
	 w(299)=0.1172710929040113203047476128842699d-02
	 w(300)=0.1176133058766473877026862851676919d-02
	 w(301)=0.1179550271597968764646448581834867d-02
	 w(302)=0.1182962553248734457881639596795540d-02
	 w(303)=0.1186369889453522645818781812960729d-02
	 w(304)=0.1189772265967759751381958632209200d-02
	 w(305)=0.1193169668567606481627828819192628d-02
	 w(306)=0.1196562083050017291359550647007338d-02
	 w(307)=0.1199949495232799759811200929466009d-02
	 w(308)=0.1203331890954673880154460977825873d-02
	 w(309)=0.1206709256075331261579705976996576d-02
	 w(310)=0.1210081576475494243703999769432363d-02
	 w(311)=0.1213448838056974923058863562592309d-02
	 w(312)=0.1216811026742734091411054636501098d-02
	 w(313)=0.1220168128476940085669959720033205d-02
	 w(314)=0.1223520129225027549135577326539692d-02
	 w(315)=0.1226867014973756103841433989801313d-02
	 w(316)=0.1230208771731268933747151018479773d-02
	 w(317)=0.1233545385527151278535751089702820d-02
	 w(318)=0.1236876842412488837771167728606143d-02
	 w(319)=0.1240203128459926085171795469008869d-02
	 w(320)=0.1243524229763724492756294259357117d-02
	 w(321)=0.1246840132439820664618238466066855d-02
	 w(322)=0.1250150822625884380086578631861240d-02
	 w(323)=0.1253456286481376546029262968056070d-02
	 w(324)=0.1256756510187607058057745395416920d-02
	 w(325)=0.1260051479947792570390487796613015d-02
	 w(326)=0.1263341181987114174133946002835363d-02
	 w(327)=0.1266625602552774983739911906237681d-02
	 w(328)=0.1269904727914057631398467966902516d-02
	 w(329)=0.1273178544362381669126195266428875d-02
	 w(330)=0.1276447038211360878309662148377155d-02
	 w(331)=0.1279710195796860486464607377081543d-02
	 w(332)=0.1282968003477054290971619639134522d-02
	 w(333)=0.1286220447632481689549504104546568d-02
	 w(334)=0.1289467514666104617227916655559397d-02
	 w(335)=0.1292709191003364389581237278720142d-02
	 w(336)=0.1295945463092238451986045998473160d-02
	 w(337)=0.1299176317403297034664957606559251d-02
	 w(338)=0.1302401740429759713279965309291569d-02
	 w(339)=0.1305621718687551874838838272656258d-02
	 w(340)=0.1308836238715361088678513891517445d-02
	 w(341)=0.1312045287074693382289822442337142d-02
	 w(342)=0.1315248850349929421748279597094015d-02
	 w(343)=0.1318446915148380596516081077839979d-02
	 w(344)=0.1321639468100345008380833514905032d-02
	 w(345)=0.1324826495859163364296956335478095d-02
	 w(346)=0.1328007985101274772896091251483875d-02
	 w(347)=0.1331183922526272444433258634661227d-02
	 w(348)=0.1334354294856959293935903760849376d-02
	 w(349)=0.1337519088839403447323380573015566d-02
	 w(350)=0.1340678291242993650264826251822567d-02
	 w(351)=0.1343831888860494579543786491841785d-02
	 w(352)=0.1346979868508102056698358959170437d-02
	 w(353)=0.1350122217025498163706030950505470d-02
	 w(354)=0.1353258921275906260482796782955385d-02
	 w(355)=0.1356389968146145903966550916327258d-02
	 w(356)=0.1359515344546687668555164243588528d-02
	 w(357)=0.1362635037411707867670063378960048d-02
	 w(358)=0.1365749033699143176216546124918974d-02
	 w(359)=0.1368857320390745153712480607555909d-02
	 w(360)=0.1371959884492134667857450832505386d-02
	 w(361)=0.1375056713032856218314827629318891d-02
	 w(362)=0.1378147793066432160479661118936797d-02
	 w(363)=0.1381233111670416829005708955094944d-02
	 w(364)=0.1384312655946450560865333654328416d-02
	 w(365)=0.1387386413020313617716422338955382d-02
	 w(366)=0.1390454370041980007350903171286285d-02
	 w(367)=0.1393516514185671203999854653545644d-02
	 w(368)=0.1396572832649909767270626804854858d-02
	 w(369)=0.1399623312657572859491817002336475d-02
	 w(370)=0.1402667941455945661242367986192741d-02
	 w(371)=0.1405706706316774684841481176709664d-02
	 w(372)=0.1408739594536320985576465032762437d-02
	 w(373)=0.1411766593435413270446065694767341d-02
	 w(374)=0.1414787690359500904197255598352033d-02
	 w(375)=0.1417802872678706812433885116509500d-02
	 w(376)=0.1420812127787880281576032585867542d-02
	 w(377)=0.1423815443106649655449319295146203d-02
	 w(378)=0.1426812806079474928283888159088532d-02
	 w(379)=0.1429804204175700233903177867328688d-02
	 w(380)=0.1432789624889606230883058282996614d-02
	 w(381)=0.1435769055740462383462327768535815d-02
	 w(382)=0.1438742484272579137986008934412913d-02
	 w(383)=0.1441709898055359994663316038302632d-02
	 w(384)=0.1444671284683353474422604906115582d-02
	 w(385)=0.1447626631776304980646054800067249d-02
	 w(386)=0.1450575926979208555567271121035189d-02
	 w(387)=0.1453519157962358531115438200876899d-02
	 w(388)=0.1456456312421401073990092713347055d-02
	 w(389)=0.1459387378077385624751030407910422d-02
	 w(390)=0.1462312342676816230708301947252577d-02
	 w(391)=0.1465231193991702772397697604788352d-02
	 w(392)=0.1468143919819612083427565451105788d-02
	 w(393)=0.1471050507983718963483253426197311d-02
	 w(394)=0.1473950946332857084275912355662271d-02
	 w(395)=0.1476845222741569788222844521943303d-02
	 w(396)=0.1479733325110160779647030844216672d-02
	 w(397)=0.1482615241364744708283919050916627d-02
	 w(398)=0.1485490959457297644884005445159810d-02
	 w(399)=0.1488360467365707448700193963664422d-02
	 w(400)=0.1491223753093824026649368212244727d-02
	 w(401)=0.1494080804671509483938065023717253d-02
	 w(402)=0.1496931610154688165942591825184460d-02
	 w(403)=0.1499776157625396591134384719271045d-02
	 w(404)=0.1502614435191833274841859676075196d-02
	 w(405)=0.1505446430988408443640465597459235d-02
	 w(406)=0.1508272133175793640163105250932701d-02
	 w(407)=0.1511091529940971218123548174864960d-02
	 w(408)=0.1513904609497283727345918628189376d-02
	 w(409)=0.1516711360084483188593801494207904d-02
	 w(410)=0.1519511769968780257992969747651486d-02
	 w(411)=0.1522305827442893280842198654872679d-02
	 w(412)=0.1525093520826097234607094297011995d-02
	 w(413)=0.1527874838464272560892327283256509d-02
	 w(414)=0.1530649768729953886188126653961954d-02
	 w(415)=0.1533418300022378631187353959495566d-02
	 w(416)=0.1536180420767535508469943338236018d-02
	 w(417)=0.1538936119418212908351960104288340d-02
	 w(418)=0.1541685384454047172696997890186700d-02
	 w(419)=0.1544428204381570756488102770211807d-02
	 w(420)=0.1547164567734260276958882013983389d-02
	 w(421)=0.1549894463072584450082925185742123d-02
	 w(422)=0.1552617878984051914221136210242232d-02
	 w(423)=0.1555334804083258940727046769469427d-02
	 w(424)=0.1558045227011937031310653973505587d-02
	 w(425)=0.1560749136439000401961798661805666d-02
	 w(426)=0.1563446521060593353234574935956745d-02
	 w(427)=0.1566137369600137526694736599669236d-02
	 w(428)=0.1568821670808379047332542084320576d-02
	 w(429)=0.1571499413463435551743956166844685d-02
	 w(430)=0.1574170586370843101883604339139676d-02
	 w(431)=0.1576835178363602984193354062459572d-02
	 w(432)=0.1579493178302228393910876334458945d-02
	 w(433)=0.1582144575074791004363021008673666d-02
	 w(434)=0.1584789357596967421049320134235763d-02
	 w(435)=0.1587427514812085520321415225528622d-02
	 w(436)=0.1590059035691170672464686825275288d-02
	 w(437)=0.1592683909232991848988847988201110d-02
	 w(438)=0.1595302124464107613934747383902524d-02
	 w(439)=0.1597913670438911999005112594862182d-02
	 w(440)=0.1600518536239680262327449866650152d-02
	 w(441)=0.1603116710976614530657803050211391d-02
	 w(442)=0.1605708183787889324834561758727661d-02
	 w(443)=0.1608292943839696968291996841818648d-02
	 w(444)=0.1610870980326292878443690155775586d-02
	 w(445)=0.1613442282470040740746515278052937d-02
	 w(446)=0.1616006839521457565256316275336169d-02
	 w(447)=0.1618564640759258625486922885104974d-02
	 w(448)=0.1621115675490402279384632508665851d-02
	 w(449)=0.1623659933050134672230782237083192d-02
	 w(450)=0.1626197402802034321285527738228220d-02
	 w(451)=0.1628728074138056581986440221231976d-02
	 w(452)=0.1631251936478577995516027861901260d-02
	 w(453)=0.1633768979272440517552784017068063d-02
	 w(454)=0.1636279191996995628020861275317231d-02
	 w(455)=0.1638782564158148321653967884000969d-02
	 w(456)=0.1641279085290400979189581355818809d-02
	 w(457)=0.1643768744956897119010073090438078d-02
	 w(458)=0.1646251532749465029047837645564064d-02
	 w(459)=0.1648727438288661278772020855453641d-02
	 w(460)=0.1651196451223814111074942321006623d-02
	 w(461)=0.1653658561233066713876809882170632d-02
	 w(462)=0.1656113758023420371267826528360697d-02
	 w(463)=0.1658562031330777494007293803816709d-02
	 w(464)=0.1661003370919984529199820120202078d-02
	 w(465)=0.1663437766584874748969247496170765d-02
	 w(466)=0.1665865208148310917951416100989481d-02
	 w(467)=0.1668285685462227839427392584481022d-02
	 w(468)=0.1670699188407674779919295526436030d-02
	 w(469)=0.1673105706894857772071359433102325d-02
	 w(470)=0.1675505230863181795639387544280429d-02
	 w(471)=0.1677897750281292836412253289802816d-02
	 w(472)=0.1680283255147119822889620546622531d-02
	 w(473)=0.1682661735487916440540563895251551d-02
	 w(474)=0.1685033181360302823468281854733268d-02
	 w(475)=0.1687397582850307123306608586567758d-02
	 w(476)=0.1689754930073406955174542797890583d-02
	 w(477)=0.1692105213174570720515526540589780d-02
	 w(478)=0.1694448422328298806648721293782582d-02
	 w(479)=0.1696784547738664662860044130013229d-02
	 w(480)=0.1699113579639355752861242898518213d-02
	 w(481)=0.1701435508293714383445806209781665d-02
	 w(482)=0.1703750323994778409171021572206456d-02
	 w(483)=0.1706058017065321812896013311896252d-02
	 w(484)=0.1708358577857895162006110898110231d-02
	 w(485)=0.1710651996754865940154417997748101d-02
	 w(486)=0.1712938264168458754351972990074240d-02
	 w(487)=0.1715217370540795417238412785622411d-02
	 w(488)=0.1717489306343934904365573608655245d-02
	 w(489)=0.1719754062079913186326984918507905d-02
	 w(490)=0.1722011628280782935566735859434619d-02
	 w(491)=0.1724261995508653107701717539016569d-02
	 w(492)=0.1726505154355728397191769039585783d-02
	 w(493)=0.1728741095444348567192780363283365d-02
	 w(494)=0.1730969809427027653428331497102048d-02
	 w(495)=0.1733191286986493041915973457365156d-02
	 w(496)=0.1735405518835724420384784531365434d-02
	 w(497)=0.1737612495717992603221362975122549d-02
	 w(498)=0.1739812208406898229781946148209176d-02
	 w(499)=0.1742004647706410335908875467133547d-02
	 w(500)=0.1744189804450904798490156635637408d-02
	 w(501)=0.1746367669505202652901395361256027d-02
	 w(502)=0.1748538233764608283169920190372576d-02
	 w(503)=0.1750701488154947484701436186562292d-02
	 w(504)=0.1752857423632605399410085937036175d-02
	 w(505)=0.1755006031184564323093327797233566d-02
	 w(506)=0.1757147301828441384893575371847417d-02
	 w(507)=0.1759281226612526098689076979563028d-02
	 w(508)=0.1761407796615817786257049256313142d-02
	 w(509)=0.1763527002948062872052615115665346d-02
	 w(510)=0.1765638836749792049447633002814433d-02
	 w(511)=0.1767743289192357318274041748316735d-02
	 w(512)=0.1769840351477968893516883346920205d-02
	 w(513)=0.1771930014839731985002704653373357d-02
	 w(514)=0.1774012270541683447929578298676875d-02
	 w(515)=0.1776087109878828304085523084626996d-02
	 w(516)=0.1778154524177176133602644709427361d-02
	 w(517)=0.1780214504793777337094858910357465d-02
	 w(518)=0.1782267043116759268027600978720069d-02
	 w(519)=0.1784312130565362235168468105277846d-02
	 w(520)=0.1786349758589975374968284148868014d-02
	 w(521)=0.1788379918672172393722620184577334d-02
	 w(522)=0.1790402602324747179364348578499769d-02
	 w(523)=0.1792417801091749282738353351405101d-02
	 w(524)=0.1794425506548519268210065231345582d-02
	 w(525)=0.1796425710301723933460036053031936d-02
	 w(526)=0.1798418403989391398317314037442608d-02
	 w(527)=0.1800403579280946062484928976299755d-02
	 w(528)=0.1802381227877243432011344450465369d-02
	 w(529)=0.1804351341510604814362282926690518d-02
	 w(530)=0.1806313911944851881947878901194399d-02
	 w(531)=0.1808268930975341103960665188961717d-02
	 w(532)=0.1810216390428998046380447992128883d-02
	 w(533)=0.1812156282164351540002677517078410d-02
	 w(534)=0.1814088598071567716347472645574391d-02
	 w(535)=0.1816013330072483911307010498141724d-02
	 w(536)=0.1817930470120642436389544655610724d-02
	 w(537)=0.1819840010201324217418869325002852d-02
	 w(538)=0.1821741942331582300548600846411911d-02
	 w(539)=0.1823636258560275225451202635917782d-02
	 w(540)=0.1825522950968100265542234943541272d-02
	 w(541)=0.1827402011667626535100866672484307d-02
	 w(542)=0.1829273432803327963148242954079257d-02
	 w(543)=0.1831137206551616133945859199664264d-02
	 w(544)=0.1832993325120872993976649953682984d-02
	 w(545)=0.1834841780751483425272059049344440d-02
	 w(546)=0.1836682565715867684948916316837070d-02
	 w(547)=0.1838515672318513710820505412038108d-02
	 w(548)=0.1840341092896009292946767218552701d-02
	 w(549)=0.1842158819817074110989143725419481d-02
	 w(550)=0.1843968845482591637236128294585339d-02
	 w(551)=0.1845771162325640905166149803938342d-02
	 w(552)=0.1847565762811528143414980280944819d-02
	 w(553)=0.1849352639437818275015418326415531d-02
	 w(554)=0.1851131784734366281777563865274419d-02
	 w(555)=0.1852903191263348433678563549066369d-02
	 w(556)=0.1854666851619293383131270470960305d-02
	 w(557)=0.1856422758429113124001826735824274d-02
	 w(558)=0.1858170904352133815246742853203905d-02
	 w(559)=0.1859911282080126469040613887365770d-02
	 w(560)=0.1861643884337337503266178803602566d-02
	 w(561)=0.1863368703880519158238996491373054d-02
	 w(562)=0.1865085733498959777539579520191797d-02
	 w(563)=0.1866794966014513952826394791123807d-02
	 w(564)=0.1868496394281632532503708882897184d-02
	 w(565)=0.1870190011187392494118825054650051d-02
	 w(566)=0.1871875809651526680363828554794211d-02
	 w(567)=0.1873553782626453398557527095026170d-02
	 w(568)=0.1875223923097305883483844077763507d-02
	 w(569)=0.1876886224081961623463493411843657d-02
	 w(570)=0.1878540678631071549536336512805852d-02
	 w(571)=0.1880187279828089087632394358095221d-02
	 w(572)=0.1881826020789299073610060251686529d-02
	 w(573)=0.1883456894663846531040632244532738d-02
	 w(574)=0.1885079894633765311618857954498957d-02
	 w(575)=0.1886695013914006598079758829652356d-02
	 w(576)=0.1888302245752467269502575699536585d-02
	 w(577)=0.1889901583430018128883252757966458d-02
	 w(578)=0.1891493020260531992857452915528869d-02
	 w(579)=0.1893076549590911643456673747959865d-02
	 w(580)=0.1894652164801117641780610045479736d-02
	 w(581)=0.1896219859304196003469486235594847d-02
	 w(582)=0.1897779626546305735860659705401891d-02
	 w(583)=0.1899331460006746236714374286645785d-02
	 w(584)=0.1900875353197984554394121885265951d-02
	 w(585)=0.1902411299665682509387649434498547d-02
	 w(586)=0.1903939292988723677055228024364241d-02
	 w(587)=0.1905459326779240231492381208137924d-02
	 w(588)=0.1906971394682639650394850105742257d-02
	 w(589)=0.1908475490377631280814154012506866d-02
	 w(590)=0.1909971607576252765692686776954755d-02
	 w(591)=0.1911459740023896331067871230792906d-02
	 w(592)=0.1912939881499334933835476435654537d-02
	 w(593)=0.1914412025814748269962785451937457d-02
	 w(594)=0.1915876166815748643042884732863558d-02
	 w(595)=0.1917332298381406693081930099216195d-02
	 w(596)=0.1918780414424276985411828554647171d-02
	 w(597)=0.1920220508890423459621359955546445d-02
	 w(598)=0.1921652575759444738399347750787823d-02
	 w(599)=0.1923076609044499296184073652756486d-02
	 w(600)=0.1924492602792330487513717189482872d-02
	 w(601)=0.1925900551083291434973187615998508d-02
	 w(602)=0.1927300448031369776633302628745368d-02
	 w(603)=0.1928692287784212272878855727552719d-02
	 w(604)=0.1930076064523149272522701902892523d-02
	 w(605)=0.1931451772463219038103579589375682d-02
	 w(606)=0.1932819405853191930265975517300065d-02
	 w(607)=0.1934178958975594451120928210044197d-02
	 w(608)=0.1935530426146733146487255413756010d-02
	 w(609)=0.1936873801716718366913280704649846d-02
	 w(610)=0.1938209080069487887379724895829078d-02
	 w(611)=0.1939536255622830385585018657430167d-02
	 w(612)=0.1940855322828408778714883968565364d-02
	 w(613)=0.1942166276171783418598623634554048d-02
	 w(614)=0.1943469110172435145155150125805611d-02
	 w(615)=0.1944763819383788198032377422971544d-02
	 w(616)=0.1946050398393232986344192384146986d-02
	 w(617)=0.1947328841822148716409815381493074d-02
	 w(618)=0.1948599144325925877400953584191357d-02
	 w(619)=0.1949861300593988584802744289647145d-02
	 w(620)=0.1951115305349816781595080122849042d-02
	 w(621)=0.1952361153350968297061502732279984d-02
	 w(622)=0.1953598839389100763133446807275383d-02
	 w(623)=0.1954828358289993388178211823747384d-02
	 w(624)=0.1956049704913568588139634890251980d-02
	 w(625)=0.1957262874153913474941034411975036d-02
	 w(626)=0.1958467860939301202060591013862529d-02
	 w(627)=0.1959664660232212167189929263323186d-02
	 w(628)=0.1960853267029355071887261205193875d-02
	 w(629)=0.1962033676361687838137050564480275d-02
	 w(630)=0.1963205883294438381728754683268572d-02
	 w(631)=0.1964369882927125242367799834648505d-02
	 w(632)=0.1965525670393578070432544495989579d-02
	 w(633)=0.1966673240861957970291584463968250d-02
	 w(634)=0.1967812589534777700096353351849296d-02
	 w(635)=0.1968943711648921727964572023170538d-02
	 w(636)=0.1970066602475666144470700882660315d-02
	 w(637)=0.1971181257320698431360149662420634d-02
	 w(638)=0.1972287671524137086404600406624684d-02
	 w(639)=0.1973385840460551104316400768692580d-02
	 w(640)=0.1974475759538979313640586488609243d-02
	 w(641)=0.1975557424202949569543694012218483d-02
	 w(642)=0.1976630829930497802419126646449634d-02
	 w(643)=0.1977695972234186922229440411989101d-02
	 w(644)=0.1978752846661125578506518855379119d-02
	 w(645)=0.1979801448792986775931209513389852d-02
	 w(646)=0.1980841774246026345414598481242898d-02
	 w(647)=0.1981873818671101270603703620343831d-02
	 w(648)=0.1982897577753687869734971348080633d-02
	 w(649)=0.1983913047213899832759566679438713d-02
	 w(650)=0.1984920222806506113665051235141629d-02
	 w(651)=0.1985919100320948677918649291222424d-02
	 w(652)=0.1986909675581360104957907617831551d-02
	 w(653)=0.1987891944446581045655160838161206d-02
	 w(654)=0.1988865902810177534682820329081476d-02
	 w(655)=0.1989831546600458157707111280902422d-02
	 w(656)=0.1990788871780491073338489432067752d-02
	 w(657)=0.1991737874348120889767576193005963d-02
	 w(658)=0.1992678550335985396016058369282596d-02
	 w(659)=0.1993610895811532147732606485066891d-02
	 w(660)=0.1994534906877034907464473791210147d-02
	 w(661)=0.1995450579669609939336046415389035d-02
	 w(662)=0.1996357910361232158066223772248704d-02
	 w(663)=0.1997256895158751132257117296749535d-02
	 w(664)=0.1998147530303906941887164791423023d-02
	 w(665)=0.1999029812073345889942367185439622d-02
	 w(666)=0.1999903736778636068119964287727679d-02
	 w(667)=0.2000769300766282776539476175314710d-02
	 w(668)=0.2001626500417743797396647189035665d-02
	 w(669)=0.2002475332149444522496440109219769d-02
	 w(670)=0.2003315792412792934601838951371206d-02
	 w(671)=0.2004147877694194442535829953649612d-02
	 w(672)=0.2004971584515066569974541721575638d-02
	 w(673)=0.2005786909431853497870137148281921d-02
	 w(674)=0.2006593849036040460442661638241036d-02
	 w(675)=0.2007392399954167994680664326172850d-02
	 w(676)=0.2008182558847846043291021398204129d-02
	 w(677)=0.2008964322413767911039003286763906d-02
	 w(678)=0.2009737687383724074420240421586588d-02
	 w(679)=0.2010502650524615844606855373999600d-02
	 w(680)=0.2011259208638468883610642627830024d-02
	 w(681)=0.2012007358562446573606790845210674d-02
	 w(682)=0.2012747097168863239362256366734896d-02
	 w(683)=0.2013478421365197223713510790235656d-02
	 w(684)=0.2014201328094103816038999808379240d-02
	 w(685)=0.2014915814333428033672265049700622d-02
	 w(686)=0.2015621877096217256202295458096320d-02
	 w(687)=0.2016319513430733712608289759559902d-02
	 w(688)=0.2017008720420466821176626799527191d-02
	 w(689)=0.2017689495184145382148455987017273d-02
	 w(690)=0.2018361834875749623046935750240721d-02
	 w(691)=0.2019025736684523096633763789924202d-02
	 w(692)=0.2019681197834984431445259008694051d-02
	 w(693)=0.2020328215586938934858871294898562d-02
	 w(694)=0.2020966787235490048641611844648633d-02
	 w(695)=0.2021596910111050656932513414045861d-02
	 w(696)=0.2022218581579354246611846801965386d-02
	 w(697)=0.2022831799041465920010436969789063d-02
	 w(698)=0.2023436559933793259913039505565282d-02
	 w(699)=0.2024032861728097046810355633619654d-02
	 w(700)=0.2024620701931501828354881654080483d-02
	 w(701)=0.2025200078086506340976406567525871d-02
	 w(702)=0.2025770987770993783613589695428391d-02
	 w(703)=0.2026333428598241943518668344678876d-02
	 w(704)=0.2026887398216933174092963981633564d-02
	 w(705)=0.2027432894311164224711473975260332d-02
	 w(706)=0.2027969914600455922495454737474012d-02
	 w(707)=0.2028498456839762705992521029062659d-02
	 w(708)=0.2029018518819482010724405309127151d-02
	 w(709)=0.2029530098365463506563140282098282d-02
	 w(710)=0.2030033193339018186897047236568137d-02
	 w(711)=0.2030527801636927309548532371789838d-02
	 w(712)=0.2031013921191451189406313068166095d-02
	 w(713)=0.2031491549970337842735315974777759d-02
	 w(714)=0.2031960685976831483128108857400375d-02
	 w(715)=0.2032421327249680869062348371934891d-02
	 w(716)=0.2032873471863147503029346298137051d-02
	 w(717)=0.2033317117927013682199477284382566d-02
	 w(718)=0.2033752263586590400590771813353144d-02
	 w(719)=0.2034178907022725102707658898378498d-02
	 w(720)=0.2034597046451809288617443958126552d-02
	 w(721)=0.2035006680125785970432728390803325d-02
	 w(722)=0.2035407806332156980168598575404474d-02
	 w(723)=0.2035800423393990128944033364263347d-02
	 w(724)=0.2036184529669926217497600595558203d-02
	 w(725)=0.2036560123554185897988134743985688d-02
	 w(726)=0.2036927203476576387051709539871092d-02
	 w(727)=0.2037285767902498030086841218977452d-02
	 w(728)=0.2037635815332950716740480014589067d-02
	 w(729)=0.2037977344304540147567969567484225d-02
	 w(730)=0.2038310353389483951840776105574949d-02
	 w(731)=0.2038634841195617656476411530677539d-02
	 w(732)=0.2038950806366400506065596942484819d-02
	 w(733)=0.2039258247580921133972335626737922d-02
	 w(734)=0.2039557163553903084483187133239263d-02
	 w(735)=0.2039847553035710185982656767106860d-02
	 w(736)=0.2040129414812351775132237610938522d-02
	 w(737)=0.2040402747705487772031265083732248d-02
	 w(738)=0.2040667550572433606338367021888240d-02
	 w(739)=0.2040923822306164994332915335797499d-02
	 w(740)=0.2041171561835322566896508449795025d-02
	 w(741)=0.2041410768124216348395136971018206d-02
	 w(742)=0.2041641440172830086443308351356821d-02
	 w(743)=0.2041863577016825432532029703604752d-02
	 w(744)=0.2042077177727545973503171405519100d-02
	 w(745)=0.2042282241412021113853357671153313d-02
	 w(746)=0.2042478767212969808851153884948573d-02
	 w(747)=0.2042666754308804148451944178039979d-02
	 w(748)=0.2042846201913632791995516475446958d-02
	 w(749)=0.2043017109277264253671996054668283d-02
	 w(750)=0.2043179475685210038742392528081747d-02
	 w(751)=0.2043333300458687630500649090847900d-02
	 w(752)=0.2043478582954623327964706860127915d-02
	 w(753)=0.2043615322565654934284721167742577d-02
	 w(754)=0.2043743518720134295857190754306256d-02
	 w(755)=0.2043863170882129692134384945766957d-02
	 w(756)=0.2043974278551428076119078070553722d-02
	 w(757)=0.2044076841263537165535224594570860d-02
	 w(758)=0.2044170858589687384665832709474214d-02
	 w(759)=0.2044256330136833656849918404408315d-02
	 w(760)=0.2044333255547657047631046380062331d-02
	 w(761)=0.2044401634500566258550588523913147d-02
	 w(762)=0.2044461466709698971579455054248134d-02
	 w(763)=0.2044512751924923044182677855393758d-02
	 w(764)=0.2044555489931837555011849964905532d-02
	 w(765)=0.2044589680551773700221049632691517d-02
	 w(766)=0.2044615323641795540402501849531627d-02
	 w(767)=0.2044632419094700598138854735613504d-02
	 w(768)=0.2044640966839020306169572685914276d-02
c
	 X1(  1)=0.9999987751809603899785742468541436d+00
	 X1(  2)=0.9999935465070646777485727571466477d+00
	 X1(  3)=0.9999841397443226880267420764283078d+00
	 X1(  4)=0.9999705528427260376722561075471150d+00
	 X1(  5)=0.9999527856031551597830102724983022d+00
	 X1(  6)=0.9999308380387709598643117631376301d+00
	 X1(  7)=0.9999047102209633149676913285317450d+00
	 X1(  8)=0.9998744022506416737002235754981086d+00
	 X1(  9)=0.9998399142506053104694209753913120d+00
	 X1( 10)=0.9998012463630057165406049312938967d+00
	 X1( 11)=0.9997583987483593112386259721101727d+00
	 X1( 12)=0.9997113715851157245414839686993341d+00
	 X1( 13)=0.9996601650694509926115168643979200d+00
	 X1( 14)=0.9996047794151608771814692598375738d+00
	 X1( 15)=0.9995452148536022832206941391298502d+00
	 X1( 16)=0.9994814716336592680862944449344103d+00
	 X1( 17)=0.9994135500217222909268162052713949d+00
	 X1( 18)=0.9993414503016749055204667230968016d+00
	 X1( 19)=0.9992651727748847913692946272645427d+00
	 X1( 20)=0.9991847177601973895479985589514096d+00
	 X1( 21)=0.9991000855939311399951419211024310d+00
	 X1( 22)=0.9990112766298737207791142567244880d+00
	 X1( 23)=0.9989182912392789208944857855116902d+00
	 X1( 24)=0.9988211298108639143363785356927242d+00
	 X1( 25)=0.9987197927508067856802791690549741d+00
	 X1( 26)=0.9986142804827442085738844160963738d+00
	 X1( 27)=0.9985045934477692110115212022625303d+00
	 X1( 28)=0.9983907321044289822713474936368031d+00
	 X1( 29)=0.9982726969287226902440775187314026d+00
	 X1( 30)=0.9981504884140992871653405510432951d+00
	 X1( 31)=0.9980241070714552880842730159361264d+00
	 X1( 32)=0.9978935534291325107664296875769969d+00
	 X1( 33)=0.9977588280329157687847211136486779d+00
	 X1( 34)=0.9976199314460305117175354126435264d+00
	 X1( 35)=0.9974768642491404079256519202678464d+00
	 X1( 36)=0.9973296270403448665046288046012697d+00
	 X1( 37)=0.9971782204351764958330037199918124d+00
	 X1( 38)=0.9970226450665984967454075087698327d+00
	 X1( 39)=0.9968629015850019888136760676238474d+00
	 X1( 40)=0.9966989906582032685604836203119653d+00
	 X1( 41)=0.9965309129714409986888771137886647d+00
	 X1( 42)=0.9963586692273733276088347181208512d+00
	 X1( 43)=0.9961822601460749386941226296746120d+00
	 X1( 44)=0.9960016864650340288206018377029133d+00
	 X1( 45)=0.9958169489391492158290612132290808d+00
	 X1( 46)=0.9956280483407263746277901639783863d+00
	 X1( 47)=0.9954349854594754017070544340358479d+00
	 X1( 48)=0.9952377611025069078828607293134795d+00
	 X1( 49)=0.9950363760943288391235063943814146d+00
	 X1( 50)=0.9948308312768430253414018549992920d+00
	 X1( 51)=0.9946211275093416570560524435768072d+00
	 X1( 52)=0.9944072656685036898530687986110986d+00
	 X1( 53)=0.9941892466483911765795527228510410d+00
	 X1( 54)=0.9939670713604455272288882296456670d+00
	 X1( 55)=0.9937407407334836964784114607892021d+00
	 X1( 56)=0.9935102557136942988520756809222647d+00
	 X1( 57)=0.9932756172646336514874144856574001d+00
	 X1( 58)=0.9930368263672217444921123574802191d+00
	 X1( 59)=0.9927938840197381388805358373050146d+00
	 X1( 60)=0.9925467912378177920848363231043380d+00
	 X1( 61)=0.9922955490544468110388479799464071d+00
	 X1( 62)=0.9920401585199581328360853831728576d+00
	 X1( 63)=0.9917806207020271329657876520252450d+00
	 X1( 64)=0.9915169366856671611332340645897933d+00
	 X1( 65)=0.9912491075732250046725317735505000d+00
	 X1( 66)=0.9909771344843762795617994668742164d+00
	 X1( 67)=0.9907010185561207490521828865870820d+00
	 X1( 68)=0.9904207609427775699234730322225553d+00
	 X1( 69)=0.9901363628159804663802837447695951d+00
	 X1( 70)=0.9898478253646728316038055039643836d+00
	 X1( 71)=0.9895551497951027569751060641725293d+00
	 X1( 72)=0.9892583373308179889868121709102356d+00
	 X1( 73)=0.9889573892126608138607935625550349d+00
	 X1( 74)=0.9886523066987628698901921069177935d+00
	 X1( 75)=0.9883430910645398875248047797940522d+00
	 X1( 76)=0.9880297436026863572194472928368640d+00
	 X1( 77)=0.9877122656231701250655023044709179d+00
	 X1( 78)=0.9873906584532269162263980434805204d+00
	 X1( 79)=0.9870649234373547861982747142579444d+00
	 X1( 80)=0.9867350619373084999175813807096420d+00
	 X1( 81)=0.9864010753320938387378086755502523d+00
	 X1( 82)=0.9860629650179618352980056707927271d+00
	 X1( 83)=0.9857207324084029363061551563070085d+00
	 X1( 84)=0.9853743789341410932608926211163314d+00
	 X1( 85)=0.9850239060431277811354523195615976d+00
	 X1( 86)=0.9846693152005359450481105696879284d+00
	 X1( 87)=0.9843106078887538749437732869007356d+00
	 X1( 88)=0.9839477856073790083117229226082069d+00
	 X1( 89)=0.9835808498732116609649005115045622d+00
	 X1( 90)=0.9832098022202486859064523480761476d+00
	 X1( 91)=0.9828346441996770603096187080028285d+00
	 X1( 92)=0.9824553773798674006373846949851531d+00
	 X1( 93)=0.9820720033463674059286513306243504d+00
	 X1( 94)=0.9816845237018952292780189399448886d+00
	 X1( 95)=0.9812929400663327775366051772464947d+00
	 X1( 96)=0.9808972540767189392616470880738106d+00
	 X1( 97)=0.9804974673872427409429607653383842d+00
	 X1( 98)=0.9800935816692364315346537400102815d+00
	 X1( 99)=0.9796855986111684953208045208059486d+00
	 X1(100)=0.9792735199186365931441409017872366d+00
	 X1(101)=0.9788573473143604320270640021568952d+00
	 X1(102)=0.9784370825381745632146786743585367d+00
	 X1(103)=0.9780127273470211086698030787272129d+00
	 X1(104)=0.9775842835149424160502410202526432d+00
	 X1(105)=0.9771517528330736421989102037186030d+00
	 X1(106)=0.9767151371096352651777280011596507d+00
	 X1(107)=0.9762744381699255248764637409750838d+00
	 X1(108)=0.9758296578563127922280730106011497d+00
	 X1(109)=0.9753807980282278670623350938488912d+00
	 X1(110)=0.9749278605621562046299195105647439d+00
	 X1(111)=0.9744708473516300708293117531040425d+00
	 X1(112)=0.9740097603072206261693317773040224d+00
	 X1(113)=0.9735446013565299385002816552102191d+00
	 X1(114)=0.9730753724441829245470610773742715d+00
	 X1(115)=0.9726020755318192202778911439532260d+00
	 X1(116)=0.9721247125980849801425881417109586d+00
	 X1(117)=0.9716432856386246052146298001911286d+00
	 X1(118)=0.9711577966660724002715568832818741d+00
	 X1(119)=0.9706682477100441598485529276187876d+00
	 X1(120)=0.9701746408171286833003445096011061d+00
	 X1(121)=0.9696769780508792189068636286720802d+00
	 X1(122)=0.9691752614918048370584126542558899d+00
	 X1(123)=0.9686694932373617325563708137647329d+00
	 X1(124)=0.9681596754019444560657794141010573d+00
	 X1(125)=0.9676458101168770747564409022644522d+00
	 X1(126)=0.9671278995304042621694644938435135d+00
	 X1(127)=0.9666059458076823173464884419132956d+00
	 X1(128)=0.9660799511307701132591060926466749d+00
	 X1(129)=0.9655499176986199745763196862670818d+00
	 X1(130)=0.9650158477270684848081424204098305d+00
	 X1(131)=0.9644777434488272228637656043028431d+00
	 X1(132)=0.9639356071134734290630038024813664d+00
	 X1(133)=0.9633894409874406006400267014158684d+00
	 X1(134)=0.9628392473540090167786820362660755d+00
	 X1(135)=0.9622850285132961932190092922501467d+00
	 X1(136)=0.9617267867822472664748390496252675d+00
	 X1(137)=0.9611645244946253077026677763652669d+00
	 X1(138)=0.9605982440010015662622925912460506d+00
	 X1(139)=0.9600279476687456430099850248047403d+00
	 X1(140)=0.9594536378820155933652770987910500d+00
	 X1(141)=0.9588753170417479601927271282509574d+00
	 X1(142)=0.9582929875656477365403265259753464d+00
	 X1(143)=0.9577066518881782582765025581675585d+00
	 X1(144)=0.9571163124605510266679654640686974d+00
	 X1(145)=0.9565219717507154609409416119605453d+00
	 X1(146)=0.9559236322433485808686274202910553d+00
	 X1(147)=0.9553212964398446194279916263164171d+00
	 X1(148)=0.9547149668583045655693461361532114d+00
	 X1(149)=0.9541046460335256371423981398708105d+00
	 X1(150)=0.9534903365169906840227884234868648d+00
	 X1(151)=0.9528720408768575214834128565997290d+00
	 X1(152)=0.9522497616979481938551158799365962d+00
	 X1(153)=0.9516235015817381685216364612530562d+00
	 X1(154)=0.9509932631463454602939784306376027d+00
	 X1(155)=0.9503590490265196862096683471204499d+00
	 X1(156)=0.9497208618736310508026550872522129d+00
	 X1(157)=0.9490787043556592618898961826283169d+00
	 X1(158)=0.9484325791571823769209665667503064d+00
	 X1(159)=0.9477824889793655799373158216387615d+00
	 X1(160)=0.9471284365399498891880902406945801d+00
	 X1(161)=0.9464704245732407954497260458476781d+00
	 X1(162)=0.9458084558300968310968099133917886d+00
	 X1(163)=0.9451425330779180699719925733971006d+00
	 X1(164)=0.9444726591006345581030306514973852d+00
	 X1(165)=0.9437988366986946753153211184087878d+00
	 X1(166)=0.9431210686890534277885817009655009d+00
	 X1(167)=0.9424393579051606716066193879356641d+00
	 X1(168)=0.9417537071969492673494177335646600d+00
	 X1(169)=0.9410641194308231657770620208132124d+00
	 X1(170)=0.9403705974896454246553094937222548d+00
	 X1(171)=0.9396731442727261567728998033321887d+00
	 X1(172)=0.9389717626958104092009885331783148d+00
	 X1(173)=0.9382664556910659738453741776264004d+00
	 X1(174)=0.9375572262070711293424762382355614d+00
	 X1(175)=0.9368440772088023143503091789577718d+00
	 X1(176)=0.9361270116776217322859838393081807d+00
	 X1(177)=0.9354060326112648875615545446590075d+00
	 X1(178)=0.9346811430238280533703165735011428d+00
	 X1(179)=0.9339523459457556710759448418498367d+00
	 X1(180)=0.9332196444238276812571506439020171d+00
	 X1(181)=0.9324830415211467864608190445315669d+00
	 X1(182)=0.9317425403171256457168750521754557d+00
	 X1(183)=0.9309981439074740008684120090502303d+00
	 X1(184)=0.9302498554041857347709007183700220d+00
	 X1(185)=0.9294976779355258614145826842322290d+00
	 X1(186)=0.9287416146460174480244354680075940d+00
	 X1(187)=0.9279816686964284691923825643239849d+00
	 X1(188)=0.9272178432637585930967043689692009d+00
	 X1(189)=0.9264501415412258998638907491540868d+00
	 X1(190)=0.9256785667382535321284594324653934d+00
	 X1(191)=0.9249031220804562778465479033863629d+00
	 X1(192)=0.9241238108096270854193697343564923d+00
	 X1(193)=0.9233406361837235111829092808615662d+00
	 X1(194)=0.9225536014768540993205114358689888d+00
	 X1(195)=0.9217627099792646942553056669270866d+00
	 X1(196)=0.9209679649973246855796858483032568d+00
	 X1(197)=0.9201693698535131855793494495151831d+00
	 X1(198)=0.9193669278864051394096814493803517d+00
	 X1(199)=0.9185606424506573679825499101384147d+00
	 X1(200)=0.9177505169169945436218614681535783d+00
	 X1(201)=0.9169365546721950985465060750436164d+00
	 X1(202)=0.9161187591190770662396011546704967d+00
	 X1(203)=0.9152971336764838557632259261259745d+00
	 X1(204)=0.9144716817792699590781169795138123d+00
	 X1(205)=0.9136424068782865914280762788276277d+00
	 X1(206)=0.9128093124403672648491226034080491d+00
	 X1(207)=0.9119724019483132948636970251925879d+00
	 X1(208)=0.9111316789008792404205123521034984d+00
	 X1(209)=0.9102871468127582771409155473096542d+00
	 X1(210)=0.9094388092145675039329109586044503d+00
	 X1(211)=0.9085866696528331830342707606190966d+00
	 X1(212)=0.9077307316899759135464373238953446d+00
	 X1(213)=0.9068709989042957385212002778295999d+00
	 X1(214)=0.9060074748899571856624088280273472d+00
	 X1(215)=0.9051401632569742417052574215288731d+00
	 X1(216)=0.9042690676311952605359601245404007d+00
	 X1(217)=0.9033941916542878051149060855853008d+00
	 X1(218)=0.9025155389837234232666652012344782d+00
	 X1(219)=0.9016331132927623574004895806402115d+00
	 X1(220)=0.9007469182704381882252326178408180d+00
	 X1(221)=0.8998569576215424125228834260821160d+00
	 X1(222)=0.8989632350666089550451900650738300d+00
	 X1(223)=0.8980657543418986145981203990233418d+00
	 X1(224)=0.8971645191993834443791845593246417d+00
	 X1(225)=0.8962595334067310666329178497866236d+00
	 X1(226)=0.8953508007472889216900975231222584d+00
	 X1(227)=0.8944383250200684514565411739496232d+00
	 X1(228)=0.8935221100397292174176085346388073d+00
	 X1(229)=0.8926021596365629532248022248377310d+00
	 X1(230)=0.8916784776564775519311364922879042d+00
	 X1(231)=0.8907510679609809879422161904619759d+00
	 X1(232)=0.8898199344271651737502411664831820d+00
	 X1(233)=0.8888850809476897515184238795879435d+00
	 X1(234)=0.8879465114307658195835804349331787d+00
	 X1(235)=0.8870042298001395939449272986965321d+00
	 X1(236)=0.8860582399950760048073877570387366d+00
	 X1(237)=0.8851085459703422282479836924616677d+00
	 X1(238)=0.8841551516961911530741594752732210d+00
	 X1(239)=0.8831980611583447829431557041317660d+00
	 X1(240)=0.8822372783579775738118211768604389d+00
	 X1(241)=0.8812728073116997067865218297676936d+00
	 X1(242)=0.8803046520515402964430754494591699d+00
	 X1(243)=0.8793328166249305346869107344518368d+00
	 X1(244)=0.8783573050946867702239187636804127d+00
	 X1(245)=0.8773781215389935237127341140951222d+00
	 X1(246)=0.8763952700513864386694517588669015d+00
	 X1(247)=0.8754087547407351681960544701202408d+00
	 X1(248)=0.8744185797312261976040937444850705d+00
	 X1(249)=0.8734247491623456030054352649787088d+00
	 X1(250)=0.8724272671888617459421476076792042d+00
	 X1(251)=0.8714261379808079041278802952159259d+00
	 X1(252)=0.8704213657234648383733443901665557d+00
	 X1(253)=0.8694129546173432957687756088974920d+00
	 X1(254)=0.8684009088781664491965264191039866d+00
	 X1(255)=0.8673852327368522732471997611851909d+00
	 X1(256)=0.8663659304394958566130029035169365d+00
	 X1(257)=0.8653430062473516510322655035519796d+00
	 X1(258)=0.8643164644368156568593311993752385d+00
	 X1(259)=0.8632863092994075453342969987630893d+00
	 X1(260)=0.8622525451417527176273393638351652d+00
	 X1(261)=0.8612151762855643007326302079395711d+00
	 X1(262)=0.8601742070676250802871100263745508d+00
	 X1(263)=0.8591296418397693703896490728191800d+00
	 X1(264)=0.8580814849688648204963908678216131d+00
	 X1(265)=0.8570297408367941594683353832760802d+00
	 X1(266)=0.8559744138404368768474819864107451d+00
	 X1(267)=0.8549155083916508414381146473105093d+00
	 X1(268)=0.8538530289172538572700740143159626d+00
	 X1(269)=0.8527869798590051570211227406773513d+00
	 X1(270)=0.8517173656735868329757719025072863d+00
	 X1(271)=0.8506441908325852055981974812759687d+00
	 X1(272)=0.8495674598224721297971366927372120d+00
	 X1(273)=0.8484871771445862389609144271731191d+00
	 X1(274)=0.8474033473151141268410102221118164d+00
	 X1(275)=0.8463159748650714673628360171194274d+00
	 X1(276)=0.8452250643402840724426544398090493d+00
	 X1(277)=0.8441306203013688878898265417618015d+00
	 X1(278)=0.8430326473237149274738367415353069d+00
	 X1(279)=0.8419311499974641452358012382618769d+00
	 X1(280)=0.8408261329274922461244243324321553d+00
	 X1(281)=0.8397176007333894350366249292413860d+00
	 X1(282)=0.8386055580494411043433130032675271d+00
	 X1(283)=0.8374900095246084599810529701771937d+00
	 X1(284)=0.8363709598225090861906077405423787d+00
	 X1(285)=0.8352484136213974489836137216249057d+00
	 X1(286)=0.8341223756141453384188931840748596d+00
	 X1(287)=0.8329928505082222497701662208234851d+00
	 X1(288)=0.8318598430256757036671799939612613d+00
	 X1(289)=0.8307233579031115052925280910105371d+00
	 X1(290)=0.8295833998916739427166875936630998d+00
	 X1(291)=0.8284399737570259244540558986916603d+00
	 X1(292)=0.8272930842793290563230234212971716d+00
	 X1(293)=0.8261427362532236576933720545593181d+00
	 X1(294)=0.8249889344878087172045426538549913d+00
	 X1(295)=0.8238316838066217880385678610401411d+00
	 X1(296)=0.8226709890476188228317192787965116d+00
	 X1(297)=0.8215068550631539483091703497703864d+00
	 X1(298)=0.8203392867199591797272282869209166d+00
	 X1(299)=0.8191682888991240752079400397977687d+00
	 X1(300)=0.8179938664960753300511285652299640d+00
	 X1(301)=0.8168160244205563111091665990797069d+00
	 X1(302)=0.8156347675966065313100456972479202d+00
	 X1(303)=0.8144501009625410644145485279648848d+00
	 X1(304)=0.8132620294709299000935822525140186d+00
	 X1(305)=0.8120705580885772394119803268751218d+00
	 X1(306)=0.8108756917965007308053291912925798d+00
	 X1(307)=0.8096774355899106466366250874327355d+00
	 X1(308)=0.8084757944781890004198146525531250d+00
	 X1(309)=0.8072707734848686047975209859262304d+00
	 X1(310)=0.8060623776476120703605045636050587d+00
	 X1(311)=0.8048506120181907453966556924519606d+00
	 X1(312)=0.8036354816624635966575621421419049d+00
	 X1(313)=0.8024169916603560312309421735647452d+00
	 X1(314)=0.8011951471058386596074793926571446d+00
	 X1(315)=0.7999699531069060000308416991646715d+00
	 X1(316)=0.7987414147855551242199120691404860d+00
	 X1(317)=0.7975095372777642445525040071030279d+00
	 X1(318)=0.7962743257334712428000792276766179d+00
	 X1(319)=0.7950357853165521405032294762030085d+00
	 X1(320)=0.7937939212047995110779283722172453d+00
	 X1(321)=0.7925487385899008337428027578079258d+00
	 X1(322)=0.7913002426774167893579162537117982d+00
	 X1(323)=0.7900484386867594982658005685089635d+00
	 X1(324)=0.7887933318511707002257125694726154d+00
	 X1(325)=0.7875349274176998765323372064727547d+00
	 X1(326)=0.7862732306471823144103980818247153d+00
	 X1(327)=0.7850082468142171137768787781002905d+00
	 X1(328)=0.7837399812071451364627989916729902d+00
	 X1(329)=0.7824684391280268979867300711423664d+00
	 X1(330)=0.7811936258926204019724747257698304d+00
	 X1(331)=0.7799155468303589173035754486560460d+00
	 X1(332)=0.7786342072843286981075555915954629d+00
	 X1(333)=0.7773496126112466466630360323561655d+00
	 X1(334)=0.7760617681814379193231089895538357d+00
	 X1(335)=0.7747706793788134755485887643198275d+00
	 X1(336)=0.7734763516008475701449970206092672d+00
	 X1(337)=0.7721787902585551887973776562615335d+00
	 X1(338)=0.7708780007764694269972733638198711d+00
	 X1(339)=0.7695739885926188124564326326483188d+00
	 X1(340)=0.7682667591585045711020522010633268d+00
	 X1(341)=0.7669563179390778367485958280367500d+00
	 X1(342)=0.7656426704127168045414657175403922d+00
	 X1(343)=0.7643258220712038282680379938055916d+00
	 X1(344)=0.7630057784197024616318082916818712d+00
	 X1(345)=0.7616825449767344435856277919154117d+00
	 X1(346)=0.7603561272741566278202438955517527d+00
	 X1(347)=0.7590265308571378565045931938200451d+00
	 X1(348)=0.7576937612841357783745274488024096d+00
	 X1(349)=0.7563578241268736112668859549572351d+00
	 X1(350)=0.7550187249703168491960599011769941d+00
	 X1(351)=0.7536764694126499140704261965484647d+00
	 X1(352)=0.7523310630652527521462596593769395d+00
	 X1(353)=0.7509825115526773753169634973685841d+00
	 X1(354)=0.7496308205126243473356886261707741d+00
	 X1(355)=0.7482759955959192150696425827849343d+00
	 X1(356)=0.7469180424664888848846185887277319d+00
	 X1(357)=0.7455569668013379442585047042638192d+00
	 X1(358)=0.7441927742905249287227619886077080d+00
	 X1(359)=0.7428254706371385342310891407366301d+00
	 X1(360)=0.7414550615572737750547192404148795d+00
	 X1(361)=0.7400815527800080873040219382494242d+00
	 X1(362)=0.7387049500473773781763117561244980d+00
	 X1(363)=0.7373252591143520210299900542491864d+00
	 X1(364)=0.7359424857488127963853746973481615d+00
	 X1(365)=0.7345566357315267789527975092849556d+00
	 X1(366)=0.7331677148561231707887752416844516d+00
	 X1(367)=0.7317757289290690806812849969733599d+00
	 X1(368)=0.7303806837696452498653998387428480d+00
	 X1(369)=0.7289825852099217241707646915164436d+00
	 X1(370)=0.7275814390947334727026165769408959d+00
	 X1(371)=0.7261772512816559531582767531716452d+00
	 X1(372)=0.7247700276409806238812654178635033d+00
	 X1(373)=0.7233597740556904027554123017684622d+00
	 X1(374)=0.7219464964214350730415587185553606d+00
	 X1(375)=0.7205302006465066362596684461713907d+00
	 X1(376)=0.7191108926518146122193861949358512d+00
	 X1(377)=0.7176885783708612863023033666666771d+00
	 X1(378)=0.7162632637497169040994113265664175d+00
	 X1(379)=0.7148349547469948135074424944146277d+00
	 X1(380)=0.7134036573338265543880192129080367d+00
	 X1(381)=0.7119693774938368958937495678400773d+00
	 X1(382)=0.7105321212231188215656281163010108d+00
	 X1(383)=0.7090918945302084623063178242943087d+00
	 X1(384)=0.7076487034360599773341074231913966d+00
	 X1(385)=0.7062025539740203832225558643743390d+00
	 X1(386)=0.7047534521898043311310525823354338d+00
	 X1(387)=0.7033014041414688323317388675066752d+00
	 X1(388)=0.7018464158993879321384518002751961d+00
	 X1(389)=0.7003884935462273323435679060994750d+00
	 X1(390)=0.6989276431769189622688389574735394d+00
	 X1(391)=0.6974638708986354985365271707926864d+00
	 X1(392)=0.6959971828307648336673614240563050d+00
	 X1(393)=0.6945275851048844936120500539050175d+00
	 X1(394)=0.6930550838647360043232992768362174d+00
	 X1(395)=0.6915796852661992074754993186816711d+00
	 X1(396)=0.6901013954772665254394529276725737d+00
	 X1(397)=0.6886202206780171756197330887725531d+00
	 X1(398)=0.6871361670605913342624684495407454d+00
	 X1(399)=0.6856492408291642498415662097102201d+00
	 X1(400)=0.6841594481999203061315930170486368d+00
	 X1(401)=0.6826667954010270350757447500268304d+00
	 X1(402)=0.6811712886726090795575459524777380d+00
	 X1(403)=0.6796729342667221061851291159052760d+00
	 X1(404)=0.6781717384473266681971529805249038d+00
	 X1(405)=0.6766677074902620185996275456110643d+00
	 X1(406)=0.6751608476832198736431215424197372d+00
	 X1(407)=0.6736511653257181267500357279771668d+00
	 X1(408)=0.6721386667290745130018325045105461d+00
	 X1(409)=0.6706233582163802242963190563770583d+00
	 X1(410)=0.6691052461224734752852874231601348d+00
	 X1(411)=0.6675843367939130202030206932839633d+00
	 X1(412)=0.6660606365889516206963798061889401d+00
	 X1(413)=0.6645341518775094647673902919536430d+00
	 X1(414)=0.6630048890411475369394526543866593d+00
	 X1(415)=0.6614728544730409397585040161898218d+00
	 X1(416)=0.6599380545779521667406620919606461d+00
	 X1(417)=0.6584004957722043268780855357054653d+00
	 X1(418)=0.6568601844836543208149872233274080d+00
	 X1(419)=0.6553171271516659688059390763882940d+00
	 X1(420)=0.6537713302270830905688086104761089d+00
	 X1(421)=0.6522228001722025371448684988972412d+00
	 X1(422)=0.6506715434607471748788210793143943d+00
	 X1(423)=0.6491175665778388216316798965282540d+00
	 X1(424)=0.6475608760199711353396500680167891d+00
	 X1(425)=0.6460014782949824550323484792656766d+00
	 X1(426)=0.6444393799220285944239035625139135d+00
	 X1(427)=0.6428745874315555881906726844693389d+00
	 X1(428)=0.6413071073652723910495129649906668d+00
	 X1(429)=0.6397369462761235297507386688588106d+00
	 X1(430)=0.6381641107282617081000951557458790d+00
	 X1(431)=0.6365886072970203651242757385123066d+00
	 X1(432)=0.6350104425688861864947036862002637d+00
	 X1(433)=0.6334296231414715693244970147258372d+00
	 X1(434)=0.6318461556234870404537286344865879d+00
	 X1(435)=0.6302600466347136283382888690800875d+00
	 X1(436)=0.6286713028059751886578513222600004d+00
	 X1(437)=0.6270799307791106837586365503283725d+00
	 X1(438)=0.6254859372069464160468609935671309d+00
	 X1(439)=0.6238893287532682154489511322415854d+00
	 X1(440)=0.6222901120927935810547948593591273d+00
	 X1(441)=0.6206882939111437770604936029345808d+00
	 X1(442)=0.6190838809048158831272697841990273d+00
	 X1(443)=0.6174768797811547992733747641926851d+00
	 X1(444)=0.6158672972583252054160325087079179d+00
	 X1(445)=0.6142551400652834756806437898010350d+00
	 X1(446)=0.6126404149417495475946648402786482d+00
	 X1(447)=0.6110231286381787462837629848954464d+00
	 X1(448)=0.6094032879157335637880398876868290d+00
	 X1(449)=0.6077808995462553936163006781155084d+00
	 X1(450)=0.6061559703122362206565343487515556d+00
	 X1(451)=0.6045285070067902665609574532481346d+00
	 X1(452)=0.6028985164336255907241592746400900d+00
	 X1(453)=0.6012660054070156469730722797014716d+00
	 X1(454)=0.5996309807517707960876768244748635d+00
	 X1(455)=0.5979934493032097742715337283560492d+00
	 X1(456)=0.5963534179071311176914224885100950d+00
	 X1(457)=0.5947108934197845432055465621395493d+00
	 X1(458)=0.5930658827078422853999503004543120d+00
	 X1(459)=0.5914183926483703900529747743401673d+00
	 X1(460)=0.5897684301287999641477618869253647d+00
	 X1(461)=0.5881160020468983825529978217404836d+00
	 X1(462)=0.5864611153107404514922680261967066d+00
	 X1(463)=0.5848037768386795289225765779142506d+00
	 X1(464)=0.5831439935593186019427629252608458d+00
	 X1(465)=0.5814817724114813213527286325567270d+00
	 X1(466)=0.5798171203441829934845658940165926d+00
	 X1(467)=0.5781500443166015294268582078815478d+00
	 X1(468)=0.5764805512980483517636017225983896d+00
	 X1(469)=0.5748086482679392589493733795850954d+00
	 X1(470)=0.5731343422157652474425490813375015d+00
	 X1(471)=0.5714576401410632917185517086425328d+00
	 X1(472)=0.5697785490533870822852848957296811d+00
	 X1(473)=0.5680970759722777218230840465786206d+00
	 X1(474)=0.5664132279272343795716911385728695d+00
	 X1(475)=0.5647270119576849040869344105154157d+00
	 X1(476)=0.5630384351129563944899680699725004d+00
	 X1(477)=0.5613475044522457303321006792587067d+00
	 X1(478)=0.5596542270445900601984138893947254d+00
	 X1(479)=0.5579586099688372491735456863353713d+00
	 X1(480)=0.5562606603136162852931842929584807d+00
	 X1(481)=0.5545603851773076451049903330062754d+00
	 X1(482)=0.5528577916680136184628358086628034d+00
	 X1(483)=0.5511528869035285926784188710195326d+00
	 X1(484)=0.5494456780113092961544832716135874d+00
	 X1(485)=0.5477361721284450016240407728091701d+00
	 X1(486)=0.5460243764016276891201636643242563d+00
	 X1(487)=0.5443102979871221688010828819757290d+00
	 X1(488)=0.5425939440507361637554950520229228d+00
	 X1(489)=0.5408753217677903529131490896304894d+00
	 X1(490)=0.5391544383230883741859497622471157d+00
	 X1(491)=0.5374313009108867879649818874097211d+00
	 X1(492)=0.5357059167348650010990245689382284d+00
	 X1(493)=0.5339782930080951514802900849910344d+00
	 X1(494)=0.5322484369530119533632867253151869d+00
	 X1(495)=0.5305163558013825035428690325595642d+00
	 X1(496)=0.5287820567942760485177025330377757d+00
	 X1(497)=0.5270455471820337127655331451459750d+00
	 X1(498)=0.5253068342242381882568140280774674d+00
	 X1(499)=0.5235659251896833853334046788512237d+00
	 X1(500)=0.5218228273563440450792186013077598d+00
	 X1(501)=0.5200775480113453133098568559480706d+00
	 X1(502)=0.5183300944509322763084252536264206d+00
	 X1(503)=0.5165804739804394584348928784849862d+00
	 X1(504)=0.5148286939142602817365090154690791d+00
	 X1(505)=0.5130747615758164876869544146195853d+00
	 X1(506)=0.5113186842975275211820611474402504d+00
	 X1(507)=0.5095604694207798769200930993196245d+00
	 X1(508)=0.5078001242958964082947363955912282d+00
	 X1(509)=0.5060376562821055989291057766833824d+00
	 X1(510)=0.5042730727475107969792291192870175d+00
	 X1(511)=0.5025063810690594123356279449026898d+00
	 X1(512)=0.5007375886325120768517668638663300d+00
	 X1(513)=0.4989667028324117677282994713484729d+00
	 X1(514)=0.4971937310720528941821922412277083d+00
	 X1(515)=0.4954186807634503475299614535119082d+00
	 X1(516)=0.4936415593273085148144111404787578d+00
	 X1(517)=0.4918623741929902561044124452908303d+00
	 X1(518)=0.4900811327984858455974166538727473d+00
	 X1(519)=0.4882978425903818766545454856839726d+00
	 X1(520)=0.4865125110238301308982530110478619d+00
	 X1(521)=0.4847251455625164115027038012753980d+00
	 X1(522)=0.4829357536786293408071616123224906d+00
	 X1(523)=0.4811443428528291223828320525170031d+00
	 X1(524)=0.4793509205742162676837512893622618d+00
	 X1(525)=0.4775554943403002874124609089464581d+00
	 X1(526)=0.4757580716569683477313565534430898d+00
	 X1(527)=0.4739586600384538914507449269597470d+00
	 X1(528)=0.4721572670073052243247901769667067d+00
	 X1(529)=0.4703539000943540665866765271008774d+00
	 X1(530)=0.4685485668386840698544593566848071d+00
	 X1(531)=0.4667412747875992995392216922172194d+00
	 X1(532)=0.4649320314965926828872972957756562d+00
	 X1(533)=0.4631208445293144227884652041203613d+00
	 X1(534)=0.4613077214575403774821636897009264d+00
	 X1(535)=0.4594926698611404062939141801451119d+00
	 X1(536)=0.4576756973280466815342876855566769d+00
	 X1(537)=0.4558568114542219666928877424722777d+00
	 X1(538)=0.4540360198436278610599647890351738d+00
	 X1(539)=0.4522133301081930109084172372465247d+00
	 X1(540)=0.4503887498677812873690743044667459d+00
	 X1(541)=0.4485622867501599311322949070752703d+00
	 X1(542)=0.4467339483909676641090556037750441d+00
	 X1(543)=0.4449037424336827681848387038683880d+00
	 X1(544)=0.4430716765295911311997692263560197d+00
	 X1(545)=0.4412377583377542602885864083459550d+00
	 X1(546)=0.4394019955249772627141719154308714d+00
	 X1(547)=0.4375643957657767943284928018309007d+00
	 X1(548)=0.4357249667423489757949526036355685d+00
	 X1(549)=0.4338837161445372767062787238479769d+00
	 X1(550)=0.4320406516698003677322084825728543d+00
	 X1(551)=0.4301957810231799409313698590367193d+00
	 X1(552)=0.4283491119172684983618860436243531d+00
	 X1(553)=0.4265006520721771091253654472046570d+00
	 X1(554)=0.4246504092155031349790707811467503d+00
	 X1(555)=0.4227983910822979246511922240421373d+00
	 X1(556)=0.4209446054150344769942805297016007d+00
	 X1(557)=0.4190890599635750731120262049392087d+00
	 X1(558)=0.4172317624851388775947005944457932d+00
	 X1(559)=0.4153727207442695089987038531482651d+00
	 X1(560)=0.4135119425128025797057933633092326d+00
	 X1(561)=0.4116494355698332052976941637061146d+00
	 X1(562)=0.4097852077016834835819204010050816d+00
	 X1(563)=0.4079192667018699434047636883799841d+00
	 X1(564)=0.4060516203710709633875305629894868d+00
	 X1(565)=0.4041822765170941607222369715889275d+00
	 X1(566)=0.4023112429548437501630928817912644d+00
	 X1(567)=0.4004385275062878733502347147806378d+00
	 X1(568)=0.3985641380004258986022873231016178d+00
	 X1(569)=0.3966880822732556913144606939786689d+00
	 X1(570)=0.3948103681677408550990094439475882d+00
	 X1(571)=0.3929310035337779438050054838899604d+00
	 X1(572)=0.3910499962281636445544959743416655d+00
	 X1(573)=0.3891673541145619319322398586882311d+00
	 X1(574)=0.3872830850634711934663368560568117d+00
	 X1(575)=0.3853971969521913265371828158626119d+00
	 X1(576)=0.3835096976647908068523047815649438d+00
	 X1(577)=0.3816205950920737286247479817354316d+00
	 X1(578)=0.3797298971315468165928052615412321d+00
	 X1(579)=0.3778376116873864100189971867049694d+00
	 X1(580)=0.3759437466704054188063281944282680d+00
	 X1(581)=0.3740483099980202518699607311675823d+00
	 X1(582)=0.3721513095942177179025653050421557d+00
	 X1(583)=0.3702527533895218986717197905495046d+00
	 X1(584)=0.3683526493209609949878461546813378d+00
	 X1(585)=0.3664510053320341454812870259922189d+00
	 X1(586)=0.3645478293726782183272382011969375d+00
	 X1(587)=0.3626431293992345760573662769852384d+00
	 X1(588)=0.3607369133744158135970531074714846d+00
	 X1(589)=0.3588291892672724696673207195716376d+00
	 X1(590)=0.3569199650531597116906016691526747d+00
	 X1(591)=0.3550092487137039943396305895647351d+00
	 X1(592)=0.3530970482367696918688428706808119d+00
	 X1(593)=0.3511833716164257043677760103720428d+00
	 X1(594)=0.3492682268529120380760782009804016d+00
	 X1(595)=0.3473516219526063598998371503590736d+00
	 X1(596)=0.3454335649279905262690499899807053d+00
	 X1(597)=0.3435140637976170864761623910141371d+00
	 X1(598)=0.3415931265860757606357116926927508d+00
	 X1(599)=0.3396707613239598924052149452961573d+00
	 X1(600)=0.3377469760478328766075482821981015d+00
	 X1(601)=0.3358217788001945618951689612560379d+00
	 X1(602)=0.3338951776294476285966357548932121d+00
	 X1(603)=0.3319671805898639418859871201160524d+00
	 X1(604)=0.3300377957415508804156397439845057d+00
	 X1(605)=0.3281070311504176405535726362793490d+00
	 X1(606)=0.3261748948881415163656639288596523d+00
	 X1(607)=0.3242413950321341554841489400490923d+00
	 X1(608)=0.3223065396655077910032688719078143d+00
	 X1(609)=0.3203703368770414495432797280155725d+00
	 X1(610)=0.3184327947611471356240906689929719d+00
	 X1(611)=0.3164939214178359924899000620042856d+00
	 X1(612)=0.3145537249526844395262959285034753d+00
	 X1(613)=0.3126122134768002864113853510931579d+00
	 X1(614)=0.3106693951067888241426146651552137d+00
	 X1(615)=0.3087252779647188930810389334749248d+00
	 X1(616)=0.3067798701780889281548952820138015d+00
	 X1(617)=0.3048331798797929813644301618879482d+00
	 X1(618)=0.3028852152080867217300254960797362d+00
	 X1(619)=0.3009359843065534128257629690539885d+00
	 X1(620)=0.2989854953240698680406594228716146d+00
	 X1(621)=0.2970337564147723837098994342019233d+00
	 X1(622)=0.2950807757380226502584836624404203d+00
	 X1(623)=0.2931265614583736414998034795550205d+00
	 X1(624)=0.2911711217455354822317437169259673d+00
	 X1(625)=0.2892144647743412942730060929315851d+00
	 X1(626)=0.2872565987247130210824360169839911d+00
	 X1(627)=0.2852975317816272311042250007592567d+00
	 X1(628)=0.2833372721350808999819498451209343d+00
	 X1(629)=0.2813758279800571717844981113325398d+00
	 X1(630)=0.2794132075164910993870171272241951d+00
	 X1(631)=0.2774494189492353641501109226543741d+00
	 X1(632)=0.2754844704880259750405960335252577d+00
	 X1(633)=0.2735183703474479473372130594077740d+00
	 X1(634)=0.2715511267469009610647762061506806d+00
	 X1(635)=0.2695827479105649993003277913302489d+00
	 X1(636)=0.2676132420673659664949488366889286d+00
	 X1(637)=0.2656426174509412869549604174610522d+00
	 X1(638)=0.2636708822996054836263333833418901d+00
	 X1(639)=0.2616980448563157373262064094764490d+00
	 X1(640)=0.2597241133686374265654940778820708d+00
	 X1(641)=0.2577490960887096481066478298324059d+00
	 X1(642)=0.2557730012732107184007131675804983d+00
	 X1(643)=0.2537958371833236560479064190488342d+00
	 X1(644)=0.2518176120847016454260137114300942d+00
	 X1(645)=0.2498383342474334816309935286923654d+00
	 X1(646)=0.2478580119460089968742423534375718d+00
	 X1(647)=0.2458766534592844684810604150952732d+00
	 X1(648)=0.2438942670704480086349314837216712d+00
	 X1(649)=0.2419108610669849360123069613940910d+00
	 X1(650)=0.2399264437406431294526602310253572d+00
	 X1(651)=0.2379410233873983638086523250538340d+00
	 X1(652)=0.2359546083074196281213244735795025d+00
	 X1(653)=0.2339672068050344262653069828028430d+00
	 X1(654)=0.2319788271886940602091071797725888d+00
	 X1(655)=0.2299894777709388960356118381542271d+00
	 X1(656)=0.2279991668683636128680115716896269d+00
	 X1(657)=0.2260079028015824348464261469280488d+00
	 X1(658)=0.2240156938951943463005805243713005d+00
	 X1(659)=0.2220225484777482902639516870946875d+00
	 X1(660)=0.2200284748817083504748759578870220d+00
	 X1(661)=0.2180334814434189170101755397061260d+00
	 X1(662)=0.2160375765030698356969314394824066d+00
	 X1(663)=0.2140407684046615414480977517361399d+00
	 X1(664)=0.2120430654959701756677194858205592d+00
	 X1(665)=0.2100444761285126878715827185817389d+00
	 X1(666)=0.2080450086575119216691918425591849d+00
	 X1(667)=0.2060446714418616852530340582621130d+00
	 X1(668)=0.2040434728440918065411560272723028d+00
	 X1(669)=0.2020414212303331731191417606743883d+00
	 X1(670)=0.2000385249702827571276443643302709d+00
	 X1(671)=0.1980347924371686252416871985303556d+00
	 X1(672)=0.1960302320077149338880123343074088d+00
	 X1(673)=0.1940248520621069098468159019284625d+00
	 X1(674)=0.1920186609839558163842710285283503d+00
	 X1(675)=0.1900116671602639050622995512597168d+00
	 X1(676)=0.1880038789813893533721135694558153d+00
	 X1(677)=0.1859953048410111883381071638837755d+00
	 X1(678)=0.1839859531360941962387372629595299d+00
	 X1(679)=0.1819758322668538185910906745560314d+00
	 X1(680)=0.1799649506367210345458917275211459d+00
	 X1(681)=0.1779533166523072298397617789905883d+00
	 X1(682)=0.1759409387233690524515980417969826d+00
	 X1(683)=0.1739278252627732551099947705036304d+00
	 X1(684)=0.1719139846864615247986848145984888d+00
	 X1(685)=0.1698994254134152994070339029403865d+00
	 X1(686)=0.1678841558656205716726737644284074d+00
	 X1(687)=0.1658681844680326805634133158419795d+00
	 X1(688)=0.1638515196485410902456196586514237d+00
	 X1(689)=0.1618341698379341567863125221070195d+00
	 X1(690)=0.1598161434698638827362670698620819d+00
	 X1(691)=0.1577974489808106597414706515577344d+00
	 X1(692)=0.1557780948100479993303291289821992d+00
	 X1(693)=0.1537580893996072520240678384062910d+00
	 X1(694)=0.1517374411942423149178210662827233d+00
	 X1(695)=0.1497161586413943278799521144757558d+00
	 X1(696)=0.1476942501911563585171936133581176d+00
	 X1(697)=0.1456717242962380760532447062751243d+00
	 X1(698)=0.1436485894119304142685080768351149d+00
	 X1(699)=0.1416248539960702236486955210470177d+00
	 X1(700)=0.1396005265090049128900758792988254d+00
	 X1(701)=0.1375756154135570799091836383664273d+00
	 X1(702)=0.1355501291749891325048503908745884d+00
	 X1(703)=0.1335240762609678988204645987174602d+00
	 X1(704)=0.1314974651415292277544077477039973d+00
	 X1(705)=0.1294703042890425794666570029457857d+00
	 X1(706)=0.1274426021781756061295858780753799d+00
	 X1(707)=0.1254143672858587230710352160992063d+00
	 X1(708)=0.1233856080912496704577669453799248d+00
	 X1(709)=0.1213563330756980656674526207409213d+00
	 X1(710)=0.1193265507227099464973876868250164d+00
	 X1(711)=0.1172962695179123053581607084577754d+00
	 X1(712)=0.1152654979490176146005445007028470d+00
	 X1(713)=0.1132342445057883431239131593949868d+00
	 X1(714)=0.1112025176800014644145254410408680d+00
	 X1(715)=0.1091703259654129561620507689360793d+00
	 X1(716)=0.1071376778577222916027493500090678d+00
	 X1(717)=0.1051045818545369227377524741220237d+00
	 X1(718)=0.1030710464553367555749230341902282d+00
	 X1(719)=0.1010370801614386175428096513831894d+00
	 X1(720)=0.9900269147596071722524051470358054d-01
	 X1(721)=0.9696788890378709656513514826675190d-01
	 X1(722)=0.9493268095153207568614380249020351d-01
	 X1(723)=0.9289707612750469048075502701726465d-01
	 X1(724)=0.9086108294167312311354222341330658d-01
	 X1(725)=0.8882470990562912558824959435969957d-01
	 X1(726)=0.8678796553255243652744690310658812d-01
	 X1(727)=0.8475085833717519131351083220917696d-01
	 X1(728)=0.8271339683574632573971848394497071d-01
	 X1(729)=0.8067558954599597332026569617526723d-01
	 X1(730)=0.7863744498709985640804935665974287d-01
	 X1(731)=0.7659897167964367126907878584736624d-01
	 X1(732)=0.7456017814558746726240652284191131d-01
	 X1(733)=0.7252107290823002027449349147039094d-01
	 X1(734)=0.7048166449217320055694754306478040d-01
	 X1(735)=0.6844196142328633511659776960103060d-01
	 X1(736)=0.6640197222867056480688975512677027d-01
	 X1(737)=0.6436170543662319626960908485505640d-01
	 X1(738)=0.6232116957660204887596195981344609d-01
	 X1(739)=0.6028037317918979681606267042553125d-01
	 X1(740)=0.5823932477605830648589796477850164d-01
	 X1(741)=0.5619803289993296932085800651049154d-01
	 X1(742)=0.5415650608455703022494265315323374d-01
	 X1(743)=0.5211475286465591174477019830931028d-01
	 X1(744)=0.5007278177590153413753351015207956d-01
	 X1(745)=0.4803060135487663148206566433580284d-01
	 X1(746)=0.4598822013903906398219371142183110d-01
	 X1(747)=0.4394564666668612661157513729478800d-01
	 X1(748)=0.4190288947691885424922686969395052d-01
	 X1(749)=0.3985995710960632345497135485554474d-01
	 X1(750)=0.3781685810534995103403827529010772d-01
	 X1(751)=0.3577360100544778954007390284669306d-01
	 X1(752)=0.3373019435185881986582288038640217d-01
	 X1(753)=0.3168664668716724107075940054799593d-01
	 X1(754)=0.2964296655454675759495630118728566d-01
	 X1(755)=0.2759916249772486400849152406132605d-01
	 X1(756)=0.2555524306094712744570168616245627d-01
	 X1(757)=0.2351121678894146787360219174288694d-01
	 X1(758)=0.2146709222688243634380236746734701d-01
	 X1(759)=0.1942287792035549137725253325143549d-01
	 X1(760)=0.1737858241532127363116772715152311d-01
	 X1(761)=0.1533421425807987899747998413568821d-01
	 X1(762)=0.1328978199523513028217762565422492d-01
	 X1(763)=0.1124529417365884761489594961531918d-01
	 X1(764)=0.9200759340455117738129018631746654d-02
	 X1(765)=0.7156186042924562325436928215714351d-02
	 X1(766)=0.5111582828528605478026995941665003d-02
	 X1(767)=0.3066958244853740549090747454018110d-02
	 X1(768)=0.1022320839575796445281385554075653d-02

Cã==================================================================¬

56       CONTINUE
C
        RETURN
        END
Cã==================================================================¬

c***************************************************************
c           subroutine Olmon(wavel,Rm,ref,eps,KD,KS) Dec 2019
c	2019 New Bulk Constants
c	Olmon et al. PHYSICAL REVIEW B 86, 235147 (2012)
c***************************************************************
C     Calculations of ref=real(n)+i*imag(n) and eps=real(eps)+i*imag(eps)
c     for  GOLD according to Olmon PRB 2012 spline for Bulk Single Cristalyne sample (SC data)

c     with account for size-limiting effect.
c     wavel - is the input wavelength nanometers,
c     Rm    - is the size-correction parameters Rm=aLeff/As
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

c	 Deps=(wavel/wp)**2*[ 1.d0/(1.d0+aim*gamb_om)-1.d0/(1.d0+aim*gama_om)]
c
c	if KD=1 approximation a9-11 is used
c      Deps = Deps2*[B+ i]
c      Deps2=(5.d0/Rm)*(wavel/520.d0)**3.d0
c	 B=(gama+gamb)/2/omega
c
c  	 if KD=2 B=0 is used
c	KS is the parameter for sample choice
c	if KS=0 TS sample constants are used
c	if KS=1 CS crystalline sample constants are used
c	for details see Olmon et al PRB 2012	 
C--------------------------------------------------------

       subroutine Olmon(wavel,Rm,ref,eps,KD,KS)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MN,MK
      DIMENSION XN(1000),YN(1000),MN(1000),XK(1000),YK(1000),MK(1000)

	complex*16 ref, eps,epsb,Deps,aim
C--------------------------------------------------------

	aim=(0.d0,1.d0)

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
		  mn(	1)=-0.1524273019E-03
		  mn(	2)=-0.1524273019E-03
		  mn(	3)=-0.2399223450E-03
		  mn(	4)=-0.2065229207E-03
		  mn(	5)=-0.2038618210E-03
		  mn(	6)=-0.1021415376E-03
		  mn(	7)=	5.391957312E-005
		  mn(	8)=	0.000133552102
		  mn(	9)=	0.0001923957612
		  mn(	10)=	5.940299086E-005
		  mn(	11)=-0.7057478066E-04
		  mn(	12)=-0.1957361553E-03
		  mn(	13)=-0.1053251955E-03
		  mn(	14)=-0.3207085087E-04
		  mn(	15)=	2.803498725E-005
		  mn(	16)=-0.1695365121E-03
		  mn(	17)=-0.4333706503E-03
		  mn(	18)=-0.5433802081E-03
		  mn(	19)=-0.4206389112E-03
		  mn(	20)=-0.1725724564E-03
		  mn(	21)=	0.0003905884493
		  mn(	22)=	0.0004110716963
		  mn(	23)=	0.0003738358622
		  mn(	24)=	0.0002440987556
		  mn(	25)=	0.0001827119006
		  mn(	26)=	9.272259691E-005
		  mn(	27)=	8.644361086E-005
		  mn(	28)=	3.9107112E-005
		  mn(	29)=	5.582399758E-005
		  mn(	30)=	2.045619589E-005
		  mn(	31)=	2.168515224E-005
		  mn(	32)=	2.897205493E-005
		  mn(	33)=	1.639084117E-005
		  mn(	34)=	1.305458033E-005
		  mn(	35)=	1.418166918E-005
		  mn(	36)=	6.436813542E-006
		  mn(	37)=	1.189836085E-005
		  mn(	38)=	1.217334038E-005
		  mn(	39)=	5.625120807E-007
		  mn(	40)=	1.055297079E-005
		  mn(	41)=	5.13266774E-006
		  mn(	42)=	2.666896676E-006
		  mn(	43)=	6.703240244E-006
		  mn(	44)=	3.204760292E-006
		  mn(	45)=	2.832329986E-006
		  mn(	46)=	2.669423561E-006
		  mn(	47)=	5.554927862E-006
		  mn(	48)=-0.1367072964E-05
		  mn(	49)=	6.093788393E-006
		  mn(	50)=-0.1970564037E-06
		  mn(	51)=	2.486507343E-006
		  mn(	52)=	1.788494009E-006
		  mn(	53)=	1.076536017E-006
		  mn(	54)=	1.547939399E-006
		  mn(	55)=	3.118225136E-006
		  mn(	56)=	7.136710083E-007
		  mn(	57)=-0.1863466392E-05
		  mn(	58)=	4.629005982E-006
		  mn(	59)=-0.1250770117E-05
		  mn(	60)=	4.001654718E-006
		  mn(	61)=-0.3891119482E-05
		  mn(	62)=	7.338061616E-006
		  mn(	63)=-0.5853453691E-05
		  mn(	64)=	4.554228978E-006
		  mn(	65)=-0.3989695329E-06
		  mn(	66)=	2.357217756E-007
		  mn(	67)=	2.290409003E-006
		  mn(	68)=-0.1842857399E-05
		  mn(	69)=	2.075303662E-006
		  mn(	70)=	2.287708556E-006
		  mn(	71)=-0.1733206262E-05
		  mn(	72)=	1.428907365E-006
		  mn(	73)=-0.1018213069E-05
		  mn(	74)=-0.4949059562E-06
		  mn(	75)=	5.506115788E-006
		  mn(	76)=-0.4152897421E-05
		  mn(	77)=	1.536423846E-006
		  mn(	78)=	7.660805399E-007
		  mn(	79)=-0.2006094547E-05
		  mn(	80)=	4.006454642E-006
		  mn(	81)=-0.2775141494E-06
		  mn(	82)=	2.723930442E-008
		  mn(	83)=-0.3177846095E-05
		  mn(	84)=	9.31762294E-006
		  mn(	85)=-0.1071569811E-04
		  mn(	86)=	9.689506723E-006
		  mn(	87)=-0.5235758405E-05
		  mn(	88)=	1.585956179E-006
		  mn(	89)=	1.560154939E-006
		  mn(	90)=-0.5261896434E-05
		  mn(	91)=	1.044322901E-005
		  mn(	92)=-0.9596132540E-05
		  mn(	93)=	6.512041612E-006
		  mn(	94)=-0.1693813576E-05
		  mn(	95)=	4.757966665E-007
		  mn(	96)=-0.9859204715E-06
		  mn(	97)=	2.741953236E-006
		  mn(	98)=-0.1878797056E-05
		  mn(	99)=	1.387793512E-006
		  mn(	100)=-0.1393813877E-05
		  mn(	101)=	3.735767225E-006
		  mn(	102)=-0.7960032629E-05
		  mn(	103)=	4.917638956E-007
		  mn(	104)=	2.994570167E-006
		  mn(	105)=-0.3830659164E-05
		  mn(	106)=	3.402913423E-006
		  mn(	107)=-0.2534655876E-05
		  mn(	108)=	1.220741021E-006
		  mn(	109)=-0.1713654575E-05
		  mn(	110)=	2.728168942E-006
		  mn(	111)=-0.2569732393E-05
		  mn(	112)=	2.372683077E-006
		  mn(	113)=	6.643309589E-009
		  mn(	114)=-0.1384538121E-05
		  mn(	115)=	2.628337047E-006
		  mn(	116)=-0.2049675089E-05
		  mn(	117)=	7.557014498E-007
		  mn(	118)=	2.131368532E-006
		  mn(	119)=-0.1998186968E-05
		  mn(	120)=	1.311094605E-006
		  mn(	121)=-0.1517119164E-06
		  mn(	122)=	8.809150761E-007
		  mn(	123)=-0.1635932008E-05
		  mn(	124)=	2.780304358E-006
		  mn(	125)=-0.1634469609E-05
		  mn(	126)=	8.728103937E-007
		  mn(	127)=	1.510507322E-007
		  mn(	128)=	7.090971599E-007
		  mn(	129)=	1.415964489E-007
		  mn(	130)=	1.099441215E-006
		  mn(	131)=-0.1458696617E-05
		  mn(	132)=	1.290647555E-006
		  mn(	133)=-0.6166590866E-06
		  mn(	134)=	3.917620158E-006
		  mn(	135)=-0.5957273842E-05
		  mn(	136)=	5.177586804E-006
		  mn(	137)=	2.600721755E-007
		  mn(	138)=-0.3023698269E-05
		  mn(	139)=	2.904323436E-006
		  mn(	140)=	8.095759409E-007
		  mn(	141)=-0.3067922983E-05
		  mn(	142)=	2.659760748E-006
		  mn(	143)=	1.932717772E-006
		  mn(	144)=-0.3443067892E-06
		  mn(	145)=-0.3546215489E-05
		  mn(	146)=	5.606799317E-006
		  mn(	147)=-0.2705077369E-05
		  mn(	148)=	2.214035954E-006
		  mn(	149)=-0.1538359169E-05
		  mn(	150)=	9.435752797E-007
		  mn(	151)=	7.117802797E-007
		  mn(	152)=	1.123949355E-006
		  mn(	153)=-0.2203542743E-05
		  mn(	154)=	4.649228595E-006
		  mn(	155)=-0.5164061826E-05
		  mn(	156)=	6.955902358E-006
		  mn(	157)=-0.4788222946E-05
		  mn(	158)=	3.009281177E-006
		  mn(	159)=	1.589432214E-006
		  mn(	160)=-0.3099240585E-05
		  mn(	161)=	1.58476775E-006
		  mn(	162)=	5.577203083E-006
		  mn(	163)=-0.2686300413E-04
		  mn(	164)=	6.120818588E-005
		  mn(	165)=-0.3482661527E-04
		  mn(	166)=	6.696582969E-006
		  mn(	167)=	4.0624334E-006
		  mn(	168)=-0.1364679681E-05
		  mn(	169)=	2.53474697E-006
		  mn(	170)=-0.2143997068E-05
		  mn(	171)=	2.787550172E-006
		  mn(	172)=-0.6255001610E-05
		  mn(	173)=	7.936717321E-006
		  mn(	174)=-0.3868442197E-05
		  mn(	175)=	2.392981403E-006
		  mn(	176)=-0.1251098839E-05
		  mn(	177)=-0.2210865803E-05
		  mn(	178)=	4.437512573E-006
		  mn(	179)=-0.2414404935E-05
		  mn(	180)=	1.901494611E-006
		  mn(	181)=-0.1080520262E-05
		  mn(	182)=-0.2268088327E-05
		  mn(	183)=	4.123621892E-006
		  mn(	184)=-0.1768926988E-05
		  mn(	185)=-0.3550874595E-06
		  mn(	186)=	2.154980054E-006
		  mn(	187)=-0.1326924485E-05
		  mn(	188)=-0.1709332219E-06
		  mn(	189)=-0.1287072675E-05
		  mn(	190)=	2.732089752E-006
		  mn(	191)=-0.1827711411E-05
		  mn(	192)=	1.287024571E-006
		  mn(	193)=	2.920029567E-006
		  mn(	194)=-0.7175356279E-05
		  mn(	195)=	9.058908873E-006
		  mn(	196)=-0.9032067325E-06
		  mn(	197)=	1.688217796E-007
		  mn(	198)=-0.3610011570E-06
		  mn(	199)=	1.160119306E-006
		  mn(	200)=	8.23396094E-007
		  mn(	201)=	5.110984171E-007
		  mn(	202)=-0.6156611359E-05
		  mn(	203)=	6.224635636E-006
		  mn(	204)=	1.897806813E-007
		  mn(	205)=-0.3036211723E-06
		  mn(	206)=	3.011880417E-007
		  mn(	207)=	3.261762821E-006
		  mn(	208)=-0.9264942739E-05
		  mn(	209)=	8.868951069E-006
		  mn(	210)=-0.1065054013E-05
		  mn(	211)=-0.8473766156E-06
		  mn(	212)=	1.150155764E-006
		  mn(	213)=-0.2177408122E-06
		  mn(	214)=	3.134353779E-006
		  mn(	215)=-0.8979854359E-05
		  mn(	216)=	8.461929633E-006
		  mn(	217)=-0.1021665487E-05
		  mn(	218)=-0.1317566356E-05
		  mn(	219)=	2.995484704E-006
		  mn(	220)=-0.1709952767E-05
		  mn(	221)=	5.545309251E-007
		  mn(	222)=	2.125854162E-006
		  mn(	223)=-0.6481266053E-05
		  mn(	224)=	5.640608911E-006
		  mn(	225)=	9.952886822E-007
		  mn(	226)=-0.1748468704E-05
		  mn(	227)=	2.701846371E-006
		  mn(	228)=-0.1575401221E-05
		  mn(	229)=	3.122886446E-007
		  mn(	230)=	2.250004897E-006
		  mn(	231)=-0.7427713728E-05
		  mn(	232)=	9.293522365E-006
		  mn(	233)=-0.8782556495E-05
		  mn(	234)=	1.767860055E-005
		  mn(	235)=-0.1724386760E-04
		  mn(	236)=	7.874787961E-007
		  mn(	237)=	1.087440861E-005
		  mn(	238)=-0.1838195127E-05
		  mn(	239)=-0.6850361185E-05
		  mn(	240)=	6.354277015E-006
		  mn(	241)=-0.2677163186E-05
		  mn(	242)=	1.102921471E-006
		  mn(	243)=-0.4956912182E-05
		  mn(	244)=	1.547213116E-005
		  mn(	245)=-0.1966874760E-04
		  mn(	246)=	2.007226097E-005
		  mn(	247)=-0.2469556866E-04
		  mn(	248)=	2.122942164E-005
		  mn(	249)=-0.9614708173E-05
		  mn(	250)=	1.398569444E-005
		  mn(	251)=-0.1325532962E-04
		  mn(	252)=-0.5950502285E-07
		  mn(	253)=	1.024096979E-005
		  mn(	254)=-0.9572533393E-05
		  mn(	255)=	6.693742982E-006
		  mn(	256)=-0.2445001490E-05
		  mn(	257)=-0.1595892729E-06
		  mn(	258)=-0.1595892729E-06
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
	    mk(	1)=-0.2056974797E-03
	    mk(	2)=-0.2056974797E-03
	    mk(	3)=-0.1068065588E-03
	    mk(	4)=	2.688397435E-005
	    mk(	5)=	0.0001183978268
	    mk(	6)=	0.0001703125192
	    mk(	7)=	0.000213886397
	    mk(	8)=	0.0001172744509
	    mk(	9)=-0.7917934408E-04
	    mk(	10)=-0.1701008492E-03
	    mk(	11)=-0.1914573904E-03
	    mk(	12)=-0.8722156273E-04
	    mk(	13)=	6.038513654E-005
	    mk(	14)=	2.832088618E-005
	    mk(	15)=-0.1230924566E-03
	    mk(	16)=-0.2418463267E-03
	    mk(	17)=-0.1676036327E-03
	    mk(	18)=	0.0001676511814
	    mk(	19)=	0.0004789119822
	    mk(	20)=	0.000722128429
	    mk(	21)=	0.0005407566278
	    mk(	22)=	0.0002561071177
	    mk(	23)=	5.128294627E-005
	    mk(	24)=-0.5852862237E-04
	    mk(	25)=-0.1062370393E-03
	    mk(	26)=-0.5809402905E-04
	    mk(	27)=-0.9587203859E-04
	    mk(	28)=-0.2927795607E-05
	    mk(	29)=-0.1187081651E-03
	    mk(	30)=-0.1380908200E-04
	    mk(	31)=	3.011215267E-005
	    mk(	32)=-0.1097369471E-03
	    mk(	33)=	4.058991964E-006
	    mk(	34)=-0.3597852373E-04
	    mk(	35)=-0.3887003001E-04
	    mk(	36)=	8.135027812E-005
	    mk(	37)=-0.1060730336E-03
	    mk(	38)=-0.2741568969E-04
	    mk(	39)=	3.854510182E-005
	    mk(	40)=-0.4791596583E-04
	    mk(	41)=-0.1335478074E-04
	    mk(	42)=	4.216835303E-005
	    mk(	43)=-0.8320096898E-04
	    mk(	44)=	1.740474506E-005
	    mk(	45)=	3.575885323E-005
	    mk(	46)=-0.6391073138E-04
	    mk(	47)=	1.494621512E-005
	    mk(	48)=	5.342520799E-005
	    mk(	49)=-0.1642259652E-03
	    mk(	50)=	0.000153261479
	    mk(	51)=-0.5989278296E-04
	    mk(	52)=-0.4435751202E-04
	    mk(	53)=	4.905884956E-005
	    mk(	54)=-0.6075869421E-04
	    mk(	55)=	7.080490163E-005
	    mk(	56)=-0.5006182511E-04
	    mk(	57)=-0.5082416871E-04
	    mk(	58)=	7.85017941E-005
	    mk(	59)=-0.7959787112E-04
	    mk(	60)=	0.0001260874676
	    mk(	61)=-0.2125139603E-03
	    mk(	62)=	0.000239937889
	    mk(	63)=-0.1979182266E-03
	    mk(	64)=	0.0001031185805
	    mk(	65)=-0.4357142927E-04
	    mk(	66)=-0.3104021235E-04
	    mk(	67)=	6.618865489E-005
	    mk(	68)=-0.3827516211E-04
	    mk(	69)=-0.1255056031E-04
	    mk(	70)=	4.887492889E-005
	    mk(	71)=-0.1031202277E-04
	    mk(	72)=-0.4127333641E-04
	    mk(	73)=	1.773812258E-005
	    mk(	74)=-0.6318548379E-04
	    mk(	75)=	8.710424826E-005
	    mk(	76)=-0.2624160969E-04
	    mk(	77)=-0.1024047049E-04
	    mk(	78)=-0.2261527568E-04
	    mk(	79)=-0.4390308411E-04
	    mk(	80)=	0.0001131857357
	    mk(	81)=-0.5364249211E-04
	    mk(	82)=	1.510185882E-005
	    mk(	83)=-0.9128112172E-04
	    mk(	84)=	0.0002076084199
	    mk(	85)=-0.2548356662E-03
	    mk(	86)=	0.000175775122
	    mk(	87)=-0.4321478085E-04
	    mk(	88)=-0.2045299000E-04
	    mk(	89)=-0.1609070508E-04
	    mk(	90)=-0.5074319542E-04
	    mk(	91)=	0.0001989441456
	    mk(	92)=-0.2746986030E-03
	    mk(	93)=	0.0002247755234
	    mk(	94)=-0.1117551322E-03
	    mk(	95)=	2.78139445E-005
	    mk(	96)=-0.1650973092E-04
	    mk(	97)=	2.255585476E-005
	    mk(	98)=-0.2901001452E-04
	    mk(	99)=	2.179430865E-005
	    mk(	100)=-0.1001771327E-04
	    mk(	101)=	8.80581103E-006
	    mk(	102)=-0.3339589332E-04
	    mk(	103)=	1.557357936E-008
	    mk(	104)=	2.618454061E-005
	    mk(	105)=-0.4469140921E-04
	    mk(	106)=	2.661110254E-005
	    mk(	107)=-0.3871237640E-04
	    mk(	108)=	0.0002236356442
	    mk(	109)=-0.2365730980E-03
	    mk(	110)=	6.226602689E-005
	    mk(	111)=	1.679391785E-006
	    mk(	112)=-0.5348922556E-04
	    mk(	113)=	0.0002305097825
	    mk(	114)=-0.2400842347E-03
	    mk(	115)=	6.842699935E-005
	    mk(	116)=-0.1182175862E-04
	    mk(	117)=	2.240200245E-006
	    mk(	118)=-0.5174987641E-04
	    mk(	119)=	0.0002304669954
	    mk(	120)=-0.2385694852E-03
	    mk(	121)=	6.840085604E-005
	    mk(	122)=-0.4986672689E-05
	    mk(	123)=-0.1669400647E-04
	    mk(	124)=	1.873465661E-005
	    mk(	125)=-0.2433643417E-04
	    mk(	126)=	2.657102842E-005
	    mk(	127)=-0.4614250040E-04
	    mk(	128)=	0.0001980460975
	    mk(	129)=-0.1975890376E-03
	    mk(	130)=	3.485954556E-005
	    mk(	131)=	6.32749399E-006
	    mk(	132)=-0.1648188110E-04
	    mk(	133)=	9.663241093E-006
	    mk(	134)=	2.461637308E-005
	    mk(	135)=-0.5769647750E-04
	    mk(	136)=	5.818428873E-005
	    mk(	137)=-0.2391347114E-04
	    mk(	138)=-0.1119725147E-04
	    mk(	139)=	2.022202635E-005
	    mk(	140)=-0.1403739152E-04
	    mk(	141)=-0.1165710249E-04
	    mk(	142)=	1.355163327E-005
	    mk(	143)=	1.52492754E-005
	    mk(	144)=-0.1169590712E-04
	    mk(	145)=-0.1528974344E-04
	    mk(	146)=	2.642432002E-005
	    mk(	147)=-0.2495404711E-04
	    mk(	148)=	2.746008004E-005
	    mk(	149)=-0.1475548351E-04
	    mk(	150)=-0.1390843013E-04
	    mk(	151)=	2.521323699E-005
	    mk(	152)=-0.1376991226E-04
	    mk(	153)=-0.1470449367E-04
	    mk(	154)=	2.82899601E-005
	    mk(	155)=-0.2187864294E-04
	    mk(	156)=	1.533923347E-005
	    mk(	157)=	4.361476917E-005
	    mk(	158)=-0.2324228631E-03
	    mk(	159)=	0.000246954308
	    mk(	160)=-0.6895262338E-04
	    mk(	161)=-0.1305734226E-04
	    mk(	162)=	7.785949066E-005
	    mk(	163)=-0.3386467637E-03
	    mk(	164)=	0.0007748249087
	    mk(	165)=-0.4090410153E-03
	    mk(	166)=-0.1353338039E-03
	    mk(	167)=	0.0004721961759
	    mk(	168)=-0.4939666099E-03
	    mk(	169)=	0.000469479856
	    mk(	170)=-0.4648155946E-03
	    mk(	171)=	0.0003573443869
	    mk(	172)=-0.1699676120E-03
	    mk(	173)=-0.7975126030E-04
	    mk(	174)=	0.0001329450659
	    mk(	175)=-0.3595832188E-04
	    mk(	176)=	1.038513955E-005
	    mk(	177)=-0.5230534190E-04
	    mk(	178)=	4.267573689E-005
	    mk(	179)=-0.5638982736E-04
	    mk(	180)=	0.0001429472027
	    mk(	181)=-0.1355180711E-03
	    mk(	182)=-0.1134252412E-03
	    mk(	183)=	0.0003031247883
	    mk(	184)=-0.4175554897E-03
	    mk(	185)=	0.0002303931365
	    mk(	186)=-0.2426309411E-03
	    mk(	187)=	0.0004081644099
	    mk(	188)=-0.3900467590E-03
	    mk(	189)=	9.016395449E-005
	    mk(	190)=	6.935937211E-005
	    mk(	191)=-0.3375169639E-04
	    mk(	192)=	2.543160308E-005
	    mk(	193)=-0.1076663665E-03
	    mk(	194)=	4.762030851E-005
	    mk(	195)=	4.960001629E-005
	    mk(	196)=-0.2099677794E-04
	    mk(	197)=-0.3301101301E-05
	    mk(	198)=-0.7651052891E-05
	    mk(	199)=-0.2711406005E-05
	    mk(	200)=-0.1775324660E-04
	    mk(	201)=	3.768255664E-005
	    mk(	202)=-0.1679410543E-03
	    mk(	203)=	0.0003553848668
	    mk(	204)=-0.2890638089E-03
	    mk(	205)=	0.0002737233476
	    mk(	206)=-0.3561114657E-03
	    mk(	207)=	0.0003673192324
	    mk(	208)=-0.4124951753E-03
	    mk(	209)=	0.0003310844259
	    mk(	210)=-0.3310944402E-04
	    mk(	211)=-0.2335613916E-03
	    mk(	212)=	0.0002396543562
	    mk(	213)=-0.7701317672E-04
	    mk(	214)=	3.338422991E-005
	    mk(	215)=-0.9185894828E-04
	    mk(	216)=	7.647561007E-005
	    mk(	217)=-0.2759584847E-04
	    mk(	218)=	3.501512179E-008
	    mk(	219)=-0.6400008719E-05
	    mk(	220)=-0.7786476221E-05
	    mk(	221)=	4.509850835E-006
	    mk(	222)=-0.4272364763E-04
	    mk(	223)=	0.0001336504311
	    mk(	224)=-0.1848693440E-03
	    mk(	225)=	0.0002447843083
	    mk(	226)=-0.2707702816E-03
	    mk(	227)=	0.0002577210848
	    mk(	228)=-0.2537834158E-03
	    mk(	229)=	0.0001952617005
	    mk(	230)=-0.3808913961E-04
	    mk(	231)=-0.7523675594E-04
	    mk(	232)=	0.0001057710681
	    mk(	233)=-0.1885282923E-03
	    mk(	234)=	0.0001297305257
	    mk(	235)=	0.0001188394073
	    mk(	236)=-0.1648188147E-03
	    mk(	237)=	4.317438174E-005
	    mk(	238)=-0.3904107566E-04
	    mk(	239)=	8.257237207E-005
	    mk(	240)=-0.8411750069E-04
	    mk(	241)=-0.1123723625E-04
	    mk(	242)=	9.918423324E-005
	    mk(	243)=	3.335442452E-006
	    mk(	244)=-0.1433911623E-03
	    mk(	245)=	0.0001346720585
	    mk(	246)=-0.2687995515E-04
	    mk(	247)=-0.5700657294E-04
	    mk(	248)=	5.433887629E-005
	    mk(	249)=-0.2314161276E-04
	    mk(	250)=	8.946712047E-006
	    mk(	251)=-0.4165567041E-04
	    mk(	252)=	0.0001284365424
	    mk(	253)=-0.1497577798E-03
	    mk(	254)=	9.574634766E-005
	    mk(	255)=-0.9787508130E-04
	    mk(	256)=	0.0001057708793
	    mk(	257)=-0.2690996341E-04
	    mk(	258)=-0.2690996341E-04

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

			  mn(	1)=-0.1245561126E-03
			  mn(	2)=-0.1245561126E-03
			  mn(	3)=-0.2012568965E-03
			  mn(	4)=-0.2086908594E-03
			  mn(	5)=-0.1732610296E-03
			  mn(	6)=-0.1049726901E-03
			  mn(	7)=	3.716833379E-005
			  mn(	8)=	0.0001424316554
			  mn(	9)=	0.0001738580517
			  mn(	10)=	6.387356337E-005
			  mn(	11)=-0.6859806460E-04
			  mn(	12)=-0.2093654918E-03
			  mn(	13)=-0.1120377185E-03
			  mn(	14)=	7.95959396E-006
			  mn(	15)=	5.737959546E-005
			  mn(	16)=-0.1511173547E-03
			  mn(	17)=-0.4153923482E-03
			  mn(	18)=-0.5691245305E-03
			  mn(	19)=-0.3998182133E-03
			  mn(	20)=-0.1732426697E-03
			  mn(	21)=	0.0003705188885
			  mn(	22)=	0.00039697138
			  mn(	23)=	0.0003506137853
			  mn(	24)=	0.0002274908354
			  mn(	25)=	0.0001715019514
			  mn(	26)=	8.020038773E-005
			  mn(	27)=	8.264325156E-005
			  mn(	28)=	3.649838617E-005
			  mn(	29)=	5.183298811E-005
			  mn(	30)=	2.181534156E-005
			  mn(	31)=	1.963400408E-005
			  mn(	32)=	2.875700484E-005
			  mn(	33)=	1.447066802E-005
			  mn(	34)=	1.51095061E-005
			  mn(	35)=	1.379696645E-005
			  mn(	36)=	5.865069427E-006
			  mn(	37)=	1.293621391E-005
			  mn(	38)=	9.913243222E-006
			  mn(	39)=	4.011450864E-006
			  mn(	40)=	9.04153167E-006
			  mn(	41)=	2.923914511E-006
			  mn(	42)=	6.670597012E-006
			  mn(	43)=	4.452439979E-006
			  mn(	44)=	3.379751789E-006
			  mn(	45)=	3.794329437E-006
			  mn(	46)=	4.036785639E-006
			  mn(	47)=	1.69392801E-007
			  mn(	48)=	6.062427409E-006
			  mn(	49)=-0.4585804443E-06
			  mn(	50)=	3.782796168E-006
			  mn(	51)=	2.122178375E-006
			  mn(	52)=	2.514431197E-007
			  mn(	53)=	2.58461856E-006
			  mn(	54)=	1.900263696E-006
			  mn(	55)=	1.110794106E-006
			  mn(	56)=	1.016002781E-006
			  mn(	57)=-0.1110017853E-06
			  mn(	58)=	4.06069729E-006
			  mn(	59)=-0.2109946649E-05
			  mn(	60)=	2.891307593E-006
			  mn(	61)=-0.3103421283E-06
			  mn(	62)=	3.419512263E-006
			  mn(	63)=-0.4220785397E-05
			  mn(	64)=	5.600058803E-006
			  mn(	65)=-0.2227528568E-05
			  mn(	66)=	1.308567299E-006
			  mn(	67)=	6.924733636E-007
			  mn(	68)=	1.311887598E-006
			  mn(	69)=-0.1946027768E-05
			  mn(	70)=	4.254207103E-006
			  mn(	71)=-0.2087087851E-05
			  mn(	72)=	1.761883482E-006
			  mn(	73)=-0.1155178119E-05
			  mn(	74)=	5.31829998E-007
			  mn(	75)=	2.360817671E-006
			  mn(	76)=-0.1552544322E-05
			  mn(	77)=	1.403638175E-006
			  mn(	78)=-0.4714864649E-06
			  mn(	79)=-0.1916647019E-05
			  mn(	80)=	5.716862821E-006
			  mn(	81)=-0.4872569539E-05
			  mn(	82)=	4.969801675E-006
			  mn(	83)=-0.5353010929E-05
			  mn(	84)=	8.064836021E-006
			  mn(	85)=-0.7727012389E-05
			  mn(	86)=	5.147519974E-006
			  mn(	87)=-0.1971880920E-06
			  mn(	88)=-0.6137218928E-06
			  mn(	89)=	6.737655991E-008
			  mn(	90)=-0.2159772240E-05
			  mn(	91)=	6.066784999E-006
			  mn(	92)=-0.6698422523E-05
			  mn(	93)=	6.240665985E-006
			  mn(	94)=-0.2893421642E-05
			  mn(	95)=	8.199605866E-007
			  mn(	96)=-0.9658626330E-06
			  mn(	97)=	2.504224748E-006
			  mn(	98)=-0.1476178525E-05
			  mn(	99)=	8.663161872E-007
			  mn(	100)=-0.2757184830E-06
			  mn(	101)=-0.9937820269E-07
			  mn(	102)=	3.809138601E-007
			  mn(	103)=	6.302161968E-007
			  mn(	104)=	6.177404837E-007
			  mn(	105)=-0.8269604926E-06
			  mn(	106)=	1.838722039E-007
			  mn(	107)=	1.953304203E-007
			  mn(	108)=-0.1314023050E-05
			  mn(	109)=	5.593663053E-006
			  mn(	110)=	3.420028262E-007
			  mn(	111)=-0.3802321191E-06
			  mn(	112)=	1.945082501E-006
			  mn(	113)=-0.5034963908E-06
			  mn(	114)=	1.061914818E-006
			  mn(	115)=-0.5658039066E-06
			  mn(	116)=	2.295890398E-006
			  mn(	117)=-0.1353052415E-05
			  mn(	118)=	2.294584904E-007
			  mn(	119)=	1.767192934E-006
			  mn(	120)=	2.059537394E-007
			  mn(	121)=	4.469433179E-007
			  mn(	122)=-0.3673122810E-06
			  mn(	123)=	2.818782999E-006
			  mn(	124)=-0.1897982404E-05
			  mn(	125)=	7.442075011E-007
			  mn(	126)=	1.847841678E-006
			  mn(	127)=-0.7544378103E-07
			  mn(	128)=	8.041819273E-007
			  mn(	129)=-0.2313260706E-06
			  mn(	130)=	2.681357973E-006
			  mn(	131)=-0.1669101473E-05
			  mn(	132)=	7.906170673E-007
			  mn(	133)=	1.295970625E-006
			  mn(	134)=	3.039897132E-006
			  mn(	135)=-0.3995633004E-05
			  mn(	136)=	3.0460976E-006
			  mn(	137)=	2.025860059E-006
			  mn(	138)=-0.2231455894E-05
			  mn(	139)=	3.515536593E-006
			  mn(	140)=-0.1893974097E-05
			  mn(	141)=	6.647143614E-007
			  mn(	142)=	1.700775702E-006
			  mn(	143)=	2.673521752E-006
			  mn(	144)=-0.1630506933E-05
			  mn(	145)=	3.035713364E-007
			  mn(	146)=	2.810269876E-006
			  mn(	147)=-0.4704328427E-06
			  mn(	148)=	1.444382724E-006
			  mn(	149)=	3.030136079E-007
			  mn(	150)=-0.2593088346E-06
			  mn(	151)=	2.973159615E-006
			  mn(	152)=	4.745381979E-007
			  mn(	153)=-0.2620766273E-05
			  mn(	154)=	6.20024768E-006
			  mn(	155)=-0.3617900137E-05
			  mn(	156)=	4.387847274E-006
			  mn(	157)=-0.4930832898E-06
			  mn(	158)=-0.2698844084E-06
			  mn(	159)=	3.548302401E-006
			  mn(	160)=	2.532251367E-007
			  mn(	161)=-0.2533038688E-05
			  mn(	162)=	1.169675679E-005
			  mn(	163)=-0.4229311773E-04
			  mn(	164)=	7.109217955E-005
			  mn(	165)=-0.5917844237E-04
			  mn(	166)=	1.528557243E-005
			  mn(	167)=	2.113164455E-006
			  mn(	168)=-0.5251909602E-06
			  mn(	169)=	1.667892348E-006
			  mn(	170)=-0.3440146113E-05
			  mn(	171)=	5.485588319E-006
			  mn(	172)=-0.6532299769E-05
			  mn(	173)=	6.005809968E-006
			  mn(	174)=-0.1735594453E-05
			  mn(	175)=	3.769913715E-007
			  mn(	176)=	6.307864669E-007
			  mn(	177)=-0.2820630331E-05
			  mn(	178)=	2.831958311E-006
			  mn(	179)=	2.020014881E-006
			  mn(	180)=-0.1982616624E-05
			  mn(	181)=	2.590724496E-006
			  mn(	182)=-0.4616440967E-05
			  mn(	183)=	4.361113175E-006
			  mn(	184)=	1.044879778E-006
			  mn(	185)=-0.6978464396E-06
			  mn(	186)=-0.7717524888E-06
			  mn(	187)=	1.377384647E-006
			  mn(	188)=	2.552012673E-006
			  mn(	189)=-0.4444099871E-05
			  mn(	190)=	4.035500645E-006
			  mn(	191)=	1.037337116E-006
			  mn(	192)=-0.1607164126E-05
			  mn(	193)=	2.257709291E-006
			  mn(	194)=-0.5938416448E-05
			  mn(	195)=	9.120221636E-006
			  mn(	196)=-0.3604777080E-05
			  mn(	197)=	2.230133555E-006
			  mn(	198)=-0.6767697440E-06
			  mn(	199)=-0.7632395082E-06
			  mn(	200)=	6.406535801E-007
			  mn(	201)=	3.387399765E-006
			  mn(	202)=-0.9098115769E-05
			  mn(	203)=	9.044961048E-006
			  mn(	204)=-0.1753262229E-05
			  mn(	205)=	3.900952182E-007
			  mn(	206)=	3.474521734E-007
			  mn(	207)=	2.609460355E-006
			  mn(	208)=-0.6493106182E-05
			  mn(	209)=	6.784664862E-006
			  mn(	210)=-0.3429290194E-05
			  mn(	211)=	3.874486809E-006
			  mn(	212)=-0.1310150847E-05
			  mn(	213)=-0.1716342439E-05
			  mn(	214)=	5.082704026E-006
			  mn(	215)=-0.8436910855E-05
			  mn(	216)=	5.938629742E-006
			  mn(	217)=	1.283718812E-006
			  mn(	218)=-0.1490539180E-05
			  mn(	219)=	1.589383751E-006
			  mn(	220)=-0.1819922668E-05
			  mn(	221)=	2.622593725E-006
			  mn(	222)=	1.138230663E-007
			  mn(	223)=-0.6127826509E-05
			  mn(	224)=	7.73851902E-006
			  mn(	225)=-0.3525366195E-05
			  mn(	226)=	3.317356911E-006
			  mn(	227)=-0.1868776662E-05
			  mn(	228)=	1.089536181E-006
			  mn(	229)=-0.2546114578E-06
			  mn(	230)=	2.075693316E-006
			  mn(	231)=-0.5962883145E-05
			  mn(	232)=	5.157775461E-006
			  mn(	233)=	6.370254464E-007
			  mn(	234)=-0.9733084086E-06
			  mn(	235)=	1.967705243E-007
			  mn(	236)=	1.855634162E-006
			  mn(	237)=-0.1383811813E-05
			  mn(	238)=	6.212574018E-007
			  mn(	239)=	3.40426093E-007
			  mn(	240)=-0.1931148351E-05
			  mn(	241)=-0.2854066021E-05
			  mn(	242)=	1.038482945E-005
			  mn(	243)=	1.914624707E-007
			  mn(	244)=-0.1423563482E-04
			  mn(	245)=	1.319591991E-005
			  mn(	246)=-0.1710125306E-05
			  mn(	247)=-0.9326241537E-05
			  mn(	248)=	1.859140527E-005
			  mn(	249)=-0.1486778801E-04
			  mn(	250)=	7.450749045E-007
			  mn(	251)=	8.872967398E-006
			  mn(	252)=-0.3440007895E-05
			  mn(	253)=	1.771570592E-006
			  mn(	254)=-0.6736671156E-05
			  mn(	255)=	5.549128839E-006
			  mn(	256)=-0.2082301635E-05
			  mn(	257)=-0.1730124334E-06
			  mn(	258)=-0.1730124334E-06
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
			  mk(	1)=-0.2184725378E-03
			  mk(	2)=-0.2184725378E-03
			  mk(	3)=-0.1009746558E-03
			  mk(	4)=	1.745580891E-005
			  mk(	5)=	9.013217706E-005
			  mk(	6)=	0.0001708878976
			  mk(	7)=	0.0001811027744
			  mk(	8)=	0.0001272354121
			  mk(	9)=-0.8629999345E-04
			  mk(	10)=-0.1512370398E-03
			  mk(	11)=-0.1999760698E-03
			  mk(	12)=-0.7199018306E-04
			  mk(	13)=	6.786485065E-005
			  mk(	14)=	4.239827116E-005
			  mk(	15)=-0.1247740365E-03
			  mk(	16)=-0.2515828531E-03
			  mk(	17)=-0.1859123197E-03
			  mk(	18)=	0.0001305039993
			  mk(	19)=	0.0004654615829
			  mk(	20)=	0.0006949400228
			  mk(	21)=	0.000535880421
			  mk(	22)=	0.0002512260722
			  mk(	23)=	1.532711276E-005
			  mk(	24)=-0.3195247352E-04
			  mk(	25)=-0.1210770415E-03
			  mk(	26)=-0.1730528312E-04
			  mk(	27)=-0.1310171955E-03
			  mk(	28)=	4.201119387E-005
			  mk(	29)=-0.1507348555E-03
			  mk(	30)=	2.536350854E-005
			  mk(	31)=	1.439024925E-005
			  mk(	32)=-0.9593952551E-04
			  mk(	33)=-0.2329699757E-04
			  mk(	34)=	1.137385227E-006
			  mk(	35)=-0.4181745294E-04
			  mk(	36)=	5.559796599E-005
			  mk(	37)=-0.7399659220E-04
			  mk(	38)=-0.5768552671E-04
			  mk(	39)=	7.871886153E-005
			  mk(	40)=-0.7287336315E-04
			  mk(	41)=-0.5588204713E-05
			  mk(	42)=	3.306056447E-005
			  mk(	43)=-0.5828671159E-04
			  mk(	44)=-0.5901823112E-05
			  mk(	45)=	4.118325351E-005
			  mk(	46)=-0.6747119137E-04
			  mk(	47)=	3.162579566E-005
			  mk(	48)=	4.459903188E-005
			  mk(	49)=-0.1514763427E-03
			  mk(	50)=	0.0001293496349
			  mk(	51)=-0.5661381330E-04
			  mk(	52)=-0.2642564237E-04
			  mk(	53)=	4.023087914E-005
			  mk(	54)=-0.5055937273E-04
			  mk(	55)=	4.482617587E-005
			  mk(	56)=-0.2580115079E-04
			  mk(	57)=-0.5454037546E-04
			  mk(	58)=	7.494068938E-005
			  mk(	59)=-0.7026625383E-04
			  mk(	60)=	9.805835428E-005
			  mk(	61)=-0.1822182243E-03
			  mk(	62)=	0.0002276274668
			  mk(	63)=-0.2079764961E-03
			  mk(	64)=	0.0001180701164
			  mk(	65)=-0.4692161956E-04
			  mk(	66)=-0.2735382143E-04
			  mk(	67)=	5.978241831E-005
			  mk(	68)=-0.2886741570E-04
			  mk(	69)=-0.3856270031E-04
			  mk(	70)=	8.952860388E-005
			  mk(	71)=-0.3872764839E-04
			  mk(	72)=-0.2773837181E-04
			  mk(	73)=-0.3240625925E-05
			  mk(	74)=-0.4682839908E-04
			  mk(	75)=	0.0001043358955
			  mk(	76)=-0.6790318180E-04
			  mk(	77)=	1.825016957E-005
			  mk(	78)=-0.2968524077E-04
			  mk(	79)=-0.3946653633E-04
			  mk(	80)=	0.0001067973708
			  mk(	81)=-0.5082207956E-04
			  mk(	82)=	1.462548411E-005
			  mk(	83)=-0.8785005922E-04
			  mk(	84)=	0.0001986725731
			  mk(	85)=-0.2443364409E-03
			  mk(	86)=	0.0001721776797
			  mk(	87)=-0.6112958120E-04
			  mk(	88)=-0.3524929943E-05
			  mk(	89)=-0.6914745129E-06
			  mk(	90)=-0.6667146362E-04
			  mk(	91)=	0.000194064391
			  mk(	92)=-0.2628679965E-03
			  mk(	93)=	0.0002147593076
			  mk(	94)=-0.1072421937E-03
			  mk(	95)=	2.671339631E-005
			  mk(	96)=-0.1572761982E-04
			  mk(	97)=	2.135330274E-005
			  mk(	98)=-0.2733358525E-04
			  mk(	99)=	2.006873707E-005
			  mk(	100)=-0.7315279057E-05
			  mk(	101)=	2.279410522E-007
			  mk(	102)=-0.1381154859E-05
			  mk(	103)=-0.8485344029E-06
			  mk(	104)=-0.1107589334E-06
			  mk(	105)=-0.5554785554E-06
			  mk(	106)=-0.1643503994E-05
			  mk(	107)=-0.3638884690E-06
			  mk(	108)=	3.646326051E-005
			  mk(	109)=-0.1328374414E-03
			  mk(	110)=	1.954190362E-005
			  mk(	111)=	8.078567103E-006
			  mk(	112)=	2.214381249E-005
			  mk(	113)=-0.1410858014E-03
			  mk(	114)=	0.0002561579763
			  mk(	115)=-0.2115479561E-03
			  mk(	116)=	1.169992928E-005
			  mk(	117)=	0.0001860548754
			  mk(	118)=-0.2061448555E-03
			  mk(	119)=	6.632145529E-005
			  mk(	120)=-0.3492763313E-04
			  mk(	121)=	2.752903188E-005
			  mk(	122)=-0.5018932075E-04
			  mk(	123)=	0.000201387471
			  mk(	124)=-0.2024121914E-03
			  mk(	125)=	4.07096111E-005
			  mk(	126)=-0.6213046265E-05
			  mk(	127)=	1.480154733E-005
			  mk(	128)=-0.1941022025E-04
			  mk(	129)=	1.882109066E-005
			  mk(	130)=-0.2028348040E-04
			  mk(	131)=	1.894512592E-005
			  mk(	132)=-0.1780677743E-04
			  mk(	133)=	9.552070267E-006
			  mk(	134)=	1.971480014E-005
			  mk(	135)=-0.4516569525E-04
			  mk(	136)=	3.435260853E-005
			  mk(	137)=	3.827277087E-005
			  mk(	138)=-0.2279975229E-03
			  mk(	139)=	0.0002359730529
			  mk(	140)=-0.7134831561E-04
			  mk(	141)=	9.762602466E-006
			  mk(	142)=-0.8310592001E-05
			  mk(	143)=	7.330437885E-005
			  mk(	144)=-0.2325381260E-03
			  mk(	145)=	0.0002107370061
			  mk(	146)=-0.5260881132E-04
			  mk(	147)=	5.464398353E-005
			  mk(	148)=-0.2045028466E-03
			  mk(	149)=	0.0002156634871
			  mk(	150)=-0.8708167100E-04
			  mk(	151)=	9.477879859E-005
			  mk(	152)=-0.2305662506E-03
			  mk(	153)=	0.0001890601274
			  mk(	154)=	2.836059556E-005
			  mk(	155)=-0.2407659185E-03
			  mk(	156)=	0.000299378404
			  mk(	157)=-0.2876279655E-03
			  mk(	158)=	0.0002067370639
			  mk(	159)=	2.304922047E-005
			  mk(	160)=-0.2289332745E-03
			  mk(	161)=	0.0002471709519
			  mk(	162)=-0.1952186896E-03
			  mk(	163)=-0.8363428677E-04
			  mk(	164)=	0.0006025169935
			  mk(	165)=-0.4270598614E-03
			  mk(	166)=	0.0003665920567
			  mk(	167)=-0.2758913618E-03
			  mk(	168)=	5.227445494E-005
			  mk(	169)=	0.0001796300758
			  mk(	170)=-0.2294777649E-03
			  mk(	171)=	0.0002434994648
			  mk(	172)=-0.2935776163E-03
			  mk(	173)=	0.0002703563238
			  mk(	174)=-0.2265061466E-03
			  mk(	175)=	0.0001019952717
			  mk(	176)=	0.000103684503
			  mk(	177)=-0.2448010275E-03
			  mk(	178)=	0.0002748521575
			  mk(	179)=-0.4767323190E-03
			  mk(	180)=	0.0003917978795
			  mk(	181)=-0.3724579994E-03
			  mk(	182)=	0.0002369276749
			  mk(	183)=-0.6676965816E-04
			  mk(	184)=	9.718158396E-006
			  mk(	185)=-0.9210702721E-05
			  mk(	186)=-0.4492788871E-05
			  mk(	187)=-0.8518499809E-05
			  mk(	188)=	3.098439266E-006
			  mk(	189)=-0.3885215752E-04
			  mk(	190)=	2.581202969E-005
			  mk(	191)=-0.2889695366E-05
			  mk(	192)=-0.4763846950E-04
			  mk(	193)=	0.0001595863175
			  mk(	194)=-0.2458283724E-03
			  mk(	195)=	0.000416801228
			  mk(	196)=-0.3587257086E-03
			  mk(	197)=	9.477398324E-005
			  mk(	198)=	8.707481517E-005
			  mk(	199)=-0.3333149163E-03
			  mk(	200)=	0.000371584864
			  mk(	201)=-0.3575582205E-03
			  mk(	202)=	0.0002103005432
			  mk(	203)=	2.250781923E-005
			  mk(	204)=-0.1095159532E-03
			  mk(	205)=	9.828610219E-005
			  mk(	206)=-0.3847922307E-04
			  mk(	207)=	2.268860268E-005
			  mk(	208)=-0.8495666013E-04
			  mk(	209)=	7.210329775E-005
			  mk(	210)=-0.2473547495E-04
			  mk(	211)=-0.4620578240E-05
			  mk(	212)=	1.168938078E-005
			  mk(	213)=-0.7291458205E-04
			  mk(	214)=	0.0002483833527
			  mk(	215)=-0.2908480887E-03
			  mk(	216)=	6.344881669E-005
			  mk(	217)=	0.0002064579096
			  mk(	218)=-0.2881794984E-03
			  mk(	219)=	0.000292800329
			  mk(	220)=-0.3017963686E-03
			  mk(	221)=	0.0002811063902
			  mk(	222)=-0.2597568199E-03
			  mk(	223)=	0.0001444815459
			  mk(	224)=	9.067872248E-006
			  mk(	225)=-0.2214577019E-04
			  mk(	226)=	4.828257986E-005
			  mk(	227)=-0.2009741102E-03
			  mk(	228)=	0.0001879765492
			  mk(	229)=-0.4908984277E-04
			  mk(	230)=-0.2038473680E-04
			  mk(	231)=	0.0001007653021
			  mk(	232)=-0.1148190884E-03
			  mk(	233)=	3.707988179E-005
			  mk(	234)=-0.6189515432E-04
			  mk(	235)=	0.0001820710515
			  mk(	236)=-0.2227488097E-03
			  mk(	237)=	0.0002150869394
			  mk(	238)=-0.2094660488E-03
			  mk(	239)=	0.000145411185
			  mk(	240)=-0.1432291149E-03
			  mk(	241)=	0.0001468414227
			  mk(	242)=-0.4518427879E-04
			  mk(	243)=	5.769432629E-006
			  mk(	244)=-0.6422902007E-05
			  mk(	245)=-0.8016803387E-05
			  mk(	246)=	1.077815432E-005
			  mk(	247)=-0.6213939901E-04
			  mk(	248)=	6.018391543E-005
			  mk(	249)=-0.5716523743E-04
			  mk(	250)=	0.0001419106992
			  mk(	251)=-0.1729332579E-03
			  mk(	252)=	0.0001645241643
			  mk(	253)=-0.1602238964E-03
			  mk(	254)=	0.0001047501143
			  mk(	255)=-0.1040425212E-03
			  mk(	256)=	0.0001080151449
			  mk(	257)=-0.2676916600E-04
			  mk(	258)=-0.2676916600E-04
c
		  end if

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
		   ref=cdsqrt(eps)

c---------------------------------------------------------------------


		   RETURN
		   END
c***************************************************************

c***************************************************************
c           subroutine Ag2021(wavel,Rm,ref,eps) 2021 July
c	Corrected a bag in version of 2015 (lambda > or < 600 nm)
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
c        deps1=(WL/WLP)**2*(1.D0/(1.D0+(GB*WL/WLP)**2)-
c     +        1.D0/(1.D0+(GR*WL/WLP)**2))
c        deps2=(WL/WLP)**2*((GR*WL/WLP)/(1.D0+(GR*WL/WLP)**2)-
c     +        (GB*WL/WLP)/(1.D0+(GB*WL/WLP)**2))
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
c     refre=dsqrt((dsqrt(eps1**2.d0+eps2**2.d0)+eps1)/2.d0)
c     refim=dsqrt((dsqrt(eps1**2.d0+eps2**2.d0)-eps1)/2.d0)
c------------------------------------------------------------------
      subroutine Ag2021(wavel,Rm,ref,eps)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MN,MK
      DIMENSION XN(100),YN(100),MN(100),XK(100),YK(100),MK(100)
	complex*16 ref,eps,im
c     
      xx=wavel
	if(xx.gt.600) go to 1 

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
	if(xx.gt.600.) then
		refre=0.0755 - 0.000155*xx + 2.15E-7*xx*xx
		go to 29
	end if

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


29        J=2
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
       eps2=2.d0*refre*refim
C==================================================================
c      Step 2: Calculation of coorection terms
C    Ïîïðàâêà è ðåàëüíîé è ìíèìîé ÷àñòè EPSILON=EPS1+iEPS2
c      WL=WAVEL
c      WLP=136.1D0
c      GB=0.0019D0
c      pi=4.d0*datan(1.d0)
c      GR=GB+0.0047*136.1/(2.0*pi*Rm)
c        deps1=(WL/WLP)**2*(1.D0/(1.D0+(GB*WL/WLP)**2)-
c     +        1.D0/(1.D0+(GR*WL/WLP)**2))
c        deps2=(WL/WLP)**2*((GR*WL/WLP)/(1.D0+(GR*WL/WLP)**2)-
c     +        (GB*WL/WLP)/(1.D0+(GB*WL/WLP)**2))
c*********************************************
C        DEPS1=(OMPB/OM)**2*(1.D0/(1.D0+(GAMB/OM)**2)-
C     +        1.D0/(1.D0+(GAMJ/OM)**2))
C        DEPS2=(OMPB/OM)**2*((GAMJ/OM)/(1.D0+(GAMJ/OM)**2)-
C     +        (GAMB/OM)/(1.D0+(GAMB/OM)**2))

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
       deps1=0.d0
	 As=1.d0
	 deps2=2.5d0*(As/Rm)*(wavel/400.d0)**3.d0 
       
       eps1=eps1  + deps1
       eps2=eps2  + deps2
c---------------------------------------------------------------------
c      Step 3: Recalculation of tefre and refim: Bohren-Huffman, p.281
       refre=dsqrt((dsqrt(eps1**2.d0+eps2**2.d0)+eps1)/2.d0)
       refim=dsqrt((dsqrt(eps1**2.d0+eps2**2.d0)-eps1)/2.d0)

	im=(0.d0,1.d0)
	eps=eps1+im*eps2
	ref=refre+im*refim
c---------------------------------------------------------------------
       RETURN
       END
c***************************************


C--------------------------------------------------------

       subroutine g34JCh3(wavel,Rm,ref,eps,KD)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MN,MK
      DIMENSION XN(100),YN(100),MN(100),XK(100),YK(100),MK(100)

	complex*16 ref,eps,epsb,Deps,aim
C--------------------------------------------------------
c	Rm=aLeff/As !!!!!

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