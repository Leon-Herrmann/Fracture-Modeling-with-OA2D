C==============================User element CPTC8
C==============================Start of user element subroutine
C==============================modified from
C==http://madamcode.blogspot.com/2017/02/simple-uel-for-abaqus-plane-stress-2d.html
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3 NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4 PERIOD)
	 

C  
C==============================Allocation of variables/arrays and initialization
      IMPLICIT REAL*8 (A-H,O-Z) 
C
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1 SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2 U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3 PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5 JPROPS(*)
C	 

	 
      REAL*8 rg(9),sg(9),wrg(9),wsg(9),
     1  DNDR_DS(2,8),JAC(2,2),JINV(2,2),DNDX_DY(2,8),D(6,6),G,T(6,6),
     2  B(6,NDOFEL),P(NDOFEL,NDOFEL),C1,S1,G12,
     3  STRESS(6),STRAIN(6),JDET,tm,CONST,GP,atemp,btemp,TRIANGLE,GAMMA130,GAMMA230,
     4  c_angle,T_sigma(6,6),T_epsilon(6,6),E,NU,E1,NU23,G23,E2,NU13,G13,E3,NU12,E110,E330,GAMMA120
C
	  INTEGER ig,i
      include 'koutvar.inc.f' 
C	  
	  data elemshift /0/  ! Set variable equal to 0 as initial value
C	  	  	  
C==============================RHS = Load Vector
C==============================AMATRX = STIFFNESS MATRIX
C==============================SVARS = Solution dependent variables like stress and strain
C
C==============================import material properties from inputfile
	  E1    = PROPS(1)
	  NU23  = PROPS(2)
	  G23   = PROPS(3)	
	  E2    = PROPS(4)
	  NU13  = PROPS(5)
	  G13   = PROPS(6)
	  E3    = PROPS(7)
	  NU12  = PROPS(8)
	  G12   = PROPS(9) 
	  tm    = PROPS(10)
	  tB	= PROPS(11)
C
C==============================initializing AMATRX and RHS to 0
       RHS(:,1) = 0.d0
       AMATRX = 0.d0
C	
C==============================Define D the constitutive matrix (for plane stress analysis)
C    
	D=0.d0
C	
	CONST=-2.0d0*E2*E3*NU12*NU13*NU23-E1*E3*NU23*NU23-E2*E2*NU12*NU12-E2*E3*NU13*NU13+E2*E1
C	
	D(1,1)=E1*E1*(-E3*NU23*NU23+E2)/CONST
	D(1,2)=E1*E2*(E3*NU13*NU23+NU12*E2)/CONST
	D(1,6)=E1*(NU12*NU23+NU13)*E3*E2/CONST
	D(2,1)=E1*E2*(E3*NU13*NU23+NU12*E2)/CONST
	D(2,2)=(-NU13*NU13*E3+E1)*E2*E2/CONST
	D(2,6)=(NU12*NU13*E2+NU23*E1)*E3*E2/CONST
	D(6,1)=E1*(NU12*NU23+NU13)*E3*E2/CONST
	D(6,2)=(NU12*NU13*E2+NU23*E1)*E3*E2/CONST
	D(6,6)=E3*(-NU12*NU12*E2+E1)*E2/CONST
	D(3,3)=G12
	D(4,4)=G13
	D(5,5)=G23
C	
	C1=cos(tm)
	S1=sin(tm)
C    
C==============================Transformation matrix MatrixInverse(T_t_2_sigma)
	T_sigma=0.d0
	
	T_sigma(1,1)=C1*C1
	T_sigma(1,4)=2.d0*S1*C1
	T_sigma(1,6)=S1*S1
C	
	T_sigma(2,2)=1.d0
C	
	T_sigma(3,3)=C1
	T_sigma(3,5)=S1
C	
	T_sigma(4,1)=-S1*C1
	T_sigma(4,4)=2.d0*C1*C1-1.d0
	T_sigma(4,6)=S1*C1
C	
	T_sigma(5,3)=-S1
	T_sigma(5,5)=C1
C	
	T_sigma(6,1)=S1*S1
	T_sigma(6,4)=-2.d0*S1*C1
	T_sigma(6,6)=C1*C1
C
C==============================Transformation matrix T_t_2_epsilon
	T_epsilon=0.d0
C	
	T_epsilon(1,1)=C1*C1
	T_epsilon(1,4)=-S1*C1
	T_epsilon(1,6)=S1*S1
C	
	T_epsilon(2,2)=1.d0
C	
	T_epsilon(3,3)=C1
	T_epsilon(3,5)=-S1
C	
	T_epsilon(4,1)=2.d0*S1*C1
	T_epsilon(4,4)=C1*C1-S1*S1
	T_epsilon(4,6)=-2.d0*S1*C1
C	
	T_epsilon(5,3)=S1
	T_epsilon(5,5)=C1
C	
	T_epsilon(6,1)=S1*S1
	T_epsilon(6,4)=S1*C1
	T_epsilon(6,6)=C1*C1
C	
C==============================Transform constitutive matrix	
	D=matmul(T_sigma,matmul(D,T_epsilon))
C		
		IF(NDOFEL==8*3+3) THEN
			TRIANGLE=1
		ELSE IF(NDOFEL==6*3+3) THEN
			TRIANGLE=0
		ENDIF
C==================================================================================Quadrilateral element
		IF(TRIANGLE==1) THEN
C==============================Define Gauss points and weighting factors using Abaqus order
		  GP=0.7745966692
		  atemp=0.555555555555555555556
		  btemp=0.888888888888888888889
		  rg(:)=(/-GP,0.0d0,GP,-GP,0.0d0,GP,-GP,0.0d0,GP/)					!xi coordinates
		  wrg(:)=(/atemp,btemp,atemp,atemp,btemp,atemp,atemp,btemp,atemp/)	!weighting factors xi
		  sg(:)=(/-GP,-GP,-GP,0.0d0,0.0d0,0.0d0,GP,GP,GP/)					!eta coordinates
		  wsg(:)=(/atemp,atemp,atemp,btemp,btemp,btemp,atemp,atemp,atemp/)	!weighting factors eta
C      
C==============================reset output
			stressavg11out(1,jelem)=0.
			stressavg22out(1,jelem)=0.
			stressavg12out(1,jelem)=0.
			stressavg13out(1,jelem)=0.
			stressavg23out(1,jelem)=0.
			stressavg33out(1,jelem)=0.
C==============================Loop over Gauss points
		ig=1
		  DO ig=1,9  
C         
C==============================Derivative of shape functions in respect to dr=xi, ds=eta coordinates 
C==============================and evaluated at the current Gauss point
C		 
			 dndr_ds(1,1)=-(sg(ig)-1)*(sg(ig)+2*rg(ig))/4.0d0
			 dndr_ds(1,2)=(sg(ig)-2*rg(ig))*(sg(ig)-1)/4.0d0
			 dndr_ds(1,3)=(sg(ig)+2*rg(ig))*(1+sg(ig))/4.0d0
			 dndr_ds(1,4)=-(sg(ig)-2*rg(ig))*(1+sg(ig))/4.0d0
			 dndr_ds(1,5)=(sg(ig)-1)*rg(ig)
			 dndr_ds(1,6)=-sg(ig)*sg(ig)*0.5d0+0.5d0
			 dndr_ds(1,7)=-(1+sg(ig))*rg(ig)
			 dndr_ds(1,8)=sg(ig)*sg(ig)*0.5d0-0.5d0
C		 
			 dndr_ds(2,1)=-(rg(ig)-1)*(rg(ig)+2*sg(ig))/4.0d0
			 dndr_ds(2,2)=-(rg(ig)+1)*(rg(ig)-2*sg(ig))/4.0d0
			 dndr_ds(2,3)=(rg(ig)+1)*(rg(ig)+2*sg(ig))/4.0d0
			 dndr_ds(2,4)=(rg(ig)-1)*(rg(ig)-2*sg(ig))/4.0d0
			 dndr_ds(2,5)=rg(ig)*rg(ig)*0.5d0-0.5d0
			 dndr_ds(2,6)=-(rg(ig)+1)*sg(ig)
			 dndr_ds(2,7)=-rg(ig)*rg(ig)*0.5d0+0.5d0
			 dndr_ds(2,8)=(rg(ig)-1)*sg(ig)
C         
C==============================Calculate the Jacobian Matrix = JAC
			 JAC = 0.d0
			  DO i=1,8
				JAC(1,1) = JAC(1,1) + dndr_ds(1,i)*COORDS(1,i)
				JAC(1,2) = JAC(1,2) + dndr_ds(1,i)*COORDS(2,i)
				JAC(2,1) = JAC(2,1) + dndr_ds(2,i)*COORDS(1,i)
				JAC(2,2) = JAC(2,2) + dndr_ds(2,i)*COORDS(2,i)
			  END DO
C==============================Determinant of the Jacobian Matrix = JDET
			JDET = JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)
C 
C==============================Inverse of Jacobian = JINV
			 JINV(1,1) =  JAC(2,2)/JDET
			 JINV(1,2) = -JAC(1,2)/JDET
			 JINV(2,1) = -JAC(2,1)/JDET
			 JINV(2,2) =  JAC(1,1)/JDET
C		 
C==============================Derivative of shape functions with respect to x, y coordinates
			DNDX_DY=matmul(JINV,DNDR_DS)
C		 
C==============================Define strain displacement matrix = B		 
			 B=0.d0
			 c_angle=cos(tB)/sin(tB)
			 DO ib=0,7
			 	B(1,1+3*ib)=DNDX_DY(1,1+ib);	B(3,1+3*ib)=DNDX_DY(2,1+ib);	B(4,1+3*ib)=-c_angle*DNDX_DY(1,1+ib);
				B(2,2+3*ib)=DNDX_DY(2,1+ib);	B(3,2+3*ib)=DNDX_DY(1,1+ib);	B(5,2+3*ib)=-c_angle*DNDX_DY(1,1+ib);
				B(4,3+3*ib)=DNDX_DY(1,1+ib);	B(5,3+3*ib)=DNDX_DY(2,1+ib);	B(6,3+3*ib)=-c_angle*DNDX_DY(1,1+ib);
			 END DO
C	
			B(4,NDOFEL-2)=c_angle;
			B(5,NDOFEL-1)=0.5d0*c_angle;  
			B(6,NDOFEL)=1;
C		 	
C==============================Calculate the integrand = P
			P=matmul(matmul(transpose(B),D),B)*JDET
C		 
C==============================Calculate the stiffness matrix = AMATRX
			AMATRX = AMATRX +P*wrg(ig)*wsg(ig)
C
C==============================Calculate STRAIN and STRESS
			STRAIN=MATMUL(B,U)
			STRESS=MATMUL(D,STRAIN)
C
C==============================Save output to common block
			strain11out(ig,jelem)=STRAIN(1)
			strain22out(ig,jelem)=STRAIN(2)
			strain12out(ig,jelem)=STRAIN(3)
			strain13out(ig,jelem)=STRAIN(4)
			strain23out(ig,jelem)=STRAIN(5)
			strain33out(ig,jelem)=STRAIN(6)
			stress11out(ig,jelem)=STRESS(1)
			stress22out(ig,jelem)=STRESS(2)
			stress12out(ig,jelem)=STRESS(3)
			stress13out(ig,jelem)=STRESS(4)
			stress23out(ig,jelem)=STRESS(5)
			stress33out(ig,jelem)=STRESS(6)
C
			stressavg11out(1,jelem)=stressavg11out(1,jelem)+STRESS(1)*wrg(ig)*wsg(ig)*JDET
			stressavg22out(1,jelem)=stressavg22out(1,jelem)+STRESS(2)*wrg(ig)*wsg(ig)*JDET
			stressavg12out(1,jelem)=stressavg12out(1,jelem)+STRESS(3)*wrg(ig)*wsg(ig)*JDET
			stressavg13out(1,jelem)=stressavg13out(1,jelem)+STRESS(4)*wrg(ig)*wsg(ig)*JDET
			stressavg23out(1,jelem)=stressavg23out(1,jelem)+STRESS(5)*wrg(ig)*wsg(ig)*JDET
			stressavg33out(1,jelem)=stressavg33out(1,jelem)+STRESS(6)*wrg(ig)*wsg(ig)*JDET
C
			elemshift=max(elemshift,jelem)
			elemout(jelem)=jelem
C		
			END DO
C==================================================================================Triangle element
		ELSE IF(TRIANGLE==0) THEN
C==============================Define Gauss points and weighting factors
		  GP=(1.0d0/(3.0d0))				!Gauss point constant
		  rg(1:3)=(/2*GP,0.5*GP,0.5*GP/)				!xi coordinates
		  wrg(1:3)=(/(1.d0/6.d0),(1.d0/6.d0),(1.d0/6.d0)/)	!weighting factors xi
		  sg(1:3)=(/0.5*GP,2*GP,0.5*GP/)				!eta coordinates
		  wsg(1:3)=(/(1.d0/6.d0),(1.d0/6.d0),(1.d0/6.d0)/)	!weighting factors eta
C      
C==============================reset output
			stressavg11out(1,jelem)=0.
			stressavg22out(1,jelem)=0.
			stressavg12out(1,jelem)=0.
			stressavg13out(1,jelem)=0.
			stressavg23out(1,jelem)=0.
			stressavg33out(1,jelem)=0.
C==============================Loop over Gauss points
		  DO ig=1,3  
C         
C==============================Derivative of shape functions in respect to dr=xi, ds=eta coordinates 
C==============================and evaluated at the current Gauss point
			 dndr_ds(1,1) =  4.D0*rg(ig)-1.D0					!N1_dr
			 dndr_ds(1,2) =  0.D0								!N2_dr
			 dndr_ds(1,3) =  -3.D0+4.D0*sg(ig)+4.D0*rg(ig)		!N3_dr
			 dndr_ds(1,4) =  4.D0*sg(ig)
			 dndr_ds(1,5) =  -4.D0*sg(ig)
			 dndr_ds(1,6) =  -8.D0*rg(ig)+4.D0-4.D0*sg(ig)
C												
			 dndr_ds(2,1) =  0.D0								!N1_ds
			 dndr_ds(2,2) =  4.D0*sg(ig)-1.D0					!N2_ds
			 dndr_ds(2,3) =  -3.D0+4.D0*sg(ig)+4.D0*rg(ig)		!N3_ds
			 dndr_ds(2,4) =  4.D0*rg(ig)
			 dndr_ds(2,5) =  -8.D0*sg(ig)+4.D0-4.D0*rg(ig)
			 dndr_ds(2,6) =  -4.D0*rg(ig)
C         
C==============================Calculate the Jacobian Matrix = JAC
			 JAC = 0.d0
			  DO i=1,6
				JAC(1,1) = JAC(1,1) + dndr_ds(1,i)*COORDS(1,i)
				JAC(1,2) = JAC(1,2) + dndr_ds(1,i)*COORDS(2,i)
				JAC(2,1) = JAC(2,1) + dndr_ds(2,i)*COORDS(1,i)
				JAC(2,2) = JAC(2,2) + dndr_ds(2,i)*COORDS(2,i)
			  END DO
C==============================Determinant of the Jacobian Matrix = JDET			 
			 JDET = JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)
C 
C==============================Inverse of Jacobian = JINV  
			 JINV(1,1) =  JAC(2,2)/JDET
			 JINV(1,2) = -JAC(1,2)/JDET
			 JINV(2,1) = -JAC(2,1)/JDET
			 JINV(2,2) =  JAC(1,1)/JDET
C		 
C==============================Derivative of shape functions with respect to x, y coordinates
			 DNDX_DY=matmul(JINV,DNDR_DS)
C		 
C==============================Define strain displacement matrix = B		 
			 B=0.d0
			 c_angle=cos(tB)/sin(tB)
			 DO ib=0,5
				B(1,1+3*ib)=DNDX_DY(1,1+ib);	B(3,1+3*ib)=DNDX_DY(2,1+ib);	B(4,1+3*ib)=-c_angle*DNDX_DY(1,1+ib);
				B(2,2+3*ib)=DNDX_DY(2,1+ib);	B(3,2+3*ib)=DNDX_DY(1,1+ib);	B(5,2+3*ib)=-c_angle*DNDX_DY(1,1+ib);
				B(4,3+3*ib)=DNDX_DY(1,1+ib);	B(5,3+3*ib)=DNDX_DY(2,1+ib);	B(6,3+3*ib)=-c_angle*DNDX_DY(1,1+ib);
			 END DO
C		 
			B(4,NDOFEL-2)=c_angle;
			B(5,NDOFEL-1)=0.5d0*c_angle;  
			B(6,NDOFEL)=1;	 
C==============================Calculate the integrand = P
			 P=matmul(matmul(transpose(B),D),B)*JDET
C		 
C==============================Calculate the stiffness matrix = AMATRX
			AMATRX = AMATRX+P*wrg(ig)
C
C==============================Calculate STRAIN and STRESS
			STRAIN=MATMUL(B,U)
			STRESS=MATMUL(D,STRAIN)
C
C==============================Save output to common block
			strain11out(ig,jelem)=STRAIN(1)
			strain22out(ig,jelem)=STRAIN(2)
			strain12out(ig,jelem)=STRAIN(3)
			strain13out(ig,jelem)=STRAIN(4)
			strain23out(ig,jelem)=STRAIN(5)
			strain33out(ig,jelem)=STRAIN(6)
			stress11out(ig,jelem)=STRESS(1)
			stress22out(ig,jelem)=STRESS(2)
			stress12out(ig,jelem)=STRESS(3)
			stress13out(ig,jelem)=STRESS(4)
			stress23out(ig,jelem)=STRESS(5)
			stress33out(ig,jelem)=STRESS(6)
C
			stressavg11out(1,jelem)=stressavg11out(1,jelem)+STRESS(1)*wrg(ig)*JDET
			stressavg22out(1,jelem)=stressavg22out(1,jelem)+STRESS(2)*wrg(ig)*JDET
			stressavg12out(1,jelem)=stressavg12out(1,jelem)+STRESS(3)*wrg(ig)*JDET
			stressavg13out(1,jelem)=stressavg13out(1,jelem)+STRESS(4)*wrg(ig)*JDET
			stressavg23out(1,jelem)=stressavg23out(1,jelem)+STRESS(5)*wrg(ig)*JDET
			stressavg33out(1,jelem)=stressavg33out(1,jelem)+STRESS(6)*wrg(ig)*JDET
C
			elemshift=max(elemshift,jelem)
			elemout(jelem)=jelem
C		
		   END DO
C	
		END IF
C	  
C==============================Calculate the right hand side = RHS
		RHS(:,1) = -MATMUL((AMATRX),U)
		write(6,*)'OBSRHS',RHS
		write(6,*)'OBSU',U
C		
	END
C
C==============================End of user element subroutine
C
C==============================Start of user material subroutine
C==============================modified from Lars Mikkelsen
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC) 
C
      INCLUDE 'ABA_PARAM.INC' 
      include 'koutvar.inc.f' 
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3)
C	 
C==============================Jacobian matrix of constitutive matrix is set to 0
C==============================Stiffness matrix is 0
      do k1=1,ntens
       do k2=1,ntens
        DDSDDE(k1,k2)=0.0d0
       enddo
      enddo
C      
C==============================Assigning stresses and strains to statevariables
	    STATEV(1)=strain11out(npt,noel-elemshift)
	    STATEV(2)=strain22out(npt,noel-elemshift)
	    STATEV(3)=strain12out(npt,noel-elemshift)
		STATEV(4)=strain13out(npt,noel-elemshift)
	    STATEV(5)=strain23out(npt,noel-elemshift)
	    STATEV(6)=strain33out(npt,noel-elemshift)
	    STATEV(7)=stress11out(npt,noel-elemshift)
	    STATEV(8)=stress22out(npt,noel-elemshift)
	    STATEV(9)=stress12out(npt,noel-elemshift)
		STATEV(10)=stress13out(npt,noel-elemshift)
	    STATEV(11)=stress23out(npt,noel-elemshift)
	    STATEV(12)=stress33out(npt,noel-elemshift)
		STATEV(13)=stressavg11out(npt,noel-elemshift)
	    STATEV(14)=stressavg22out(npt,noel-elemshift)
	    STATEV(15)=stressavg12out(npt,noel-elemshift)
		STATEV(16)=stressavg13out(npt,noel-elemshift)
	    STATEV(17)=stressavg23out(npt,noel-elemshift)
	    STATEV(18)=stressavg33out(npt,noel-elemshift)
		
C	
C==============================Control check
	   if (elemout(noel-elemshift).ne.noel-elemshift) then 
	        call xit
	   end if
C            
      RETURN
      END
C	  
C==============================End of user material subroutine
C
C
C==============================Start of DISP subroutine for shear case in stiffness calculation
C
      SUBROUTINE DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2),COORDS(3)
C		
		U(1)=U(1)*COORDS(2)
C
		RETURN
		END
C==============================End of DISP subroutine
		
		
		