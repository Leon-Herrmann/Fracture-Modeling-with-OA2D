C==============================Start of commonblock
C==============================modified from Lars Mikkelsen
	  integer nig   ! Number of integration points
      parameter(nig=9)   
C
      integer maxelem,elemshift
      parameter(maxelem=80)
      real*8 Emodout(nig,maxelem)
      real*8 strain11out(nig,maxelem)
	  real*8 strain22out(nig,maxelem)
	  real*8 strain12out(nig,maxelem)
	  real*8 strain13out(nig,maxelem)
	  real*8 strain23out(nig,maxelem)
	  real*8 strain33out(nig,maxelem)
	  real*8 stress11out(nig,maxelem)
	  real*8 stress22out(nig,maxelem)
	  real*8 stress12out(nig,maxelem)
	  real*8 stress13out(nig,maxelem)
	  real*8 stress23out(nig,maxelem)
	  real*8 stress33out(nig,maxelem)
	  real*8 stressavg11out(nig,maxelem)
	  real*8 stressavg22out(nig,maxelem)
	  real*8 stressavg12out(nig,maxelem)
	  real*8 stressavg13out(nig,maxelem)
	  real*8 stressavg23out(nig,maxelem)
	  real*8 stressavg33out(nig,maxelem)
	  real*8 elemout(maxelem)
      common/kkoutint/elemshift
      common/kkoutvar/Emodout,elemout,
     &       stress11out,stress22out,stress12out,stress13out,stress23out,stress33out,
     &       strain11out,strain22out,strain12out,strain13out,strain23out,strain33out,
     &       stressavg11out,stressavg22out,stressavg12out,stressavg13out,stressavg23out,stressavg33out 
C==============================End of commonblock
