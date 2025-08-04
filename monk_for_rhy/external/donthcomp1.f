      SUBROUTINE donthcomp1(Ear,Ne,Param,Ifl,Photar,Photer)
      implicit none
      INTEGER Ne, Ifl
      REAL Param(*), Ear(0:Ne) , Photar(Ne), Photer(Ne), prim(0:Ne)
c
c     driver for the Comptonization code solving Kompaneets equation
c     seed photons - (disc) blackbody
c     reflection + Fe line with smearing
c
c     Version optimized for a number of data files but with the same values
c     of parameters:
c
c     number of model parameters: 16
c     1: photon spectral index
c     2: plasma temperature in keV
c     3: (disc)blackbody temperature in keV
c     4: type of seed spectrum (0-blackbody, 1-diskbb)
c     5: redshift
c     6: whether to output total spectrum or scattered spectrum
c           (0-total; 1-scattered)

      INTEGER n , ierr, np
      REAL xninv , normfac , normlum, SPP, pa0(6)

      integer nth
      REAL xth(900), spt(900), sca(900), spout(900)

      LOGICAL recalc

      real    zfactor
      real outopt
      integer i,j, jl

      SAVE pa0,normfac,normlum,xth,nth,spt,sca,spout

      DATA pa0/6*9999./

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         Photer(i) = 0.0
      ENDDO
c     xtot is the energy array (units m_e c^2)
c     spnth is the nonthermal spectrum alone (E F_E)
c     sptot is the total spectrum array (E F_E), = spref if no reflection

      ierr = 0
      zfactor = 1.0 + param(5)
      outopt = param(6)

c  calculate internal source spectrum if input parameters have changed

      np = 6
      recalc = .FALSE.
      DO n = 1 , np
         IF ( Param(n).NE.pa0(n) ) recalc = .TRUE.
      END DO
      IF ( recalc ) THEN

c       if (param(4) .lt. 0.5) then
      call thCompton(param(3)/511.,param(2)/511.,param(1),
     >                     xth,nth,spt,sca)
c         else
c            call thCompton(param(3)/511.,param(2)/511.,param(1),
c     >           xth,nth,spt)
c         end if

         xninv = 511.0/zfactor
         normfac = 1/SPP(xninv,xth,nth,spt)
         normfac = 1.

c Calculate luminosity normalization  (used in another model!)

         normlum = 0.0
         do i = 2, nth - 1
            normlum = normlum +
     >          0.5 * (spt(i)/xth(i) + spt(i-1)/xth(i-1))
     >          * (xth(i) - xth(i-1))
         enddo
         normlum = normlum * normfac

      end if

c     zero arrays

      do i=1,ne
         photar(i) = 0
         prim(i) = 0
      end do
      prim(0) = 0

c     put primary into final array only if scale >= 0.
      if (outopt < 0.5) then
        do i=1,900
            spout(i) = spt(i)
        enddo
      else
        do i=1,900
            spout(i) = sca(i)
        enddo
      endif

      j = 1
      do i=0,ne
         do while (j .le. nth .and. 511*xth(j) .lt. ear(i)*zfactor)
            j = j + 1
         end do
         if (j. le. nth) then
            if (j .gt. 1) then
               jl = j - 1
               prim(i) = spout(jl)+(ear(i)/511*zfactor - xth(jl))*
     >                 (spout(jl+1)-spout(jl))/
     >                 (xth(jl+1)-xth(jl))
            else
               prim(i) = spout(1)
            endif
         endif
      end do
      do i=1,ne
         photar(i) = 0.5*(prim(i)/ear(i)**2+prim(i-1)/ear(i-1)**2)
     >           *(ear(i)-ear(i-1))*normfac
      end do


      RETURN
      END


      subroutine thCompton(tempbb,theta,gamma,x,jmax,sptot,spsca)
c  version: January 96
c
c
c     Thermal Comptonization; solves Kompaneets eq. with some
c     relativistic corrections. See Lightman & Zdziarski (1987), ApJ
c     The seed spectrum is blackbody.
      IMPLICIT NONE
      real delta,xmin,xmax,deltal,xnr,xr,taukn,arg,flz,planck,pi
      real sptot(900),x(900),dphesc(900),dphdot(900)
     >  ,rel(900),bet(900),c2(900),spsca(900)
      integer j,jmax,jmaxth,jnr,jrel
      DOUBLE PRECISION w,w1,z1,z2,z3,z4,z5,z6
      real tautom
c  input parameters:
      real tempbb,theta,gamma

c use internally Thomson optical depth
      tautom=sqrt(2.25+3/(theta*((gamma+.5)**2-2.25)))-1.5
c
      pi=3.14159
c clear arrays (important for repeated calls)
      do 31 j=1,900
        dphesc(j)=0
        dphdot(j)=0
        rel(j)=0
        bet(j)=0
        c2(j)=0
        sptot(j)=0
  31  continue
C
C JMAX - # OF PHOTON ENERGIES
C
c delta is the 10-log interval of the photon array.
      delta=0.02
      deltal=delta*log(10.)
      xmin=1e-4*tempbb
      xmax=40*theta
      jmax=min(899,int(log10(xmax/xmin)/delta)+1)
C
C X - ARRAY FOR PHOTON ENERGIES
C
      do 4 j=1,jmax+1
         x(j)=xmin*10.**((j-1)*delta)
  4   continue
c
c compute c2(x), and rel(x) arrays
c c2(x) is the relativistic correction to Kompaneets equation
c rel(x) is the Klein-Nishina cross section divided by the Thomson crossection
      do 500 j=1,jmax
        w=x(j)
c c2 is the Cooper's coefficient calculated at w1
c w1 is x(j+1/2) (x(i) defined up to jmax+1)
        w1=sqrt(x(j)*x(j+1))
        c2(j)=SNGL(w1**4/(1+4.6*w1+1.1*w1*w1))
        if (w.le.0.05) then
c use asymptotic limit for rel(x) for x less than 0.05
          rel(j)=SNGL(1-2*w+26*w*w/5)
        else
          z1=(1+w)/w**3
          z2=1+2*w
          z3=log(z2)
          z4=2*w*(1+w)/z2
          z5=z3/2/w
          z6=(1+3*w)/z2/z2
          rel(j)=SNGL(0.75*(z1*(z4-z3)+z5-z6))
        end if
 500  continue

c the thermal emision spectrum
        jmaxth=min(900,int(log10(50*tempbb/xmin)/delta))
        if (jmaxth .gt. jmax) then
c           print *,'thcomp: ',jmaxth,jmax
           jmaxth = jmax
        end if
        planck=15/(pi*tempbb)**4
        do 5 j=1,jmaxth
          dphdot(j)=planck*x(j)**2/(exp(x(j)/tempbb)-1)
  5     continue

c
c compute beta array, the probability of escape per Thomson time.
c bet evaluated for spherical geometry and nearly uniform sources.
c Between x=0.1 and 1.0, a function flz modifies beta to allow
c the increasingly large energy change per scattering to gradually
c eliminate spatial diffusion
        jnr=INT(log10(0.1/xmin)/delta+1)
        jnr=min(jnr,jmax-1)
        jrel=INT(log10(1/xmin)/delta+1)
        jrel=min(jrel,jmax)
        xnr=x(jnr)
        xr=x(jrel)
        do  501 j=1, jnr-1
          taukn=tautom*rel(j)
          bet(j)= 1/tautom/(1+taukn/3)
  501   continue
        do 600 j=jnr,jrel
          taukn=tautom*rel(j)
          arg=(x(j)-xnr)/(xr-xnr)
          flz=1-arg
          bet(j)=1/tautom/(1+taukn/3*flz)
  600   continue
        do 601 j=jrel+1,jmax
          bet(j)=1/tautom
  601   continue
c
      call thermlc
     >  (tautom,theta,deltal,x,jmax,dphesc,dphdot,bet,c2)
c
c     the spectrum in E F_E
      do 497 j=1,jmax-1
          sptot(j)=dphesc(j)*x(j)**2
          spsca(j)=(dphesc(j) - dphdot(j) * exp(-1. * tautom))
     >     * (x(j) ** 2)
c          write(1,*) x(j), sptot(j)
 497  continue
c     the input spectrum
c      do 498 j=1,jmaxth
c         write(2,*) x(j), dphdot(j)*x(j)**2
c 498  continue
c
      return
      end

      subroutine thermlc
     >  (tautom,theta,deltal,x,jmax,dphesc,dphdot,bet,c2)
c This program computes the effects of Comptonization by
c nonrelativistic thermal electrons in a sphere including escape, and
c relativistic corrections up to photon energies of 1 MeV.
c the dimensionless photon energy is x=hv/(m*c*c)
c
c The input parameters and functions are:
c dphdot(x), the photon production rate
c tautom, the Thomson scattering depth
c theta, the temperature in units of m*c*c
c c2(x), and bet(x), the coefficients in the K-equation and the
c   probability of photon escape per Thomson time, respectively,
c   including Klein-Nishina corrections
c The output parameters and functions are:
c dphesc(x), the escaping photon density
      implicit none
      integer jmax
      real tautom,theta,deltal
      real x(900),dphesc(900), dphdot(900),bet(900),c2(900)
      integer j,jj
      real a(900),b(900),c(900),c20,w1,w2,t1,t2,t3,x32,aa
      real d(900),alp(900),g(900),gam(900),u(900)
c u(x) is the dimensionless photon occupation number

      c20=tautom/deltal
c
c determine u
c define coefficients going into equation
c a(j)*u(j+1)+b(j)*u(j)+c(j)*u(j-1)=d(j)
      do 1 j=2,jmax-1
        w1=sqrt(x(j)*x(j+1))
        w2=sqrt(x(j-1)*x(j))
c  w1 is x(j+1/2)
c  w2 is x(j-1/2)
        a(j)=-c20*c2(j)*(theta/deltal/w1+0.5)
        t1=-c20*c2(j)*(0.5-theta/deltal/w1)
        t2=c20*c2(j-1)*(theta/deltal/w2+0.5)
        t3=x(j)**3*(tautom*bet(j))
        b(j)=t1+t2+t3
        c(j)=c20*c2(j-1)*(0.5-theta/deltal/w2)
        d(j)=x(j)*dphdot(j)
   1  continue

c define constants going into boundary terms
c u(1)=aa*u(2) (zero flux at lowest energy)
c u(jx2) given from region 2 above
      x32=sqrt(x(1)*x(2))
      aa=(theta/deltal/x32+0.5)/(theta/deltal/x32-0.5)
c
c zero flux at the highest energy
      u(jmax)=0
c
c invert tridiagonal matrix
      alp(2)=b(2)+c(2)*aa
      gam(2)=a(2)/alp(2)
      do 2 j=3,jmax-1
        alp(j)=b(j)-c(j)*gam(j-1)
        gam(j)=a(j)/alp(j)
  2   continue
      g(2)=d(2)/alp(2)
      do 3 j=3,jmax-2
        g(j)=(d(j)-c(j)*g(j-1))/alp(j)
  3   continue
      g(jmax-1)=
     > (d(jmax-1)-a(jmax-1)*u(jmax)-c(jmax-1)*g(jmax-2))/alp(jmax-1)
      u(jmax-1)=g(jmax-1)
      do 4 j=3,jmax-1
        jj=jmax+1-j
        u(jj)=g(jj)-gam(jj)*u(jj+1)
  4   continue
      u(1)=aa*u(2)

c compute new value of dph(x) and new value of dphesc(x)
      do 9 j=1,jmax
        dphesc(j)=x(j)*x(j)*u(j)*bet(j)*tautom
  9   continue
      return
      end


c Currently we need only blackbody Comptonization; not from disk blackbody
c#if 0
c      subroutine thdsCompton(tempbb,theta,gamma,x,jmax,sptot)
cc  version: January 96
cc
cc
cc     Thermal Comptonization; solves Kompaneets eq. with some
cc     relativistic corrections. See Lightman & Zdziarski (1987), ApJ
cc     The seed spectrum is DISK blackbody.
c      IMPLICIT NONE
c      real delta,xmin,xmax,deltal,xnr,xr,taukn,arg,flz,pi
c      real sptot(900),x(900),dphesc(900),dphdot(900)
c     >  ,rel(900),bet(900),c2(900)
c      integer j,jmax,jmaxth,jnr,jrel
c      DOUBLE PRECISION w,w1,z1,z2,z3,z4,z5,z6
c      real tautom
cc  input parameters:
c      real tempbb,theta,gamma
c
c      real    ear(0:5000),photar(5000),photer(5000), parth(10)
c      integer ifl,ne
c
cc use internally Thomson optical depth
c      tautom=sqrt(2.25+3/(theta*((gamma+.5)**2-2.25)))-1.5
cc
c      pi=3.14159
cc clear arrays (important for repeated calls)
c      do 31 j=1,900
c        dphesc(j)=0
c        dphdot(j)=0
c        rel(j)=0
c        bet(j)=0
c        c2(j)=0
c        sptot(j)=0
c  31  continue
cC
cC JMAX - # OF PHOTON ENERGIES
cC
cc delta is the 10-log interval of the photon array.
c      delta=0.02
c      deltal=delta*log(10.)
c      xmin=1e-4*tempbb
c      xmax=40*theta
c      jmax=min(899,int(log10(xmax/xmin)/delta)+1)
cC
cC X - ARRAY FOR PHOTON ENERGIES
cC
c      do 4 j=1,jmax+1
c         x(j)=xmin*10.**((j-1)*delta)
c  4   continue
cc
cc compute c2(x), and rel(x) arrays
cc c2(x) is the relativistic correction to Kompaneets equation
cc rel(x) is the Klein-Nishina cross section divided by the Thomson crossection
c      do 500 j=1,jmax
c        w=x(j)
cc c2 is the Cooper's coefficient calculated at w1
cc w1 is x(j+1/2) (x(i) defined up to jmax+1)
c        w1=sqrt(x(j)*x(j+1))
c        c2(j)=SNGL(w1**4/(1+4.6*w1+1.1*w1*w1))
c        if (w.le.0.05) then
cc use asymptotic limit for rel(x) for x less than 0.05
c          rel(j)=SNGL(1-2*w+26*w*w/5)
c        else
c          z1=(1+w)/w**3
c          z2=1+2*w
c          z3=log(z2)
c          z4=2*w*(1+w)/z2
c          z5=z3/2/w
c          z6=(1+3*w)/z2/z2
c          rel(j)=SNGL(0.75*(z1*(z4-z3)+z5-z6))
c        end if
c 500  continue
c
cc the thermal emision spectrum
c        jmaxth=min(900,int(log10(50*tempbb/xmin)/delta))
c        if (jmaxth .gt. jmax) then
c           print *,'thcomp: ',jmaxth,jmax
c           jmaxth = jmax
c        end if
cc        planck=15/(pi*tempbb)**4
cc        do 5 j=1,jmaxth
cc          dphdot(j)=planck*x(j)**2/(exp(x(j)/tempbb)-1)
cc  5     continue
c
c        do j=1,jmaxth-1
c           ear(j-1) = 511*sqrt(x(j)*x(j+1))
c        end do
c        parth(1) = tempbb*511
c        ne = jmaxth-2
c        call XSDSKB(ear,ne,parth,ifl,photar,photer)
c
c        do j=1,ne
c           dphdot(j+1) = 511*photar(j)/(ear(j)-ear(j-1))
c        end do
c        jmaxth = ne+1
c        dphdot(1) = dphdot(2)
c
cc
cc compute beta array, the probability of escape per Thomson time.
cc bet evaluated for spherical geometry and nearly uniform sources.
cc Between x=0.1 and 1.0, a function flz modifies beta to allow
cc the increasingly large energy change per scattering to gradually
cc eliminate spatial diffusion
c        jnr=INT(log10(0.1/xmin)/delta+1)
c        jnr=min(jnr,jmax-1)
c        jrel=INT(log10(1/xmin)/delta+1)
c        jrel=min(jrel,jmax)
c        xnr=x(jnr)
c        xr=x(jrel)
c        do  501 j=1, jnr-1
c          taukn=tautom*rel(j)
c          bet(j)= 1/tautom/(1+taukn/3)
c  501   continue
c        do 600 j=jnr,jrel
c          taukn=tautom*rel(j)
c          arg=(x(j)-xnr)/(xr-xnr)
c          flz=1-arg
c          bet(j)=1/tautom/(1+taukn/3*flz)
c  600   continue
c        do 601 j=jrel+1,jmax
c          bet(j)=1/tautom
c  601   continue
cc
c      call thermlc
c     >  (tautom,theta,deltal,x,jmax,dphesc,dphdot,bet,c2)
cc
cc     the spectrum in E F_E
c      do 497 j=1,jmax-1
c          sptot(j)=dphesc(j)*x(j)**2
cc          write(1,*) x(j), sptot(j)
c 497  continue
c
cc      print *,'jmax: ',jmax,jmaxth
cc      open(33,file='spec.dat')
ccc     the input spectrum
cc      do 498 j=1,min(jmaxth,jmax-1)
cc         write(33,*) 511*x(j), dphdot(j)*x(j), dphesc(j)*x(j)
cc 498  continue
cc      close(33)
c
c      return
c      end
c
c#endif

c ------------------------------------------------------------------ c

      REAL FUNCTION SPP(Y,Xnonth,Nnonth,Spnth)

      IMPLICIT NONE
      INTEGER Nnonth
      REAL Y , Xnonth(Nnonth) , Spnth(Nnonth)
      INTEGER il , ih
      REAL xx
      SAVE ih
      DATA ih/2/

      xx = 1/Y
      IF ( xx.LT.Xnonth(ih) ) ih = 2
      DO WHILE ( ih.LT.Nnonth .AND. xx.GT.Xnonth(ih) )
         ih = ih + 1
      ENDDO

      il = ih - 1
      SPP = Spnth(il) + (Spnth(ih)-Spnth(il))*(xx-Xnonth(il))
     &      /(Xnonth(ih)-Xnonth(il))

      RETURN
      END
