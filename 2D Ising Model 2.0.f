C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c-------Posar Implicit None
c-------Tenir en compte quina precisió volem
c-------Debug:  gfortran -g -fbounds-check -Wall -fbacktrace

C	------------Pràctiques amb l'ordinador------------

C	PROGRAMA FET PER: 	DAVID FERNANDEZ BONET
C	DATA:				30/04/2020
C	CONTACTE:			davferdz@gmail.com

C     --------------------------------------------
	Program P1ex3
	IMPLICIT NONE

c	Definim variables i paràmetres importants
	INTEGER*4 SEED,NSEED,SEED0,I,J,L, II
	PARAMETER (L=30)
	
c--
	INTEGER*2 S(1:L,1:L)
	INTEGER*4 PBC(0:L+1)
	DOUBLE PRECISION genrand_real2, MAG, MAGNE, ENE, ENERG, TEMP(200)

	
c--
	INTEGER DE, IMC, IPAS, suma, N
	DOUBLE PRECISION DELTA, x,y, w(-8:8), time1, time2
	CHARACTER*3 NOM
	CHARACTER*30 date

	INTEGER MCINI,MCD, MCTOT
	PARAMETER (MCTOT=40000)
	PARAMETER (MCINI=2000)
	PARAMETER (MCD=20)
	DOUBLE PRECISION ssum,sume,sume2,summ,summ2,sumam,vare,varm
	DOUBLE PRECISION erre, errm, sus, cal





	N= L*L
	NOM= "A120"
      nseed=200
      seed0=23456




	CALL cpu_time(time1)

      OPEN(unit=12,FILE=nom//'.dat')
      write(12,*) '# l,temp,ssum,sume,sume2,vare,summ,sumam,summ2',
     &             'varm,erre,errv,cal,sus'

	PBC(0)=L
	PBC(L+1)=1
	DO I=1,L
		PBC(I)=I
	END DO

C	Definim les 200 temperatures en un vector

      TEMP(1)=2.2550
      
      DO II=2,200
            TEMP(II)=TEMP(1)+(II-1)*0.005d0
      END DO
      
C	INICI LOOP TEMPERATURES 

C	Afegim el loop més exterior que correrà per les 200 temperatures
C	De manera que TEMP(II) ens doni la temperatura desitjada
	DO II=1,200

C	Millorem el rendiment del programa amb l funció w
 		DO DE=-8,8
            	w(de)=dexp(-dfloat(de)/temp(II))
      	END DO

C	Posem els contadors a 0
      ssum=0.d0
      sume=0.d0
      sume2=0.d0
      summ=0.d0
      summ2=0.d0
      sumam=0.d0

C	Seed loop: computem diferents llavors
C	per generar valors aleatoris
      DO seed=seed0,seed0+nseed-1,1
c            write(*,*) 'llavor= ',seed

	CALL init_genrand(SEED)
C	Computa la "lattice" inicial amb valor d'spin 1 o -1

	DO I=1,L
		DO J=1,L
			IF (genrand_real2().lt.0.5d0) THEN
				S(I,J)=1
			ELSE
				S(I,J)=-1
			ENDIF
		END DO
	END DO

	IMC=0
	MAG=MAGNE(S,L)
	ENE=ENERG(S,L,PBC)



C	Comença el bucle per tots els passos MC. El codi genera 2 valors
C	aleatoris per triar un punt de la lattice, i canvia l'spin
C	en funció de l'interacció amb els seus veïns.

 	DO IMC=1,MCTOT
            DO IPAS=1,N
 			x=genrand_real2()
                  I=INT(L*x)+1
                  y=genrand_real2()
                  J=INT(L*y)+1

			suma=S(I,PBC(J+1))+S(I,PBC(J-1))+
     $		S(PBC(I+1),J)+S(PBC(I-1),J)
			DE = 2*S(I,J)*suma

			IF (DE .LE. 0) THEN
				S(I,J) = -S(I,J)
				ENE= ENE+DE
			END IF
			IF (DE .GT. 0) THEN
				DELTA=genrand_real2()
				IF (DELTA .LT. W(DE)) THEN
					S(I,J) = -S(I,J)
					ENE= ENE+DE
				END IF
			END IF

		END DO

            IF ((imc.GT.mcini).AND.(mcd*(imc/mcd).EQ.imc)) THEN
                  mag=magne(s,l)

                  ssum=ssum+1.d0
                  sume=sume+ene
                  sume2=sume2+ene*ene

                  summ=summ+mag
                  sumam=sumam+dabs(mag)
                  summ2=summ2+mag*mag
            ENDIF	
C	Acabem el loop de MC
	END DO

C	Acabem el loop de les seed
	END DO

C	Resultats: normalitzem promitgjos i calculem variàncies
C	--

C	Promig energia
      sume=sume/ssum
      sume2=sume2/ssum
C	Promig magnetització/ABS(magnetització)
      summ=summ/ssum
      sumam=sumam/ssum
      summ2=summ2/ssum
C	Variància energia i magnetització
      vare=sume2-sume*sume
      varm=summ2-sumam*sumam 



	erre=dsqrt(vare)/(n*dsqrt(ssum))
      errm=dsqrt(varm)/(n*dsqrt(ssum))
C	Capacitat calorífica  
      cal=(vare)/(n*temp(ii)*temp(ii))
C	Susceptibilitat
      sus=(varm)/(n*temp(ii))
      write(12,*)l,temp(ii),ssum,sume/n,sume2/(n*n),vare,summ/n,
     &       sumam/n,dsqrt(summ2/(n*n)),varm,erre,errm,cal,sus
      write(12,*)
      write(12,*)
C	(Separació entre temperatures pel fitxer de dades)


      call cpu_time(time2)
      call fdate(date)
      write(*,*) date
      write(*,*) 'cputime= ',time2-time1

	END DO
C	Acabem el loop de les 200 temperatures

	


	END





C	----------SUBROUTINES I FUNCIONS------------


C	FUNCIÓ ENERGIA

	REAL*8 FUNCTION ENERG(S,L,PBC)
	INTEGER*2 S(1:L,1:L)
	INTEGER*4 I,J,L
	INTEGER*4 PBC(0:L+1)
	REAL*8 ENE
	ENE=0.0d0

	DO I=1,L
		DO J=1,L
			ENE= ENE-S(I,J)*S(PBC(I+1),J) -S(I,J)*S(I, PBC(J+1))
		END DO
	END DO
	ENERG=ENE
	RETURN
	END


C	FUNCIÓ MAGNE
      real*8 function magne(s,l)
      integer*2 s(1:L,1:L)
	INTEGER*4 i,j,l
      real*8 mag
      mag=0.d0
      do i=1,l
            do j=1,l
                  mag=mag+s(i,j)
            enddo
      enddo
      magne=mag
      return
      end

C	SUBRUTINA WRITECONFIG

	SUBROUTINE WRITECONFIG(S,L)
	INTEGER I,J,L
	INTEGER*2 S(1:L,1:L)
	

	
	OPEN (14, FILE="P1-configuration.dat")
	DO I=1,L
		DO J=1,L
			IF (S(I,J) .EQ. 1) THEN
				WRITE(14,*) I,J
			END IF
			
		END DO
	END DO

	CLOSE (14)
	END
	
	


c     initialize mt(0:N-1) with a seed
c-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end
c-----------------------------------------------------------------------
c     initialize by an array with array-length
c     init_key is the array for initializing keys
c     key_length is its length
c-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
c
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)
     &           +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
c
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0xffffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0x7fffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1]-real-interval
c-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1) with 53-bit resolution
c-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
c-----------------------------------------------------------------------
c     initialize large number (over 32-bit constant number)
c-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
