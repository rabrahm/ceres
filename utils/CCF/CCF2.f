      DOUBLE PRECISION FUNCTION CCF(M_L, M_H, WAV, SPEC,
     +     WEIGHT, SN, V_R, SNW, N, M)
c     This function calculates the cross-correlation function
c     of a spectrum (SPE, at wavlenths WAV) with a mask. M_L
c     and M_H contain the beginning and end of regions where the mask
c     equals 1, ordered.

      IMPLICIT NONE

c     Speed of light, km/s
      DOUBLE PRECISION C
      PARAMETER (C=2.99792458D5)

      INTEGER N, M, COND
      DOUBLE PRECISION M_L(N), M_H(N), WEIGHT(N)
      DOUBLE PRECISION WAV(M), SPEC(M), SN(M)
      DOUBLE PRECISION V_R, GAMMA, SNW
      DOUBLE PRECISION FRACTION, PIX_INIT, PIX_END
      DOUBLE PRECISION M_LLOC(N), M_HLOC(N)

      INTEGER I, J

C     DOPPLER FACTOR
      GAMMA = DSQRT(1. + V_R / C) / DSQRT(1. - V_R / C)

C     DOPPLER SHIFT MASK
      DO I=1,N
         M_LLOC(I) = M_L(I) * GAMMA
         M_HLOC(I) = M_H(I) * GAMMA
      END DO
         

C     I marks where we are in terms of masks
      I = 1
      CCF = 0.0D0
      
      SNW = 0.0D0
      COND = 0
      DO J=2,(M-1)
         PIX_INIT = 0.5*(WAV(J-1) + WAV(J))
         PIX_END  = 0.5*(WAV(J) + WAV(J+1))
         DO WHILE ((M_HLOC(I) < PIX_INIT) .and. (COND .eq. 0))
            IF (I .eq. N) THEN
               COND = 1
            END IF
            IF (COND .eq. 0) THEN
               I = I + 1
            END IF
         END DO
         IF ((PIX_END < M_HLOC(I)) .AND. 
     +      (PIX_INIT > M_LLOC(I))) THEN
            CCF = CCF + SPEC(J) * WEIGHT(I) * SN(J)
            SNW = SNW + SN(J)*WEIGHT(I)
         ELSE IF ((PIX_END < M_HLOC(I)) .AND.
     +           (PIX_INIT < M_LLOC(I)) .AND.
     +           (PIX_END > M_LLOC(I))) THEN
            FRACTION = (PIX_END - M_LLOC(I)) / (PIX_END - PIX_INIT)
            CCF = CCF + SPEC(J) * WEIGHT(I) * FRACTION * SN(J)
            SNW = SNW + FRACTION*SN(J)*WEIGHT(I)
          ELSE IF ((PIX_END > M_HLOC(I)) .AND.
     +           (PIX_INIT > M_LLOC(I)) .AND.
     +           (PIX_INIT < M_HLOC(I))) THEN
            FRACTION = (M_HLOC(I) - PIX_INIT) / (PIX_END - PIX_INIT)
            CCF = CCF + SPEC(J) * WEIGHT(I) * FRACTION * SN(J)
            SNW = SNW + FRACTION*SN(J)*WEIGHT(I)
         ELSE IF ((PIX_END > M_HLOC(I)) .AND.
     +           (PIX_INIT < M_LLOC(I))) THEN
            FRACTION = (M_HLOC(I) - M_LLOC(I)) / (PIX_END - PIX_INIT)
            CCF = CCF + SPEC(J) * WEIGHT(I) * FRACTION * SN(J)
            SNW = SNW + FRACTION*SN(J)*WEIGHT(I)
         END IF
      END DO

C      CCF = CCF / SNW

      END

      DOUBLE PRECISION FUNCTION CCFPIX(M_L, M_H, X, THAR,
     +     DELTA, N, M)
c     This function calculates the cross-correlation function

      IMPLICIT NONE

      INTEGER N, M, COND
      DOUBLE PRECISION M_L(N), M_H(N)
      DOUBLE PRECISION X(M), THAR(M)
      DOUBLE PRECISION DELTA
      DOUBLE PRECISION FRACTION, PIX_INIT, PIX_END

      DOUBLE PRECISION M_LLOC(N), M_HLOC(N)

      INTEGER I, J


C     SHIFT MASK BY DELTA
      DO I=1,N
         M_LLOC(I) = M_L(I) + DELTA
         M_HLOC(I) = M_H(I) + DELTA
      END DO
         
C     I marks where we are in terms of masks
      I = 1
      CCFPIX = 0.0D0

      

      COND = 0
      DO J=2,(M-1)
         PIX_INIT = 0.5*(X(J-1) + X(J))
         PIX_END  = 0.5*(X(J) + X(J+1))
         DO WHILE ((M_HLOC(I) < PIX_INIT) .and. (COND .eq. 0))
            IF (I .eq. N) THEN
               COND = 1
            END IF
            IF (COND .eq. 0) THEN
               I = I + 1
            END IF
         END DO
         IF ((PIX_END < M_HLOC(I)) .AND. 
     +      (PIX_INIT > M_LLOC(I))) THEN
            CCFPIX = CCFPIX + THAR(J) 
         ELSE IF ((PIX_END < M_HLOC(I)) .AND.
     +           (PIX_INIT < M_LLOC(I)) .AND.
     +           (PIX_END > M_LLOC(I))) THEN
            FRACTION = (PIX_END - M_LLOC(I)) / (PIX_END - PIX_INIT)
            CCFPIX = CCFPIX + THAR(J) * FRACTION
         ELSE IF ((PIX_END > M_HLOC(I)) .AND.
     +           (PIX_INIT > M_LLOC(I)) .AND.
     +           (PIX_INIT < M_HLOC(I))) THEN
            FRACTION = (M_HLOC(I) - PIX_INIT) / (PIX_END - PIX_INIT)
            CCFPIX = CCFPIX + THAR(J) * FRACTION
         ELSE IF ((PIX_END > M_HLOC(I)) .AND.
     +           (PIX_INIT < M_LLOC(I))) THEN
            FRACTION = (M_HLOC(I) - M_LLOC(I)) / (PIX_END - PIX_INIT)
            CCFPIX = CCFPIX + THAR(J) * FRACTION
         END IF
      END DO

      END


      DOUBLE PRECISION FUNCTION CCFCOS(M_L, M_H, WAV, SPEC,
     +     WEIGHT, SN, V_R, SNW, N, M)
c     This function calculates the cross-correlation function
c     of a spectrum (SPE, at wavelengths WAV) with a mask. M_L
c     and M_H contain the beginning and end of regions where the mask
c     equals 1, ordered.
c     This version uses a sine function instead of a box --> better lobes
c     The function used is sin((lambda-lambda_init)*pi/(2*half_width))

      IMPLICIT NONE

c     Speed of light, km/s
      DOUBLE PRECISION C
      PARAMETER (C=2.99792458D5)

      INTEGER N, M, COND
      DOUBLE PRECISION M_L(N), M_H(N), WEIGHT(N)
      DOUBLE PRECISION WAV(M), SPEC(M), SN(M)
      DOUBLE PRECISION V_R, GAMMA, SNW
      DOUBLE PRECISION PIX_INIT, PIX_END
      DOUBLE PRECISION M_LLOC(N), M_HLOC(N), HW(N)
      DOUBLE PRECISION WD, NORM, ARG

      INTEGER I, J

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C     DOPPLER FACTOR
      GAMMA = DSQRT(1. + V_R / C) / DSQRT(1. - V_R / C)

C     DOPPLER SHIFT MASK
      DO I=1,N
         M_LLOC(I) = M_L(I) * GAMMA
         M_HLOC(I) = M_H(I) * GAMMA
         HW(I)     = 0.5D0*( M_HLOC(I) - M_LLOC(I) )
      END DO
         

C     I marks where we are in terms of masks
      I = 1
      CCFCOS = 0.0D0     
      SNW = 0.0D0
      COND = 0
      DO J=2,(M-1)
         PIX_INIT = 0.5*(WAV(J-1) + WAV(J))
         PIX_END  = 0.5*(WAV(J) + WAV(J+1))
         DO WHILE ((M_HLOC(I) < PIX_INIT) .and. (COND .eq. 0))
            IF (I .eq. N) THEN
               COND = 1
            END IF
            IF (COND .eq. 0) THEN
               I = I + 1
            END IF
         END DO
         NORM = HW(I) 
         ARG  = PI / (2.0D0*HW(I))
         IF ((PIX_END < M_HLOC(I)) .AND. 
     +      (PIX_INIT > M_LLOC(I))) THEN
            WD = NORM*( DCOS((PIX_INIT-M_LLOC(I))*ARG) 
     +           - DCOS((PIX_END-M_LLOC(I))*ARG) )
            WD = WD / (PIX_END-PIX_INIT)
            CCFCOS = CCFCOS + WD*SPEC(J)
            SNW = SNW + WD
         ELSE IF ((PIX_END < M_HLOC(I)) .AND.
     +           (PIX_INIT < M_LLOC(I)) .AND.
     +           (PIX_END > M_LLOC(I))) THEN
            WD = NORM*( 1.0D0 
     +           - DCOS((PIX_END-M_LLOC(I))*ARG) )
            WD = WD / (PIX_END-PIX_INIT)
            CCFCOS = CCFCOS + SPEC(J) * WD
            SNW = SNW + WD
          ELSE IF ((PIX_END > M_HLOC(I)) .AND.
     +           (PIX_INIT > M_LLOC(I)) .AND.
     +           (PIX_INIT < M_HLOC(I))) THEN
            WD = NORM*( DCOS((PIX_INIT-M_LLOC(I))*ARG)  
     +           + 1.0D0)
            WD = WD / (PIX_END-PIX_INIT)
            CCFCOS = CCFCOS + SPEC(J) * WD
            SNW = SNW + WD
         ELSE IF ((PIX_END > M_HLOC(I)) .AND.
     +           (PIX_INIT < M_LLOC(I))) THEN
            WD = 2 * NORM
            WD = WD / (PIX_END-PIX_INIT)
            CCFCOS = CCFCOS + SPEC(J) * WD
            SNW = SNW + WD
         END IF
      END DO

C      CCF = CCF / SNW

      END

