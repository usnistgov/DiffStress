CFEE
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE VOIGT   ---->   VERSION OF 09/02/98
C
C     TRANSFORMS 6X1 MATRIX T1 INTO SECOND ORDER TENSOR T2 IF IOPT=1
C     AND VICEVERSA IF IOPT=2.
C     TRANSFORMS 6X6 MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=3
C     AND VICEVERSA IF IOPT=4.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE VOIGT(T1,T2,C2,C4,IOPT)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
      DIMENSION IJV(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

      IF(IOPT.EQ.1) THEN
      DO 30 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      T2(I1,I2)=T1(I)
   30 T2(I2,I1)=T1(I)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 40 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
   40 T1(I)=T2(I1,I2)
      ENDIF
C
      IF (IOPT.EQ.3) THEN
      DO 10 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 10 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
   10 C4(I2,I1,J2,J1)=C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.4) THEN
      DO 20 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 20 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
   20 C2(I,J)=C4(I1,I2,J1,J2)
      ENDIF
C
      RETURN
      END
