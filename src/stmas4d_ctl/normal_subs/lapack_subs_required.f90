MODULE LAPACK_subs
  IMPLICIT NONE
  CONTAINS



!==============================================================================
!========================    The Following Code is   ==========================
!======================== from LAPACK and is minimal ==========================
!========================  to run the routine SGESV  ==========================
!==============================================================================

      SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
      !EXTERNAL           SGETRF, SGETRS, XERBLA
      INTRINSIC          MAX

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESV ', -INFO )
         RETURN
      END IF

      CALL SGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN

         CALL SGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,INFO )
      END IF
      RETURN

      END SubRoutine SGESV

      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            I, IINFO, J, JB, NB
      !EXTERNAL           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
      INTEGER            ILAENV
      !EXTERNAL           ILAENV
      INTRINSIC          MAX, MIN

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRF', -INFO )
         RETURN
      END IF
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
         CALL SGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
            CALL SGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
            CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
            IF( J+JB.LE.N ) THEN
               CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,IPIV, 1 )
               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),LDA )
               IF( J+JB.LE.M ) THEN
                  CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1,N-J-JB+1, JB, -ONE, &
                              A( J+JB, J ), LDA,A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
      END Subroutine SGETRF


      SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      LOGICAL            NOTRAN
      LOGICAL            LSAME
      !EXTERNAL           LSAME
      !EXTERNAL           SLASWP, STRSM, XERBLA
      INTRINSIC          MAX

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRS', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
      IF( NOTRAN ) THEN
         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
         CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,ONE, A, LDA, B, LDB )
         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,NRHS, ONE, A, LDA, B, LDB )
      ELSE
         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,ONE, A, LDA, B, LDB )
         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,A, LDA, B, LDB )
         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
      RETURN
      END SubRoutine SGETRS


      SUBROUTINE XERBLA( SRNAME, INFO )
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
      INTRINSIC          LEN_TRIM
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
      STOP
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ','an illegal value' )
      END SubRoutine XERBLA


      SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      REAL ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      REAL A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      !EXTERNAL LSAME
      !EXTERNAL XERBLA
      INTRINSIC MAX
      REAL TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')).AND.(.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')).AND.(.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('SGEMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0).OR.(N.EQ.0).OR.(((ALPHA.EQ.ZERO).OR.(K.EQ.0)).AND.(BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (NOTB) THEN
          IF (NOTA) THEN
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      IF (B(L,J).NE.ZERO) THEN
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      END IF
   80             CONTINUE
   90         CONTINUE
          ELSE
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*B(J,L)
                          DO 150 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  150                     CONTINUE
                      END IF
  160             CONTINUE
  170         CONTINUE
          ELSE
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
      RETURN
      END SubRoutine SGEMM



      SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      REAL               SFMIN
      INTEGER            I, J, JP
      REAL               SLAMCH
      INTEGER            ISAMAX
      !EXTERNAL           SLAMCH, ISAMAX
      !EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
      INTRINSIC          MAX, MIN
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETF2', -INFO )
         RETURN
      END IF
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
      SFMIN = SLAMCH('S')
      DO 10 J = 1, MIN( M, N )
         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
            IF( JP.NE.J ) CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
            IF( J.LT.M ) THEN 
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN 
                  CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 ) 
               ELSE 
                 DO 20 I = 1, M-J 
                    A( J+I, J ) = A( J+I, J ) / A( J, J ) 
   20            CONTINUE 
               END IF 
            END IF 
         ELSE IF( INFO.EQ.0 ) THEN
            INFO = J
         END IF
         IF( J.LT.MIN( M, N ) ) THEN
            CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
      END SubRoutine SGETF2


      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
      INTEGER            INCX, K1, K2, LDA, N
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      REAL               TEMP
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
      RETURN
      END SubRoutine SLASWP



      SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      REAL ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      REAL A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      !EXTERNAL LSAME
      !EXTERNAL XERBLA
      INTRINSIC MAX
      REAL TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')).AND.(.NOT.LSAME(TRANSA,'T')).AND.(.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('STRSM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END SubRoutine STRSM


      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
      INTEGER            IEEECK, IPARMQ
      !EXTERNAL           IEEECK, IPARMQ
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
      ILAENV = -1
      RETURN
   10 CONTINUE
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
         IF( ( IC.GE.129 .AND. IC.LE.137 ).OR.(IC.GE.145 .AND. IC.LE.153 ).OR.( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ).OR.( IC.GE.145 .AND. IC.LE.153 ) &
                  .OR.( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
      GO TO ( 50, 60, 70 )ISPEC
   50 CONTINUE
      NB = 1
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF'.OR.C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
   60 CONTINUE
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
   70 CONTINUE
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
   80 CONTINUE
      ILAENV = 6
      RETURN
   90 CONTINUE
      ILAENV = 2
      RETURN
  100 CONTINUE
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
  110 CONTINUE
      ILAENV = 1
      RETURN
  120 CONTINUE
      ILAENV = 50
      RETURN
  130 CONTINUE
      ILAENV = 25
      RETURN
  140 CONTINUE
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
  150 CONTINUE
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
  160 CONTINUE
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
      END Function ILAENV


      LOGICAL FUNCTION LSAME(CA,CB)
      CHARACTER CA,CB

      INTRINSIC ICHAR
      INTEGER INTA,INTB,ZCODE
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
      ZCODE = ICHAR('Z')
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.&
             INTA.GE.145 .AND. INTA.LE.153 .OR.&
             INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.&
             INTB.GE.145 .AND. INTB.LE.153 .OR.&
             INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB


      END Function LSAME


      REAL             FUNCTION SLAMCH( CMACH )
      CHARACTER          CMACH
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,RND, SFMIN, SMALL, T
      LOGICAL            LSAME
      !EXTERNAL           LSAME
      !EXTERNAL           SLAMC2
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,EMAX, RMAX, PREC
      DATA               FIRST / .TRUE. /
      IF( FIRST ) THEN
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
      SLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
      END Function SLAMCH



      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
      REAL               SLAMC3
      !EXTERNAL           SLAMC3
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
      DATA               FIRST / .TRUE. /
      IF( FIRST ) THEN
         ONE = 1
         A = 1
         C = 1
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
         B = 1
         C = SLAMC3( A, B )
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) ) LRND = .FALSE.
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
         LT = 0
         A = 1
         C = 1
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
      END IF
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
      END SubRoutine SLAMC1


      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,SIXTH, SMALL, THIRD, TWO, ZERO
      REAL               SLAMC3
      !EXTERNAL           SLAMC3
      !EXTERNAL           SLAMC1, SLAMC4, SLAMC5
      INTRINSIC          ABS, MAX, MIN
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,LRMIN, LT
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
      IF( FIRST ) THEN
         ZERO = 0
         ONE = 1
         TWO = 2
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
         B = LBETA
         A = B**( -LT )
         LEPS = A
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS ) B = LEPS
         LEPS = 1
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
         IF( A.LT.LEPS ) LEPS = A
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
               IWARN = .TRUE.
            END IF
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
               IWARN = .TRUE.
            END IF
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
               IWARN = .TRUE.
            END IF
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
            IWARN = .TRUE.
         END IF
         FIRST = .FALSE.
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
         IEEE = IEEE .OR. LIEEE1
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
      RETURN
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
           '  EMIN = ', I8, /&
           ' If, after inspection, the value EMIN looks',&
           ' acceptable please comment out ',&
           / ' the IF block as marked within the code of routine',&
           ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
      END SubRoutine SLAMC2


      REAL             FUNCTION SLAMC3( A, B )
      REAL               A, B
      SLAMC3 = A + B
      RETURN
      END Function SLAMC3



      SUBROUTINE SLAMC4( EMIN, START, BASE )
      INTEGER            BASE
      INTEGER            EMIN
      REAL               START
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
      REAL               SLAMC3
      !EXTERNAL           SLAMC3
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ).AND.( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
      RETURN
      END SubRoutine SLAMC4

      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
      REAL               SLAMC3
      !EXTERNAL           SLAMC3
      INTRINSIC          MOD
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
         EMAX = EMAX - 1
      END IF
      IF( IEEE ) THEN
         EMAX = EMAX - 1
      END IF
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE ) OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE ) Y = OLDY
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
      RMAX = Y
      RETURN
      END SubRoutine SLAMC5


      INTEGER FUNCTION ISAMAX(N,SX,INCX)
      INTEGER INCX,N
      REAL SX(*)
      REAL SMAX
      INTEGER I,IX
      INTRINSIC ABS
      ISAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ISAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (ABS(SX(IX)).LE.SMAX) GO TO 5
          ISAMAX = I
          SMAX = ABS(SX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
          IF (ABS(SX(I)).LE.SMAX) GO TO 30
          ISAMAX = I
          SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END Function ISAMAX


      SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      REAL ALPHA
      INTEGER INCX,INCY,LDA,M,N
      REAL A(LDA,*),X(*),Y(*)
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      REAL TEMP
      INTEGER I,INFO,IX,J,JY,KX
      !EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('SGER  ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
      RETURN
      END SubRoutine SGER


      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      REAL STEMP
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          STEMP = SX(IX)
          SX(IX) = SY(IY)
          SY(IY) = STEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
          STEMP = SX(I+1)
          SX(I+1) = SY(I+1)
          SY(I+1) = STEMP
          STEMP = SX(I+2)
          SX(I+2) = SY(I+2)
          SY(I+2) = STEMP
   50 CONTINUE
      RETURN
      END Subroutine SSWAP


      SUBROUTINE SSCAL(N,SA,SX,INCX)
      REAL SA
      INTEGER INCX,N
      REAL SX(*)
      INTEGER I,M,MP1,NINCX
      INTRINSIC MOD
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END SubRoutine SSCAL



      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
      INTEGER            NH, NS
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ).OR.( ISPEC.EQ.IACC22 ) ) THEN
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) NS = 4
         IF( NH.GE.60 ) NS = 10
         IF( NH.GE.150 ) NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 ) NS = 64
         IF( NH.GE.3000 ) NS = 128
         IF( NH.GE.6000 ) NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
      IF( ISPEC.EQ.INMIN ) THEN
         IPARMQ = NMIN
      ELSE IF( ISPEC.EQ.INIBL ) THEN
         IPARMQ = NIBBLE
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
         IPARMQ = NS
      ELSE IF( ISPEC.EQ.INWIN ) THEN
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
         IPARMQ = 0
         IF( NS.GE.KACMIN ) IPARMQ = 1
         IF( NS.GE.K22MIN ) IPARMQ = 2
      ELSE
         IPARMQ = -1
      END IF
      END Function IPARMQ



      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
      INTEGER            ISPEC
      REAL               ONE, ZERO
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,NEGZRO, NEWZRO, POSINF
      IEEECK = 1
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( ISPEC.EQ.0 ) RETURN
      NAN1 = POSINF + NEGINF
      NAN2 = POSINF / NEGINF
      NAN3 = POSINF / POSINF
      NAN4 = POSINF*ZERO
      NAN5 = NEGINF*NEGZRO
      NAN6 = NAN5*0.0
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
      RETURN
      END  FUNCTION IEEECK

END MODULE LAPACK_subs
