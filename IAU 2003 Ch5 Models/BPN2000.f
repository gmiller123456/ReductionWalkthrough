      SUBROUTINE BPN2000 ( X, Y, S, RBPN )
*+
*  - - - - - - - -
*   B P N 2 0 0 0
*  - - - - - - - -
*
*  CEO-based bias-precession-nutation matrix.
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     X,Y           d      CIP coordinates
*     S             d      the quantity s (radians)
*
*  Returned:
*     RBPN        d(3,3)   intermediate-to-celestial matrix ("Q")
*
*  Calls SOFA routines iau_IR, iau_RZ, iau_RXR
*
*  This revision:  2002 November 26
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION X, Y, S, RBPN(3,3)

      DOUBLE PRECISION X2, Y2, R2, R, Z, A, AXY, RR(3,3), RL(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Prepare to evaluate expression (10).
      X2 = X*X
      Y2 = Y*Y
      R2 = X2 + Y2
      R = SQRT ( R2 )
      Z = SQRT ( 1D0 - R2 )
      A = 1D0 / ( 1D0 + Z )
      AXY = A*X*Y

*  Right-hand matrix.
      CALL iau_IR ( RR )
      CALL iau_RZ ( S, RR )

*  Left-hand matrix.
      RL(1,1) = 1D0-A*X2
      RL(1,2) = -AXY
      RL(1,3) = X
      RL(2,1) = -AXY
      RL(2,2) = 1D0-A*Y2
      RL(2,3) = Y
      RL(3,1) = -X
      RL(3,2) = -Y
      RL(3,3) = 1D0-A*R2

*  The result is the product of the two matrices.
      CALL iau_RXR ( RL, RR, RBPN )

      END
