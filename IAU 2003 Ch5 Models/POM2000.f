      SUBROUTINE POM2000 ( XP, YP, SP, RPOM )
*+
*  - - - - - - - -
*   P O M 2 0 0 0
*  - - - - - - - -
*
*  Form the matrix of polar motion, IAU 2000.
*
*  Annex to IERS Conventions 2000, Chapter 5
*
*  Given:
*     XP,YP      d      coordinates of the pole (radians)
*     SP         d      the quantity s' (radians)
*
*  Returned:
*     RPOM     d(3,3)   polar-motion matrix
*
*  The returned rotation matrix is the first to be applied when
*  transforming a TRS vector into a CRS vector.
*
*  Calls the SOFA routines iau_IR, iau_RX, iau_RY, iau_RZ.
*
*  This revision:  2002 November 25
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION XP, YP, SP, RPOM(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Construct the matrix.
      CALL iau_IR ( RPOM )
      CALL iau_RX ( YP, RPOM )
      CALL iau_RY ( XP, RPOM )
      CALL iau_RZ ( -SP, RPOM )

      END
