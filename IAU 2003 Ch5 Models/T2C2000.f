      SUBROUTINE T2C2000 ( RPOM, THETA, RBPN, RT2C )
*+
*  - - - - - - - -
*   T 2 C 2 0 0 0
*  - - - - - - - -
*
*  Form the TRS-to-CRS matrix, IAU 2000, from components.
*
*  Annex to IERS Conventions 2000, Chapter 5
*
*  Given:
*     RPOM     d(3,3)   polar motion matrix (W)
*
*  followed by either (for the CEO-based transformation):
*     THETA      d      Earth Rotation Angle (radians, giving matrix R)
*     RBPN     d(3,3)   intermediate-to-celestial matrix (Q)
*
*  or alternatively (for the classical, equinox-based, transformation):
*     THETA      d      Greenwich Sidereal Time (radians)
*     RBPN     d(3,3)   true-to-celestial matrix
*
*  Returned:
*     RT2C     d(3,3)   terrestrial-to-celestial matrix
*
*  Calls the SOFA routines iau_CR, iau_RZ, iau_RXR.
*
*  This revision:  2002 November 25
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RPOM(3,3), THETA, RBPN(3,3), RT2C(3,3)

      DOUBLE PRECISION R(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Polar motion.
      CALL iau_CR ( RPOM, R )

*  Earth rotation.
      CALL iau_RZ ( -THETA, R )

*  CIP motion.
      CALL iau_RXR ( RBPN, R, RT2C )

      END
