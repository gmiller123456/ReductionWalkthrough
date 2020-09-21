      SUBROUTINE CBPN2000 ( DATE1, DATE2, DPSI, DEPS, RBPNC )
*+
*  - - - - - - - - -
*   C B P N 2 0 0 0
*  - - - - - - - - -
*
*  Classical bias-precession-nutation matrix.
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     DATE1,DATE2   d      TT date (JD = DATE1+DATE2)
*     DPSI,DEPS     d      nutation (luni-solar + planetary, radians)
*
*  Returned:
*     RBPNC       d(3,3)   true-to-celestial matrix
*
*  Calls SOFA routines iau_IR, iau_RX, iau_RY, iau_RZ
*
*  This revision:  2002 November 26
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS, RBPNC(3,3)

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  J2000 obliquity (Lieske et al. 1977)
      DOUBLE PRECISION EPS0
      PARAMETER ( EPS0 = 84381.448D0 * DAS2R )

*  The ICRS RA of the J2000 equinox (Chapront et al., 2002)
      DOUBLE PRECISION DRA0
      PARAMETER ( DRA0 = -0.0146D0 * DAS2R )

      DOUBLE PRECISION T

*  The precession and obliquity corrections (radians per century)
      DOUBLE PRECISION PRECOR, OBLCOR
      PARAMETER ( PRECOR = -0.29965D0 * DAS2R,
     :            OBLCOR = -0.02524D0 * DAS2R )

*  The frame bias corrections in longitude and obliquity
      DOUBLE PRECISION DPSIBI, DEPSBI
      PARAMETER ( DPSIBI = -0.041775D0 * DAS2R,
     :            DEPSBI = -0.0068192D0 * DAS2R )

      DOUBLE PRECISION DPSIPR, DEPSPR, EPSA80, PSIA77, OMA77,
     :                 CHIA, PSIA, OMA, EPSA

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  Precession rate contributions with respect to IAU 1976/80.
      DPSIPR = PRECOR * T
      DEPSPR = OBLCOR * T

*  IAU 1980 mean obliquity of date.
      EPSA80 = EPS0 + (  -46.8150D0 +
     :                (   -0.00059D0 +
     :                (    0.001813D0 ) * T ) * T ) * T * DAS2R

*  Precession angles (Lieske et al. 1977)
      PSIA77 =        ( 5038.7784D0 +
     :                (   -1.07259D0 +
     :                (   -0.001147D0 ) * T ) * T ) * T * DAS2R
      OMA77  = EPS0 + (
     :                (    0.05127D0 +
     :                (   -0.007726D0 ) * T ) * T ) * T * DAS2R
      CHIA   =        (   10.5526D0 +
     :                (   -2.38064D0 +
     :                (   -0.001125D0 ) * T ) * T ) * T * DAS2R

*  Apply IAU 2000A precession corrections.
      PSIA = PSIA77 + DPSIPR
      OMA  = OMA77  + DEPSPR
      EPSA = EPSA80 + DEPSPR

*  Initialize the true-to-celestial matrix.
      CALL iau_IR ( RBPNC )

*  Remove IAU 2000A nutation (pure: luni-solar and planetary).
      CALL iau_RX ( EPSA+DEPS, RBPNC )
      CALL iau_RZ ( DPSI, RBPNC )
      CALL iau_RX ( -EPSA, RBPNC )

*  Remove precession (Lieske et al. 1977 plus corrections).
      CALL iau_RZ ( -CHIA, RBPNC )
      CALL iau_RX ( OMA, RBPNC )
      CALL iau_RZ ( PSIA, RBPNC )
      CALL iau_RX ( -EPS0, RBPNC )

*  Remove frame bias.
      CALL iau_RX ( DEPSBI, RBPNC )
      CALL iau_RY ( -DPSIBI*SIN(EPS0), RBPNC )
      CALL iau_RZ ( -DRA0, RBPNC )

      END
