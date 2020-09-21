      DOUBLE PRECISION FUNCTION GST2000 ( UTA, UTB, TTA, TTB, DPSI )
*+
*  - - - - - - - -
*   G S T 2 0 0 0
*  - - - - - - - -
*
*  Greenwich Sidereal Time (model consistent with IAU 2000 resolutions).
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     UTA, UTB     d      UT1 date (JD = UTA+UTB)
*     TTA, TTB     d      TT date (JD = TTA+TTB)
*     DPSI         d      nutation in longitude (radians)
*
*  The result is the Greenwich Mean (apparent) Sidereal Time (radians),
*  in the range 0 to 2pi.
*
*  Calls SOFA routine iau_ANP and IERS routines GMST2000 and EE2000.
*
*  This revision:  2002 December 9
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION UTA, UTB, TTA, TTB, DPSI

      DOUBLE PRECISION iau_ANP, GMST2000, EE2000

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Greenwich Sidereal Time, IAU 2000.
      GST2000 = iau_ANP ( GMST2000 ( UTA, UTB, TTA, TTB ) +
     :                    EE2000 ( TTA, TTB, DPSI ) )

      END
