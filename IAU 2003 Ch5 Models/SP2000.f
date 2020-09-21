      DOUBLE PRECISION FUNCTION SP2000 ( DATE1, DATE2 )
*+
*  - - - - - - -
*   S P 2 0 0 0
*  - - - - - - -
*
*  The quantity s', positioning the Terrestrial Ephemeris Origin on the
*  equator of the Celestial Intermediate Pole.
*
*  Annex to IERS Conventions 2000, Chapter 5
*
*  Given:
*     DATE1,DATE2    d      TT date (JD = DATE1+DATE2)
*
*  Returned:
*     SP2000         d      the quantity s' in radians
*
*  This revision:  2002 November 25
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Time since J2000, in Julian centuries
      DOUBLE PRECISION T

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  Approximate S'.
      SP2000 = -47D-6 * T * DAS2R

      END
