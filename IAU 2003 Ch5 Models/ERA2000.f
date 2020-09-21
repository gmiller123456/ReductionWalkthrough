      DOUBLE PRECISION FUNCTION ERA2000 ( DJ1, DJ2 )
*+
*  - - - - - - - -
*   E R A 2 0 0 0
*  - - - - - - - -
*
*  Earth rotation angle (IAU 2000 model).
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     DJ1,DJ2     d      UT1 date (JD = DJ1+DJ2)
*
*  The result is the Earth Rotation Angle (radians), in the range 0 to
*  2pi.
*
*  Calls the SOFA routine iau_ANP.
*
*  This revision:  2003 May 4
*
*-----------------------------------------------------------------------

      DOUBLE PRECISION DJ1, DJ2

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

      DOUBLE PRECISION D1, D2, T, F

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Days since fundamental epoch.
      IF ( DJ1 .LT. DJ2 ) THEN
         D1 = DJ1
         D2 = DJ2
      ELSE
         D1 = DJ2
         D2 = DJ1
      END IF
      T = D1 + ( D2-DJ0 )

*  Fractional part of T (days).
      F = MOD ( D1, 1D0 ) + MOD ( D2, 1D0 )

*  Earth rotation angle at this UT1.
      ERA2000 = iau_ANP ( D2PI * ( F + 0.7790572732640D0
     :                               + 0.00273781191135448D0 * T ) )

*  Finished.

      END
