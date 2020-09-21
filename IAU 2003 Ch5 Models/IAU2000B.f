CTITLE MHB_2000_SHORT
      subroutine MHB_2000_SHORT( jd, dpsi_ls, deps_ls,
     .                    dpsi_prec, deps_prec,
     .                    dpsi_tot , deps_tot    )
 
*     Subroutine to compute the complete MHB_2000 nutation series
*     with the associated corrections for planetary nutations,
*     the precession constant change and a rate of change of oblquity.)
 
* USAGE:
*     call  MHB_2000_SHORT( jd, dpsi_ls, deps_ls, 
*    .        dPsi_prec, deps_prec, dpsi_tot , deps_tot )
 
* where:
*     <jd>    is the full julian date including fractional part of
*             of the day (REAL*8 input)
*     <dpsi_ls> and <deps_ls> are the luni-solar nutation in
*             longitude and oblquity (mas) (REAL*8 OUTPUT)
*     <dpsi_prec> and <deps_prec> are the contributions to the
*             nutations in longitude and obliquity due changes in
*             the precession constant and rate of change of
*             obliquity (mas) (REAL*8 OUTPUT).
*     <dpsi_tot> and <deps_tot> are the total nutations in longitude
*             and obliquity including the correction for the precession
*             constant (when precession is computed using the IAU 1976
*             precession constant), and are obtained by summing all
*             of the above corrections (mas) (REAL*8 OUTPUT).
 
* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning
*               message will be printed.  See individual modules
*               for further restrictions.
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
*     dpsi_ls and deps_ls      - luni-solar nutation in
*             longitude and oblquity (mas) (REAL* OUTPUT)
*     dpsi_prec and deps_prec  - contributions to the
*             nutations in longitude and obliquity due changes in
*             the precession constant and rate of change of
*             obliquity (mas) (REAL* OUTPUT).
*     dpsi_tot and deps_tot    - total nutations in longitude
*             and obliquity including the correction for the precession
*             constant (when precession is computed using the IAU 1976
*             precession constant), and are obtained by summing all
*             of the above corrections (mas) (REAL* OUTPUT).
 
 
      real*8 jd, dpsi_ls, deps_ls, 
     .    dpsi_prec, deps_prec,
     .    dpsi_bias, deps_bias,
     .    dpsi_bias_nut, deps_bias_nut,
     .    dpsi_tot, deps_tot
 
*---------------------------------------------------------------
 
*     Call each of the routines needed for each contribution.
 
*     Luni-solar nutation
      call ls_nut_short( jd, dpsi_ls, deps_ls )
 
*     Precession and obliquity rate contributions (NOTE: IAU-1976
*     precession constant assumed to be used in the basic calculation
*     of precession).
 
      call prec_nut_short( jd, dpsi_prec, deps_prec )
 
*     Set up bias corrections

      dpsi_bias = -41.775d0
      deps_bias =  -6.819d0
      dpsi_bias_nut = -1.584d0
      deps_bias_nut = +1.634d0

*     Now add up all of the terms to get the total nutation angles
 
      dpsi_tot = dpsi_ls + dpsi_prec + dpsi_bias + dpsi_bias_nut
      deps_tot = deps_ls + deps_prec + deps_bias + deps_bias_nut
 
*     Thats all
      return
      end
 
CTITLE LS_NUT_SHORT
 
      subroutine ls_nut_short( jd, dpsi_ls, deps_ls )
 
*     Routine to compute the MHB_2000 luni-solar contributions
*     the nutations in longitude and obliquity.  The MHB_2000 is
*     based on:
 
*     (1) The Souchay and Kinoshita Rigid Earth nutation series
*     SKRE1997.   There are many duplicate argument terms in this
*     series and all of these have been compacted in single argument
*     terms.
 
*     (2) Value of the Retrograde FCN resonance factors from
*     the Mathews et al., 2000 , nutation formulation (full complex
*     estimates, and scaling parameter R from the same
*     theory.  

*     (3) The effects of annual modulation of geodetic precession.
*     The correction applied is
*      0  1  0  0  0 -0.150 (correction to in-phase nutation in
*                         longitude).

*     (4) A prograde annual nutation has been estimated along with the
*      resonance coefficients.  This probably reflects the influence of
*      S1 atmospheric tide.

*     (6) The new Simons et al., fundamental arguments are used in this
*      version.  (The largest change from KSV_1995_1 was 0.007 mas for
*      semiannual nutation.  All other changes, including the 18.6 year
*      nutation were 0.001-0.002 mas.)
 
*     REFERENCES:
* NEW Version based on: Corrections and new developments in rigid Earth
*     nutation theory: Lunisolar influence including indirect planetary
*     effects, J. Souchay and H. Kinioshita, Astron. and Astrophys., 1995.
* (Version here based on data files: SKRE1997.DPSI and SKRE1997.DEPS
*  and generated with ks_plan.f)
*     Souchay, J., and H. Kinoshita, Corrections and new developments in
*         rigid Earth nutation theory: I. Lunisolar influence including indirect
*         planetary effects, Astron. Astrophys.,  312, 1017--1030, 1996.
*     Souchay, J., and H. Kinoshita, Corrections and new developments in
*         rigid Earth nutation theory: II. Influence of second-order geopotential
*         and direct planetray effect, Astron. Astrophys., 318, 639--652, 1997.
*     Souchay, J., B. Loysel, H, Kinoshita, and M. Folgueira, Corrections
*         and new developments in rigid Earth nutation theory: III. final tables
*         REN-2000 including crossed-nutation and spin-orbit coupling effects,
*         Astron. Astrophys. Suppl., 135, 111-131, 1999.    
*     Kinoshita, H., and J. Souchay, The theory of the nutations for
*         the rigid Earth at the second order, Celes. Mech. and Dynam.
*         Astron., 48, 187--266, 1990.
*     Mathews, P. M., B. A. Buffett, and T. A. Herring, Modeling of Nutation-Precession:
*         Insights into the Earth's Interior and New Nonrigid Earth Nutation Series
*         to be submitted, J. Geophys. Res, 2000.
*     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*         Francou, G., Laskar, J., 1994, "Numerical Expressions for 
*         Precession Formulae and Mean Elements for the Moon and
*         Planets," Astron. Astrophys., 282, pp. 663-683.

 
 
* USAGE:
*     call ls_nut_short( jd, dpsi_ls, deps_ls )
*     where <jd>    is a full julian date with fractional part
*                   of the day added (REAL*8 INPUT)
*     and <dpsi_ls> and <deps_ls> are the nutations
*                   in longitude and obliquity in milliarcsec.
*                   (REAL*8 OUTPUT)
 
* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning
*               message will be printed.
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
* dpsi_ls  - The nutation in longitude (mas).
* deps_ls  - The nutation in obliquity (mas).
 
 
 
      real*8 jd, dpsi_ls, deps_ls
 
* LOCAL VARIABLES
 
*   epoch       - Julian date (jd passed in unless the JD
*                 appears to be an MJD in which case it is
*                 converted to JD (2 400 000.5d0 added)
*   ls_arg(5)   - The arguments for the Luni-solar nutations.
*                 (l, l', F, D and Omega).  All in Radians.
 
 
      real*8 epoch, ls_arg(5)
 
***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(*,100) jd
 100      format('**WARNING** MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          epoch = jd + 2 400 000.5d0
      else
          epoch = jd
      end if
 
***** Get the fundamental arguments at this epoch
 
      call ls_angles_short( epoch, ls_arg)
 
*     Now compute the luni-solare nutations by summing over all
*     terms in the series.
 
      call eval_ls_nut_short( epoch, ls_arg, dpsi_ls, deps_ls )
 
*     Thats all
      return
      end
 
 
CTITLE LS_ANGLES_SHORT
 
      subroutine ls_angles_short( epoch, ls_arg )
 
*     Routine to compute the value of the fundamental argument
*     for Brown's arguments.  Arguments based on the IERS
*     standards.

* MOD TAH 960206: Changed arguments to use Simons et al., 1994 
* values:
*     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*          Francou, G., Laskar, J., 1994, "Numerical Expressions for 
*          Precession Formulae and Mean Elements for the Moon and
*          Planets," Astron. Astrophys., 282, pp. 663-683.


* PHYSICAL CONSTANTS
 
*   pi          - Define here to full precision
*   rad_to_deg  - Conversion from radians to degs.
*   DJ2000      - Julian date of J2000
*   sec360      - number of seconds in 360 degreees.
 
 
      real*8 pi, rad_to_deg, DJ2000, sec360
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( sec360        = 1296000.d0           )
 
*     Computed quanities
      parameter ( rad_to_deg    = 180.d0   /pi         )
 
*-------------------------------------------------------------------
 
* PASSED VARIABLES
 
* INPUT
* epoch  - Julian date for arguments (jd + fraction of day, REAL*8)
 
* OUTPUT
* ls_arg(5) -  Brown's arguments (radians, REAL*8)
 
 
      real*8 epoch, ls_arg(5)
 
* LOCAL VARIABLES
*      cent             - Julian centuries to DJ2000.
*      el               - Mean longitude of moon minus mean
*                       - longitude of moon's perigee (arcsec)
*      elc(2)           - Coefficients for computing el
*      elp              - Mean longitude of the sun minus mean
*                       - longitude of sun perigee (arcsec)
*      elpc(2)          - Coeffiecents for computing elp
*      f                - Moon's mean longitude minus omega (sec)
*      fc(2)            - Coefficients for computing f
*      d                - Mean elongation of the moon from the
*                       - sun (arcsec)
*      dc(2)            - coefficients for computing d
*      om               - longitude of the ascending node of the
*                       - moon's mean orbit on the elliptic
*                       - measured from the mean equinox of date
*      omc(2)           - Coefficients for computing om.
 
 
      real*8 cent, el, elc(2), elp, elpc(2),
     .    f, fc(2), d, dc(2), om, omc(2)
 
****  DATA statements for the fundamental arguments.
*     Simons et al., 1994 values
 
      data elc    / 1717915923.2178d0,   485868.249036d0/
      data elpc   / 129596581.0481d0,   1287104.79305d0/
      data fc     / 1739527262.8478d0,    335779.526232d0/
      data dc     / 1602961601.2090d0,   1072260.70369d0/
* MOD TAH MHB_2000: 960606: Replaced <Om> with expression from b.3 of 
*     Simon et al., 1994 since b.3 is computed with new precession constant
*     (Only the rate changes).   
      data omc    / -6962890.5431d0,     450160.398036d0/

 
****  Get the number of centuries to current time
 
      cent = (epoch-dj2000) / 36525.d0
 
****  Compute angular arguments and their time derivatives
* New formulas adding in the higher order term.

      el = elc(1) * cent + elc(2)
      el = mod( el, sec360 )
c
      elp = elpc(1) * cent + elpc(2)
      elp = mod( elp, sec360 )
c
      f = fc(1) * cent + fc(2)
      f = mod( f, sec360 )
c
      d = dc(1) * cent + dc(2)
      d = mod( d, sec360 )
c
      om = omc(1) * cent + omc(2)
      om = mod( om, sec360 )
c
c
 
****  Now save the values.  Convert values from arcseconds to radians
 
      ls_arg(1) = el / (3600.d0*rad_to_deg)
      ls_arg(2) = elp/ (3600.d0*rad_to_deg)
      ls_arg(3) = f  / (3600.d0*rad_to_deg)
      ls_arg(4) = d  / (3600.d0*rad_to_deg)
      ls_arg(5) = om / (3600.d0*rad_to_deg)
 
***** Thats all
      return
      end
 
CTITLE EVAL_LS_NUT_SHORT
 
      subroutine eval_ls_nut_short( epoch, ls_arg, dpsi_ls, deps_ls )
 
*     Routine to compute the nutations in longitude and obliquity
*     by summing over all terms in the nutations series.
 
* NOTE: ls_angles must be called before routine.
 
* PARAMETERS:
 
* num_ls  - Number of terms in the nutations series
 
      integer*4 num_ls
 
      parameter ( num_ls      =  78)
 
*   DJ2000      - Julian date of J2000
*   pi          - Pi.
 
 
      real*8 pi, DJ2000
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
 
* PASSED PARAMETERS:
 
* INPUT:
* epoch      - Julian date at which nutation angles are needed.
* ls_arg(5)  - Five arguments for the nutatoions (l,l',F,D and Om)
*              computed at the epoch that the nutations need to be
*              evaluated (rad) (REAL*8)
 
* OUTPUT:
* dpsi_ls, deps_ls   - nutations in longitude and obliquity (mas)
*              (REAL*8)
 
 
      real*8 epoch, ls_arg(5), dpsi_ls, deps_ls
 
* LOCAL VARIABLES:
 
*  i and j   - Counters to loop over the coeffients and the argumemts
 
 
      integer*4 i,j
 
*  arg       - Final summed argumemt for the nutations
*              contributions (rads)
*  cent      - Number of centuries since J2000.
*  dpsi_lsu and deps_lsu - Nutations in longitude and oblquity
*              in micro-arc-sec (units that the data statements
*              are in)
*  carg, sarg -- Cosine and sine of arguments.
 
 
      real*8 arg, cent, dpsi_lsu, deps_lsu, carg, sarg
 
*     RFCN Freq.  -1.00231810920 cyc per sidreal day, Period   430.2082 solar days
*  Units now in 0.1 microarcsec
*  IX01-IX10(11,10)  -- Invidual declarations of the coefficents
*                           of the nutation series so no data statement has more than 10 lines
*                           The first 5 values are the arguments for l lp F D Om
*                           The remaining elements are:
*                            6 - Nutation in longitude psi (sin, uas)
*                            7 - dpsi/dt (uasec/cent)
*                            8 - Nutation in oblquity eps (cos, uas)
*                            9 - deps/dt (uas/cent)
*                           10 - Out-of-phase longitude (cos, uas)
*                           11 - Out-of-phase obliquity (sin, uas)

      integer*4 IX01(11,10), IX02(11,10), IX03(11,10), IX04(11,10), 
     .          IX05(11,10), IX06(11,10), IX07(11,10), IX08(11,8)

      integer*4 nutc_int(11,78)
      equivalence (nutc_int(1,  1),IX01(1,1))
      equivalence (nutc_int(1, 11),IX02(1,1))
      equivalence (nutc_int(1, 21),IX03(1,1))
      equivalence (nutc_int(1, 31),IX04(1,1))
      equivalence (nutc_int(1, 41),IX05(1,1))
      equivalence (nutc_int(1, 51),IX06(1,1))
      equivalence (nutc_int(1, 61),IX07(1,1))
      equivalence (nutc_int(1, 71),IX08(1,1))

      data IX01/   0,   0,   0,   0,   1,-172064161,-174666,
     .                    92052331,  9086, 33386, 15377,
     .            0,   0,   2,  -2,   2, -13170906,   -1675,
     .                     5730336, -3015,-13696, -4587,
     .            0,   0,   2,   0,   2,  -2276413,    -234,
     .                      978459,  -485,  2796,  1374,
     .            0,   0,   0,   0,   2,   2074554,     207,
     .                     -897492,   470,  -698,  -291,
     .            0,   1,   0,   0,   0,   1475877,   -3633,
     .                       73871,  -184, 11817, -1924,
     .            0,   1,   2,  -2,   2,   -516821,    1226,
     .                      224386,  -677,  -524,  -174,
     .            1,   0,   0,   0,   0,    711159,      73,
     .                       -6750,     0,  -872,   358,
     .            0,   0,   2,   0,   1,   -387298,    -367,
     .                      200728,    18,   380,   318,
     .            1,   0,   2,   0,   2,   -301461,     -36,
     .                      129025,   -63,   816,   367,
     .            0,  -1,   2,  -2,   2,    215829,    -494,
     .                      -95929,   299,   111,   132 /
      data IX02/   0,   0,   2,  -2,   1,    128227,    137,
     .                      -68982,    -9,   181,    39,
     .           -1,   0,   2,   0,   2,    123457,      11,
     .                      -53311,    32,    19,    -4,
     .           -1,   0,   0,   2,   0,    156994,      10,
     .                       -1235,     0,  -168,    82,
     .            1,   0,   0,   0,   1,     63110,      63,
     .                      -33228,     0,    27,    -9,
     .           -1,   0,   0,   0,   1,    -57976,     -63,
     .                       31429,     0,  -189,   -75,
     .           -1,   0,   2,   2,   2,    -59641,     -11,
     .                       25543,   -11,   149,    66,
     .            1,   0,   2,   0,   1,    -51613,     -42,
     .                       26366,     0,   129,    78,
     .           -2,   0,   2,   0,   1,     45893,      50,
     .                      -24236,   -10,    31,    20,
     .            0,   0,   0,   2,   0,     63384,      11,
     .                       -1220,     0,  -150,    29,
     .            0,   0,   2,   2,   2,    -38571,      -1,
     .                       16452,   -11,   158,    68 /
      data IX03/  -2,   0,   0,   2,   0,    -47722,       0,
     .                         477,     0,   -18,   -25,
     .            2,   0,   2,   0,   2,    -31046,      -1,
     .                       13238,   -11,   131,    59,
     .            1,   0,   2,  -2,   2,     28593,       0,
     .                      -12338,    10,    -1,    -3,
     .           -1,   0,   2,   0,   1,     20441,      21,
     .                      -10758,     0,    10,    -3,
     .            2,   0,   0,   0,   0,     29243,       0,
     .                        -609,     0,   -74,    13,
     .            0,   0,   2,   0,   0,     25887,       0,
     .                        -550,     0,   -66,    11,
     .            0,   1,   0,   0,   1,    -14053,     -25,
     .                        8551,    -2,    79,   -45,
     .           -1,   0,   0,   2,   1,     15164,      10,
     .                       -8001,     0,    11,    -1,
     .            0,   2,   2,  -2,   2,    -15794,      72,
     .                        6850,   -42,   -16,    -5,
     .            0,   0,  -2,   2,   0,     21783,      0,
     .                        -167,     0,    13,    13 /
      data IX04/   1,   0,   0,  -2,   1,    -12873,     -10,
     .                        6953,     0,   -37,   -14,
     .            0,  -1,   0,   0,   1,    -12654,      11,
     .                        6415,     0,    63,    26,
     .           -1,   0,   2,   2,   1,    -10204,       0,
     .                        5222,     0,    25,    15,
     .            0,   2,   0,   0,   0,     16707,     -85,
     .                         168,    -1,   -10,    10,
     .            1,   0,   2,   2,   2,     -7691,       0,
     .                        3268,     0,    44,    19,
     .           -2,   0,   2,   0,   0,    -11024,       0,
     .                         104,     0,   -14,     2,
     .            0,   1,   2,   0,   2,      7566,     -21,
     .                       -3250,     0,   -11,    -5,
     .            0,   0,   2,   2,   1,     -6637,     -11,
     .                        3353,     0,    25,    14,
     .            0,  -1,   2,   0,   2,     -7141,      21,
     .                        3070,     0,     8,     4,
     .            0,   0,   0,   2,   1,     -6302,    -11,
     .                        3272,     0,     2,     4 /
      data IX05/   1,   0,   2,  -2,   1,      5800,      10,
     .                       -3045,     0,     2,    -1,
     .            2,   0,   2,  -2,   2,      6443,       0,
     .                       -2768,     0,    -7,    -4,
     .           -2,   0,   0,   2,   1,     -5774,     -11,
     .                        3041,     0,   -15,    -5,
     .            2,   0,   2,   0,   1,     -5350,       0,
     .                        2695,     0,    21,    12,
     .            0,  -1,   2,  -2,   1,     -4752,     -11,
     .                        2719,     0,    -3,    -3,
     .            0,   0,   0,  -2,   1,     -4940,     -11,
     .                        2720,     0,   -21,    -9,
     .           -1,  -1,   0,   2,   0,      7350,       0,
     .                         -51,     0,    -8,     4,
     .            2,   0,   0,  -2,   1,      4065,       0,
     .                       -2206,     0,     6,     1,
     .            1,   0,   0,   2,   0,      6579,       0,
     .                        -199,     0,   -24,     2,
     .            0,   1,   2,  -2,   1,      3579,      0,
     .                       -1900,     0,     5,     1 /
      data IX06/   1,  -1,   0,   0,   0,      4725,       0,
     .                         -41,     0,    -6,     3,
     .           -2,   0,   2,   0,   2,     -3075,       0,
     .                        1313,     0,    -2,    -1,
     .            3,   0,   2,   0,   2,     -2904,       0,
     .                        1233,     0,    15,     7,
     .            0,  -1,   0,   2,   0,      4348,       0,
     .                         -81,     0,   -10,     2,
     .            1,  -1,   2,   0,   2,     -2878,       0,
     .                        1232,     0,     8,     4,
     .            0,   0,   0,   1,   0,     -4230,       0,
     .                         -20,     0,     5,    -2,
     .           -1,  -1,   2,   2,   2,     -2819,       0,
     .                        1207,     0,     7,     3,
     .           -1,   0,   2,   0,   0,     -4056,       0,
     .                          40,     0,     5,    -2,
     .            0,  -1,   2,   2,   2,     -2647,       0,
     .                        1129,     0,    11,     5,
     .           -2,   0,   0,   0,   1,     -2294,      0,
     .                        1266,     0,   -10,    -4 /
      data IX07/   1,   1,   2,   0,   2,      2481,       0,
     .                       -1062,     0,    -7,    -3,
     .            2,   0,   0,   0,   1,      2179,       0,
     .                       -1129,     0,    -2,    -2,
     .           -1,   1,   0,   1,   0,      3276,       0,
     .                          -9,     0,     1,     0,
     .            1,   1,   0,   0,   0,     -3389,       0,
     .                          35,     0,     5,    -2,
     .            1,   0,   2,   0,   0,      3339,       0,
     .                        -107,     0,   -13,     1,
     .           -1,   0,   2,  -2,   1,     -1987,       0,
     .                        1073,     0,    -6,    -2,
     .            1,   0,   0,   0,   2,     -1981,       0,
     .                         854,     0,     0,     0,
     .           -1,   0,   0,   1,   0,      4026,       0,
     .                        -553,     0,  -353,  -139,
     .            0,   0,   2,   1,   2,      1660,       0,
     .                        -710,     0,    -5,    -2,
     .           -1,   0,   2,   4,   2,     -1521,      0,
     .                         647,     0,     9,     4 /
      data IX08/  -1,   1,   0,   1,   1,      1314,       0,
     .                        -700,     0,     0,     0,
     .            0,  -2,   2,  -2,   1,     -1283,       0,
     .                         672,     0,     0,     0,
     .            1,   0,   2,   2,   1,     -1331,       0,
     .                         663,     0,     8,     4,
     .           -2,   0,   2,   2,   2,      1383,       0,
     .                        -594,     0,    -2,    -2,
     .           -1,   0,   0,   0,   2,      1405,       0,
     .                        -610,     0,     4,     2,
     .            1,   1,   2,  -2,   2,      1290,       0,
     .                        -556,     0,     0,     0,
     .           -2,   0,   2,   4,   2,     -1214,       0,
     .                         518,     0,     5,     2,
     .           -1,   0,   4,   0,   2,      1146,       0,
     .                        -490,     0,    -3,    -1 /
 
****  Initialize the values and sum over the series
 
      dpsi_lsu = 0.0d0
      deps_lsu = 0.0d0
 
      cent = (epoch-DJ2000) / 36525.d0
 
      do 200 i = num_ls, 1, -1
 
*         Sum the mulitpliers by the arguments to the argument of
*         nutation
          arg = 0.d0
          do 150 j = 1,5
 
*            Sum into the argument for nutation.
             arg = arg + nutc_int(j,i)*ls_arg(j)
 150      continue
 
          arg = mod(arg, 2.d0*pi)
          carg = cos(arg)
          sarg = sin(arg)
 
****      Now add contributions to dpsi and deps
          dpsi_lsu = dpsi_lsu +
     .               (nutc_int( 6,i)+ nutc_int(7,i)*cent)*sarg +
     .                nutc_int(10,i)*carg
          deps_lsu = deps_lsu +
     .               (nutc_int( 8,i)+ nutc_int(9,i)*cent)*carg +
     .                nutc_int(11,i)*sarg
 
 200  continue
 
*     Convert values from 0.1 micro-arc-sec to mill-arc-second
      dpsi_ls = dpsi_lsu * 1.d-4
      deps_ls = deps_lsu * 1.d-4
 
****  Thats all
      return
      end
 
CTITLE PREC_NUT_SHORT
 
      subroutine prec_nut_short( jd, dpsi_prec, deps_prec )
 
*     Routine to evaluate the corrections to the nutations in longitude
*     and obliquity due to the corrections to the IAU-1976 Luni-solar
*     precession constant and the secular rate of change of the obliquity
*     of the ecliptic.
 
* PARAMETERS:
 
*   DJ2000      - Julian date of J2000
 
 
      real*8 DJ2000
 
      parameter ( DJ2000        = 2451545.d0           )
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
* dpsi_prec   - Contribution to the nutation in longitude (mas).  Should
*          be added to standard nutation in longitude values. Value
*          valid only when the IAU-1976 precession constant used to
*          compute the transformation to mean system.
* deps_prec   - Contribution to the nutation in obliquity (mas).  Should
*          be added to standard nutation in obliquity values.
 
 
      real*8 jd, dpsi_prec, deps_prec
 
* LOCAL VARIABLES
 
*   epoch       - Julian date (jd passed in unless the JD
*                 appears to be an MJD in which case it is
*                 converted to JD (2 400 000.5d0 added)
*   cent        - Number of Julian centuries since J2000.0
*   DpsiDt      - Correction to precession constant as a
*                 linear rate of change of nutation in
*                 longitude. (arc-second/century)
*   DepsDt      - Correction to rate of change of oblquity
*                 (arc-second/century)
 
 
 
      real*8 epoch, cent, DpsiDt,  DepsDt
*
*     Theoretical estimaate of adjustment to precession constant is
*     -0.29965 "/cent; estimate is -0.2997 +- 0.0008 "/cent.
  
      data  DpsiDt  /  -0.2997d0  /
      data  DepsDt  /  -0.0252d0  /
 
***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(*,100) jd
 100      format('**WARNING** MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          epoch = jd + 2 400 000.5d0
      else
          epoch = jd
      end if
 
****  Compute the number of centuries
 
      cent = (epoch - DJ2000)/36525.d0
 
      dpsi_prec = DpsiDt*cent*1000.d0
      deps_prec = DepsDt*cent*1000.d0
 
*     Thats all
      return
      end
