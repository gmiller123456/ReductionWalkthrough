
const DJC = (36525.0);

/* Days per Julian millennium */
const DJM = (365250.0);

/* Reference epoch (J2000.0), Julian Date */
const DJ00 = (2451545.0);

const D2PI = (6.283185307179586476925287);
const DAS2R =(4.848136811095359935899141e-6);
const TURNAS =(1296000.0);


function fmod(a,m){
    return a%m;
}

class nutation{
    static iauGst00b(uta, utb){
        let gmst00, ee00b, gst;
     
     
        gmst00 = this.iauGmst00(uta, utb, uta, utb);
        ee00b = this.iauEe00b(uta, utb);
        gst = this.iauAnp(gmst00 + ee00b);
     
        return gst;
     }

     static iauGmst00(uta, utb, tta, ttb){
        let t, gmst;
     
     
     /* TT Julian centuries since J2000.0. */
        t = ((tta - DJ00) + ttb) / DJC;
     
     /* Greenwich Mean Sidereal Time, IAU 2000. */
        gmst = this.iauAnp(this.iauEra00(uta, utb) +
                        (     0.014506   +
                        (  4612.15739966 +
                        (     1.39667721 +
                        (    -0.00009344 +
                        (     0.00001882 )
               * t) * t) * t) * t) * DAS2R);
     
        return gmst;
     }

     static iauEra00(dj1, dj2){
        let d1, d2, t, f, theta;
     
     /* Days since fundamental epoch. */
        if (dj1 < dj2) {
           d1 = dj1;
           d2 = dj2;
        } else {
           d1 = dj2;
           d2 = dj1;
        }
        t = d1 + (d2- DJ00);
     
     /* Fractional part of T (days). */
        f = fmod(d1, 1.0) + fmod(d2, 1.0);
     
     /* Earth rotation angle at this UT1. */
        theta = this.iauAnp(D2PI * (f + 0.7790572732640
                                 + 0.00273781191135448 * t));
     
        return theta;
     }

     static iauEe00b(date1, date2){
        let  epsa, ee;
     
     
     /* IAU 2000 precession-rate adjustments. */
        const pr=this.iauPr00(date1, date2);
        const dpsipr=pr[0];
        const depspr=pr[1];
     
     /* Mean obliquity, consistent with IAU 2000 precession-nutation. */
        epsa = this.iauObl80(date1, date2) + depspr;
     
     /* Nutation in longitude. */
        const t = ((date1 - DJ00) + date2) / DJC;
        const nut=this.nutation(t);
        const dpsi=nut[0];
        const deps=nut[1];
     
     /* Equation of the equinoxes. */
        ee = this.iauEe00(date1, date2, epsa, dpsi);
     
        return ee;
     }

     static iauObl80(date1, date2){
        let t, eps0;
     
     
     /* Interval between fundamental epoch J2000.0 and given date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
     
     /* Mean obliquity of date. */
        eps0 = DAS2R * (84381.448  +
                       (-46.8150   +
                       (-0.00059   +
                       ( 0.001813) * t) * t) * t);
     
        return eps0;
     }

     static iauPr00(date1, date2){
        let t;
     
     /* Precession and obliquity corrections (radians per century) */
        const PRECOR = -0.29965 * DAS2R,
                            OBLCOR = -0.02524 * DAS2R;
     
     
     /* Interval between fundamental epoch J2000.0 and given date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
     
     /* Precession rate contributions with respect to IAU 1976/80. */
        let ret=new Array();
        ret[0] = PRECOR * t;
        ret[1] = OBLCOR * t;
     
        return ret;
     }
     
    static iauEe00(date1, date2, epsa, dpsi){
        let ee;
     /* Equation of the equinoxes. */
        ee = dpsi * Math.cos(epsa) + this.iauEect00(date1, date2);
     
        return ee;
     }

     static iauAnp(a){
        let w;
     
     
        w = fmod(a, D2PI);
        if (w < 0) w += D2PI;
     
        return w;
     }
     
    static iauEect00( date1,  date2){
    /* Time since J2000.0, in Julian centuries */
    let t;

    /* Miscellaneous */
        let i, j;
    let a, s0, s1;

    /* Fundamental arguments */
    let fa=new Array();

    /* Returned value. */
    let eect;

    /* ----------------------------------------- */
    /* The series for the EE complementary terms */
    /* ----------------------------------------- */

    /* Terms of order t^0 */
    const e0 = [

    /* 1-10 */
        [[ 0,  0,  0,  0,  1,  0,  0,  0], 2640.96e-6, -0.39e-6 ],
        [[ 0,  0,  0,  0,  2,  0,  0,  0],   63.52e-6, -0.02e-6 ],
        [[ 0,  0,  2, -2,  3,  0,  0,  0],   11.75e-6,  0.01e-6 ],
        [[ 0,  0,  2, -2,  1,  0,  0,  0],   11.21e-6,  0.01e-6 ],
        [[ 0,  0,  2, -2,  2,  0,  0,  0],   -4.55e-6,  0.00e-6 ],
        [[ 0,  0,  2,  0,  3,  0,  0,  0],    2.02e-6,  0.00e-6 ],
        [[ 0,  0,  2,  0,  1,  0,  0,  0],    1.98e-6,  0.00e-6 ],
        [[ 0,  0,  0,  0,  3,  0,  0,  0],   -1.72e-6,  0.00e-6 ],
        [[ 0,  1,  0,  0,  1,  0,  0,  0],   -1.41e-6, -0.01e-6 ],
        [[ 0,  1,  0,  0, -1,  0,  0,  0],   -1.26e-6, -0.01e-6 ],

    /* 11-20 */
        [[ 1,  0,  0,  0, -1,  0,  0,  0],   -0.63e-6,  0.00e-6 ],
        [[ 1,  0,  0,  0,  1,  0,  0,  0],   -0.63e-6,  0.00e-6 ],
        [[ 0,  1,  2, -2,  3,  0,  0,  0],    0.46e-6,  0.00e-6 ],
        [[ 0,  1,  2, -2,  1,  0,  0,  0],    0.45e-6,  0.00e-6 ],
        [[ 0,  0,  4, -4,  4,  0,  0,  0],    0.36e-6,  0.00e-6 ],
        [[ 0,  0,  1, -1,  1, -8, 12,  0],   -0.24e-6, -0.12e-6 ],
        [[ 0,  0,  2,  0,  0,  0,  0,  0],    0.32e-6,  0.00e-6 ],
        [[ 0,  0,  2,  0,  2,  0,  0,  0],    0.28e-6,  0.00e-6 ],
        [[ 1,  0,  2,  0,  3,  0,  0,  0],    0.27e-6,  0.00e-6 ],
        [[ 1,  0,  2,  0,  1,  0,  0,  0],    0.26e-6,  0.00e-6 ],

    /* 21-30 */
        [[ 0,  0,  2, -2,  0,  0,  0,  0],   -0.21e-6,  0.00e-6 ],
        [[ 0,  1, -2,  2, -3,  0,  0,  0],    0.19e-6,  0.00e-6 ],
        [[ 0,  1, -2,  2, -1,  0,  0,  0],    0.18e-6,  0.00e-6 ],
        [[ 0,  0,  0,  0,  0,  8,-13, -1],   -0.10e-6,  0.05e-6 ],
        [[ 0,  0,  0,  2,  0,  0,  0,  0],    0.15e-6,  0.00e-6 ],
        [[ 2,  0, -2,  0, -1,  0,  0,  0],   -0.14e-6,  0.00e-6 ],
        [[ 1,  0,  0, -2,  1,  0,  0,  0],    0.14e-6,  0.00e-6 ],
        [[ 0,  1,  2, -2,  2,  0,  0,  0],   -0.14e-6,  0.00e-6 ],
        [[ 1,  0,  0, -2, -1,  0,  0,  0],    0.14e-6,  0.00e-6 ],
        [[ 0,  0,  4, -2,  4,  0,  0,  0],    0.13e-6,  0.00e-6 ],

    /* 31-33 */
        [[ 0,  0,  2, -2,  4,  0,  0,  0],   -0.11e-6,  0.00e-6 ],
        [[ 1,  0, -2,  0, -3,  0,  0,  0],    0.11e-6,  0.00e-6 ],
        [[ 1,  0, -2,  0, -1,  0,  0,  0],    0.11e-6,  0.00e-6 ]
    ];

    /* Terms of order t^1 */
    const e1 = [
        [[ 0,  0,  0,  0,  1,  0,  0,  0],    -0.87e-6,  0.00e-6 ]
    ];

    /*--------------------------------------------------------------------*/

    /* Interval between fundamental epoch J2000.0 and current date (JC). */
    t = ((date1 - DJ00) + date2) / DJC;

    /* Fundamental Arguments (from IERS Conventions 2003) */

    /* Mean anomaly of the Moon. */
    fa[0] = this.iauFal03(t);

    /* Mean anomaly of the Sun. */
    fa[1] = this.iauFalp03(t);

    /* Mean longitude of the Moon minus that of the ascending node. */
    fa[2] = this.iauFaf03(t);

    /* Mean elongation of the Moon from the Sun. */
    fa[3] = this.iauFad03(t);

    /* Mean longitude of the ascending node of the Moon. */
    fa[4] = this.iauFaom03(t);

    /* Mean longitude of Venus. */
    fa[5] = this.iauFave03(t);

    /* Mean longitude of Earth. */
    fa[6] = this.iauFae03(t);

    /* General precession in longitude. */
    fa[7] = this.iauFapa03(t);

    /* Evaluate the EE complementary terms. */
    s0 = 0.0;
    s1 = 0.0;

    for (i = e0.length-1; i >= 0; i--) {
        a = 0.0;
        for (j = 0; j < 8; j++) {
            a += (e0[i][0][j]) * fa[j];
        }
        s0 += e0[i][1] * Math.sin(a) + e0[i][2] * Math.cos(a);
    }

    for (i = e1.length-1; i >= 0; i--) {
        a = 0.0;
        for (j = 0; j < 8; j++) {
            a += (e1[i][0][j]) * fa[j];
        }
        s1 += e1[i][1] * Math.sin(a) + e1[i][2] * Math.cos(a);
    }

    eect = (s0 + s1 * t ) * DAS2R;

    return eect;

    }
    static iauFapa03(t)
    {
       let a;
    
    
    /* General accumulated precession in longitude. */
       a = (0.024381750 + 0.00000538691 * t) * t;
    
       return a;
    
    }

    static iauFae03(t){
    let a;


    /* Mean longitude of Earth (IERS Conventions 2003). */
    a = fmod(1.753470314 + 628.3075849991 * t, D2PI);

    return a;

    }


    static iauFave03(t){
        let a;
    
    
    /* Mean longitude of Venus (IERS Conventions 2003). */
        a = fmod(3.176146697 + 1021.3285546211 * t, D2PI);
    
        return a;
    }

    static iauFaom03(t){
        let a;
     
     
     /* Mean longitude of the Moon's ascending node */
     /* (IERS Conventions 2003).                    */
        a = fmod(          450160.398036 +
                  t * ( - 6962890.5431 +
                  t * (         7.4722 +
                  t * (         0.007702 +
                  t * (       - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;
     
        return a;
     
     }
     
     static iauFad03(t){
        let a;
     
     
     /* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
        a = fmod(          1072260.703692 +
                  t * ( 1602961601.2090 +
                  t * (        - 6.3706 +
                  t * (          0.006593 +
                  t * (        - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;
     
        return a;
     
     }

     static iauFaf03(t){
        let a;
     
     
     /* Mean longitude of the Moon minus that of the ascending node */
     /* (IERS Conventions 2003).                                    */
        a = fmod(           335779.526232 +
                  t * ( 1739527262.8478 +
                  t * (       - 12.7512 +
                  t * (        - 0.001037 +
                  t * (          0.00000417 ) ) ) ), TURNAS ) * DAS2R;
     
        return a;
     }
     
     static iauFalp03(t){
        let a;
     
     
     /* Mean anomaly of the Sun (IERS Conventions 2003). */
        a = fmod(         1287104.793048 +
                  t * ( 129596581.0481 +
                  t * (       - 0.5532 +
                  t * (         0.000136 +
                  t * (       - 0.00001149 ) ) ) ), TURNAS ) * DAS2R;
     
        return a;
     }

     static iauFal03( t){
        let a;
     
     
     /* Mean anomaly of the Moon (IERS Conventions 2003). */
        a = fmod(           485868.249036  +
                  t * ( 1717915923.2178 +
                  t * (         31.8792 +
                  t * (          0.051635 +
                  t * (        - 0.00024470 ) ) ) ), TURNAS ) * DAS2R;
     
        return a;
     
     }
     static nutation(T){
		//IAU 2000B Nutation.
		//Transformed from reference implementation NU2000B.f
		//Annexe to IERS Conventions 2000, Chapter 5

		const NALS=[
			[  0,    0,    0,    0,    1,],
			[  0,    0,    2,   -2,    2,],
			[  0,    0,    2,    0,    2,],
			[  0,    0,    0,    0,    2,],
			[  0,    1,    0,    0,    0,],
			[  0,    1,    2,   -2,    2,],
			[  1,    0,    0,    0,    0,],
			[  0,    0,    2,    0,    1,],
			[  1,    0,    2,    0,    2,],
			[  0,   -1,    2,   -2,    2,],
			[  0,    0,    2,   -2,    1,],
			[ -1,    0,    2,    0,    2,],
			[ -1,    0,    0,    2,    0,],
			[  1,    0,    0,    0,    1,],
			[ -1,    0,    0,    0,    1,],
			[ -1,    0,    2,    2,    2,],
			[  1,    0,    2,    0,    1,],
			[ -2,    0,    2,    0,    1,],
			[  0,    0,    0,    2,    0,],
			[  0,    0,    2,    2,    2,],
			[  0,   -2,    2,   -2,    2,],
			[ -2,    0,    0,    2,    0,],
			[  2,    0,    2,    0,    2,],
			[  1,    0,    2,   -2,    2,],
			[ -1,    0,    2,    0,    1,],
			[  2,    0,    0,    0,    0,],
			[  0,    0,    2,    0,    0,],
			[  0,    1,    0,    0,    1,],
			[ -1,    0,    0,    2,    1,],
			[  0,    2,    2,   -2,    2,],
			[  0,    0,   -2,    2,    0,],
			[  1,    0,    0,   -2,    1,],
			[  0,   -1,    0,    0,    1,],
			[ -1,    0,    2,    2,    1,],
			[  0,    2,    0,    0,    0,],
			[  1,    0,    2,    2,    2,],
			[ -2,    0,    2,    0,    0,],
			[  0,    1,    2,    0,    2,],
			[  0,    0,    2,    2,    1,],
			[  0,   -1,    2,    0,    2,],
			[  0,    0,    0,    2,    1,],
			[  1,    0,    2,   -2,    1,],
			[  2,    0,    2,   -2,    2,],
			[ -2,    0,    0,    2,    1,],
			[  2,    0,    2,    0,    1,],
			[  0,   -1,    2,   -2,    1,],
			[  0,    0,    0,   -2,    1,],
			[ -1,   -1,    0,    2,    0,],
			[  2,    0,    0,   -2,    1,],
			[  1,    0,    0,    2,    0,],
			[  0,    1,    2,   -2,    1,],
			[  1,   -1,    0,    0,    0,],
			[ -2,    0,    2,    0,    2,],
			[  3,    0,    2,    0,    2,],
			[  0,   -1,    0,    2,    0,],
			[  1,   -1,    2,    0,    2,],
			[  0,    0,    0,    1,    0,],
			[ -1,   -1,    2,    2,    2,],
			[ -1,    0,    2,    0,    0,],
			[  0,   -1,    2,    2,    2,],
			[ -2,    0,    0,    0,    1,],
			[  1,    1,    2,    0,    2,],
			[  2,    0,    0,    0,    1,],
			[ -1,    1,    0,    1,    0,],
			[  1,    1,    0,    0,    0,],
			[  1,    0,    2,    0,    0,],
			[ -1,    0,    2,   -2,    1,],
			[  1,    0,    0,    0,    2,],
			[ -1,    0,    0,    1,    0,],
			[  0,    0,    2,    1,    2,],
			[ -1,    0,    2,    4,    2,],
			[ -1,    1,    0,    1,    1,],
			[  0,   -2,    2,   -2,    1,],
			[  1,    0,    2,    2,    1,],
			[ -2,    0,    2,    2,    2,],
			[ -1,    0,    0,    0,    2,],
			[  1,    1,    2,   -2,    2]
		];

		const CLS=[
			[ -172064161, -174666,  33386, 92052331,  9086, 15377,],
			[  -13170906,   -1675, -13696,  5730336, -3015, -4587,],
			[   -2276413,    -234,   2796,   978459,  -485,  1374,],
			[    2074554,     207,   -698,  -897492,   470,  -291,],
			[    1475877,   -3633,  11817,    73871,  -184, -1924,],
			[    -516821,    1226,   -524,   224386,  -677,  -174,],
			[     711159,      73,   -872,    -6750,     0,   358,],
			[    -387298,    -367,    380,   200728,    18,   318,],
			[    -301461,     -36,    816,   129025,   -63,   367,],
			[     215829,    -494,    111,   -95929,   299,   132,],
			[     128227,     137,    181,   -68982,    -9,    39,],
			[     123457,      11,     19,   -53311,    32,    -4,],
			[     156994,      10,   -168,    -1235,     0,    82,],
			[      63110,      63,     27,   -33228,     0,    -9,],
			[     -57976,     -63,   -189,    31429,     0,   -75,],
			[     -59641,     -11,    149,    25543,   -11,    66,],
			[     -51613,     -42,    129,    26366,     0,    78,],
			[      45893,      50,     31,   -24236,   -10,    20,],
			[      63384,      11,   -150,    -1220,     0,    29,],
			[     -38571,      -1,    158,    16452,   -11,    68,],
			[      32481,       0,      0,   -13870,     0,     0,],
			[     -47722,       0,    -18,      477,     0,   -25,],
			[     -31046,      -1,    131,    13238,   -11,    59,],
			[      28593,       0,     -1,   -12338,    10,    -3,],
			[      20441,      21,     10,   -10758,     0,    -3,],
			[      29243,       0,    -74,     -609,     0,    13,],
			[      25887,       0,    -66,     -550,     0,    11,],
			[     -14053,     -25,     79,     8551,    -2,   -45,],
			[      15164,      10,     11,    -8001,     0,    -1,],
			[     -15794,      72,    -16,     6850,   -42,    -5,],
			[      21783,       0,     13,     -167,     0,    13,],
			[     -12873,     -10,    -37,     6953,     0,   -14,],
			[     -12654,      11,     63,     6415,     0,    26,],
			[     -10204,       0,     25,     5222,     0,    15,],
			[      16707,     -85,    -10,      168,    -1,    10,],
			[      -7691,       0,     44,     3268,     0,    19,],
			[     -11024,       0,    -14,      104,     0,     2,],
			[       7566,     -21,    -11,    -3250,     0,    -5,],
			[      -6637,     -11,     25,     3353,     0,    14,],
			[      -7141,      21,      8,     3070,     0,     4,],
			[      -6302,     -11,      2,     3272,     0,     4,],
			[       5800,      10,      2,    -3045,     0,    -1,],
			[       6443,       0,     -7,    -2768,     0,    -4,],
			[      -5774,     -11,    -15,     3041,     0,    -5,],
			[      -5350,       0,     21,     2695,     0,    12,],
			[      -4752,     -11,     -3,     2719,     0,    -3,],
			[      -4940,     -11,    -21,     2720,     0,    -9,],
			[       7350,       0,     -8,      -51,     0,     4,],
			[       4065,       0,      6,    -2206,     0,     1,],
			[       6579,       0,    -24,     -199,     0,     2,],
			[       3579,       0,      5,    -1900,     0,     1,],
			[       4725,       0,     -6,      -41,     0,     3,],
			[      -3075,       0,     -2,     1313,     0,    -1,],
			[      -2904,       0,     15,     1233,     0,     7,],
			[       4348,       0,    -10,      -81,     0,     2,],
			[      -2878,       0,      8,     1232,     0,     4,],
			[      -4230,       0,      5,      -20,     0,    -2,],
			[      -2819,       0,      7,     1207,     0,     3,],
			[      -4056,       0,      5,       40,     0,    -2,],
			[      -2647,       0,     11,     1129,     0,     5,],
			[      -2294,       0,    -10,     1266,     0,    -4,],
			[       2481,       0,     -7,    -1062,     0,    -3,],
			[       2179,       0,     -2,    -1129,     0,    -2,],
			[       3276,       0,      1,       -9,     0,     0,],
			[      -3389,       0,      5,       35,     0,    -2,],
			[       3339,       0,    -13,     -107,     0,     1,],
			[      -1987,       0,     -6,     1073,     0,    -2,],
			[      -1981,       0,      0,      854,     0,     0,],
			[       4026,       0,   -353,     -553,     0,  -139,],
			[       1660,       0,     -5,     -710,     0,    -2,],
			[      -1521,       0,      9,      647,     0,     4,],
			[       1314,       0,      0,     -700,     0,     0,],
			[      -1283,       0,      0,      672,     0,     0,],
			[      -1331,       0,      8,      663,     0,     4,],
			[       1383,       0,     -2,     -594,     0,    -2,],
			[       1405,       0,      4,     -610,     0,     2,],
			[ 	    1290,       0,      0,     -556,     0,     0]
		];

		const TURNAS = 1296000;
		const DAS2R = 4.848136811095359935899141E-6;
		const D2PI = 6.283185307179586476925287;
		const U2R = DAS2R/1E7;
		const DMAS2R = DAS2R / 1E3;

		const DPPLAN = - 0.135 * DMAS2R;
		const DEPLAN = + 0.388 * DMAS2R;

		//const T=(jd_tdb-2451545.5)/36525.0;

		//  Mean anomaly of the Moon.
		const EL  = ((485868.249036 + T * 1717915923.2178) % TURNAS) * DAS2R;
  
  		//  Mean anomaly of the Sun.
		const ELP = ((1287104.79305 + T *  129596581.0481) % TURNAS ) * DAS2R;
  
  		//  Mean argument of the latitude of the Moon.
		const F   = ((335779.526232 + T * 1739527262.8478) % TURNAS ) * DAS2R;
  
  		//  Mean elongation of the Moon from the Sun.
		const D   = ((1072260.70369 + T * 1602961601.2090) % TURNAS ) * DAS2R;
  
  		//  Mean longitude of the ascending node of the Moon.
		const OM  = ((450160.398036 - T *    6962890.5431) % TURNAS ) * DAS2R;
  
  		//  Initialize the nutation values.
		let DP = 0;
		let DE = 0;
  
		//  Summation of luni-solar nutation series (in reverse order).
		for(let I=NALS.length-1;I>=0;I--){
  
  			//     Argument and functions.
		   const ARG = (  NALS[I][0]  * EL  +
							NALS[I][1]  * ELP +
							NALS[I][2]  * F   +
							NALS[I][3]  * D   +
							NALS[I][4]  * OM) % D2PI ;
			
		   const SARG = Math.sin(ARG);
		   const CARG = Math.cos(ARG);
  
  			//     Term.
		   DP = DP + ( CLS[I][0] + CLS[I][1] * T ) * SARG +   CLS[I][2] * CARG
		   DE = DE + ( CLS[I][3] + CLS[I][4] * T ) * CARG +   CLS[I][5] * SARG
  
		}
  
  		//  Convert from 0.1 microarcsec units to radians.
		const DPSILS = DP * U2R
		const DEPSLS = DE * U2R
  
		//  ------------------
		//  PLANETARY NUTATION
		//  ------------------
		
		//  Fixed terms to allow for long-period nutation.
		const DPSIPL = DPPLAN
		const DEPSPL = DEPLAN
  
		//  -----
		//  TOTAL
		//  -----
		
		//  Add planetary and luni-solar components.
		const DPSI = DPSIPL + DPSILS
		const DEPS = DEPSPL + DEPSLS
  
		return [DPSI,DEPS];
	}
     
     
}