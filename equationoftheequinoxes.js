
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

class EquationOfTheEquinoxes{
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
        const nut=Reduce7.nutation(t);
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
     
     
}