#include <stdio.h>
#include "eph_manager.h"
#include "novas.h"

void test(double jd_utc, short body, on_surface* geo_loc, short int leap_secs, double ut1_utc)
{
    double jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
    double jd_ut1 = jd_utc + ut1_utc / 86400.0;
    double delta_t = 32.184 + leap_secs - ut1_utc;
    //delta_t=69.402165200000;
    printf("JD: %15.12f\r\n",jd_utc);
    printf("TT: %15.12f\r\n",jd_tt);
    printf("Delta T: %15.12f\r\n",delta_t);

    in_space dummy;
    observer location;
    cat_entry dummy_star;
    object bodyObject;

    make_cat_entry("DUMMY", "xxx", 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &dummy_star);
    make_object(0, body, "", &dummy_star, &bodyObject);
    make_observer (1,geo_loc,&dummy, &location);

    short int coord_sys = 2;        /* True equator and equinox of date */

    sky_pos radec;
    place (jd_tt,&bodyObject,&location,delta_t,coord_sys,0, &radec);

    double zd, az, rar, decr;
    equ2hor(jd_ut1, delta_t, 0, 0.379409,0.033159, geo_loc, radec.ra, radec.dec, 0, &zd, &az, &rar, &decr);
    //equ2hor(jd_ut1, delta_t, 0, 0.033159, 0.379409, geo_loc, radec.ra, radec.dec, 0, &zd, &az, &rar, &decr);

    printf("RA:%15.10f  DEC:%15.10f  R:%15.12f  ALT:%15.10f  AZ:%15.10f\n", rar*15, decr, radec.dis,90-zd,az);
}

int main(int argc, char *args[]){
    double jd_utc = julian_date(2020, 3, 14, 0);
    //double jd_utc = julian_date(2020, 3, 14, 16.72277778);
    const short int leap_secs = 37;
    const double ut1_utc = -0.2181652;

    on_surface geo_loc;
    make_on_surface(38.2463888888, -85.763611111111, 0.0 , 10.0, 1010.0, &geo_loc);

    short int de_num = 0;
    double jd_beg, jd_end;
    ephem_open("JPLEPH.421", &jd_beg, &jd_end, &de_num);

    //test(jd_utc,1,&geo_loc,leap_secs,ut1_utc);
    printf("\r\n");
    test(jd_utc,11,&geo_loc,leap_secs,ut1_utc);
    test(jd_utc, 4,&geo_loc,leap_secs,ut1_utc);

    ephem_close();
}
