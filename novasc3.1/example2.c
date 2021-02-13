#include <stdio.h>
#include "eph_manager.h"
#include "novas.h"

void test(double jd_utc, short body, on_surface* geo_loc, short int leap_secs, double ut1_utc)
{
    cat_entry dummy_star;
    object bodyObject;

    make_cat_entry("DUMMY", "xxx", 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &dummy_star);
    make_object(0, body, "", &dummy_star, &bodyObject);

    double jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
    double jd_ut1 = jd_utc + ut1_utc / 86400.0;
    double delta_t = 32.184 + leap_secs - ut1_utc;

    double ra, dec, dis, zd, az, rar, decr;
    topo_planet(jd_tt, &bodyObject, delta_t, geo_loc, 0, &ra, &dec, &dis);
    equ2hor(jd_ut1, delta_t, 0, 0.033159, 0.379409, geo_loc, ra, dec, 0, &zd, &az, &rar, &decr);

    printf("RA:%15.10f  DEC:%15.10f  R:%15.12f  ALT:%15.10f  AZ:%15.10f\n", rar*15, decr, dis,90-zd,az);
}

int main(int argc, char *args[]){
    double jd_utc = julian_date(2020, 3, 14, 16.72277778);
    const short int leap_secs = 37;
    const double ut1_utc = -0.2181652;

    on_surface geo_loc;
    make_on_surface(38.2464000, -85.7636, 0.0 , 10.0, 1010.0, &geo_loc);

    short int de_num = 0;
    double jd_beg, jd_end;
    ephem_open("JPLEPH", &jd_beg, &jd_end, &de_num);

    test(jd_utc,1,&geo_loc,leap_secs,ut1_utc);
    //test(jd_utc, 4,&geo_loc,leap_secs,ut1_utc);

    ephem_close();
}
