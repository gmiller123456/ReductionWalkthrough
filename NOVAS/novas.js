/*
Partial Javascript implementation of the USNO's public domain NOVAS library
Greg Miller (gmiller@gregmiller.net) http://www.astrogreg.com
Released as public domain
*/

const au= 149597870.691;
const T0 = 2451545.00000000;
const ASEC360 = 1296000.0;
const ASEC2RAD = 4.848136811095359935899141e-6;
const TWOPI=Math.PI*2;
const DEG2RAD=Math.PI/180.0;
const ERAD = 6378136.6;
const F = 0.003352819697896;
const ANGVEL = 7.2921150e-5;
const AU_KM = 1.4959787069098932e+8;
const C_AUDAY = 173.1446326846693;
const AU = 1.4959787069098932e+11;
const RAD2DEG = 57.295779513082321;

function ter2cel ( jd_ut_high,  jd_ut_low,  delta_t,
                    method, option,  xp,  yp,  vec1,  vec2)
{
   var error = 0;
   var rs, j;

   var jd_ut1, jd_tt, dummy, secdiff, jd_tdb, gast, r_cio, theta;
   

 var v1=new Array();
 var v2=new Array();
 var v3=new Array();
 var v4=new Array();
 var x=new Array();
 var y=new Array();
 var z=new Array();

/*
   Compute the TT Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_ut1 = jd_ut_high + jd_ut_low;
   jd_tt = jd_ut1 + (delta_t / 86400.0);

/*
   Compute the TDB Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_tdb = jd_tt;
   secdiff=tdb2tt (jd_tdb);
   jd_tdb = jd_tt + secdiff / 86400.0;
   switch (method)
   {
      case (0):
      	//TODO: CIO method
         break;

      case (1):

/*
   Equinox mode.

   Apply polar motion.
*/

         if ((xp == 0.0) && (yp == 0.0))
         {
            for (j = 0; j < 3; j++)
            {
               v1[j] = vec1[j];
            }
         }
          else
            wobble (jd_tdb,0,xp,yp,vec1, v1);  //TODO: implement wobble

/*
   Apply Earth rotation.
*/

         gast=sidereal_time (jd_ut_high,delta_t,1,1);
         spin (-gast * 15.0,v1, v2);
/*
   'option' = 1 skips remaining transformations.
*/

         if (option == 1)
         {
            vec2[0] = v2[0];
            vec2[1] = v2[1];
            vec2[2] = v2[2];
         }
          else
         {

/*
   Apply precession, nutation, and frame tie.
*/

            //nutation (jd_tdb,-1,accuracy,v2, v3);
            //precession (jd_tdb,v3,T0, v4);
            //frame_tie (v4,-1, vec2);
            v3=nutation (jd_tdb,-1,accuracy,v2);
            v4=precession (jd_tdb,v3,T0);
            vec2=frame_tie (v4,-1);
         }
         break;

/*
   Invalid value of 'method'.
*/

      default:
         error = 2;
         break;
   }
   return (error);
}

function spin (angle, pos1, pos2)
{
   var xx, yx, zx, xy, yy, zy, xz, yz, zz;
   var angr, cosang, sinang;

      angr = angle * DEG2RAD;
      cosang = Math.cos (angr);
      sinang = Math.sin (angr);

/*
   Rotation matrix follows.
*/

      xx =  cosang;
      yx =  sinang;
      zx =  0.0;
      xy =  -sinang;
      yy =  cosang;
      zy =  0.0;
      xz =  0.0;
      yz =  0.0;
      zz =  1.0;


/*
   Perform rotation.
*/

   pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
   pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
   pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];

   return;
}


function tdb2tt (tdb_jd)

{
	var tt_jd, secdiff;
   var t;

   t = (tdb_jd - T0) / 36525.0;

/*
   Expression given in USNO Circular 179, eq. 2.6.
*/

   secdiff = 0.001657 * Math.sin ( 628.3076 * t + 6.2401)
            + 0.000022 * Math.sin ( 575.3385 * t + 4.2970)
            + 0.000014 * Math.sin (1256.6152 * t + 6.1969)
            + 0.000005 * Math.sin ( 606.9777 * t + 4.0212)
            + 0.000005 * Math.sin (  52.9691 * t + 0.4444)
            + 0.000002 * Math.sin (  21.3299 * t + 5.5431)
            + 0.000010 * t * Math.sin ( 628.3076 * t + 4.2490);

   tt_jd = tdb_jd - secdiff / 86400.0;

    return secdiff;
}

function equ2hor ( jd_ut1,  delta_t,  xp,  yp, location,  ra,  dec, ref_option)
{
   var j;

   var sinlat, coslat, sinlon, coslon, sindc, cosdc, sinra, cosra,
      uze, une, uwe, uz, un, uw, p, pz, pn, pw,
      proj, zd0, zd1, refr, sinzd, coszd, sinzd0, coszd0, pr;

   //Return values:
   var zd, az, rar, decr;

	uze=new Array();
	une=new Array();
	uwe=new Array();
	uz=new Array();
	un=new Array();
	uw=new Array();
	p=new Array();
	pr=new Array();

/*
   Preliminaries.
*/

   rar = ra;
   decr = dec;

   sinlat = Math.sin (location.latitude * DEG2RAD);
   coslat = Math.cos (location.latitude * DEG2RAD);
   sinlon = Math.sin (location.longitude * DEG2RAD);
   coslon = Math.cos (location.longitude * DEG2RAD);
   sindc = Math.sin (dec * DEG2RAD);
   cosdc = Math.cos (dec * DEG2RAD);
   sinra = Math.sin (ra * 15.0 * DEG2RAD);
   cosra = Math.cos (ra * 15.0 * DEG2RAD);

/*
   Set up orthonormal basis vectors in local Earth-fixed system.

   Define vector toward local zenith in Earth-fixed system (z axis).
*/
   uze[0] = coslat * coslon;
   uze[1] = coslat * sinlon;
   uze[2] = sinlat;

/*
   Define vector toward local north in Earth-fixed system (x axis).
*/

   une[0] = -sinlat * coslon;
   une[1] = -sinlat * sinlon;
   une[2] = coslat;

/*
   Define vector toward local west in Earth-fixed system (y axis).
*/

   uwe[0] = sinlon;
   uwe[1] = -coslon;
   uwe[2] = 0.0;

/*
   Obtain vectors in celestial system.

   Rotate Earth-fixed orthonormal basis vectors to celestial system
   (wrt equator and equinox of date).
*/

   ter2cel (jd_ut1,0.0,delta_t,1,1,xp,yp,uze, uz);
   ter2cel (jd_ut1,0.0,delta_t,1,1,xp,yp,une, un);
   ter2cel (jd_ut1,0.0,delta_t,1,1,xp,yp,uwe, uw);

/*
   Define unit vector 'p' toward object in celestial system
   (wrt equator and equinox of date).
*/

   p[0] = cosdc * cosra;
   p[1] = cosdc * sinra;
   p[2] = sindc;

/*
   Compute coordinates of object wrt orthonormal basis.

   Compute components of 'p' - projections of 'p' onto rotated
   Earth-fixed basis vectors.
*/

   pz = p[0] * uz[0] + p[1] * uz[1] + p[2] * uz[2];
   pn = p[0] * un[0] + p[1] * un[1] + p[2] * un[2];
   pw = p[0] * uw[0] + p[1] * uw[1] + p[2] * uw[2];

/*
   Compute azimuth and zenith distance.
*/

   proj = Math.sqrt (pn * pn + pw * pw);

   if (proj > 0.0)
      az = -Math.atan2 (pw, pn) * RAD2DEG;

   if (az < 0.0)
      az += 360.0;

   if (az >= 360.0)
      az -= 360.0;

   zd = Math.atan2 (proj, pz) * RAD2DEG;

/*
   Apply atmospheric refraction if requested.
*/
	//TODO: refraction

	var ret=new Array();
	ret[0]=zd;
	ret[1]=az;

	return ret;
}



function vector2radec (pos)
{
   let xyproj;
   let ra,dec;

   xyproj = Math.sqrt (pos[0] * pos[0] + pos[1] * pos[1]);
   if ((xyproj == 0.0) && (pos[2] == 0))
   {
      ra = 0.0;
      dec = 0.0;
      return 1;
   }
    else if (xyproj == 0.0)
   {
      ra = 0.0;
      if (pos[2] < 0.0)
         dec = -90.0;
       else
         dec = 90.0;
      return 2;
   }
    else
   {
      ra = Math.atan2 (pos[1], pos[0]) / ASEC2RAD / 54000.0;
      dec = Math.atan2 (pos[2], xyproj) / ASEC2RAD / 3600.0;

      if (ra < 0.0)
         ra += 24.0;
   }

   const temp=new Array();
   temp[0]=ra;
   temp[1]=dec;
   return temp;
}


function aberration (pos, ve, lighttime)
{
   var p1mag, vemag, beta, dot, cosd, gammai, p, q, r;

   if (lighttime == 0.0)
   {
      p1mag = Math.sqrt (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] *
         pos[2]);
      lighttime = p1mag / C_AUDAY;
   }
    else
      p1mag = lighttime * C_AUDAY;

   vemag = Math.sqrt (ve[0] * ve[0] + ve[1] * ve[1] + ve[2] * ve[2]);
   beta = vemag / C_AUDAY;
   dot = pos[0] * ve[0] + pos[1] * ve[1] + pos[2] * ve[2];

   cosd = dot / (p1mag * vemag);
   gammai = Math.sqrt (1.0 - beta * beta);

   p = beta * cosd;
   q = (1.0 + p / (1.0 + gammai)) * lighttime;
   r = 1.0 + p;

   const pos2=new Array();
   pos2[0] = (gammai * pos[0] + q * ve[0]) / r;
   pos2[1] = (gammai * pos[1] + q * ve[1]) / r;
   pos2[2] = (gammai * pos[2] + q * ve[2]) / r;

   return pos2;
}
function grav_def(jd_tdb,loc_code,pos1,pos_obs){
	//TODO:  Actually implement
	return pos1;
}

function limb_angle (pos_obj, pos_obs)
{

   var pi, halfpi, rade;
   var disobj, disobs, aprad, zdlim, coszd, zdobj;

      pi = TWOPI / 2.0;
      halfpi = pi / 2.0;
      rade = ERAD / AU;

/*
   Compute the distance to the object and the distance to the observer.
*/

   disobj = Math.sqrt (pos_obj[0] * pos_obj[0] +
                  pos_obj[1] * pos_obj[1] +
                  pos_obj[2] * pos_obj[2]);

   disobs = Math.sqrt (pos_obs[0] * pos_obs[0] +
                  pos_obs[1] * pos_obs[1] +
                  pos_obs[2] * pos_obs[2]);

/*
   Compute apparent angular radius of Earth's limb.
*/

   if (disobs >= rade)
   {
      aprad = Math.asin (rade / disobs);
   }
    else
   {
      aprad = halfpi;
   }

/*
   Compute zenith distance of Earth's limb.
*/

   zdlim = pi - aprad;

/*
   Compute zenith distance of observed object.
*/

   coszd = (pos_obj[0] * pos_obs[0] + pos_obj[1] * pos_obs[1] +
      pos_obj[2] * pos_obs[2]) / (disobj * disobs);

   if (coszd <= -1.0)
   {
      zdobj = pi;
   }
    else if (coszd >= 1.0)
   {
      zdobj = 0.0;
   }
    else
   {
      zdobj = Math.acos (coszd);
   }


   const temp=new Array();
/*
   Angle of object wrt limb is difference in zenith distances.
*/
	temp[0]=(zdlim - zdobj) * RAD2DEG;
   //*limb_ang = (zdlim - zdobj) * RAD2DEG;

/*
   Nadir angle of object as a fraction of angular radius of limb.
*/

	temp[1]=(pi - zdobj) / aprad;
   //*nadir_ang = (pi - zdobj) / aprad;

   return temp;
}


function geo_posvel(jd_tt,delta_t,obs){
	const jd_tdb=tt2tdb(jd_tt);
	const jd_ut1 = jd_tt - (delta_t / 86400.0);

	const gmst=sidereal_time (jd_ut1,delta_t,0,1);
	const eqeq=e_tilt(jd_tdb)[2];
	const gast=gmst+eqeq/3600.0;
	
	const pv1=terra(obs,gast);

   const pos2=nutation (jd_tdb,-1,pv1);
   const pos3=precession (jd_tdb,pos2,T0);
   const pos4=frame_tie (pos3,-1);

   const vel1=new Array();
   vel1[0]=pv1[3];
   vel1[1]=pv1[4];
   vel1[2]=pv1[5];

   const vel2=nutation (jd_tdb,-1,vel1);
   const vel3=precession (jd_tdb,vel2,T0);
   const vel4=frame_tie (vel3,-1);

   const temp=new Array();
   temp[0]=pos4[0];
   temp[1]=pos4[1];
   temp[2]=pos4[2];
   temp[3]=vel4[0];
   temp[4]=vel4[1];
   temp[5]=vel4[2];

   return temp;

}

function frame_tie(pos1,direction){
	const pos2=new Array();
/*
   'xi0', 'eta0', and 'da0' are ICRS frame biases in arcseconds taken
   from IERS (2003) Conventions, Chapter 5.
*/

   const xi0  = -0.0166170;
   const eta0 = -0.0068192;
   const da0  = -0.01460;
   var  xx, yx, zx, xy, yy, zy, xz, yz, zz;

/*
   Compute elements of rotation matrix to first order the first time
   this function is called.  Elements will be saved for future use and
   not recomputed.
*/

      xx =  1.0;
      yx = -da0  * ASEC2RAD;
      zx =  xi0  * ASEC2RAD;
      xy =  da0  * ASEC2RAD;
      yy =  1.0;
      zy =  eta0 * ASEC2RAD;
      xz = -xi0  * ASEC2RAD;
      yz = -eta0 * ASEC2RAD;
      zz =  1.0;

/*
   Include second-order corrections to diagonal elements.
*/

      xx = 1.0 - 0.5 * (yx * yx + zx * zx);
      yy = 1.0 - 0.5 * (yx * yx + zy * zy);
      zz = 1.0 - 0.5 * (zy * zy + zx * zx);


/*
   Perform the rotation in the sense specified by 'direction'.
*/

   if (direction < 0)
   {

/*
   Perform rotation from dynamical system to ICRS.
*/

      pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
      pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
      pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
   }
    else
   {

/*
   Perform rotation from ICRS to dynamical system.
*/

      pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
      pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
      pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
   }



	return pos2;
}

function precession(jd_tdb1, pos1, jd_tdb2){
	const pos2=new Array();

   var xx, yx, zx, xy, yy, zy, xz, yz, zz;
   var eps0 = 84381.406;
   var  t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;

/*
   Check to be sure that either 'jd_tdb1' or 'jd_tdb2' is equal to T0.
*/

   if ((jd_tdb1 != T0) && (jd_tdb2 != T0))
      return (error = 1);

/*
   't' is time in TDB centuries between the two epochs.
*/

   t = (jd_tdb2 - jd_tdb1) / 36525.0;

   if (jd_tdb2 == T0)
      t = -t;


/*
   Numerical coefficients of psi_a, omega_a, and chi_a, along with
   epsilon_0, the obliquity at J2000.0, are 4-angle formulation from
   Capitaine et al. (2003), eqs. (4), (37), & (39).
*/

      psia   = ((((-    0.0000000951  * t
                   +    0.000132851 ) * t
                   -    0.00114045  ) * t
                   -    1.0790069   ) * t
                   + 5038.481507    ) * t;

      omegaa = ((((+    0.0000003337  * t
                   -    0.000000467 ) * t
                   -    0.00772503  ) * t
                   +    0.0512623   ) * t
                   -    0.025754    ) * t + eps0;

      chia   = ((((-    0.0000000560  * t
                   +    0.000170663 ) * t
                   -    0.00121197  ) * t
                   -    2.3814292   ) * t
                   +   10.556403    ) * t;

      eps0 = eps0 * ASEC2RAD;
      psia = psia * ASEC2RAD;
      omegaa = omegaa * ASEC2RAD;
      chia = chia * ASEC2RAD;

      sa = Math.sin (eps0);
      ca = Math.cos (eps0);
      sb = Math.sin (-psia);
      cb = Math.cos (-psia);
      sc = Math.sin (-omegaa);
      cc = Math.cos (-omegaa);
      sd = Math.sin (chia);
      cd = Math.cos (chia);
/*
   Compute elements of precession rotation matrix equivalent to
   R3(chi_a) R1(-omega_a) R3(-psi_a) R1(epsilon_0).
*/

      xx =  cd * cb - sb * sd * cc;
      yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;
      zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
      xy = -sd * cb - sb * cd * cc;
      yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;
      zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
      xz =  sb * sc;
      yz = -sc * cb * ca - sa * cc;
      zz = -sc * cb * sa + cc * ca;


   if (jd_tdb2 == T0)
   {

/*
   Perform rotation from epoch to J2000.0.
*/
      pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
      pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
      pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
   }
    else
   {

/*
   Perform rotation from J2000.0 to epoch.
*/

      pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
      pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
      pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
   }

	return pos2;
}

function nutation (jd_tdb, direction, pos){
	const pos2=new Array();
   var cobm, sobm, cobt, sobt, cpsi, spsi, xx, yx, zx, xy, yy, zy,
      xz, yz, zz, oblm, oblt, eqeq, psi, eps;

/*
   Call 'e_tilt' to get the obliquity and nutation angles.
*/

   const values=e_tilt (jd_tdb);

    oblm=values[0];
    oblt=values[1];
    eqeq=values[2];
    psi=values[3];
    eps=values[4];
   cobm = Math.cos (oblm * DEG2RAD);
   sobm = Math.sin (oblm * DEG2RAD);
   cobt = Math.cos (oblt * DEG2RAD);
   sobt = Math.sin (oblt * DEG2RAD);
   cpsi = Math.cos (psi * ASEC2RAD);
   spsi = Math.sin (psi * ASEC2RAD);

/*
   Nutation rotation matrix follows.
*/

   xx = cpsi;
   yx = -spsi * cobm;
   zx = -spsi * sobm;
   xy = spsi * cobt;
   yy = cpsi * cobm * cobt + sobm * sobt;
   zy = cpsi * sobm * cobt - cobm * sobt;
   xz = spsi * sobt;
   yz = cpsi * cobm * sobt - sobm * cobt;
   zz = cpsi * sobm * sobt + cobm * cobt;

   if (!direction)
   {

/*
   Perform rotation.
*/

      pos2[0] = xx * pos[0] + yx * pos[1] + zx * pos[2];
      pos2[1] = xy * pos[0] + yy * pos[1] + zy * pos[2];
      pos2[2] = xz * pos[0] + yz * pos[1] + zz * pos[2];
   }
    else
   {

/*
   Perform inverse rotation.
*/

      pos2[0] = xx * pos[0] + xy * pos[1] + xz * pos[2];
      pos2[1] = yx * pos[0] + yy * pos[1] + yz * pos[2];
      pos2[2] = zx * pos[0] + zy * pos[1] + zz * pos[2];
   }

	return pos2;
}

function terra(location,st){
   var j;

   var erad_km, ht_km;
   var df, df2, phi, sinphi, cosphi, c, s, ach, ash, stlocl, sinst, cosst;

  erad_km = ERAD / 1000.0;

/*
   Compute parameters relating to geodetic to geocentric conversion.
*/

   df = 1.0 - F;
   df2 = df * df;

   phi = location.latitude * DEG2RAD;
   sinphi = Math.sin (phi);
   cosphi = Math.cos (phi);
   c = 1.0 / Math.sqrt (cosphi * cosphi + df2 * sinphi * sinphi);
   s = df2 * c;
   ht_km = location.height / 1000.0;
   ach = erad_km * c + ht_km;
   ash = erad_km * s + ht_km;

/*
   Compute local sidereal time factors at the observer's longitude.
*/

   stlocl = (st * 15.0 + location.longitude) * DEG2RAD;
   sinst = Math.sin (stlocl);
   cosst = Math.cos (stlocl);

/*
   Compute position vector components in kilometers.
*/
	const pv=new Array();

   pv[0] = ach * cosphi * cosst;
   pv[1] = ach * cosphi * sinst;
   pv[2] = ash * sinphi;

/*
   Compute velocity vector components in kilometers/sec.
*/

   pv[3] = -ANGVEL * ach * cosphi * sinst;
   pv[4] =  ANGVEL * ach * cosphi * cosst;
   pv[5] =  0.0;

/*
   Convert position and velocity components to AU and AU/DAY.
*/

   for (j = 0; j < 3; j++)
   {
      pv[j] /= AU_KM;
      pv[j+3] /= AU_KM;
      pv[j+3] *= 86400.0;
   }

   return pv;

}

function e_tilt(jd_tdb){
	const t = (jd_tdb - T0) / 36525.0;

	const angles=nutation_angles(t);
	const c_terms = ee_ct (jd_tdb,0.0) / ASEC2RAD;

	//TODO: define PSI_COR and EPS_COR
   const d_psi=angles[0];
   const d_eps=angles[1];
   //d_psi = dp + PSI_COR;
   //d_eps = de + EPS_COR;

	let mean_ob = mean_obliq (jd_tdb);
	let true_ob = mean_ob + d_eps;

   mean_ob /= 3600.0;
   true_ob /= 3600.0;


   let eq_eq = d_psi * Math.cos (mean_ob * DEG2RAD) + c_terms;
   eq_eq /= 15.0;

   const temp=new Array();
   temp[0] = mean_ob;
   temp[1] = true_ob;
   temp[2] = eq_eq;
   temp[3] = d_psi;
   temp[4] = d_eps;

   return temp;
}

function mean_obliq (jd_tdb){
   const t = (jd_tdb - T0) / 36525.0;

/*
   Compute the mean obliquity in arcseconds.  Use expression from the
   reference's eq. (39) with obliquity at J2000.0 taken from eq. (37)
   or Table 8.
*/

   const epsilon = (((( -  0.0000000434   * t
                  -  0.000000576  ) * t
                  +  0.00200340   ) * t
                  -  0.0001831    ) * t
                  - 46.836769     ) * t + 84381.406;

   return (epsilon);

}

function norm_ang(angle){

   let a;

   a = fmod (angle,TWOPI);
   if (a < 0.0)
         a += TWOPI;

   return (a);

}

function ee_ct(jd_high,jd_low){
   var i, j;

   var t, fa, fa2, s0, s1, a, c_terms;
   fa=new Array();
   fa2=new Array();

/*
   Argument coefficients for t^0.
*/

   const ke0_t = [
     [0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0],
     [0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1],
     [0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0],
     [1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0]];

/*
   Argument coefficients for t^1.
*/

   const ke1 =
      [0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0];

/*
   Sine and cosine coefficients for t^0.
*/

   const se0_t = [
      [+2640.96e-6,          -0.39e-6],
      [  +63.52e-6,          -0.02e-6],
      [  +11.75e-6,          +0.01e-6],
      [  +11.21e-6,          +0.01e-6],
      [   -4.55e-6,          +0.00e-6],
      [   +2.02e-6,          +0.00e-6],
      [   +1.98e-6,          +0.00e-6],
      [   -1.72e-6,          +0.00e-6],
      [   -1.41e-6,          -0.01e-6],
      [   -1.26e-6,          -0.01e-6],
      [   -0.63e-6,          +0.00e-6],
      [   -0.63e-6,          +0.00e-6],
      [   +0.46e-6,          +0.00e-6],
      [   +0.45e-6,          +0.00e-6],
      [   +0.36e-6,          +0.00e-6],
      [   -0.24e-6,          -0.12e-6],
      [   +0.32e-6,          +0.00e-6],
      [   +0.28e-6,          +0.00e-6],
      [   +0.27e-6,          +0.00e-6],
      [   +0.26e-6,          +0.00e-6],
      [   -0.21e-6,          +0.00e-6],
      [   +0.19e-6,          +0.00e-6],
      [   +0.18e-6,          +0.00e-6],
      [   -0.10e-6,          +0.05e-6],
      [   +0.15e-6,          +0.00e-6],
      [   -0.14e-6,          +0.00e-6],
      [   +0.14e-6,          +0.00e-6],
      [   -0.14e-6,          +0.00e-6],
      [   +0.14e-6,          +0.00e-6],
      [   +0.13e-6,          +0.00e-6],
      [   -0.11e-6,          +0.00e-6],
      [   +0.11e-6,          +0.00e-6],
      [   +0.11e-6,          +0.00e-6]];
/*
   Sine and cosine coefficients for t^1.
*/

   const se1 =
      [   -0.87e-6,          +0.00e-6];

/*
   Interval between fundamental epoch J2000.0 and current date.
*/

      t = ((jd_high - T0) + jd_low) / 36525.0;

/*
   High accuracy mode.
*/


/*
   Fundamental Arguments.

   Mean Anomaly of the Moon.
*/

      fa[0] = norm_ang ((485868.249036 +
                         (715923.2178 +
                         (    31.8792 +
                         (     0.051635 +
                         (    -0.00024470)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (1325.0*t, 1.0) * TWOPI);

/*
   Mean Anomaly of the Sun.
*/

      fa[1] = norm_ang ((1287104.793048 +
                         (1292581.0481 +
                         (     -0.5532 +
                         (     +0.000136 +
                         (     -0.00001149)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (99.0*t, 1.0) * TWOPI);

/*
   Mean Longitude of the Moon minus Mean Longitude of the Ascending
   Node of the Moon.
*/

      fa[2] = norm_ang (( 335779.526232 +
                         ( 295262.8478 +
                         (    -12.7512 +
                         (     -0.001037 +
                         (      0.00000417)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (1342.0*t, 1.0) * TWOPI);

/*
   Mean Elongation of the Moon from the Sun.
*/

      fa[3] = norm_ang ((1072260.703692 +
                         (1105601.2090 +
                         (     -6.3706 +
                         (      0.006593 +
                         (     -0.00003169)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (1236.0*t, 1.0) * TWOPI);

/*
   Mean Longitude of the Ascending Node of the Moon.
*/

      fa[4] = norm_ang (( 450160.398036 +
                         (-482890.5431 +
                         (      7.4722 +
                         (      0.007702 +
                         (     -0.00005939)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (-5.0*t, 1.0) * TWOPI);

      fa[ 5] = norm_ang (4.402608842 + 2608.7903141574 * t);
      fa[ 6] = norm_ang (3.176146697 + 1021.3285546211 * t);
      fa[ 7] = norm_ang (1.753470314 +  628.3075849991 * t);
      fa[ 8] = norm_ang (6.203480913 +  334.0612426700 * t);
      fa[ 9] = norm_ang (0.599546497 +   52.9690962641 * t);
      fa[10] = norm_ang (0.874016757 +   21.3299104960 * t);
      fa[11] = norm_ang (5.481293872 +    7.4781598567 * t);
      fa[12] = norm_ang (5.311886287 +    3.8133035638 * t);
      fa[13] =          (0.024381750 +    0.00000538691 * t) * t;

/*
   Evaluate the complementary terms.
*/

      s0 = 0.0;
      s1 = 0.0;

      for (i = 32; i >= 0; i--)
      {
         a = 0.0;

         for (j = 0; j < 14; j++)
         {
            a += ke0_t[i][j] * fa[j];
         }

         s0 += (se0_t[i][0] * Math.sin (a) + se0_t[i][1] * Math.cos (a));
      }

      a = 0.0;

      for (j = 0; j < 14; j++)
      {
         a +=  (ke1[j]) * fa[j];
      }

      s1 += (se1[0] * Math.sin (a) + se1[1] * Math.cos (a));

      c_terms = (s0 + s1 * t);
      return c_terms*ASEC2RAD;

}

function nutation_angles (t){
	const t1 = t * 36525.0;

	const angles=iau2000a (T0,t1);

   angles[0] /= ASEC2RAD; //dpsi
   angles[1] /= ASEC2RAD; //deps

   return angles;
}


function sidereal_time(jd_ut,delta_t,gst_type,method){
	const jd_tt = jd_ut + (delta_t / 86400.0);
	const theta = era (jd_ut,0);
	const jd_tdb=tt2tdb(jd_tt);
	const t = (jd_tdb - T0) / 36525.0;

	let eqeq=0.0;
   if (((gst_type == 0) && (method == 0)) ||       /* GMST; CIO-TIO */
       ((gst_type == 1) && (method == 1)))         /* GAST; equinox */
   {
     var ee=e_tilt (jd_tdb)[2];
      eqeq = ee * 15.0;
   }

	switch(method){
		case 0:
			break;
		case 1:
		         const st = eqeq + 0.014506 +
               (((( -    0.0000000368   * t
                    -    0.000029956  ) * t
                    -    0.00000044   ) * t
                    +    1.3915817    ) * t
                    + 4612.156534     ) * t;

         var gst = fmod ((st / 3600.0 + theta), 360.0) / 15.0;

         if (gst < 0.0)
            gst += 24.0;
               return gst;
	}

	//TODO: ...
}

function era(jd_high,jd_low){
   var theta, thet1, thet2, thet3;

   thet1 = 0.7790572732640 + 0.00273781191135448 * (jd_high - T0);
   thet2 = 0.00273781191135448 * jd_low;
   thet3 = fmod (jd_high, 1.0) + fmod (jd_low, 1.0);

   theta = fmod (thet1 + thet2 + thet3, 1.0) * 360.0;
   if (theta < 0.0)
      theta += 360.0;

   return theta;	
}

function fund_args (t,a)
{

   a[0] = fmod (485868.249036 +
             t * (1717915923.2178 +
             t * (        31.8792 +
             t * (         0.051635 +
             t * (       - 0.00024470)))), ASEC360) * ASEC2RAD;

   a[1] = fmod (1287104.79305 +
             t * ( 129596581.0481 +
             t * (       - 0.5532 +
             t * (         0.000136 +
             t * (       - 0.00001149)))), ASEC360) * ASEC2RAD;

   a[2] = fmod (335779.526232 +
             t * (1739527262.8478 +
             t * (      - 12.7512 +
             t * (      -  0.001037 +
             t * (         0.00000417)))), ASEC360) * ASEC2RAD;

   a[3] = fmod (1072260.70369 +
             t * (1602961601.2090 +
             t * (       - 6.3706 +
             t * (         0.006593 +
             t * (       - 0.00003169)))), ASEC360) * ASEC2RAD;

   a[4] = fmod (450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939)))), ASEC360) * ASEC2RAD;

   return;
}

function fmod(a,b){
	return a%b;
}

function tt2tdb(jd_tt){
	/*From Novas*/
	const t = (jd_tt - T0) / 36525.0;

	const secdiff = 0.001657 * Math.sin ( 628.3076 * t + 6.2401)
		+ 0.000022 * Math.sin ( 575.3385 * t + 4.2970)
		+ 0.000014 * Math.sin (1256.6152 * t + 6.1969)
		+ 0.000005 * Math.sin ( 606.9777 * t + 4.0212)
		+ 0.000005 * Math.sin (  52.9691 * t + 0.4444)
		+ 0.000002 * Math.sin (  21.3299 * t + 5.5431)
		+ 0.000010 * t * Math.sin ( 628.3076 * t + 4.2490);

	return jd_tt+secdiff/ 86400.0;
}

