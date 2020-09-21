class Reduce6{
	//By Greg Miller (gmiller@gregmiller.net)
	//Released as public domain

	static reduce(body,jd_utc,observer){
		const jd_tt=this.convertUTCtoTT(jd_utc);
		const jd_tdb=jd_tt;

		
		const earth=this.getBodyPV(2,jd_tdb);
		const target=this.getBodyLightAdjusted(earth,body,jd_tdb);

		const precessionMatrix=this.getPrecessionMatrix(jd_tdb);
		const nutationMatrix=this.getNutationMatrix(jd_tdb);

		const geocentricTarget=[target[0]-earth[0],target[1]-earth[1],target[2]-earth[2]];

		let observerGeocentric=this.getObserverGeocentricPosition(observer,jd_utc);
		observerGeocentric=Vec.vecMatrixMul(observerGeocentric,Vec.transpose(precessionMatrix));
		observerGeocentric=Vec.vecMatrixMul(observerGeocentric,Vec.transpose(nutationMatrix));

		const topocentricTarget=[geocentricTarget[0]-observerGeocentric[0],geocentricTarget[1]-observerGeocentric[1],geocentricTarget[2]-observerGeocentric[2]];

		const radecj2000=this.xyzToRaDec(topocentricTarget);

		const precessed=Vec.vecMatrixMul(topocentricTarget,precessionMatrix);
		const nutated=Vec.vecMatrixMul(precessed,nutationMatrix);
		
		const radec=this.xyzToRaDec(nutated);
		const altaz=this.raDecToAltAz(radec[0],radec[1],observer[0],observer[1],jd_utc);
		
		return [radecj2000[0],radecj2000[1],radec[0],radec[1],altaz[0],altaz[1]];
	}

	static getPrecessionMatrix(jd_tdb){
		//Fukushima-Williams IAU 2006
		const t=(jd_tdb-2451545.5)/36525.0;

		const gamma = (-0.052928 + 10.556378*t + 0.4932044*t*t - 0.00031238*t*t*t - 0.000002788*t*t*t*t + 0.0000000260*t*t*t*t*t) /60/60*Math.PI/180;
		const phi = (+84381.412819 - 46.811016*t + 0.0511268*t*t + 0.00053289*t*t*t - 0.000000440*t*t*t*t - 0.0000000176*t*t*t*t*t) /60/60*Math.PI/180;
		const psi = (-0.041775 + 5038.481484*t + 1.5584175*t*t - 0.00018522*t*t*t - 0.000026452*t*t*t*t - 0.0000000148*t*t*t*t*t) /60/60*Math.PI/180;
		const eps = (+84381.406 - 46.836769*t - 0.0001831*t*t + 0.00200340*t*t*t - 0.000000576*t*t*t*t - 0.0000000434*t*t*t*t*t) /60/60*Math.PI/180;

		const a=Vec.getZRotationMatrix(gamma);
		const b=Vec.dot(Vec.getXRotationMatrix(phi),a);
		const c=Vec.dot(Vec.getZRotationMatrix(-psi),b);
		const d=Vec.dot(Vec.getXRotationMatrix(-eps),c);

		return d;
	}

	static getLeapSeconds(jd){
		//Source IERS Resolution B1 and http://maia.usno.navy.mil/ser7/tai-utc.dat
		//This function must be updated any time a new leap second is introduced

		if(jd > 2457754.5) return 37.0;
		if(jd > 2457204.5) return 36.0;
		if(jd > 2456109.5) return 35.0;
		if(jd > 2454832.5) return 34.0;
		if(jd > 2453736.5) return 33.0;
		if(jd > 2451179.5) return 32.0;
		if(jd > 2450630.5) return 31.0;
		if(jd > 2450083.5) return 30.0;
		if(jd > 2449534.5) return 29.0;
		if(jd > 2449169.5) return 28.0;
		if(jd > 2448804.5) return 27.0;
		if(jd > 2448257.5) return 26.0;
		if(jd > 2447892.5) return 25.0;
		if(jd > 2447161.5) return 24.0;
		if(jd > 2446247.5) return 23.0;
		if(jd > 2445516.5) return 22.0;
		if(jd > 2445151.5) return 21.0;
		if(jd > 2444786.5) return 20.0;
		if(jd > 2444239.5) return 19.0;
		if(jd > 2443874.5) return 18.0;
		if(jd > 2443509.5) return 17.0;
		if(jd > 2443144.5) return 16.0;
		if(jd > 2442778.5) return 15.0;
		if(jd > 2442413.5) return 14.0;
		if(jd > 2442048.5) return 13.0;
		if(jd > 2441683.5) return 12.0;
		if(jd > 2441499.5) return 11.0;
		if(jd > 2441317.5) return 10.0;
		if(jd > 2439887.5) return 4.2131700 + (jd - 2439126.5) * 0.002592;
		if(jd > 2439126.5) return 4.3131700 + (jd - 2439126.5) * 0.002592;
		if(jd > 2439004.5) return 3.8401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438942.5) return 3.7401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438820.5) return 3.6401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438761.5) return 3.5401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438639.5) return 3.4401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438486.5) return 3.3401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438395.5) return 3.2401300 + (jd - 2438761.5) * 0.001296;
		if(jd > 2438334.5) return 1.9458580 + (jd - 2437665.5) * 0.0011232;
		if(jd > 2437665.5) return 1.8458580 + (jd - 2437665.5) * 0.0011232;
		if(jd > 2437512.5) return 1.3728180 + (jd - 2437300.5) * 0.001296;
		if(jd > 2437300.5) return 1.4228180 + (jd - 2437300.5) * 0.001296;
		return 0.0;
	}

	static convertUTCtoTT(jd_utc){
		const leapSeconds=this.getLeapSeconds(jd_utc);
		const TAI=this.getLeapSeconds(jd_utc)/60/60/24 + jd_utc;
		const TT=TAI + 32.184/60/60/24;

		return TT;
	}

	static getBodyLightAdjusted(origin,body,jd){
		
		let jd_light=jd;
		let b;
		for(let i=0;i<3;i++){
			b=this.getBodyPV(body,jd_light);
			const r=Math.sqrt((origin[0]-b[0])*(origin[0]-b[0])+(origin[1]-b[1])*(origin[1]-b[1])+(origin[2]-b[2])*(origin[2]-b[2]));
			const lightTime=r/(c/au*60*60*24);
			jd_light=jd-lightTime;
		}
		return b;
	}

	static getObserverGeocentricPosition(observer,jd_ut){
		const itrf=this.convertGeodedicLatLonToITRFXYZ(observer[0],observer[1],observer[2]);
		const gmst=this.greenwichMeanSiderealTime(jd_ut);

		const m=Vec.getZRotationMatrix(-gmst);
		const gcrs=Vec.vecMatrixMul(itrf,m);

		return gcrs;
	}

	//Convert Geodedic Lat Lon to geocentric XYZ position vector in ITRF coordinates
	static convertGeodedicLatLonToITRFXYZ(lat,lon,height){
		//Algorithm from Explanatory Supplement to the Astronomical Almanac 3rd ed. P294
		const a=6378136.6/au;
		const f=1/298.25642;

		const C=Math.sqrt(((Math.cos(lat)*Math.cos(lat)) + (1.0-f)*(1.0-f) * (Math.sin(lat)*Math.sin(lat))));

		const S=(1-f)*(1-f)*C;
		
		const h=height;

		let r=new Array();
		r[0]=(a*C+h) * Math.cos(lat) * Math.cos(lon);
		r[1]=(a*C+h) * Math.cos(lat) * Math.sin(lon);
		r[2]=(a*S+h) * Math.sin(lat);
		
		return r;
	}

	static xyzToRaDec(target){
		const x=target[0];
		const y=target[1];
		const z=target[2];

		//Convert from Cartesian to polar coordinates 
		const r=Math.sqrt(x*x+y*y+z*z);
		let l=Math.atan2(y,x);
		let t=Math.acos(z/r);

		//Make sure RA is positive, and Dec is in range +/-90
		if(l<0){l+=2*Math.PI;}
		t=.5*Math.PI-t;

		return [l,t,r];
	}

	static raDecToAltAz(ra,dec,lat,lon,jd_ut){
		//based on Explanatory supplement eq 7.16
		const gmst=this.greenwichMeanSiderealTime(jd_ut);
		const localSiderealTime=gmst+lon;
		 
		const H=localSiderealTime - ra;

		const a=Math.asin(Math.sin(dec)*Math.sin(lat)+Math.cos(dec)*Math.cos(H)*Math.cos(lat));
		const az=Math.asin( (-(Math.cos(dec)*Math.sin(H)))/Math.cos(a)  );
		return [Math.PI-az,a];
	}

	static greenwichMeanSiderealTime(jd_ut){
		//"The New Defintion of Universal Time" S. Aoki et al --- Equations 13 and 19
		//http://adsabs.harvard.edu/full/1982A%26A...105..359A
		const du=Math.floor(jd_ut-2451545)-.5;
		const frac=(jd_ut-2451545-du)*60*60*24;
		const tu=du/36525.0;

		let gmst=24110.54841 + 8640184.812866*tu + 0.093104*tu*tu - 6.210e-6*tu*tu*tu; //eq 13
		gmst+=frac * 1.002737909350795; //eq 19 simplified
		//gmst+=frac* 1.002737909350795 + 5.900610e-11*T - 5.910e-15*T*T; //eq 19 Aoki et al.

		//gmst is in seconds, convert to radians, then mod 2pi radians to keep it in a 24hour range
		gmst=(gmst/60/60*15*Math.PI/180.0)%(2.0*Math.PI);
		if(gmst<0){gmst+=2*Math.PI;}

		return gmst;
	}

	static getBodyPV(body,jd_tdb){
		let b;
		const AU_KM = 1.4959787069098932e+8;
		
		switch (body){
			case 2: //Earth
				b=de.getEarth(jd_tdb);
				break;
			case 9: //Moon
				const e=de.getEarth(jd_tdb);
				b=de.getAllPropertiesForSeries(body,jd_tdb);
				for(let i=0;i<e.length;i++){
					b[i]=b[i]+e[i];
				}
				break;
			default:
				b=de.getAllPropertiesForSeries(body,jd_tdb);
				break;
		}
		for(let i=0;i<b.length;i++){
			b[i]=b[i]/AU_KM;
		}
		return b;

	}

	static getNutationMatrix(jd_tdb){
		const t=(jd_tdb-2451545.5)/36525.0;
		const nut=this.nutation(t);
		const dpsi=nut[0]
		const deps=nut[1];
		const eps = (+84381.406 - 46.836769*t - 0.0001831*t*t + 0.00200340*t*t*t - 0.000000576*t*t*t*t - 0.0000000434*t*t*t*t*t) /60/60*Math.PI/180;

		const a=Vec.getXRotationMatrix(eps);
		const b=Vec.dot(Vec.getZRotationMatrix(-dpsi),a);
		const c=Vec.dot(Vec.getXRotationMatrix(-(eps+deps)),b);

		return c;
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
