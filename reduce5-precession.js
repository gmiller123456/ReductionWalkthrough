class Reduce5{
	//By Greg Miller (gmiller@gregmiller.net)
	//Released as public domain

	static reduce(body,jd_utc,observer){
		const jd_tt=this.convertUTCtoTT(jd_utc);
		const jd_tdb=jd_tt;
		
		const earth=this.getBodyPV(2,jd_tdb);
		const target=this.getBodyLightAdjusted(earth,body,jd_tdb);

		const precessionMatrix=this.getPrecessionMatrix(jd_tdb);

		const geocentricTarget=[target[0]-earth[0],target[1]-earth[1],target[2]-earth[2]];

		let observerGeocentric=this.getObserverGeocentricPosition(observer,jd_tdb);
		observerGeocentric=Vec.vecMatrixMul(observerGeocentric,Vec.transpose(precessionMatrix));

		const topocentricTarget=[geocentricTarget[0]-observerGeocentric[0],geocentricTarget[1]-observerGeocentric[1],geocentricTarget[2]-observerGeocentric[2]];

		const radecj2000=this.xyzToRaDec(topocentricTarget);

		const precessed=Vec.vecMatrixMul(topocentricTarget,precessionMatrix);
		
		const radec=this.xyzToRaDec(precessed);
		const altaz=this.raDecToAltAz(radec[0],radec[1],observer[0],observer[1],jd_tdb);
		
		return [radecj2000[0],radecj2000[1],radec[0],radec[1],altaz[0],altaz[1]];
	}

	static getPrecessionMatrix(jd_tbd){
		//Fukushima-Williams IAU 2006
		const t=(jd_tbd-2451545.5)/36525.0;

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

	//Input and output in radians
	static raDecToAltAz(ra,dec,lat,lon,jd_ut){
		//based on Meeus eq 13.5,13.6
		const gmst=this.greenwichMeanSiderealTime(jd_ut);
		const localSiderealTime=gmst+lon;
		const H=localSiderealTime - ra;
		let A=Math.atan2(Math.sin(H), Math.cos(H)*Math.sin(lat) - Math.tan(dec)*Math.cos(lat))
		const h=Math.asin(Math.sin(lat)*Math.sin(dec) + Math.cos(lat)*Math.cos(dec)*Math.cos(H));
		
		A=(A+Math.PI)%(2*Math.PI);
		
		return [A,h];
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
}
