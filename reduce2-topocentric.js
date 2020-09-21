class Reduce2{
	//By Greg Miller (gmiller@gregmiller.net)
	//Released as public domain

	static reduce(body,jd_tdb,observer){
		const earth=this.getBodyPV(2,jd_tdb);
		const target=this.getBodyPV(body,jd_tdb);
		const observerGeocentric=this.getObserverGeocentricPosition(observer,jd_tdb);
		
		const geocentricTarget=[target[0]-earth[0],target[1]-earth[1],target[2]-earth[2]];
		const topocentricTarget=[geocentricTarget[0]-observerGeocentric[0],geocentricTarget[1]-observerGeocentric[1],geocentricTarget[2]-observerGeocentric[2]];
		
		const radec=this.xyzToRaDec(topocentricTarget);
		const altaz=this.raDecToAltAz(radec[0],radec[1],observer[0],observer[1],jd_tdb);
		
		return [radec[0],radec[1],radec[0],radec[1],altaz[0],altaz[1]];
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
}
