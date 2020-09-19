class Reduce1{
	static reduce(body){
		const jd_tdb=2451537.5;
		const earth=Reduce1.getBodyPV(2,jd_tdb);
		console.log("Earth:",earth);
		const target=Reduce1.getBodyLightAdjusted(earth,body,jd_tdb);
		console.log("Target: ",target);
	}

	static getBodyLightAdjusted(origin,body,jd){
		
		let jd_light=jd;
		let b;
		for(let i=0;i<3;i++){
			b=Reduce1.getBodyPV(body,jd_light);
			const r=Math.sqrt((origin[0]-b[0])*(origin[0]-b[0])+(origin[1]-b[1])*(origin[1]-b[1])+(origin[2]-b[2])*(origin[2]-b[2]));
			const lightTime=r/(c/au*60*60*24);
			jd_light=jd_light-lightTime;
		}
		return b;
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
