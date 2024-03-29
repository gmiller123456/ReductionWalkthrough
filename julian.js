 
//Special "Math.floor()" function used by dateToJulianDate()
function INT(d){
	if(d>0){
		return Math.floor(d);
	}
	return Math.floor(d)-1;
}

function gregorianDateToJulianDate(year, month, day, hour, min, sec){
	let isGregorian=true;
	if(year<1582 || (year == 1582 && (month < 10 || (month==10 && day < 5)))){
		isGregorian=false;
	}

	if (month < 3){
		year = year - 1;
		month = month + 12;
	}

	let b = 0;
	if (isGregorian){
	let a = INT(year / 100.0);
		b = 2 - a + INT(a / 4.0);
	}

	let jd=INT(365.25 * (year + 4716)) + INT(30.6001 * (month + 1)) + day + b - 1524.5;
	jd+=hour/24.0;
	jd+=min/24.0/60.0;
	jd+=sec/24.0/60.0/60.0;
	return jd;
}	

//From Meeus, CH7
function julainDateToGregorian(jd){
	let temp=jd+.5;
	let Z=Math.trunc(temp);
	let F=temp-Z;
	let A=Z;
	if(Z>=2299161){
		let alpha=INT((Z-1867216.25)/36524.25);
		A=Z+1+alpha-INT(alpha/4);
	}
	let B=A+1524;
	let C=INT((B-122.1)/365.25);
	let D=INT(365.25*C);
	let E=INT((B-D)/30.6001);

	let day=B-D-INT(30.6001*E)+F;
	let month=E-1;
	if(E>13){
		month=E-13;
	}
	let year=C-4716;
	if(month<3){
		year=C-4715;
	}

	return [year,month,day,0,0,0];

}
