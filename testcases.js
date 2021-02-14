const de=new DE405();
const au=149597870691; // meters
const c=299792458; // m/s
const toRad=Math.PI/180;  //Multiply by this to convert degrees to radians

const observer=[38.2464000*toRad,-85.7636*toRad,0]; //Louisville, Kentucky -85° 45' 48.96" 

function runTests(className){
	//2458923.196782407 Horizons value for 2020-03-14 16:43:22
	const jd=gregorianDateToJulianDate(2020,3,14,16,43,22);

	//Test data from JPL Horizons, date chosen at random
	/*                    JDUT R.A.___(ICRF)___DEC   R.A._(a-appar)_DEC.   Azi_(a-appr)_Elev   L_Ap_Sid_Time   1-way_down_LT       TDB-UT*/

	test(className,"Sun" , jd,[354.76044,  -2.26890, 355.01011,  -2.16102, 154.5692,  46.5898, 22.5222447106,    8.26968075,    69.185580]);
	test(className,"Mer" , jd,[330.91329, -11.42432, 331.17363, -11.32938, 188.5413,  40.0342, 22.5222447106,    6.39316934,    69.185580]);
	test(className,"Ven" , jd,[ 36.87090,  16.93354,  37.14179,  17.02072,  95.0213,  34.3682, 22.5222447106,    6.54428376,    69.185580]);
	test(className,"Mar" , jd,[290.15547, -22.79504, 290.45091, -22.75695, 224.5099,  14.5311, 22.5222447106,   13.30932578,    69.185580]);
	test(className,"Jup" , jd,[293.38974, -21.72807, 293.68175, -21.68434, 222.6401,  17.1504, 22.5222447106,   46.42294086,    69.185580]);
	test(className,"Sat" , jd,[301.28547, -20.36694, 301.57148, -20.30997, 216.8290,  22.2722, 22.5222447106,   87.83578976,    69.185580]);
	test(className,"Ura" , jd,[ 31.92196,  12.42367,  32.18530,  12.51588, 102.9224,  35.5201, 22.5222447106,  171.02647122,    69.185580]);
	test(className,"Nep" , jd,[349.69439,  -5.54002, 349.94553,  -5.43342, 162.8510,  44.8948, 22.5222447106,  257.15710804,    69.185580]);
	test(className,"Plu" , jd,[296.33960, -22.01640, 296.63097, -21.96754, 220.0897,  18.4463, 22.5222447106,  286.80349947,    69.185580]);
	test(className,"Moon", jd,[240.64213, -18.70797, 240.93143, -18.76252, 259.0324, -16.7669, 22.5222447106,    0.02082157,    69.185580]);

	
}

function test(className,name,jd,testData){
	let planet=0;
	switch (name){
		case "Sun": planet=10; break;
		case "Mer": planet=0; break;
		case "Ven": planet=1; break;
		case "Mar": planet=3; break;
		case "Jup": planet=4; break;
		case "Sat": planet=5; break;
		case "Ura": planet=6; break;
		case "Nep": planet=7; break;
		case "Plu": planet=8; break;
		case "Moon": planet=9; break;
	}
	const result=className.reduce(planet,jd,observer);
	displayResult("resultsTable1",name,result,testData);
}

function displayResult(tableName,name,r,t){
	const table=document.getElementById(tableName);
	const rows=table.rows.length;
	const row1=table.insertRow(rows);
	const row2=table.insertRow(rows+1);
	const row3=table.insertRow(rows+2);
	table.insertRow(rows+3).insertCell(0).innerHTML="&nbsp;";

	row1.insertCell(0).innerHTML=name+" Horizons";
	row2.insertCell(0).innerHTML=name+" Computed";
	row3.insertCell(0).innerHTML=name+" Error";

	for(let i=1;i<7;i++){
		const cell1=row1.insertCell(i);
		const cell2=row2.insertCell(i);
		const cell3=row3.insertCell(i);

		cell1.innerHTML=t[i-1];
		cell2.innerHTML=r[i-1]/toRad;
		cell3.innerHTML=Math.abs(t[i-1]-r[i-1]/toRad)*60*60+"\"";
	}
}
