const au=149597870691; // meters
const c=299792458; // m/s
const toRad=Math.PI/180;  //Multiply by this to convert degrees to radians

const observer=[38.2463888888*toRad,-85.763611111111*toRad,0]; //Louisville, Kentucky -85° 45' 49"  38 14' 47"

function runTests(className){
	//2458923.196782407 Horizons value for 2020-03-14 16:43:22
	const jd=gregorianDateToJulianDate(2020,3,14,16,43,22);

	//Test data from JPL Horizons, date chosen at random
	/*                    JDUT R.A.___(ICRF)___DEC   R.A._(a-appar)_DEC.   Azi_(a-appr)_Elev   L_Ap_Sid_Time   1-way_down_LT       TDB-UT*/
	test(className,"Sun" , jd,[354.760436383,  -2.268895824, 355.010107006,  -2.161023679, 154.5692,  46.5898,    8.26968075,    69.185580]);
	test(className,"Mer" , jd,[330.913290028, -11.424315216, 331.173632891, -11.329384129, 188.5413,  40.0342,    6.39316934,    69.185580]);
	test(className,"Ven" , jd,[36.870895652 ,  16.933544344,  37.141793123,  17.020722349,  95.0213,  34.3682,    6.54428376,    69.185580]);
	test(className,"Mar" , jd,[290.155465895, -22.795043107, 290.450914110, -22.756951643, 224.5099,  14.5311,   13.30932578,    69.185580]);
	test(className,"Jup" , jd,[293.389744621, -21.728065044, 293.681748056, -21.684335302, 222.6401,  17.1504,   46.42294086,    69.185580]);
	test(className,"Sat" , jd,[301.285465229, -20.366943048, 301.571483161, -20.309971390, 216.8290,  22.2722,   87.83578976,    69.185580]);
	test(className,"Ura" , jd,[31.921963455 ,  12.423672251,  32.185296307,  12.515877311, 102.9224,  35.5201,  171.02647122,    69.185580]);
	test(className,"Nep" , jd,[349.694394547,  -5.540021796, 349.945528234,  -5.433416185, 162.8510,  44.8948,  257.15710758,    69.185580]);
	test(className,"Moon", jd,[240.642130555, -18.707971308, 240.931432653, -18.762519744, 259.0324, -16.7669,    0.02082157,    69.185580]);
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
