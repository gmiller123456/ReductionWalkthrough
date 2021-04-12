function nutation2000BGeneratedTrunc(T){
    //The IAU Resolutions on Astronomical Reference Systems, Time Scales, and Earth Rotation Models Explanation and Implementation (George H. Kaplan)
    //https://arxiv.org/pdf/astro-ph/0602086.pdf page 88
    //IAU 2000B Nutation truncated to 6 terms
    let dp=0;
    let de=0;
    let arg;

    const DAS2R = 1/60/60*Math.PI/180.0;

    const T2 = T * T;
    const T3 = T * T2;
    const T4 = T * T3;

    //Fundamental Arguments p46 eq 5.17, 5.18, 5.19
    const Lp = DAS2R*(1287104.79305  + 129596581.0481*T  - 0.5532*T2 + 0.000136*T3 - 0.00001149*T4);
    const F  = DAS2R*(335779.526232 + 1739527262.8478*T - 12.7512*T2 - 0.001037*T3 + 0.00000417*T4);
    const D  = DAS2R*(1072260.70369 + 1602961601.2090*T  - 6.3706*T2 + 0.006593*T3 - 0.00003169*T4);
    const Om = DAS2R*(450160.398036    - 6962890.5431*T  + 7.4722*T2 + 0.007702*T3 - 0.00005939*T4);

    arg= Lp + 2*(F - D + Om);
    dp+=(-516821 + 1226*T)*Math.sin(arg) + -524*Math.cos(arg);
    de+=(224386 + -677*T)*Math.cos(arg) + -174*Math.sin(arg);
    arg= Lp;
    dp+=(1475877 + -3633*T)*Math.sin(arg) + 11817*Math.cos(arg);
    de+=(73871 + -184*T)*Math.cos(arg) + -1924*Math.sin(arg);
    arg= 2*Om;
    dp+=(2074554 + 207*T)*Math.sin(arg) + -698*Math.cos(arg);
    de+=(-897492 + 470*T)*Math.cos(arg) + -291*Math.sin(arg);
    arg= 2*(F + Om);
    dp+=(-2276413 + -234*T)*Math.sin(arg) + 2796*Math.cos(arg);
    de+=(978459 + -485*T)*Math.cos(arg) + 1374*Math.sin(arg);
    arg= 2*(F - D + Om);
    dp+=(-13170906 + -1675*T)*Math.sin(arg) + -13696*Math.cos(arg);
    de+=(5730336 + -3015*T)*Math.cos(arg) + -4587*Math.sin(arg);
    arg= 1*Om;
    dp+=(-172064161 + -174666*T)*Math.sin(arg) + 33386*Math.cos(arg);
    de+=(92052331 + 9086*T)*Math.cos(arg) + 15377*Math.sin(arg);

    const U2R = DAS2R/1E7;
    return [dp*U2R,de*U2R];
}

