function simplenutation(T){
    const omega=(125.04452 - 1934.136261*T + 0.0020708*T*T + (T*T*T)/450000)*Math.PI/180;
    const L=(280.4665 + 36000.7698*T)*Math.PI/180;
    const Lp=(218.3165 + 481267.8813*T)*Math.PI/180;
    const dPsi= -17.20*Math.sin(omega) - - 1.32*Math.sin(2*L) - 0.23*Math.sin(2*Lp) + 0.21*Math.sin(2*omega);
    const dEps= 9.20*Math.cos(omega) + .57*Math.cos(2*L) + .10*Math.cos(2*Lp) - .09*Math.cos(omega);

    return [dPsi/60/60*Math.PI/180,dEps/60/60*Math.PI/180];
}
