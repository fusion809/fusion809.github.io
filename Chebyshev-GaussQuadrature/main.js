var N, g, l, theta0, thetaDot0, thetaMax, thetaMaxInitial;
var integral = 0;

function addParameters() {
    g = parseFloat(document.getElementById("g").value);
    l = parseFloat(document.getElementById("l").value);
    N = parseFloat(document.getElementById("N").value);
    theta0 = parseFloat(document.getElementById("theta0").value);
    thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    thetaMaxInitial = parseFloat(document.getElementById("thetaMaxInitial").value);
}

function f(theta) {
    thetaDot2 = thetaDot0**2 + 2*g/l * (Math.sin(theta0)-Math.sin(theta));
    return thetaDot2
}

function fPrime(theta) {
    value = -2*g/l * Math.cos(theta);
    return value
}

function findZero() {
    thetaMax = thetaMaxInitial;
    adj = f(thetaMax)/fPrime(thetaMax);

    while (Math.abs(adj) >= 1e-14) {
        thetaMax -= adj;
        adj = f(thetaMax)/fPrime(thetaMax);
    }
    document.getElementById("thetaMaxDisplay").innerHTML = thetaMax;
}

function integrate() {
    var nodes = new Array(N);
    var integrand = new Array(N);
    var transformedGrid = new Array(N);
    integral = 0;

    for ( let i = 1; i < N+1; i++) {
        nodes[i-1] = Math.cos((2*i-1)*Math.PI/(2*N));
        transformedGrid[i-1] = (thetaMax-theta0)*nodes[i-1]/2+(thetaMax+theta0)/2;
        integrand[i-1] = Math.sqrt(1-nodes[i-1]**2)*Math.pow(f(transformedGrid[i-1]),-1/2);
        integral += ((thetaMax-theta0)/2) * (Math.PI/N)*integrand[i-1];
        // console.log("i is " + i);
        // console.log("nodes[i] is " + nodes[i]);
        // console.log("transformedGrid[i] is " + transformedGrid[i]);
        // console.log("integrand[i] is " + integrand[i]);
        // console.log("integral is " + integral);
    }
    document.getElementById("integralDisplay").innerHTML = 2*Math.abs(integral);
}