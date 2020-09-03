/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param g          Acceleration due to gravity in metres per second squared.
 * @param l          Length of the pendulum in metres.
 * @param t          Time (seconds).
 * @param theta      Angle of the pendulum from the positive x-axis (in
 * radians)
 * @param thetaDot   Rate of change against time of the angle from the positive x-axis.
 * @return           [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, thetaDot) {
    return [thetaDot, -g/l*Math.cos(theta)];
}

// Initialize our global variables
var solution = {
    t: [],
    theta: [],
    thetaDot: []
};
var windowInnerWidth;
var windowInnerHeight;
var epsilon, N, g, l, theta0, thetaDot0, thetaMax, thetaMaxInitial;
var integral = 0;

/**
 * Returns the value of theta dot squared.
 * 
 * @param theta    Angle from the x-axis.
 * @return         Theta dot squared value.
 */
function thetaDotSq(theta) {
    return thetaDot0**2 + 2*g/l * (Math.sin(theta0)-Math.sin(theta));
}

/**
 * The derivative with respect to theta of the theta dot squared expression
 * 
 * @param theta    Angle from the x-axis.
 * @return         Theta derivative of the theta dot squared expression.
 */
function thetaDotSqPrime(theta) {
    return -2*g/l * Math.cos(theta)
}

/**
 * Correction to theta according to Newton's method
 * 
 * @param {*} theta 
 */
function newtonsCorrection(theta) {
    return -thetaDotSq(theta)/thetaDotSqPrime(theta);
}

/**
 * Calculates the value of theta when theta dot = 0 using Newton's method
 * then uses Chebyshev-Gauss quadrature to compute the time taken to reach this period.
 * 
 * @params         None.
 * @return         Nothing. Changes the element value/innerHTML of 
 * integralDisplay and tf to T and 4*T, respectively.
 */
function periodCalc() {
    // Parameters of the problem that are necessary to calculate the period
    g = parseFloat(document.getElementById("g").value);
    l = parseFloat(document.getElementById("l").value);
    N = parseFloat(document.getElementById("N").value);
    theta0 = parseFloat(document.getElementById("theta0").value);
    thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    thetaMaxInitial = parseFloat(document.getElementById("thetaMaxInitial").value);
    var errorTol = 1e-14;

    // Take our initial guess for thetaMax
    thetaMax = thetaMaxInitial;
    // thetaMax correction
    adj = newtonsCorrection(thetaMax);

    // Calculate when thetaDot = 0 next, which will be halfway through the problem's period
    while (Math.abs(adj) >= errorTol) {
        thetaMax += adj;
        adj = newtonsCorrection(thetaMax);
    }

    // Calculate thetaMin, which is when thetaDot = 0
    if (thetaDot0 != 0) {
        var thetaMin = theta0;
        adj = newtonsCorrection(thetaMin);
        while (Math.abs(adj) >= errorTol) {
            thetaMin += adj;
            adj = newtonsCorrection(thetaMin);
        }
    } else {
        thetaMin = theta0;
    }

    // Integrate problem from thetaMin to thetaMax to calculate the period
    // in seconds
    var nodes = new Array(N);
    var integrand = new Array(N);
    var transformedGrid = new Array(N);
    integral = 0;
    for ( let i = 1; i < N+1; i++) {
        nodes[i-1] = Math.cos((2*i-1)*Math.PI/(2*N));
        transformedGrid[i-1] = (thetaMax-thetaMin)*nodes[i-1]/2+(thetaMax+thetaMin)/2;
        integrand[i-1] = Math.sqrt(1-nodes[i-1]**2)*Math.pow(thetaDotSq(transformedGrid[i-1]),-1/2);
        integral += ((thetaMax-thetaMin)/2) * (Math.PI/N)*integrand[i-1];
    }
    period = 2*Math.abs(integral);
    document.getElementById("integralDisplay").innerHTML = period;
    document.getElementById("tf").value = 4*period;
}

/** 
 * Solve the problem using RK45.
 *
 * @params           None. Uses parameter values in the forum.
 * @return           Nothing. But it enters the solution values into the solution
 * object.
 */
function solveProblem() {
    // Obtain the parameters of the problem
    g = parseFloat(document.getElementById("g").value);
    l = parseFloat(document.getElementById("l").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    theta0 = parseFloat(document.getElementById("theta0").value);
    thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    theta = [theta0];
    thetaDot = [thetaDot0];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in theta and thetaDot
        // over the step
        k1 = dt*f(g, l, t[i], theta[i], thetaDot[i])[0];
        l1 = dt*f(g, l, t[i], theta[i], thetaDot[i])[1];
        k2 = dt*f(g, l, t[i]+dt/4, theta[i]+k1/4, thetaDot[i]+l1/4)[0];
        l2 = dt*f(g, l, t[i]+dt/4, theta[i]+k1/4, thetaDot[i]+l1/4)[1];
        k3 = dt*f(g, l, t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, thetaDot[i]+3*l1/32+9*l2/32)[0];
        l3 = dt*f(g, l, t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, thetaDot[i]+3*l1/32+9*l2/32)[1];
        k4 = dt*f(g, l, t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, thetaDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
        l4 = dt*f(g, l, t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, thetaDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
        k5 = dt*f(g, l, t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, thetaDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
        l5 = dt*f(g, l, t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, thetaDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
        k6 = dt*f(g, l, t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, thetaDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
        l6 = dt*f(g, l, t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, thetaDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

        // theta1 and thetaDot1 are our fourth order approximations
        theta1 = theta[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        thetaDot1 = thetaDot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        // theta2 and thetaDot2 are our fifth order approximations
        theta2 = theta[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        thetaDot2 = thetaDot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

        // The following are used to correct the step size
        Rtheta = Math.abs(theta1-theta2)/dt;
        RthetaDot = Math.abs(thetaDot1-thetaDot2)/dt;
        stheta = 0.84*Math.pow(epsilon/Rtheta, 1/4);                
        sthetaDot = 0.84*Math.pow(epsilon/RthetaDot, 1/4);
        R = Math.max(Rtheta, RthetaDot);
        s = Math.min(stheta, sthetaDot);
        if ( R <= epsilon ) {
            t.push(t[i]+dt);
            theta.push(theta1);
            thetaDot.push(thetaDot1);
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, theta and thetaDot to our solution object
    solution = {
        t: t,
        theta: theta,
        thetaDot: thetaDot,
    };
}

/**
 * Tabulates solution data.
 *
 * @params           None. Uses the entries of the solution object, however. 
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable() {
    if ( solution["t"].length == 0) {
        alert('You haven\'t solved the problem yet! Press the "Solve the problem" button before pressing the "Tabulate the solution" button again.');
        return
    }
    t = solution["t"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>&theta; (radians) </th>';
    tableContents += '<th>&theta; dot <br/>(radians &middot; s<sup>-1</sup>)</th>';
    tableContents += "</tr>";
    for (let j = 0; j < theta.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + theta[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + thetaDot[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '</tr>';
    }
    document.getElementById('tableOutputs').innerHTML = tableContents;
}

/**
 * Removes the solution table
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution table.
 */
function removeTable() {
    document.getElementById('tableOutputs').innerHTML = '';
}

/**
 * Generate phase plot of theta dot against theta
 * 
 * @params           None.
 * @return           Nothing. Just generates the relevant plot.
 */
function generatePhasePlot() {
    // Run solveProblem() if previously unrun
    if ( solution["t"].length == 0) {
        solveProblem();
    };

    // Extract solution data from solution object
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlot").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the phase plot
    var plot = {
        x: theta,
        y: thetaDot,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout = {
        title: "Phase plot of theta dot against theta",
        xaxis: {
            title: "theta (radians)"
        },
        yaxis: {
            title: "theta dot (radians per second)"
        }
    };
    data = [plot];

    // Generate plot
    Plotly.newPlot('phasePlot', data, layout);
}

/**
 * Remove phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removePhasePlot() {
    document.getElementById("phasePlot").innerHTML = '';
    document.getElementById("phasePlot").style = '';
}

/**
 * Generate plot of theta and theta dot against time
 * 
 * @params           None.
 * @return           Nothing. Just generates the relevant plot.
 */
function generateTimePlot() {
    // Run solveProblem() if previously unrun
    if ( solution["t"].length == 0) {
        solveProblem();
    };

    // Extract solution data from solution object
    t = solution["t"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the theta and theta dot against time plot
    var plotTheta = {
        x: t,
        y: theta,
        type: 'scatter',
        name: "theta (radians)"
    };
    var plotThetaDot = {
        x: t,
        y: thetaDot,
        type: 'scatter',
        name: "theta dot (radians per second)"
    };
    var layout = {
        title: 'theta and theta dot against time plots',
        xaxis: {
            title: 'Time (seconds)'
        }
    };
    data = [plotTheta, plotThetaDot];

    // Generate plots
    Plotly.newPlot('timePlot', data, layout);
}

/**
 * Remove time plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeTimePlot() {
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("timePlot").style = '';
}

/**
 * Generate two plots:
 * - one of thetaDot and theta against t; and
 * - a phase plot of thetaDot against theta.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    generateTimePlot();
    generatePhasePlot();
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("phasePlot").innerHTML = '';
    document.getElementById("timePlot").style = '';
    document.getElementById("phasePlot").style = '';
};