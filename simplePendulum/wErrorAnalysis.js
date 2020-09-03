/**
 * Right-hand side of our second-order ODE written as a system of 
 * first-order ODEs.
 *
 * @param g          Acceleration due to gravity in metres per second 
 * squared.
 * @param l          Length of the pendulum in metres.
 * @param t          Time.
 * @param theta      Angle of the pendulum from the positive x-axis.
 * @param thetaDot   Rate of change against time of the angle from the 
 * positive x-axis.
 * @return           [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, thetaDot) {
    return [thetaDot, -g/l*Math.cos(theta)];
};

// Initialize our global variables
var solution = {
    t: [],
    theta: [],
    thetaDot: [],
    errorThetaDot: [],
    logErrorThetaDot: [],
};
var windowInnerWidth;
var windowInnerHeight;
var epsilon, N, g, l, theta0, thetaDot0, adjMax, adjMin, thetaMin, thetaMax, j, k;

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
 * @param theta    Angle from the x-axis.
 * @return         Correction to theta min/max value according to Newton's
 * method.
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
    var thetaMaxInitial = parseFloat(document.getElementById("thetaMaxInitial").value);
    var xi = parseFloat(document.getElementById("xi").value);

    // Take our initial guess for thetaMax
    thetaMax = thetaMaxInitial;
    // thetaMax correction
    adjMax = newtonsCorrection(thetaMax);

    // Calculate when thetaDot = 0 next, which will be halfway through the problem's period
    j = 0; 
    k = 0;
    while ( ( Math.abs(adjMax) >= xi) && ( j < N )) {
        thetaMax += adjMax;
        adjMax = newtonsCorrection(thetaMax);
        j++;
    }
    thetaMax += adjMax;

    // Calculate thetaMin, which is when thetaDot = 0
    if (thetaDot0 != 0) {
        thetaMin = theta0;
        adjMin = newtonsCorrection(thetaMin);    
        while (( Math.abs(adjMin) >= xi) & (k < N)) {
            thetaMin += adjMin;
            adjMin = newtonsCorrection(thetaMin);
            k++;
        }
        thetaMin += adjMin;
    } else {
        thetaMin = theta0;
    }

    if ( (j >= N ) || (k >= N ) ) {
        document.getElementById("integralDisplay").innerHTML = "Theta min and/or theta max could not be calculated as the limit on Newton's method iterations was exceeded";
        document.getElementById("tf").value = 100;
        return
    }

    // Integrate problem from thetaMin to thetaMax to calculate the period
    // in seconds
    var nodes = 0;
    var integrand = 0;
    var transformedGrid = 0;
    integral = 0;
    for ( let i = 1; i < N+1; i++) {
        nodes = Math.cos((2*i-1)*Math.PI/(2*N));
        transformedGrid = (thetaMax-thetaMin)*nodes/2+(thetaMax+thetaMin)/2;
        integrand = Math.sqrt(1-nodes**2)*Math.pow(thetaDotSq(transformedGrid),-1/2);
        integral += ((thetaMax-thetaMin)/2) * (Math.PI/N)*integrand;
    }
    period = 2*Math.abs(integral);

    // Change what's displayed on the page accordingly
    document.getElementById("integralDisplay").innerHTML = period;
    document.getElementById("tf").value = 4*period;
    document.getElementById("jDisplay").innerHTML = j;
    document.getElementById("kDisplay").innerHTML = k;
    document.getElementById("adjMinDisplay").innerHTML = adjMin;
    document.getElementById("adjMaxDisplay").innerHTML = adjMax;
    document.getElementById("thetaMinDisplay").innerHTML = thetaMin;
    document.getElementById("thetaMaxDisplay").innerHTML = thetaMax;
}

/** 
 * Solve the problem using RK45.
 *
 * @params           None. Uses parameter values in the forum.
 * @return           Nothing. But it enters the solution values and error estimates
 * into the solution object.
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
    errorThetaDot = [0];
    logErrorThetaDot = [-Infinity];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in theta and 
        // thetaDot over the step
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
        errorThetaDot1 = Math.abs(thetaDot1 - Math.sign(thetaDot1)*Math.sqrt(Math.abs(thetaDotSq(theta1))));

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
            errorThetaDot.push(errorThetaDot1);
            logErrorThetaDot.push(Math.log10(errorThetaDot1));
            dt *= s;
            i++;
        } else {
            dt *= s;
        }
    }

    // Write our arrays to the solution object
    solution = {
        t: t,
        theta: theta,
        thetaDot: thetaDot,
        errorThetaDot: errorThetaDot,
        logErrorThetaDot: logErrorThetaDot,
    };
};

/**
 * Tabulates solution data.
 *
 * @params           None.
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable() {
    solveProblem();
    t = solution["t"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];
    errorThetaDot = solution["errorThetaDot"];
    rmsErrorthetaDot = Math.sqrt(math.dot(errorThetaDot, errorThetaDot)/(errorThetaDot.length-1));
    logErrorThetaDot = solution["logErrorThetaDot"];
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>&theta; (radians) </th>';
    tableContents += '<th>&theta; dot<br/>(radians &middot; s<sup>-1</sup>)</th>';
    tableContents += '<th>&theta; dot error</th>';
    tableContents += '</tr>';
    for (let i = 0; i < theta.length; i++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + theta[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + thetaDot[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + errorThetaDot[i].toExponential(5) + '</td>';
        tableContents += '</tr>';
    }
    tableContents += '<tr>';
    tableContents += '<td colspan = "6">RMS error in &theta; dot is: ';
    tableContents += rmsErrorthetaDot.toExponential(10);
    tableContents += '</td>';
    tableContents += '</tr>';
    document.getElementById('tableOutputs').innerHTML = tableContents;
};

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
    solveProblem();

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
 * Generate semilog plot of an error estimate in theta dot against time
 * 
 * @params           None.
 * @return           Nothing. Just generates the relevant plot.
 */
function generateErrorPlot() {
    solveProblem();

    // Extract solution data from solution object
    t = solution["t"];
    logErrorThetaDot = solution["logErrorThetaDot"];

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("errorPlot").style = "height: " + windowInnerHeight + "px;";
    
    // Logarithmic plot of the error in thetaDot
    var plot = {
        x: t,
        y: logErrorThetaDot,
        type: 'scatter',
        name: 'Semilog plot of error in theta dot',
    };
    data = [plot];
    var layout = {
        title: "Semilog plot of the error in theta dot against t",
        xaxis: {
            title: "Time (seconds)"
        },
        yaxis: {
            title: "Log of the error in theta dot"
        }
    };
    
    // Generate plot
    Plotly.newPlot('errorPlot', data, layout);
}

/**
 * Remove error plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeErrorPlot() {
    document.getElementById("errorPlot").innerHTML = '';
    document.getElementById("errorPlot").style = '';
}

/**
 * Generate plot of theta and theta dot against time
 * 
 * @params           None.
 * @return           Nothing. Just generates the relevant plot.
 */
function generateTimePlot() {
    solveProblem();

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
 * Generate three plots:
 * - one of thetaDot and theta against t;
 * - one of the log of the difference between thetaDot and that calculated
 * based on theta; and
 * - a phase plot of thetaDot against theta.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    generateTimePlot();
    generateErrorPlot();
    generatePhasePlot();
};

/**
 * Removes solution plots and associated element formatting
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("errorPlot").innerHTML = '';
    document.getElementById("phasePlot").innerHTML = '';
    document.getElementById("timePlot").style = ''
    document.getElementById("errorPlot").style = '';
    document.getElementById("phasePlot").style = '';
}