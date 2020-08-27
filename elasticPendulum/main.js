/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param g          Interaction parameter.
 * @param l0         Interaction parameter.
 * @param k          Elasticity parameter.
 * @param m          Interaction parameter.
 * @param t          Time (seconds).
 * @param x          x value.
 * @param xDot       xDot value.
 * @param theta      theta value.
 * @param thetaDot   theta dot value.
 * @return           [dx/dt, d2x/dt2, dtheta/dt, d2theta/dt2]
 */
function f(g, l0, k, m, t, x, xDot, theta, thetaDot) {
    return [xDot, (l0+x)*thetaDot**2 - k*x/m + g*Math.sin(theta), thetaDot, -g*Math.cos(theta)/(l0+x)-2*xDot*thetaDot/(l0+x)];
}

// Initialize our global variables
var solution = {
    t: [],
    x: [],
    xDot: [],
    theta: [],
    thetaDot: []
};
var epsilon;

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
    l0 = parseFloat(document.getElementById("l0").value);
    k = parseFloat(document.getElementById("k").value);
    m = parseFloat(document.getElementById("m").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    x0 = parseFloat(document.getElementById("x0").value);
    xDot0 = parseFloat(document.getElementById("xDot0").value);
    theta0 = parseFloat(document.getElementById("theta0").value);
    thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    x = [x0];
    xDot = [xDot0];
    theta = [theta0];
    thetaDot = [thetaDot0];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in x, y and z
        // over the step
        // 1st approx
        k1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[0];
        l1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[1];
        m1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[2];
        n1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[3];

        // 2nd approx
        k2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[0];
        l2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[1];
        m2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[2];
        n2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[3];

        // 3rd approx
        k3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[0];
        l3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[1];
        m3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[2];
        n3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[3];

        // 4th approx
        k4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[0];
        l4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[1];
        m4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[2];
        n4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[3];

        // 5th approx
        k5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[0];
        l5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[1];
        m5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[2];
        n5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[3];

        // 6th approx
        k6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[0];
        l6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[1];
        m6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[2];
        n6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[3];

        // x1, y1 and z1 are our fourth order approximations
        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        xDot1 = xDot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        theta1 = theta[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;
        thetaDot1 = thetaDot[i] + 25*n1/216+1408*n3/2565+2197*n4/4104-n5/5;

        // x2, y2 and z2 are our fifth order approximations
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        xDot2 = xDot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        theta2 = theta[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;
        thetaDot2 = thetaDot[i] + 16*n1/135+6656*n3/12825+28561*n4/56430-9*n5/50+2*n6/55;

        // The following are used to correct the step size
        Rx = Math.abs(x1-x2)/dt;
        RxDot = Math.abs(xDot1-xDot2)/dt;
        Rtheta = Math.abs(theta1-theta2)/dt;
        RthetaDot = Math.abs(thetaDot1-thetaDot2)/dt;
        sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        sxDot = 0.84*Math.pow(epsilon/RxDot, 1/4);
        stheta = 0.84*Math.pow(epsilon/Rtheta, 1/4);
        sthetaDot = 0.84*Math.pow(epsilon/RthetaDot, 1/4);
        R = Math.max(Rx, RxDot, Rtheta, RthetaDot);
        s = Math.min(sx, sxDot, stheta, sthetaDot);

        // If R is less than or equal to epsilon move onto the next step
        if ( R <= epsilon ) {
            t.push(t[i]+dt);
            x.push(x1);
            xDot.push(xDot1);
            theta.push(theta1);
            thetaDot.push(thetaDot1);
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, x, y and z to our solution object
    solution = {
        t: t,
        x: x,
        xDot: xDot,
        theta: theta,
        thetaDot: thetaDot
    };
}

/**
 * Tabulates solution data.
 *
 * @params           None. Uses the entries of the solution object, however. 
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable() {
    // Return an error if solveProblem() hasn't been run
    if ( solution["t"].length == 0) {
        solveProblem();
    }

    // Extract coordinate arrays from the solution object
    t = solution["t"];
    x = solution["x"];
    xDot = solution["xDot"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>x</th>';
    tableContents += '<th>x dot</th>';
    tableContents += '<th>&theta;</th>';
    tableContents += '<th>&theta; dot</th>';
    tableContents += "</tr>";
    for (let j = 0; j < x.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + x[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + xDot[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
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
    // Clear table content
    document.getElementById('tableOutputs').innerHTML = '';
}

/**
 * Generates a 2D phase plot of x against theta
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXThetaPhasePlot() {
    // Run solveProblem if unrun
    if ( solution["t"].length == 0) {
        solveProblem();
    }

    // Extra solution data from solution object
    x = solution["x"];
    theta = solution["theta"];

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXTheta").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXTheta = {
        x: x,
        y: theta,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    dataXTheta = [plotXTheta];

    // layout object
    var layoutXTheta = {
        title: 'Phase plot of the solution'
    };

    // Generate plot
    Plotly.newPlot('phasePlotXTheta', dataXTheta, layoutXTheta);
}

/**
 * Remove XTheta phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXThetaPhasePlot() {
    document.getElementById("phasePlotXTheta").innerHTML = '';
    document.getElementById("phasePlotXTheta").style = '';
}

/**
 * Generates a XY phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXXDotPhasePlot() {
    // Run solveProblem if unrun
    if ( solution["t"].length == 0) {
        solveProblem();
    }

    // Extra solution data from solution object
    x = solution["x"];
    xDot = solution["xDot"];

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXXDot").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXXDot = {
        x: x,
        y: xDot,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    dataXXDot = [plotXXDot];

    // layout object
    var layoutXXDot = {
        title: "x/k phase plot"
    };

    // Generate plot
    Plotly.newPlot('phasePlotXXDot', dataXXDot, layoutXXDot);
}

/**
 * Remove x/xdot phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXXDotPhasePlot() {
    document.getElementById("phasePlotXXDot").innerHTML = '';
    document.getElementById("phasePlotXXDot").style = '';
}

/**
 * Generates a theta dot against theta phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateThetaThetaDotPhasePlot() {
    // Run solveProblem if unrun
    if ( solution["t"].length == 0) {
        solveProblem();
    }
    
    // Extra solution data from solution object
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];
    
    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotThetaThetaDot").style = "height: " + windowInnerHeight + "px;";
    
    // Plot object and data object array
    var plotThetaThetaDot = {
        x: theta,
        y: thetaDot,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    dataThetaThetaDot = [plotThetaThetaDot];
    
    // layout object
    var layoutThetaThetaDot = {
        title: "theta dot against theta phase plot"
    };
    
    // Generate plot
    Plotly.newPlot('phasePlotThetaThetaDot', dataThetaThetaDot, layoutThetaThetaDot);
}

/**
 * Remove theta/theta dot phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeThetaThetaDotPhasePlot() {
    document.getElementById("phasePlotThetaThetaDot").innerHTML = '';
    document.getElementById("phasePlotThetaThetaDot").style = '';
}

/**
 * Generates a time plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateTimePlot() {
    // Run solveProblem if unrun
    if ( solution["t"].length == 0) {
        solveProblem();
    }

    // Extra solution data from solution object
    t = solution["t"];
    x = solution["x"];
    xDot = solution["xDot"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotTX = {
        x: t,
        y: x,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'x'
    };
    var plotTXDot = {
        x: t,
        y: xDot,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'x dot'
    };
    var plotTTheta = {
        x: t,
        y: theta,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'theta'
    };
    var plotTThetaDot = {
        x: t,
        y: thetaDot,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'theta dot'
    };
    dataTimePlot = [plotTX, plotTXDot, plotTTheta, plotTThetaDot];

    // layout object
    var layoutTimePlot = {
        title: "Time plots of the solution to the problem"
    };

    // Generate plot
    Plotly.newPlot('timePlot', dataTimePlot, layoutTimePlot);
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
 * Generate five plots:
 * - The first is a 3D phase plot of x, y and z.
 * - The second is a 2D phase plot of y against x.
 * - The third is a 2D phase plot of z against x.
 * - The fourth is a 2D phase plot of z against y.
 * - The fifth is a plot of x, y and z against time.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    generate3DPhasePlot();
    generateXYPhasePlot();
    generateXZPhasePlot();
    generateYZPhasePlot();
    generateTimePlot();
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    // Clear HTML and CSS of the plots
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("phasePlotXYZ").innerHTML = '';
    document.getElementById("phasePlotXY").innerHTML = '';
    document.getElementById("phasePlotXZ").innerHTML = '';
    document.getElementById("phasePlotYZ").innerHTML = '';
    document.getElementById("timePlot").style = '';
    document.getElementById("phasePlotXYZ").style = '';
    document.getElementById("phasePlotXY").style = '';
    document.getElementById("phasePlotXZ").style = '';
    document.getElementById("phasePlotYZ").style = '';
};