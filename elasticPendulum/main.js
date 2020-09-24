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
    var xDDot = (l0+x)*thetaDot**2 - k*x/m + g*Math.sin(theta);
    var thetaDDot = -g*Math.cos(theta)/(l0+x)-2*xDot*thetaDot/(l0+x);
    return [xDot, xDDot, thetaDot, thetaDDot];
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
 * Read inputs in the form and export them in an array.
 * 
 * @params              None.
 * @return              An array of parameter values extracted from the form.
 */
function readInputs() {
    var g = parseFloat(document.getElementById("g").value);
    var l0 = parseFloat(document.getElementById("l0").value);
    var k = parseFloat(document.getElementById("k").value);
    var m = parseFloat(document.getElementById("m").value);
    var t0 = parseFloat(document.getElementById("t0").value);
    var tf = parseFloat(document.getElementById("tf").value);
    var x0 = parseFloat(document.getElementById("x0").value);
    var xDot0 = parseFloat(document.getElementById("xDot0").value);
    var theta0 = parseFloat(document.getElementById("theta0").value);
    var thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    var epsilon = parseFloat(document.getElementById("epsilon").value);
    var dtInitial = parseFloat(document.getElementById("dtInitial").value);

    var objectOfInputs = {
        g: g,
        l0: l0,
        k: k,
        m: m,
        t0: t0,
        tf: tf,
        x0: x0,
        xDot0: xDot0, 
        theta0: theta0,
        thetaDot0: thetaDot0,
        epsilon: epsilon,
        dtInitial: dtInitial
    }
    return objectOfInputs;
}

/** 
 * Solve the problem using RK45.
 *
 * @params           None. Uses parameter values in the forum.
 * @return           Nothing. But it enters the solution values into the solution
 * object.
 */
function solveProblem(objectOfInputs) {
    // Obtain the parameters of the problem
    var g = objectOfInputs.g;
    var l0 = objectOfInputs.l0;
    var k = objectOfInputs.k;
    var m = objectOfInputs.m;
    var t0 = objectOfInputs.t0;
    var tf = objectOfInputs.tf;
    var x0 = objectOfInputs.x0;
    var xDot0 = objectOfInputs.xDot0;
    var theta0 = objectOfInputs.theta0;
    var thetaDot0 = objectOfInputs.thetaDot0;
    var epsilon = objectOfInputs.epsilon;
    var dtInitial = objectOfInputs.dtInitial;

    // Initialize the arrays used and loop variables
    var t = [t0];
    var x = [x0];
    var xDot = [xDot0];
    var theta = [theta0];
    var thetaDot = [thetaDot0];
    var dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in x, y and z
        // over the step
        // 1st approx
        var k1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[0];
        var l1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[1];
        var m1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[2];
        var n1 = dt*f(g, l0, k, m, t[i], x[i], xDot[i], theta[i], thetaDot[i])[3];

        // 2nd approx
        var k2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[0];
        var l2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[1];
        var m2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[2];
        var n2 = dt*f(g, l0, k, m, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4, theta[i]+m1/4, thetaDot[i]+n1/4)[3];

        // 3rd approx
        var k3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[0];
        var l3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[1];
        var m3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[2];
        var n3 = dt*f(g, l0, k, m, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, thetaDot[i]+3*n1/32+9*n2/32)[3];

        // 4th approx
        var k4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[0];
        var l4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[1];
        var m4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[2];
        var n4 = dt*f(g, l0, k, m, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, thetaDot[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[3];

        // 5th approx
        var k5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[0];
        var l5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[1];
        var m5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[2];
        var n5 = dt*f(g, l0, k, m, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, thetaDot[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[3];

        // 6th approx
        var k6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[0];
        var l6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[1];
        var m6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[2];
        var n6 = dt*f(g, l0, k, m, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, thetaDot[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[3];

        // x1, xDot1, theta1, thetaDot1 are our fourth order approximations
        var x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        var xDot1 = xDot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        var theta1 = theta[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;
        var thetaDot1 = thetaDot[i] + 25*n1/216+1408*n3/2565+2197*n4/4104-n5/5;

        // x2, xDot2, theta2 and thetaDot2 are our fifth order approximations
        var x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        var xDot2 = xDot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        var theta2 = theta[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;
        var thetaDot2 = thetaDot[i] + 16*n1/135+6656*n3/12825+28561*n4/56430-9*n5/50+2*n6/55;

        // The following are used to correct the step size
        var Rx = Math.abs(x1-x2)/dt;
        var RxDot = Math.abs(xDot1-xDot2)/dt;
        var Rtheta = Math.abs(theta1-theta2)/dt;
        var RthetaDot = Math.abs(thetaDot1-thetaDot2)/dt;
        var sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        var sxDot = 0.84*Math.pow(epsilon/RxDot, 1/4);
        var stheta = 0.84*Math.pow(epsilon/Rtheta, 1/4);
        var sthetaDot = 0.84*Math.pow(epsilon/RthetaDot, 1/4);
        var R = Math.max(Rx, RxDot, Rtheta, RthetaDot);
        var s = Math.min(sx, sxDot, stheta, sthetaDot);

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
function fillTable(objectOfInputs) {
    // Return an error if solveProblem() hasn't been run
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extract coordinate arrays from the solution object
    var t = solution.t;
    var x = solution.x;
    var xDot = solution.xDot;
    var theta = solution.theta;
    var thetaDot = solution.thetaDot;
    var epsilon = objectOfInputs.epsilon;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
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
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extract solution data from solution object
    var x = solution.x;
    var theta = solution.theta;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXTheta").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXTheta = {
        x: x,
        y: theta,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXTheta = [plotXTheta];

    // layout object
    var layoutXTheta = {
        title: 'Phase plot of theta against x'
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
 * Generates a xdot against x phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXXDotPhasePlot() {
    // Run solveProblem if unrun
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extract solution data from solution object
    var x = solution.x;
    var xDot = solution.xDot;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXXDot").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXXDot = {
        x: x,
        y: xDot,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXXDot = [plotXXDot];

    // layout object
    var layoutXXDot = {
        title: "xdot against x phase plot"
    };

    // Generate plot
    Plotly.newPlot('phasePlotXXDot', dataXXDot, layoutXXDot);
}

/**
 * Remove xdot against x phase plot
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
    if ( solution.t.length == 0) {
        solveProblem();
    }
    
    // Extract solution data from solution object
    var theta = solution.theta;
    var thetaDot = solution.thetaDot;
    
    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotThetaThetaDot").style = "height: " + windowInnerHeight + "px;";
    
    // Plot object and data object array
    var plotThetaThetaDot = {
        x: theta,
        y: thetaDot,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataThetaThetaDot = [plotThetaThetaDot];
    
    // layout object
    var layoutThetaThetaDot = {
        title: "theta dot against theta phase plot"
    };
    
    // Generate plot
    Plotly.newPlot('phasePlotThetaThetaDot', dataThetaThetaDot, layoutThetaThetaDot);
}

/**
 * Remove theta dot against theta phase plot
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
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extract solution data from solution object
    var t = solution.t;
    var x = solution.x;
    var xDot = solution.xDot;
    var theta = solution.theta;
    var thetaDot = solution.thetaDot;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
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
    var dataTimePlot = [plotTX, plotTXDot, plotTTheta, plotTThetaDot];

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
 * Generate four plots:
 * - The first is a phase plot of theta against x.
 * - The second is a phase plot of x dot against x.
 * - The third is a phase plot of theta dot against theta.
 * - The fourth is a plot of x, x dot, theta and theta dot against time.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    generateXThetaPhasePlot();
    generateXXDotPhasePlot();
    generateThetaThetaDotPhasePlot();
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
    // Time plots
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("timePlot").style = '';
    document.getElementById("phasePlotThetaThetaDot").innerHTML = '';
    document.getElementById("phasePlotThetaThetaDot").style = '';
    document.getElementById("phasePlotXXDot").innerHTML = '';
    document.getElementById("phasePlotXXDot").style = '';
    document.getElementById("phasePlotXTheta").innerHTML = '';
    document.getElementById("phasePlotXTheta").style = '';
};
