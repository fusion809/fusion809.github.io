/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param alpha      Parameter.
 * @param beta       Parameter.
 * @param gamma      Parameter.
 * @param delta      Parameter.
 * @param omega      Parameter.
 * @param t          Time (seconds).
 * @param x          x coordinate.
 * @param xDot       dx/dt
 * @return           [dx/dt, d2x/dt2]
 */
function f(alpha, beta, gamma, delta, omega, t, x, xDot) {
    return [xDot, - delta*xDot - alpha*x - beta*x**3 + gamma * Math.cos(omega*t)];
}

/**
 * Read inputs from the form and enter it into an object.
 * 
 * @params              None.
 * @return              An object containing all the form inputs.
 */
function readInputs() {
    // Obtain the parameters of the problem
    var alpha = parseFloat(document.getElementById("alpha").value);
    var beta = parseFloat(document.getElementById("beta").value);
    var gamma = parseFloat(document.getElementById("gamma").value);
    var delta = parseFloat(document.getElementById("delta").value);
    var omega = parseFloat(document.getElementById("omega").value);
    var t0 = parseFloat(document.getElementById("t0").value);
    var tf = parseFloat(document.getElementById("tf").value);
    var x0 = parseFloat(document.getElementById("x0").value);
    var xDot0 = parseFloat(document.getElementById("xDot0").value);
    var epsilon = parseFloat(document.getElementById("epsilon").value);
    var dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Object containing inputs
    var objectOfInputs = {
        alpha: alpha,
        beta: beta,
        gamma: gamma,
        delta: delta,
        omega: omega,
        t0: t0,
        tf: tf,
        x0: x0,
        xDot0: xDot0,
        epsilon: epsilon,
        dtInitial: dtInitial
    }
    return objectOfInputs;
}

/**
 * Approximate 
 * @param dt            Step size.
 * @param alpha         Problem parameter.
 * @param beta          Problem parameter.
 * @param gamma         Problem parameter.
 * @param delta         Problem parameter.
 * @param omega         Problem parameter.
 * @param t             An array of t values.
 * @param x             An array of x values.
 * @param xDot          An array of xDot values.
 * @param i             Counter variable.
 * @return              [x1, xDot1, x2, xDot2]
 */
function approximatorRKF45(dt, alpha, beta, gamma, delta, omega, t, x, xDot, i) {
    // Runge-Kutta-Fehlberg approximations of the change in x and xDot over the step
    var k1 = dt*f(alpha, beta, gamma, delta, omega, t[i], x[i], xDot[i])[0];
    var l1 = dt*f(alpha, beta, gamma, delta, omega, t[i], x[i], xDot[i])[1];
    var k2 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4)[0];
    var l2 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4)[1];
    var k3 = dt*f(alpha, beta, gamma, delta, omega, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32)[0];
    var l3 = dt*f(alpha, beta, gamma, delta, omega, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32)[1];
    var k4 = dt*f(alpha, beta, gamma, delta, omega, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
    var l4 = dt*f(alpha, beta, gamma, delta, omega, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
    var k5 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
    var l5 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
    var k6 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
    var l6 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

    // x1 and xDot1 are our fourth order approximations
    var x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    var xDot1 = xDot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    // x2 and xDot2 are our fifth order approximations
    var x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    var xDot2 = xDot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

    return [x1, xDot1, x2, xDot2];
}

/**
 * Check and correct step size
 * 
 * @param dt            Step size.
 * @param epsilon       Error tolerance.
 * @param t             An array of t values.
 * @param x             An array of x values.
 * @param xDot          An array of xDot values.
 * @param x1            4th order approx to x.
 * @param xDot1         4th order approx to xDot.
 * @param x2            5th order approx to x.
 * @param xDot2         5th order approx to xDot.
 * @param i             Counter variable.
 * @return              [dt, t, x, xDot, i]
 */
function stepSizeChecker(dt, epsilon, t, x, xDot, x1, xDot1, x2, xDot2, i) {
    // Initialize relevant variables
    var Rx = Math.abs(x1-x2)/dt;
    var RxDot = Math.abs(xDot1-xDot2)/dt;
    var sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
    var sxDot = 0.84*Math.pow(epsilon/RxDot, 1/4);
    var R = Math.max(Rx, RxDot);
    var s = Math.min(sx, sxDot);

    // If R is less than the margin of error, move on to next iteration
    if ( R <= epsilon ) {
        t.push(t[i]+dt);
        x.push(x1);
        xDot.push(xDot1);
        i++;
        dt *= s;
    } else {
        dt *= s;
    }

    // Return what variables are needed by solveProblem()
    return [dt, t, x, xDot, i];
}

/** 
 * Solve the problem using RKF45.
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               A solution object.
 */
function solveProblem(objectOfInputs) {
    // Obtain the parameters of the problem
    var {alpha, beta, gamma, delta, omega, t0, tf, x0, xDot0, epsilon, dtInitial} = objectOfInputs;

    // Initialize the arrays used and loop variables
    var t = [t0];
    var x = [x0];
    var xDot = [xDot0];
    var dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);
        var [x1, xDot1, x2, xDot2] = approximatorRKF45(dt, alpha, beta, gamma, delta, omega, t, x, xDot, i);
        var [dt, t, x, xDot, i] = stepSizeChecker(dt, epsilon, t, x, xDot, x1, xDot1, x2, xDot2, i);
    }

    // Write t, x and xDot to our solution object
    var solution = {
        t: t,
        x: x,
        xDot: xDot,
    };
    return solution;
}

/**
 * Tabulates solution data.
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(objectOfInputs);

    // Extract solution/parameter values
    var {t, x, xDot} = solution;
    var epsilon = objectOfInputs.epsilon;

    // Create table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>x (radians) </th>';
    tableContents += '<th>x dot <br/>(radians &middot; s<sup>-1</sup>)</th>';
    tableContents += "</tr>";
    for (let j = 0; j < x.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + x[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + xDot[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
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
 * Generate phase plot of x dot against x
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {x, xDot} = solution;

    // Height and width of plot
    adjustPlotHeight("phasePlot");

    // Characteristics of the phase plot
    var plot = {
        x: x,
        y: xDot,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout = {
        title: "Phase plot of x dot against x",
        xaxis: {
            title: "x (metres)"
        },
        yaxis: {
            title: "x dot (metres per second)"
        }
    };
    var data = [plot];

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
    rmPlot("phasePlot");
}

/**
 * Generate plot of x and x dot against time
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(objectOfInputs) {
    // Run solveProblem() if previously unrun
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {t, x, xDot} = solution;

    // Height and width of plots
    adjustPlotHeight("timePlot");

    // Characteristics of the x and x dot against time plot
    var plotx = {
        x: t,
        y: x,
        type: 'scatter',
        name: "x (metres)"
    };
    var plotxDot = {
        x: t,
        y: xDot,
        type: 'scatter',
        name: "x dot (metres per second)"
    };
    var layout = {
        title: 'x and x dot against time plots',
        xaxis: {
            title: 'Time (seconds)'
        }
    };
    var data = [plotx, plotxDot];

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
    rmPlot("timePlot");
}

/**
 * Generate two plots:
 * - one of xDot and x against t; and
 * - a phase plot of xDot against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generateTimePlot(objectOfInputs);
    generatePhasePlot(objectOfInputs);
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    removeTimePlot();
    removePhasePlot();
};
