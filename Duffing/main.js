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

// Initialize our global variables
var solution = {
    t: [],
    x: [],
    xDot: []
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
    alpha = parseFloat(document.getElementById("alpha").value);
    beta = parseFloat(document.getElementById("beta").value);
    gamma = parseFloat(document.getElementById("gamma").value);
    delta = parseFloat(document.getElementById("delta").value);
    omega = parseFloat(document.getElementById("omega").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    x0 = parseFloat(document.getElementById("x0").value);
    xDot0 = parseFloat(document.getElementById("xDot0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    x = [x0];
    xDot = [xDot0];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in x and xDot
        // over the step
        k1 = dt*f(alpha, beta, gamma, delta, omega, t[i], x[i], xDot[i])[0];
        l1 = dt*f(alpha, beta, gamma, delta, omega, t[i], x[i], xDot[i])[1];
        k2 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4)[0];
        l2 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/4, x[i]+k1/4, xDot[i]+l1/4)[1];
        k3 = dt*f(alpha, beta, gamma, delta, omega, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32)[0];
        l3 = dt*f(alpha, beta, gamma, delta, omega, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, xDot[i]+3*l1/32+9*l2/32)[1];
        k4 = dt*f(alpha, beta, gamma, delta, omega, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
        l4 = dt*f(alpha, beta, gamma, delta, omega, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, xDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
        k5 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
        l5 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, xDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
        k6 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
        l6 = dt*f(alpha, beta, gamma, delta, omega, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, xDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

        // x1 and xDot1 are our fourth order approximations
        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        xDot1 = xDot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        // x2 and xDot2 are our fifth order approximations
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        xDot2 = xDot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

        // The following are used to correct the step size
        Rx = Math.abs(x1-x2)/dt;
        RxDot = Math.abs(xDot1-xDot2)/dt;
        sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        sxDot = 0.84*Math.pow(epsilon/RxDot, 1/4);
        R = Math.max(Rx, RxDot);
        s = Math.min(sx, sxDot);
        if ( R <= epsilon ) {
            t.push(t[i]+dt);
            x.push(x1);
            xDot.push(xDot1);
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, x and xDot to our solution object
    solution = {
        t: t,
        x: x,
        xDot: xDot,
    };
}

/**
 * Tabulates solution data.
 *
 * @params           None. Uses the entries of the solution object, however. 
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable() {
    if ( solution.t.length == 0) {
        alert('You haven\'t solved the problem yet! Press the "Solve the problem" button before pressing the "Tabulate the solution.t button again.');
        return
    }
    t = solution.t;
    x = solution.x;
    xDot = solution.xDot;
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
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
 * @params           None.
 * @return           Nothing. Just generates the relevant plot.
 */
function generatePhasePlot() {
    // Run solveProblem() if previously unrun
    if ( solution.t.length == 0) {
        solveProblem();
    };

    // Extract solution data from solution object
    x = solution.x;
    xDot = solution.xDot;

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlot").style = "height: " + windowInnerHeight + "px;";

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
 * Generate plot of x and x dot against time
 * 
 * @params           None.
 * @return           Nothing. Just generates the relevant plot.
 */
function generateTimePlot() {
    // Run solveProblem() if previously unrun
    if ( solution.t.length == 0) {
        solveProblem();
    };

    // Extract solution data from solution object
    t = solution.t;
    x = solution.x;
    xDot = solution.xDot;

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";

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
    data = [plotx, plotxDot];

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
 * - one of xDot and x against t; and
 * - a phase plot of xDot against x.
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
