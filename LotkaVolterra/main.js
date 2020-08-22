/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param alpha      Interaction parameter.
 * @param beta       Interaction parameter.
 * @param gamma      Interaction parameter.
 * @param delta      Interaction parameter.
 * @param t          Time (seconds).
 * @param x          Prey population.
 * @param y          Predator population.
 * @return           [dx/dt, dy/dt]
 */
function f(alpha, beta, gamma, delta, t, x, y) {
    return [alpha*x -beta*x*y, delta*x*y-gamma*y];
}

// Initialize our global variables
var solution = {
    t: [],
    x: [],
    y: []
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
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    x0 = parseFloat(document.getElementById("x0").value);
    y0 = parseFloat(document.getElementById("y0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    x = [x0];
    y = [y0];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in x and y
        // over the step
        k1 = dt*f(alpha, beta, gamma, delta, t[i], x[i], y[i])[0];
        l1 = dt*f(alpha, beta, gamma, delta, t[i], x[i], y[i])[1];
        k2 = dt*f(alpha, beta, gamma, delta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4)[0];
        l2 = dt*f(alpha, beta, gamma, delta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4)[1];
        k3 = dt*f(alpha, beta, gamma, delta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32)[0];
        l3 = dt*f(alpha, beta, gamma, delta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32)[1];
        k4 = dt*f(alpha, beta, gamma, delta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
        l4 = dt*f(alpha, beta, gamma, delta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
        k5 = dt*f(alpha, beta, gamma, delta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
        l5 = dt*f(alpha, beta, gamma, delta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
        k6 = dt*f(alpha, beta, gamma, delta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
        l6 = dt*f(alpha, beta, gamma, delta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

        // x1 and y1 are our fourth order approximations
        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        // x2 and y2 are our fifth order approximations
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

        // The following are used to correct the step size
        Rx = Math.abs(x1-x2)/dt;
        Ry = Math.abs(y1-y2)/dt;
        sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        sy = 0.84*Math.pow(epsilon/Ry, 1/4);
        R = Math.max(Rx, Ry);
        s = Math.min(sx, sy);
        if ( R <= epsilon ) {
            t.push(t[i]+dt);
            x.push(x1);
            y.push(y1);
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, x and y to our solution object
    solution = {
        t: t,
        x: x,
        y: y,
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
        solveProblem();
    }

    // Extra solution data from solution object
    t = solution["t"];
    x = solution["x"];
    y = solution["y"];

    // Clear innerHTML in case the table already exists and we're replacing 
    // it
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>Prey population</th>';
    tableContents += '<th>Predator population</th>';
    tableContents += "</tr>";
    for (let j = 0; j < x.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + x[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + y[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
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
 * Generate phase plot
 * 
 * @params           None.
 * @return           Nothing. Just generates the phase plot.
 */
function generatePhasePlot() {
    // If solveProblem() hasn't been run, re-run it. 
    if ( solution["t"].length == 0) {
        solveProblem();
    }

    // Extract solution data from solution object
    t = solution["t"];
    x = solution["x"];
    y = solution["y"];

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlot").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the phase plot
    var plot = {
        x: x,
        y: y,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout = {
        title: "Phase plot of y against x",
        xaxis: {
            title: "x"
        },
        yaxis: {
            title: "y"
        }
    };
    data = [plot];

    // Generate plot
    Plotly.newPlot('phasePlot', data, layout);    
}

/**
 * Generate a plot of x and y against time
 * 
 * @params           None.
 * @return           Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot() {
    // Solve the problem if it hasn't been already
    if ( solution["t"].length == 0) {
        solveProblem();
    }

    // Extract solution data from solution object
    t = solution["t"];
    x = solution["x"];
    y = solution["y"];

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the x and y against time plot
    var plotx = {
        x: t,
        y: x,
        type: 'scatter',
        name: "x"
    };
    var ploty = {
        x: t,
        y: y,
        type: 'scatter',
        name: "y"
    };
    var layout = {
        title: 'y and x against time plots',
        xaxis: {
           title: 'Time (seconds)'
        }
    };
    data = [plotx, ploty];

    // Generate plots
    Plotly.newPlot('timePlot', data, layout);
}

/**
 * Generate two plots:
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    generatePhasePlot();
    generateTimePlot();
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