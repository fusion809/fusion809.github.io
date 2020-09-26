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

/**
 * Approximate the value of x and y at next step size to 4th and 5th order
 * @param dt            Step size.
 * @param alpha         Problem parameter.
 * @param beta          Problem parameter.
 * @param gamma         Problem parameter.
 * @param delta         Problem parameter.
 * @param t             An array of previous t values.
 * @param x             An array of x values.
 * @param y             An array of y values.
 * @param i             Counter variable.
 * @return              [x1, y1, x2, y2]
 */
function approximatorRKF45(dt, alpha, beta, gamma, delta, t, x, y, i) {
    // Runge-Kutta-Fehlberg approximations of the change in x and y over the step
    var k1 = dt*f(alpha, beta, gamma, delta, t[i], x[i], y[i])[0];
    var l1 = dt*f(alpha, beta, gamma, delta, t[i], x[i], y[i])[1];
    var k2 = dt*f(alpha, beta, gamma, delta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4)[0];
    var l2 = dt*f(alpha, beta, gamma, delta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4)[1];
    var k3 = dt*f(alpha, beta, gamma, delta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32)[0];
    var l3 = dt*f(alpha, beta, gamma, delta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32)[1];
    var k4 = dt*f(alpha, beta, gamma, delta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
    var l4 = dt*f(alpha, beta, gamma, delta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
    var k5 = dt*f(alpha, beta, gamma, delta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
    var l5 = dt*f(alpha, beta, gamma, delta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
    var k6 = dt*f(alpha, beta, gamma, delta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
    var l6 = dt*f(alpha, beta, gamma, delta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

    // x1 and y1 are our fourth order approximations
    var x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    var y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    // x2 and y2 are our fifth order approximations
    var x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    var y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

    return [x1, y1, x2, y2];
}

/**
 * Check and adjust the step size
 * 
 * @param dt            Step size.
 * @param epsilon       Error tolerance.
 * @param t             An array of t values.
 * @param x             An array of x values.
 * @param y             An array of y values.
 * @param x1            4th order approximation to next x value.
 * @param y1            4th order approximation to next y value.
 * @param x2            5th order approximation to next x value.
 * @param y2            5th order approximation to next y value.
 * @param i             Counter variable.
 * @return              [dt, t, x, y, i]
 */
function stepSizeChecker(dt, epsilon, t, x, y, x1, y1, x2, y2, i) {
    // Initialize variables
    var Rx = Math.abs(x1-x2)/dt;
    var Ry = Math.abs(y1-y2)/dt;
    var sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
    var sy = 0.84*Math.pow(epsilon/Ry, 1/4);
    var R = Math.max(Rx, Ry);
    var s = Math.min(sx, sy);

    // Check if R is within our margin of error and move on if it is
    if ( R <= epsilon ) {
        t.push(t[i]+dt);
        x.push(x1);
        y.push(y1);
        i++;
        dt *= s;
    } else {
        dt *= s;
    }

    // Return what needs to be used later in solveProblem()
    return [dt, t, x, y, i];
}

/** 
 * Solve the problem using RKF45.
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               A solution object.
 */
function solveProblem(objectOfInputs) {
    // Obtain the parameters of the problem
    var {alpha, beta, gamma, delta, t0, tf, x0, y0, epsilon, dtInitial} = objectOfInputs;

    // Initialize the arrays used and loop variables
    var t = [t0];
    var x = [x0];
    var y = [y0];
    var dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        dt = Math.min(dt, tf-t[i]);
        var [x1, y1, x2, y2] = approximatorRKF45(dt, alpha, beta, gamma, delta, t, x, y, i);
        var [dt, t, x, y, i] = stepSizeChecker(dt, epsilon, t, x, y, x1, y1, x2, y2, i);
    }

    // Write t, x and y to our solution object
    var solution = {
        t: t,
        x: x,
        y: y,
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
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {t, x, y} = solution;

    // Clear innerHTML in case the table already exists and we're replacing 
    // it
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
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
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the phase plot.
 */
function generatePhasePlot(objectOfInputs) {
    // Solve the problem 
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {t, x, y} = solution;

    // Height and width of plot
    adjustPlotHeight("phasePlot");

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
 * Generate a plot of x and y against time
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {t, x, y} = solution;

    // Height and width of plots
    adjustPlotHeight("timePlot");

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
    var data = [plotx, ploty];

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
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generatePhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    removePhasePlot();
    removeTimePlot();
};