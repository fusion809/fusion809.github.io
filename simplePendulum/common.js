/**
 * Reads inputs in the parameters table, adds them to an array and return
 * 
 * @params         None.
 * @return         Array of problem parameters.
 */
function readInputs() {
    // Read them from form
    var g = parseFloat(document.getElementById("g").value);
    var l = parseFloat(document.getElementById("l").value);
    var N = parseFloat(document.getElementById("N").value);
    var theta0 = parseFloat(document.getElementById("theta0").value);
    var thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    var epsilon = parseFloat(document.getElementById("epsilon").value);
    var dtInitial = parseFloat(document.getElementById("dtInitial").value);
    var t0 = parseFloat(document.getElementById("t0").value);
    var tf = parseFloat(document.getElementById("tf").value);

    // Add them to array
    var objectOfInputs = {
        g: g,
        l: l,
        N: N,
        theta0: theta0,
        thetaDot0: thetaDot0,
        epsilon: epsilon,
        dtInitial: dtInitial,
        t0: t0,
        tf: tf
    };
    return objectOfInputs;
}

/**
 * Right-hand side of our second-order ODE written as a system of 
 * first-order ODEs.
 *
 * @param g          Acceleration due to gravity in metres per second 
 * squared.
 * @param l          Length of the pendulum in metres.
 * @param t          Time.
 * @param theta      Angle of the pendulum from the positive x-axis.
 * @param thetaDot   Rate of change against time of the angle from the positive x-axis.
 * @return           [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, thetaDot) {
    return [thetaDot, -g/l*Math.cos(theta)];
};

/**
 * Returns the value of theta dot squared.
 * 
 * @param g        Acceleration due to gravity.
 * @param l        Length of the pendulum.
 * @param theta0   Initial angle from the positive x-axis.
 * @param thetaDot0 Initial rate of change in the angle from the positive x-axis.
 * @param theta    Angle from the positive x-axis.
 * @return         Theta dot squared value.
 */
function thetaDotSq(g, l, theta0, thetaDot0, theta) {
    return thetaDot0**2 + 2*g/l * (Math.sin(theta0)-Math.sin(theta));
}

/**
 * The derivative with respect to theta of the theta dot squared expression
 * 
 * @param g        Acceleration due to gravity.
 * @param l        Length of the pendulum.
 * @param theta    Angle from the x-axis.
 * @return         Theta derivative of the theta dot squared expression.
 */
function thetaDotSqPrime(g, l, theta) {
    return -2*g/l * Math.cos(theta);
}

/**
 * Correction to theta according to Newton's method
 * 
 * @param g        Acceleration due to gravity.
 * @param l        Length of the pendulum rod.
 * @param theta0   Initial angle from the positive x-axis.
 * @param thetaDot0 Initial rate of change against time of the angle from the positive x-axis.
 * @param theta    Angle from the x-axis.
 * @return         Correction to theta min/max value according to Newton's  method.
 */
function newtonsCorrection(g, l, theta0, thetaDot0, theta) {
    return -thetaDotSq(g, l, theta0, thetaDot0, theta)/thetaDotSqPrime(g, l, theta);
}

/**
 * Calculate the bounds of theta in our integration by using Newton's method
 * on our expression for theta dot squared. 
 * 
 * @param objectOfInputs      An object of problem parameters.
 * @return                    Same array with the bounds for integration on theta added to the end.
 */
function thetaBounds(objectOfInputs) {
    // Parameters of the problem that are necessary to calculate the period
    var g = objectOfInputs.g;
    var l = objectOfInputs.l;
    var theta0 = objectOfInputs.theta0;
    var thetaDot0 = objectOfInputs.thetaDot0;
    var t0 = objectOfInputs.t0;

    // Take our initial guess for thetaMax
    // Check if the problem satisfies the conditions for periodic behaviour
    if ( thetaDot0**2 + 2*g/l * (Math.sin(theta0)-1) > 0) {
        document.getElementById("integralDisplay").innerHTML = "Problem is not periodic.";
        // Too high of tf in these cases can freeze the solver up
        document.getElementById("tf").value = t0 + 10;
        // The default, 1e-11, can freeze the webpage up for these problems
        // e.g. with thetaDot0 > 5, this will likely be the case.
        document.getElementById("epsilon").value = 1e-9;
        // Add extra N to array, so that we can kill periodCalc() if we get 
        // this far
        objectOfInputs.periodic = false;
        // Clearing vars made irrelevant by non-periodicity
        document.getElementById("thetaMinDisplay").innerHTML = '';
        document.getElementById("thetaMaxDisplay").innerHTML = '';
        // Return objectOfInputs
        return objectOfInputs;
    }
    objectOfInputs.periodic = true;
    var thetaMin = Math.asin((thetaDot0**2+2*(g/l)*Math.sin(theta0))/(2*g/l));
    var thetaMax = - thetaMin - Math.PI;

    // Mention theta bounds and Newton's info in data table
    document.getElementById("thetaMinDisplay").innerHTML = thetaMin;
    document.getElementById("thetaMaxDisplay").innerHTML = thetaMax;

    // Add thetaMin and thetaMax to objectOfInputs
    objectOfInputs.thetaMin = thetaMin;
    objectOfInputs.thetaMax = thetaMax;

    return objectOfInputs;
}

/**
 * Calculates the value of theta when theta dot = 0 using Newton's method
 * then uses Chebyshev-Gauss quadrature to compute the time taken to reach this period.
 * 
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return                Nothing. Changes the element value/innerHTML of integralDisplay and tf to T and 4*T, respectively.
 */
function periodCalc(objectOfInputs) {
    // Parameters of the problem that are necessary to calculate the period
    var g = objectOfInputs.g;
    var l = objectOfInputs.l;
    var N = objectOfInputs.N;
    var t0 = objectOfInputs.t0;
    var theta0 = objectOfInputs.theta0;
    var thetaDot0 = objectOfInputs.thetaDot0;
    var nodes = 0;
    var integrand = 0;
    var transformedGrid = 0;
    var integral = 0;
    var period;

    // Kill the function if problem isn't periodic
    if (objectOfInputs.periodic) {
        var thetaMin = objectOfInputs.thetaMin;
        var thetaMax = objectOfInputs.thetaMax;
    } else {
        return;
    }

    // Integrate problem from thetaMin to thetaMax to calculate the period
    // in seconds
    for ( let i = 1; i < N+1; i++) {
        nodes = Math.cos((2*i-1)*Math.PI/(2*N));
        transformedGrid = (thetaMax-thetaMin)*nodes/2+(thetaMax+thetaMin)/2;
        integrand = Math.sqrt(1-nodes**2)*Math.pow(thetaDotSq(g, l, theta0, thetaDot0, transformedGrid),-1/2);
        integral += ((thetaMax-thetaMin)/2) * (Math.PI/N)*integrand;
    }
    period = 2*Math.abs(integral);

    // Change what's displayed on the page accordingly
    document.getElementById("integralDisplay").innerHTML = period;
    if ( period == Infinity ) {
        document.getElementById("tf").value = t0 + 100;
    } else {
        document.getElementById("tf").value = t0 + 4*period;
    }
    return period;
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
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(objectOfInputs) {
    // Run solveProblem() if previously unrun
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var theta = solution.theta;
    var thetaDot = solution.thetaDot;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
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
    document.getElementById("phasePlot").innerHTML = '';
    document.getElementById("phasePlot").style = '';
}

/**
 * Generate plot of theta and theta dot against time
 * 
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(objectOfInputs) {
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var t = solution.t;
    var theta = solution.theta;
    var thetaDot = solution.thetaDot;

    // Height and width of plots
    var windowInnerHeight = window.innerHeight;
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
    var data = [plotTheta, plotThetaDot];

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
 * Tabulates solution data.
 *
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs) {
    // Define all global variables
    var solution = solveProblem(objectOfInputs);
    var epsilon = objectOfInputs.epsilon;
    var t = solution.t;
    var theta = solution.theta;
    var thetaDot = solution.thetaDot;
    var j = 0;

    // Check the number of elements in solution
    for (m in solution) {
        j++;
    }

    // Add relevant vars if solution contains them
    if (j > 3) {
        var errorThetaDot = solution.errorThetaDot;
        var rmsErrorthetaDot = Math.sqrt(math.dot(errorThetaDot, errorThetaDot)/(errorThetaDot.length-1));
    }
    
    // Generate table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>&theta; (radians) </th>';
    tableContents += '<th>&theta; dot<br/>(radians &middot; s<sup>-1</sup>)</th>';
    if (j > 3) {
        tableContents += '<th>&theta; dot error</th>';
    }
    tableContents += '</tr>';
    for (let i = 0; i < theta.length; i++) {
        tableContents += '<tr>';
        tableContents += '<td>' + i + '</td>';
        tableContents += '<td>' + t[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + theta[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + thetaDot[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        if (j > 3) {
            tableContents += '<td>' + errorThetaDot[i].toExponential(5) + '</td>';
        }
        tableContents += '</tr>';
    }
    if (j > 3) {
        tableContents += '<tr>';
        tableContents += '<td colspan = "6">RMS error in &theta; dot is: ';
        tableContents += rmsErrorthetaDot.toExponential(10);
        tableContents += '</td>';
        tableContents += '</tr>';
    }
    document.getElementById('tableOutputs').innerHTML = tableContents;
};

/**
 * Generate two plots:
 * - one of thetaDot and theta against t;
 * - a phase plot of thetaDot against theta.
 * 
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generateTimePlot(objectOfInputs);
    generatePhasePlot(objectOfInputs);
    // The following if statement is to ensure that errorPlot is only generated if it's present in HTML
    if (!!document.getElementById("errorPlot")) {
        generateErrorPlot(objectOfInputs);
    }
};

/**
 * Removes solution plots and associated element formatting
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("phasePlot").innerHTML = '';
    document.getElementById("timePlot").style = ''
    document.getElementById("phasePlot").style = '';
    if (!!document.getElementById("errorPlot")) {
        document.getElementById("errorPlot").innerHTML = '';
        document.getElementById("errorPlot").style = '';
    }
}
