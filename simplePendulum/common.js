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
    var thetaMaxInitial = parseFloat(document.getElementById("thetaMaxInitial").value);
    var xi = parseFloat(document.getElementById("xi").value);
    var epsilon = parseFloat(document.getElementById("epsilon").value);
    var dtInitial = parseFloat(document.getElementById("dtInitial").value);
    var t0 = parseFloat(document.getElementById("t0").value);
    var tf = parseFloat(document.getElementById("tf").value);

    // Add them to array
    return [g, l, N, theta0, thetaDot0, thetaMaxInitial, xi, epsilon, dtInitial, t0, tf];

    // g:                 arrayOfInputs[0]
    // l:                 arrayOfInputs[1]
    // N:                 arrayOfInputs[2]
    // theta0:            arrayOfInputs[3]
    // thetaDot0:         arrayOfInputs[4]
    // thetaMaxInitial:   arrayOfInputs[5]
    // xi:                arrayOfInputs[6]
    // epsilon:           arrayOfInputs[7]
    // dtInitial:         arrayOfInputs[8]
    // t0:                arrayOfInputs[9]
    // tf:                arrayOfInputs[10]
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
 * @param thetaDot   Rate of change against time of the angle from the 
 * positive x-axis.
 * @return           [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, thetaDot) {
    return [thetaDot, -g/l*Math.cos(theta)];
};

/**
 * Returns the value of theta dot squared.
 * 
 * @param theta    Angle from the x-axis.
 * @return         Theta dot squared value.
 */
function thetaDotSq(g, l, theta0, thetaDot0, theta) {
    return thetaDot0**2 + 2*g/l * (Math.sin(theta0)-Math.sin(theta));
}

/**
 * The derivative with respect to theta of the theta dot squared expression
 * 
 * @param theta    Angle from the x-axis.
 * @return         Theta derivative of the theta dot squared expression.
 */
function thetaDotSqPrime(g, l, theta) {
    return -2*g/l * Math.cos(theta);
}

/**
 * Correction to theta according to Newton's method
 * 
 * @param theta    Angle from the x-axis.
 * @return         Correction to theta min/max value according to Newton's
 * method.
 */
function newtonsCorrection(g, l, theta0, thetaDot0, theta) {
    return -thetaDotSq(g, l, theta0, thetaDot0, theta)/thetaDotSqPrime(g, l, theta);
}

/**
 * Calculate the bounds of theta in our integration by using Newton's method
 * on our expression for theta dot squared. 
 * 
 * @param arrayOfInputs      An array of problem parameters.
 * @return                   Nothing. Just writes the relevant  
 */
function thetaBounds(arrayOfInputs) {
    // Parameters of the problem that are necessary to calculate the period
    var g = arrayOfInputs[0];
    var l = arrayOfInputs[1];
    var N = arrayOfInputs[2];
    var theta0 = arrayOfInputs[3];
    var thetaDot0 = arrayOfInputs[4];
    var thetaMaxInitial = arrayOfInputs[5];
    var xi = arrayOfInputs[6];

    // Take our initial guess for thetaMax
    var thetaMax = thetaMaxInitial;
    var thetaMin = theta0;
    // thetaMax correction
    var adjMax = newtonsCorrection(g, l, theta0, thetaDot0, thetaMax);
    var adjMin;

    // Calculate when thetaDot = 0 next, which will be halfway through the problem's period
    var j = 0; 
    var k = 0;
    while ( ( Math.abs(adjMax) >= xi) && ( j < N )) {
        thetaMax += adjMax;
        adjMax = newtonsCorrection(g, l, theta0, thetaDot0, thetaMax);
        j++;
    }
    thetaMax += adjMax;

    // Calculate thetaMin, which is when thetaDot = 0
    if (thetaDot0 != 0) {
        adjMin = newtonsCorrection(g, l, theta0, thetaDot0, thetaMin);    
        while (( Math.abs(adjMin) >= xi) & (k < N)) {
            thetaMin += adjMin;
            adjMin = newtonsCorrection(g, l, theta0, thetaDot0, thetaMin);
            k++;
        }
        thetaMin += adjMin;
    } else {
        thetaMin = theta0;
    }

    // Mention theta bounds and Newton's info in data table
    document.getElementById("jDisplay").innerHTML = j;
    document.getElementById("kDisplay").innerHTML = k;
    document.getElementById("thetaMinDisplay").innerHTML = thetaMin;
    document.getElementById("thetaMaxDisplay").innerHTML = thetaMax;
    document.getElementById("adjMinDisplay").innerHTML = adjMin;
    document.getElementById("adjMaxDisplay").innerHTML = adjMax;

    // If j and k reach N mention that in the table to let the user know
    // that the maximum iterations were reached
    if ( (j >= N ) || (k >= N ) ) {
        document.getElementById("integralDisplay").innerHTML = "Theta min and/or theta max could not be calculated as the limit on Newton's method iterations was exceeded";
        // Too high of tf in these cases can freeze the solver up
        document.getElementById("tf").value = 10;
        // The default, 1e-11, can freeze the webpage up for these problems
        // e.g. with thetaDot0 > 5, this will likely be the case.
        document.getElementById("epsilon").value = 1e-9;
        // Add j and k to array, so that we can kill periodCalc() if j/k>=N
        arrayOfInputs.push(j);
        arrayOfInputs.push(k);
    }

    // Add thetaMin and thetaMax to arrayOfInputs
    arrayOfInputs.push(thetaMin);
    arrayOfInputs.push(thetaMax);

    return arrayOfInputs;
}

/**
 * Calculates the value of theta when theta dot = 0 using Newton's method
 * then uses Chebyshev-Gauss quadrature to compute the time taken to reach this period.
 * 
 * @param arrayOfInputs  An array that contains all the problem parameters.
 * @return               Nothing. Changes the element value/innerHTML of 
 * integralDisplay and tf to T and 4*T, respectively.
 */
function periodCalc(arrayOfInputs) {
    // Parameters of the problem that are necessary to calculate the period
    var g = arrayOfInputs[0];
    var l = arrayOfInputs[1];
    var N = arrayOfInputs[2];
    var theta0 = arrayOfInputs[3];
    var thetaDot0 = arrayOfInputs[4];
    // Kill the function if j, k >= N
    if (arrayOfInputs[11] != "1e6") {
        var thetaMin = arrayOfInputs[11];
        var thetaMax = arrayOfInputs[12];
    } else {
        return;
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
        integrand = Math.sqrt(1-nodes**2)*Math.pow(thetaDotSq(g, l, theta0, thetaDot0, transformedGrid),-1/2);
        integral += ((thetaMax-thetaMin)/2) * (Math.PI/N)*integrand;
    }
    period = 2*Math.abs(integral);

    // Change what's displayed on the page accordingly
    document.getElementById("integralDisplay").innerHTML = period;
    document.getElementById("tf").value = 4*period;
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
 * @param arrayOfInputs  An array that contains all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(arrayOfInputs) {
    // Run solveProblem() if previously unrun
    var solution = solveProblem(arrayOfInputs);

    // Extract solution data from solution object
    var theta = solution["theta"];
    var thetaDot = solution["thetaDot"];

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
 * @param arrayOfInputs  An array that contains all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(arrayOfInputs) {
    var solution = solveProblem(arrayOfInputs);

    // Extract solution data from solution object
    var t = solution["t"];
    var theta = solution["theta"];
    var thetaDot = solution["thetaDot"];

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
 * Tabulates solution data.
 *
 * @param arrayOfInputs  An array that contains all the problem parameters.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(arrayOfInputs) {
    // Define all global variables
    var solution = solveProblem(arrayOfInputs);
    var epsilon = arrayOfInputs[7];
    var t = solution["t"];
    var theta = solution["theta"];
    var thetaDot = solution["thetaDot"];
    var j = 0;

    // Check the number of elements in solution
    for (m in solution) {
        j++;
    }

    // Add relevant vars if solution contains them
    if (j > 3) {
        var errorThetaDot = solution["errorThetaDot"];
        var rmsErrorthetaDot = Math.sqrt(math.dot(errorThetaDot, errorThetaDot)/(errorThetaDot.length-1));
    }
    
    // Generate table
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
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
 * @param arrayOfInputs  An array that contains all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(arrayOfInputs) {
    generateTimePlot(arrayOfInputs);
    generatePhasePlot(arrayOfInputs);
    if (!!document.getElementById("errorPlot")) {
        generateErrorPlot(arrayOfInputs);
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