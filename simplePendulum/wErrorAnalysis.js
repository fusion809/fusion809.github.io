/**
 * Right-hand side of our second-order ODE written as a system of 
 * first-order ODEs.
 *
 * @param g      Acceleration due to gravity in metres per second 
 * squared.
 * @param l      Length of the pendulum in metres.
 * @param t      Time.
 * @param theta  Angle of the pendulum from the positive x-axis.
 * @param dtheta Rate of change against time of the angle from the 
 * positive x-axis.
 * @return       [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, dtheta) {
    return [dtheta, -g/l*Math.cos(theta)];
};

// Initialize our global variables
var solution = {
    t: [],
    theta: [],
    dtheta: [],
    errorDtheta: [],
    logErrorDtheta: [],
};
var windowInnerWidth;
var windowInnerHeight;
var epsilon;

/** 
 * Solve the problem using RK45.
 *
 * @params       None. Uses parameter values in the forum.
 * @return       Nothing. But it enters the solution values and error estimates
 * into the solution object.
 */
function solveProblem() {
    // Obtain the parameters of the problem
    g = parseFloat(document.getElementById("g").value);
    l = parseFloat(document.getElementById("l").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    theta0 = parseFloat(document.getElementById("theta0").value);
    dtheta0 = parseFloat(document.getElementById("dtheta0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    theta = [theta0];
    dtheta = [dtheta0];
    errorDtheta = [0];
    logErrorDtheta = [-Infinity];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in theta and 
        // dtheta over the step
        k1 = dt*f(g, l, t[i], theta[i], dtheta[i])[0];
        l1 = dt*f(g, l, t[i], theta[i], dtheta[i])[1];
        k2 = dt*f(g, l, t[i]+dt/4, theta[i]+k1/4, dtheta[i]+l1/4)[0];
        l2 = dt*f(g, l, t[i]+dt/4, theta[i]+k1/4, dtheta[i]+l1/4)[1];
        k3 = dt*f(g, l, t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, dtheta[i]+3*l1/32+9*l2/32)[0];
        l3 = dt*f(g, l, t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, dtheta[i]+3*l1/32+9*l2/32)[1];
        k4 = dt*f(g, l, t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dtheta[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
        l4 = dt*f(g, l, t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dtheta[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
        k5 = dt*f(g, l, t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, dtheta[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
        l5 = dt*f(g, l, t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, dtheta[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
        k6 = dt*f(g, l, t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dtheta[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
        l6 = dt*f(g, l, t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dtheta[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

        // theta1 and dtheta1 are our fourth order approximations
        theta1 = theta[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        dtheta1 = dtheta[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        // theta2 and dtheta2 are our fifth order approximations
        theta2 = theta[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        dtheta2 = dtheta[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        errorDtheta1 = Math.abs(dtheta1 - Math.sign(dtheta1)*Math.sqrt(Math.abs(dtheta0**2+2*g*(Math.sin(theta0) - Math.sin(theta1))/l)));

        // The following are used to correct the step size
        Rtheta = Math.abs(theta1-theta2)/dt;
        Rdtheta = Math.abs(dtheta1-dtheta2)/dt;
        stheta = 0.84*Math.pow(epsilon/Rtheta, 1/4);                
        sdtheta = 0.84*Math.pow(epsilon/Rdtheta, 1/4);
        R = Math.max(Rtheta, Rdtheta);
        s = Math.min(stheta, sdtheta);
        if ( R <= epsilon ) {
            t.push(t[i]+dt);
            theta.push(theta1);
            dtheta.push(dtheta1);
            errorDtheta.push(errorDtheta1);
            logErrorDtheta.push(Math.log10(errorDtheta1));
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
        dtheta: dtheta,
        errorDtheta: errorDtheta,
        logErrorDtheta: logErrorDtheta,
    };
};

/**
 * Tabulates solution data.
 *
 * @params       None.
 * @return       Nothing. Just populates the table with the solution values. 
 */
function fillTable() {
    if ( solution["t"].length == 0) {
        alert('You haven\'t solved the problem yet! Press the "Solve the problem" button before pressing the "Tabulate the solution" button again.');
        return
    }
    t = solution["t"];
    theta = solution["theta"];
    dtheta = solution["dtheta"];
    errorDtheta = solution["errorDtheta"];
    rmsErrorDtheta = Math.sqrt(math.dot(errorDtheta, errorDtheta)/(errorDtheta.length-1));
    logErrorDtheta = solution["logErrorDtheta"];
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>Time (seconds)</th>';
    tableContents += '<th>&theta; (radians) </th>';
    tableContents += '<th>&theta; dot (radians &middot; s<sup>-1</sup>)</th>';
    tableContents += '<th>Error in &theta; dot</th>';
    tableContents += '</tr>';
    for (let j = 0; j < theta.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + theta[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + dtheta[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + errorDtheta[j].toExponential(10) + '</td>';
        tableContents += '</tr>';
    }
    tableContents += '<tr>';
    tableContents += '<td colspan = "6">RMS error in &theta; dot is: ';
    tableContents += rmsErrorDtheta.toExponential(10);
    tableContents += '</td>';
    tableContents += '</tr>';
    document.getElementById('tableOutputs').innerHTML = tableContents;
};

/**
 * Removes the solution table
 * 
 * @params       None.
 * @return       Nothing. Just removes the solution table.
 */
function removeTable() {
    document.getElementById('tableOutputs').innerHTML = '';
}

/**
 * Generate three plots:
 * - one of dtheta and theta against t;
 * - one of the log of the difference between dtheta and that calculated
 * based on theta; and
 * - a phase plot of dtheta against theta.
 * 
 * @params       None.
 * @return       Nothing. Just generates the plots.
 */
function generatePlots() {
    if ( solution["t"].length == 0) {
        alert('You haven\'t solved the problem yet! Before pressing the "Plot the solution" button again, press the "Solve the problem" button.');
        return
    };
    t = solution["t"];
    theta = solution["theta"];
    dtheta = solution["dtheta"];

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("errorPlot").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("phasePlot").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the theta and theta dot against time plot
    var plot1 = {
        x: t,
        y: theta,
        type: 'scatter',
        name: "theta (radians)"
    };
    var plot2 = {
        x: t,
        y: dtheta,
        type: 'scatter',
        name: "theta dot (radians per second)"
    };
    var layout1 = {
        title: 'theta and theta dot against time plots',
        xaxis: {
            title: 'Time (seconds)'
        }
    };
    data1 = [plot1, plot2];

    // Logarithmic plot of the error in dtheta
    var plot3 = {
        x: t,
        y: logErrorDtheta,
        type: 'scatter',
        name: 'Semilog plot of error in theta dot',
    };
    data2 = [plot3];
    var layout2 = {
        title: "Semilog plot of the error in theta dot against t",
        xaxis: {
            title: "Time (seconds)"
        },
        yaxis: {
            title: "Log of the error in theta dot"
        }
    };
    
    // Characteristics of the phase plot
    var plot4 = {
        x: theta,
        y: dtheta,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout3 = {
        title: "Phase plot of theta dot against theta",
        xaxis: {
            title: "theta (radians)"
        },
        yaxis: {
            title: "theta dot (radians per second)"
        }
    };
    data3 = [plot4];

    // Generate plots
    Plotly.newPlot('timePlot', data1, layout1);
    Plotly.newPlot('errorPlot', data2, layout2);
    Plotly.newPlot('phasePlot', data3, layout3);
};

/**
 * Removes solution plots and associated element formatting
 * 
 * @params       None.
 * @return       Nothing. Just removes the solution plots.
 */
function removePlots() {
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("errorPlot").innerHTML = '';
    document.getElementById("phasePlot").innerHTML = '';
    document.getElementById("timePlot").style = ''
    document.getElementById("errorPlot").style = '';
    document.getElementById("phasePlot").style = '';
}