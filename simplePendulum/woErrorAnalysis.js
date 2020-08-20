/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param g          Acceleration due to gravity in metres per second squared.
 * @param l          Length of the pendulum in metres.
 * @param t          Time (seconds).
 * @param theta      Angle of the pendulum from the positive x-axis (in
 * radians)
 * @param thetaDot   Rate of change against time of the angle from the positive x-axis.
 * @return           [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, thetaDot) {
    return [thetaDot, -g/l*Math.cos(theta)];
}

// Initialize our global variables
var solution = {
    t: [],
    theta: [],
    thetaDot: []
};
var T, epsilon;

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
    l = parseFloat(document.getElementById("l").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    theta0 = parseFloat(document.getElementById("theta0").value);
    thetaDot0 = parseFloat(document.getElementById("thetaDot0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // T is the period of the problem, including it can be used from the console to determine an appropriate tf
    if ( (theta0 == 0) && (thetaDot0 == 0) ) {
        T = Math.sqrt(l/(g*Math.PI))*(math.gamma(1/4)**2);
    }

    // Initialize the arrays used and loop variables
    t = [t0];
    theta = [theta0];
    thetaDot = [thetaDot0];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in theta and thetaDot
        // over the step
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
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, theta and thetaDot to our solution object
    solution = {
        t: t,
        theta: theta,
        thetaDot: thetaDot,
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
        alert('You haven\'t solved the problem yet! Press the "Solve the problem" button before pressing the "Tabulate the solution" button again.');
        return
    }
    t = solution["t"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>&theta; (radians) </th>';
    tableContents += '<th>&theta; dot <br/>(radians &middot; s<sup>-1</sup>)</th>';
    tableContents += "</tr>";
    for (let j = 0; j < theta.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
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
    document.getElementById('tableOutputs').innerHTML = '';
}

/**
 * Generate two plots:
 * - one of thetaDot and theta against t; and
 * - a phase plot of thetaDot against theta.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    if ( solution["t"].length == 0) {
        alert('You haven\'t solved the problem yet! Before pressing the "Plot the solution" button again, press the "Solve the problem" button.');
        return
    }
    t = solution["t"];
    theta = solution["theta"];
    thetaDot = solution["thetaDot"];

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("phasePlot").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the theta and thetaDot against time plot
    var plot1 = {
        x: t,
        y: theta,
        type: 'scatter',
        name: "theta (radians)"
    };
    var plot2 = {
        x: t,
        y: thetaDot,
        type: 'scatter',
        name: "theta dot (radians per second)"
    };
    var layout1 = {
        title: 'theta dot and theta against time plots',
        xaxis: {
           title: 'Time (seconds)'
        }
    };
    data1 = [plot1, plot2];

    // Characteristics of the phase plot
    var plot3 = {
        x: theta,
        y: thetaDot,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout2 = {
        title: "Phase plot of theta dot against theta",
        xaxis: {
            title: "theta (radians)"
        },
        yaxis: {
            title: "theta dot (radians per second)"
        }
    };
    data2 = [plot3];

    // Generate plots
    Plotly.newPlot('timePlot', data1, layout1);
    Plotly.newPlot('phasePlot', data2, layout2);
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