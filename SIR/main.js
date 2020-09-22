/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param beta       Interaction parameter.
 * @param gamma      Interaction parameter.
 * @param delta      Interaction parameter.
 * @param t          Time (seconds).
 * @param S          S value.
 * @param I          I value.
 * @param R          R value.
 * @return           [dx/dt, dy/dt, dz/dt]
 */
function f(beta, gamma, delta, t, S, I, R) {
    var N = S+I+R;
    return [-beta*S*I*(1-delta)/N, beta*S*I*(1-delta)/N-gamma*I*(1-delta), gamma*I*(1-delta)];
}

// Initialize our global variables
var solution = {
    t: [],
    S: [],
    I: [],
    R: []
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
    beta = parseFloat(document.getElementById("beta").value);
    gamma = parseFloat(document.getElementById("gamma").value);
    delta = parseFloat(document.getElementById("delta").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    S0 = parseFloat(document.getElementById("S0").value);
    I0 = parseFloat(document.getElementById("I0").value);
    R0 = parseFloat(document.getElementById("R0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    S = [S0];
    I = [I0];
    R = [R0];
    dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in S, I and R
        // over the step
        // 1st approx
        k1 = dt*f(beta, gamma, delta, t[i], S[i], I[i], R[i])[0];
        l1 = dt*f(beta, gamma, delta, t[i], S[i], I[i], R[i])[1];
        m1 = dt*f(beta, gamma, delta, t[i], S[i], I[i], R[i])[2];
        // 2nd approx
        k2 = dt*f(beta, gamma, delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4)[0];
        l2 = dt*f(beta, gamma, delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4)[1];
        m2 = dt*f(beta, gamma, delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4)[2];
        // 3rd approx
        k3 = dt*f(beta, gamma, delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32)[0];
        l3 = dt*f(beta, gamma, delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32)[1];
        m3 = dt*f(beta, gamma, delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32)[2];
        // 4th approx
        k4 = dt*f(beta, gamma, delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[0];
        l4 = dt*f(beta, gamma, delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[1];
        m4 = dt*f(beta, gamma, delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[2];
        // 5th approx
        k5 = dt*f(beta, gamma, delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[0];
        l5 = dt*f(beta, gamma, delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[1];
        m5 = dt*f(beta, gamma, delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[2];
        // 6th approx
        k6 = dt*f(beta, gamma, delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[0];
        l6 = dt*f(beta, gamma, delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[1];
        m6 = dt*f(beta, gamma, delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[2];

        // x1, y1 and z1 are our fourth order approximations
        x1 = S[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        y1 = I[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        z1 = R[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;

        // x2, y2 and z2 are our fifth order approximations
        x2 = S[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        y2 = I[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        z2 = R[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;

        // The following are used to correct the step size
        Rx = Math.abs(x1-x2)/dt;
        Ry = Math.abs(y1-y2)/dt;
        Rz = Math.abs(z1-z2)/dt;
        sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        sy = 0.84*Math.pow(epsilon/Ry, 1/4);
        sz = 0.84*Math.pow(epsilon/Rz, 1/4);
        RRKF45 = Math.max(Rx, Ry, Rz);
        s = Math.min(sx, sy, sz);

        // If R is less than or equal to epsilon move onto the next step
        if ( RRKF45 <= epsilon ) {
            t.push(t[i]+dt);
            S.push(x1);
            I.push(y1);
            R.push(z1);
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, S, I and R to our solution object
    solution = {
        t: t,
        S: S,
        I: I,
        R: R
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
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extract coordinate arrays from the solution object
    t = solution.t;
    S = solution.S;
    I = solution.I;
    R = solution.R;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>S</th>';
    tableContents += '<th>I</th>';
    tableContents += '<th>R</th>';
    tableContents += "</tr>";
    for (let j = 0; j < S.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + S[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + I[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + R[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
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
 * Generates a 3D phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generate3DPhasePlot() {
    // Run solveProblem if unrun
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extra solution data from solution object
    S = solution.S;
    I = solution.I;
    R = solution.R;

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXYZ").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXYZ = {
        x: S,
        y: I,
        z: R,
        type: 'scatter3d',
        mode: 'lines',
        opacity: 1,
        line: {
            width: 6,
            reversescale: false
        }
    };
    dataXYZ = [plotXYZ];

    // layout object
    var layoutXYZ = {
        title: 'Phase plot of the solution to the SIR equations. x = S, y = I, z = R'
    };

    // Generate plot
    Plotly.newPlot('phasePlotXYZ', dataXYZ, layoutXYZ);
}

/**
 * Remove 3D phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function remove3DPhasePlot() {
    document.getElementById("phasePlotXYZ").innerHTML = '';
    document.getElementById("phasePlotXYZ").style = '';
}

/**
 * Generates a XY phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXYPhasePlot() {
    // Run solveProblem if unrun
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extra solution data from solution object
    S = solution.S;
    I = solution.I;

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXY").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXY = {
        x: S,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    dataXY = [plotXY];

    // layout object
    var layoutXY = {
        title: "SI phase plot, x = S and y = I"
    };

    // Generate plot
    Plotly.newPlot('phasePlotXY', dataXY, layoutXY);
}

/**
 * Remove XY phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXYPhasePlot() {
    document.getElementById("phasePlotXY").innerHTML = '';
    document.getElementById("phasePlotXY").style = '';
}

/**
 * Generates a XZ phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXZPhasePlot() {
    // Run solveProblem if unrun
    if ( solution.t.length == 0) {
        solveProblem();
    }
    
    // Extra solution data from solution object
    S = solution.S;
    R = solution.R;
    
    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXZ").style = "height: " + windowInnerHeight + "px;";
    
    // Plot object and data object array
    var plotXZ = {
        x: S,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    dataXZ = [plotXZ];
    
    // layout object
    var layoutXZ = {
        title: "SR phase plot, x = S and y = R"
    };
    
    // Generate plot
    Plotly.newPlot('phasePlotXZ', dataXZ, layoutXZ);
}

/**
 * Remove XZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXZPhasePlot() {
    document.getElementById("phasePlotXZ").innerHTML = '';
    document.getElementById("phasePlotXZ").style = '';
}

/**
 * Generates a YZ phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateYZPhasePlot() {
    // Run solveProblem if unrun
    if ( solution.t.length == 0) {
        solveProblem();
    }

    // Extra solution data from solution object
    I = solution.I;
    R = solution.R;

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotYZ").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotYZ = {
        x: I,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    dataYZ = [plotYZ];

    // layout object
    var layoutYZ = {
        title: "IR phase plot, x = I and y = R"
    };

    // Generate plot
    Plotly.newPlot('phasePlotYZ', dataYZ, layoutYZ);
}

/**
 * Remove YZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeYZPhasePlot() {
    document.getElementById("phasePlotYZ").innerHTML = '';
    document.getElementById("phasePlotYZ").style = '';
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

    // Extra solution data from solution object
    t = solution.t;
    S = solution.S;
    I = solution.I;
    R = solution.R;

    // Height and width of plot
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotTX = {
        x: t,
        y: S,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'S'
    };
    var plotTY = {
        x: t,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'I'
    };
    var plotTZ = {
        x: t,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'R'
    };
    dataTimePlot = [plotTX, plotTY, plotTZ];

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
 * - The first is a 3D phase plot of S, I and R.
 * - The second is a 2D phase plot of I against S.
 * - The third is a 2D phase plot of R against S.
 * - The fourth is a 2D phase plot of R against I.
 * - The fifth is a plot of S, I and R against time.
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