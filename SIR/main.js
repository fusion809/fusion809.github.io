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
    // Determine N
    var N = S+I+R;
    // Calculate derivatives
    var dSdt = - beta * S * I * (1-delta)/N;
    var dIdt = beta * S * I * (1-delta)/N - gamma * I * (1-delta);
    var dRdt = gamma*I*(1-delta);
    // Put into return value
    return [dSdt, dIdt, dRdt];
}

/**
 * Calculates RKF45 approximations for S, I and R values for next step.
 * 
 * @param dt            Step size.
 * @param beta          Problem parameter.
 * @param gamma         Problem parameter.
 * @param delta         Problem parameter.
 * @param t             Array of t values.
 * @param S             Array of S values.
 * @param I             Array of I values.
 * @param R             Array of R values.
 * @param i             Loop counter index.
 * @return              An array of 4th and 5th order approximations to S, I and R at time next step.
 */
function approximatorRKF45(dt, beta, gamma, delta, t, S, I, R, i) {
    // Runge-Kutta-Fehlberg approximations of the change in S, I and R
    // over the step
    // 1st approx
    var k1 = dt*f(beta, gamma, delta, t[i], S[i], I[i], R[i])[0];
    var l1 = dt*f(beta, gamma, delta, t[i], S[i], I[i], R[i])[1];
    var m1 = dt*f(beta, gamma, delta, t[i], S[i], I[i], R[i])[2];
    // 2nd approx
    var k2 = dt*f(beta, gamma, delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4)[0];
    var l2 = dt*f(beta, gamma, delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4)[1];
    var m2 = dt*f(beta, gamma, delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4)[2];
    // 3rd approx
    var k3 = dt*f(beta, gamma, delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32)[0];
    var l3 = dt*f(beta, gamma, delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32)[1];
    var m3 = dt*f(beta, gamma, delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32)[2];
    // 4th approx
    var k4 = dt*f(beta, gamma, delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[0];
    var l4 = dt*f(beta, gamma, delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[1];
    var m4 = dt*f(beta, gamma, delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[2];
    // 5th approx
    var k5 = dt*f(beta, gamma, delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[0];
    var l5 = dt*f(beta, gamma, delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[1];
    var m5 = dt*f(beta, gamma, delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[2];
    // 6th approx
    var k6 = dt*f(beta, gamma, delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[0];
    var l6 = dt*f(beta, gamma, delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[1];
    var m6 = dt*f(beta, gamma, delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[2];

    // x1, y1 and z1 are our fourth order approximations
    var x1 = S[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    var y1 = I[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    var z1 = R[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;

    // x2, y2 and z2 are our fifth order approximations
    var x2 = S[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    var y2 = I[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
    var z2 = R[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;
    
    return [x1, y1, z1, x2, y2, z2];
}

/**
 * Step size checker and adapter function
 * @param dt            Current step size.
 * @param epsilon       Error tolerance.
 * @param t             An array of t values.
 * @param S             An array of previous S values against time.
 * @param I             An array of previous I values against time.
 * @param R             An array of previous R values against time.
 * @param x1            4th order approx for S at next step 
 * @param y1            4th order approx for I at next step
 * @param z1            4th order approx for R at next step
 * @param x2            5th order approx for S at next step
 * @param y2            5th order approx for I at next step
 * @param z2            5th order approx for R at next step
 * @param i             Counter variable value.
 * @return              Corrected dt and updated t, S, I, R and i.
 */
function stepSizeChecker(dt, epsilon, t, S, I, R, x1, y1, z1, x2, y2, z2, i) {
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

    return [dt, t, S, I, R, i];
}

/**
 * Runge-Kutta-Fehlberg 4/5th order integrator function
 * 
 * @param dtInitial     Initial guess for dt.
 * @param epsilon       Error tolerance.
 * @param beta          SIR parameter.
 * @param gamma         SIR parameter.
 * @param delta         SIR parameter.
 * @param t0            Day the simulation starts
 * @param tf            Day the simulation ends
 * @param S0            Initial number of susceptible persons.
 * @param I0            Initial number of infected persons.
 * @param R0            Initial number of recovered persons.
 * @return              Solution object containing arrays of time, S, I and R values.
 */
function RKF45(dtInitial, epsilon, beta, gamma, delta, t0, tf, S0, I0, R0) {
    // Initialize the arrays used and loop variables
    var t = [t0];
    var S = [S0];
    var I = [I0];
    var R = [R0];
    var dt = dtInitial;
    var i = 0;
    
    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Use approximatorRKF45 to make approximations
        var [x1, y1, z1, x2, y2, z2] = approximatorRKF45(dt, beta, gamma, delta, t, S, I, R, i);
    
        // Step size correction
        [dt, t, S, I, R, i] = stepSizeChecker(dt, epsilon, t, S, I, R, x1, y1, z1, x2, y2, z2, i);
    }
    
    // Write t, S, I and R to our solution object
    var solution = {
        t: t,
        S: S,
        I: I,
        R: R
    };

    return solution;
}

/**
 * Read inputs from the form and return them in an array
 * 
 * @params              Nothing.
 * @return              An array containing all the form inputs.
 */
function readInputs() {
    // Take parameter values from the form
    var beta = parseFloat(document.getElementById("beta").value);
    var gamma = parseFloat(document.getElementById("gamma").value);
    var delta = parseFloat(document.getElementById("delta").value);
    var t0 = parseFloat(document.getElementById("t0").value);
    var tf = parseFloat(document.getElementById("tf").value);
    var S0 = parseFloat(document.getElementById("S0").value);
    var I0 = parseFloat(document.getElementById("I0").value);
    var R0 = parseFloat(document.getElementById("R0").value);
    var epsilon = parseFloat(document.getElementById("epsilon").value);
    var dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Enter into arrayOfInputs
    return [dtInitial, epsilon, beta, gamma, delta, t0, tf, S0, I0, R0];
}

/** 
 * Solve the problem using RK45.
 *
 * @params           None. Uses parameter values in the forum.
 * @return           Nothing. But it enters the solution values into the solution
 * object.
 */
function solveProblem(arrayOfInputs) {
    // Obtain the parameters of the problem
    var dtInitial = arrayOfInputs[0];
    var epsilon = arrayOfInputs[1];
    var beta = arrayOfInputs[2];
    var gamma = arrayOfInputs[3];
    var delta = arrayOfInputs[4];
    var t0 = arrayOfInputs[5];
    var tf = arrayOfInputs[6];
    var S0 = arrayOfInputs[7];
    var I0 = arrayOfInputs[8];
    var R0 = arrayOfInputs[9];

    // Solve problem using RKF45 and return result
    var solution = RKF45(dtInitial, epsilon, beta, gamma, delta, t0, tf, S0, I0, R0);
    return solution;
}

/**
 * Tabulates solution data.
 *
 * @params           None. Uses the entries of the solution object, however. 
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable(arrayOfInputs) {
    var solution = solveProblem(arrayOfInputs);
    var epsilon = arrayOfInputs[1];

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
 * Set the dimensions of specified plot element
 * 
 * @param name          Name of the plot HTML element.
 * @return              None.
 */
function setPlotElementDims(name) {
    // Initialize variables
    var windowInnerHeight = window.innerHeight;

    // Set dimensions
    document.getElementById(name).style = "height: " + windowInnerHeight + "px;";
}

/**
 * Generates a 3D phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generate3DPhasePlot(arrayOfInputs) {
    // Run solveProblem if unrun
    var solution = solveProblem(arrayOfInputs);

    // Extra solution data from solution object
    S = solution.S;
    I = solution.I;
    R = solution.R;

    // Height and width of plot
    setPlotElementDims("phasePlotXYZ");

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
        title: 'Phase plot of the solution to the SIR equations. x = S, y = I and z = R'
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
function generateXYPhasePlot(arrayOfInputs) {
    var solution = solveProblem(arrayOfInputs);

    // Extra solution data from solution object
    S = solution.S;
    I = solution.I;

    // Height and width of plot
    setPlotElementDims("phasePlotXY");

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
function generateXZPhasePlot(arrayOfInputs) {
    var solution = solveProblem(arrayOfInputs);
    
    // Extra solution data from solution object
    S = solution.S;
    R = solution.R;
    
    // Height and width of plot
    setPlotElementDims("phasePlotXZ");
    
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
function generateYZPhasePlot(arrayOfInputs) {
    var solution = solveProblem(arrayOfInputs);

    // Extra solution data from solution object
    I = solution.I;
    R = solution.R;

    // Height and width of plot
    setPlotElementDims("phasePlotYZ");

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
function generateTimePlot(arrayOfInputs) {
    var solution = solveProblem(arrayOfInputs);

    // Extra solution data from solution object
    t = solution.t;
    S = solution.S;
    I = solution.I;
    R = solution.R;

    // Height and width of plot
    setPlotElementDims("timePlot");

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
function generatePlots(arrayOfInputs) {
    generate3DPhasePlot(arrayOfInputs);
    generateXYPhasePlot(arrayOfInputs);
    generateXZPhasePlot(arrayOfInputs);
    generateYZPhasePlot(arrayOfInputs);
    generateTimePlot(arrayOfInputs);
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