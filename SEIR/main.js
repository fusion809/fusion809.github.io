/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param a          Inverse of the incubation period.
 * @param beta       Infectivity parameter.
 * @param gamma      Parameter that measures the rate of recovery.
 * @param delta      Parameter that measures the effect of quarantine.
 * @param lambda     Birth rate measured in people per day.
 * @param mu         Death rate measured in people per day. 
 * @param t          Time (seconds).
 * @param S          Susceptible person count.
 * @param E          Exposed person count.
 * @param I          Infectious person count.
 * @param R          Recovered person count.
 * @return           [dS/dt, dE/dt, dI/dt, dR/dt]
 */
function f(a, beta, gamma, delta, lambda, mu, t, S, E, I, R) {
    // Determine N
    var N = S + E + I + R;
    // Calculate derivatives
    var exposure = (beta * S * I * (1-delta))/N;
    var dSdt = lambda*N - mu * S - exposure;
    var dEdt = exposure - (mu + a)*E;
    var dIdt = a*E - (mu + gamma) * I;
    var dRdt = gamma*I - mu*R;
    // Put into return value
    return [dSdt, dEdt, dIdt, dRdt];
}

/**
 * Calculates RKF45 approximations for S, I and R values for next step.
 * 
 * @param dt            Step size.
 * @param a             Problem parameter.
 * @param beta          Problem parameter.
 * @param gamma         Problem parameter.
 * @param delta         Problem parameter.
 * @param lambda        Birth rate.
 * @param mu            Death rate.
 * @param t             Array of t values.
 * @param S             Array of S values.
 * @param E             Array of E values.
 * @param I             Array of I values.
 * @param R             Array of R values.
 * @param i             Loop counter index.
 * @return              An array of 4th and 5th order approximations to S, I and R at time next step.
 */
function approximatorRKF45(dt, a, beta, gamma, delta, lambda, mu, t, S, E, I, R, i) {
    // Runge-Kutta-Fehlberg approximations of the change in S, I and R
    // over the step
    // 1st approx
    var k1 = dt*f(a, beta, gamma, delta, lambda, mu, t[i], S[i], E[i], I[i], R[i])[0];
    var l1 = dt*f(a, beta, gamma, delta, lambda, mu, t[i], S[i], E[i], I[i], R[i])[1];
    var m1 = dt*f(a, beta, gamma, delta, lambda, mu, t[i], S[i], E[i], I[i], R[i])[2];
    var n1 = dt*f(a, beta, gamma, delta, lambda, mu, t[i], S[i], E[i], I[i], R[i])[3];
    // 2nd approx
    var k2 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/4, S[i]+k1/4, E[i]+l1/4, I[i]+m1/4, R[i]+n1/4)[0];
    var l2 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/4, S[i]+k1/4, E[i]+l1/4, I[i]+m1/4, R[i]+n1/4)[1];
    var m2 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/4, S[i]+k1/4, E[i]+l1/4, I[i]+m1/4, R[i]+n1/4)[2];
    var n2 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/4, S[i]+k1/4, E[i]+l1/4, I[i]+m1/4, R[i]+n1/4)[3];
    // 3rd approx
    var k3 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, E[i]+3*l1/32+9*l2/32, I[i]+3*m1/32+9*m2/32, R[i]+3*n1/32+9*n2/32)[0];
    var l3 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, E[i]+3*l1/32+9*l2/32, I[i]+3*m1/32+9*m2/32, R[i]+3*n1/32+9*n2/32)[1];
    var m3 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, E[i]+3*l1/32+9*l2/32, I[i]+3*m1/32+9*m2/32, R[i]+3*n1/32+9*n2/32)[2];
    var n3 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, E[i]+3*l1/32+9*l2/32, I[i]+3*m1/32+9*m2/32, R[i]+3*n1/32+9*n2/32)[3];
    // 4th approx
    var k4 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, E[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, I[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, R[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[0];
    var l4 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, E[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, I[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, R[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[1];
    var m4 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, E[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, I[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, R[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[2];
    var n4 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, E[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, I[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, R[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197)[3];
    // 5th approx
    var k5 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, E[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, I[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, R[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[0];
    var l5 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, E[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, I[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, R[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[1];
    var m5 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, E[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, I[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, R[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[2];
    var n5 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, E[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, I[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, R[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104)[3];
    // 6th approx
    var k6 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, E[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, I[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, R[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[0];
    var l6 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, E[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, I[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, R[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[1];
    var m6 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, E[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, I[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, R[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[2];
    var n6 = dt*f(a, beta, gamma, delta, lambda, mu, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, E[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, I[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, R[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40)[3];

    // S1, E1, I1 and R1 are our fourth order approximations
    var S1 = S[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    var E1 = E[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    var I1 = I[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;
    var R1 = R[i] + 25*n1/216+1408*n3/2565+2197*n4/4104-n5/5;

    // S2, E2, I2 and R2 are our fifth order approximations
    var S2 = S[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    var E2 = E[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
    var I2 = I[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;
    var R2 = R[i] + 16*n1/135+6656*n3/12825+28561*n4/56430-9*n5/50+2*n6/55;
    
    return [S1, E1, I1, R1, S2, E2, I2, R2];
}

/**
 * Step size checker and adapter function
 * @param dt            Current step size.
 * @param epsilon       Error tolerance.
 * @param t             An array of t values.
 * @param S             An array of previous S values.
 * @param E             An array of previous E values.
 * @param I             An array of previous I values.
 * @param R             An array of previous R values.
 * @param S1            4th order approx for S at next step 
 * @param E1            4th order approx for E at next step
 * @param I1            4th order approx for I at next step
 * @param R1            4th order approx for R at next step
 * @param S2            5th order approx for S at next step
 * @param E2            5th order approx for E at next step
 * @param I2            5th order approx for I at next step
 * @param R2            5th order approx for R at next step
 * @param i             Counter variable value.
 * @return              Corrected dt and updated t, S, E, I, R and i.
 */
function stepSizeChecker(dt, epsilon, t, S, E, I, R, S1, E1, I1, R1, S2, E2, I2, R2, i) {
    // The following are used to correct the step size
    var RS = Math.abs(S1-S2)/dt;
    var RE = Math.abs(E1-E2)/dt;
    var RI = Math.abs(I1-I2)/dt;
    var RR = Math.abs(R1-R2)/dt;
    var sS = 0.84*Math.pow(epsilon/RS, 1/4);                
    var sE = 0.84*Math.pow(epsilon/RE, 1/4);
    var sI = 0.84*Math.pow(epsilon/RI, 1/4);
    var sR = 0.84*Math.pow(epsilon/RR, 1/4);
    var RRKF45 = Math.max(RS, RE, RI, RR);
    var s = Math.min(sS, sE, sI, sR);

    // If R is less than or equal to epsilon move onto the next step
    if ( RRKF45 <= epsilon ) {
        t.push(t[i]+dt);
        S.push(S1);
        E.push(E1);
        I.push(I1);
        R.push(R1);
        i++;
        dt *= s;
    } else {
        dt *= s;
    }

    return [dt, t, S, E, I, R, i];
}

/**
 * Runge-Kutta-Fehlberg 4/5th order integrator function
 * 
 * @param dtInitial     Initial guess for dt.
 * @param epsilon       Error tolerance.
 * @param a             Inverse of the incubation period in days.
 * @param beta          Infectivity parameter.
 * @param gamma         Recovery rate parameter.
 * @param delta         Quarantine effectiveness parameter.
 * @param lambda        Birth rate.
 * @param mu            Death rate.
 * @param t0            Day the simulation starts
 * @param tf            Day the simulation ends
 * @param S0            Initial number of susceptible persons.
 * @param E0            Initial number of exposed persons.
 * @param I0            Initial number of infectious persons.
 * @param R0            Initial number of recovered persons.
 * @return              Solution object containing arrays of time, S, I and R values.
 */
function RKF45(dtInitial, epsilon, a, beta, gamma, delta, lambda, mu, t0, tf, S0, E0, I0, R0) {
    // Initialize the arrays used and loop variables
    var t = [t0];
    var S = [S0];
    var E = [E0];
    var I = [I0];
    var R = [R0];
    var dt = dtInitial;
    var i = 0;
    
    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Use approximatorRKF45 to make approximations
        var [S1, E1, I1, R1, S2, E2, I2, R2] = approximatorRKF45(dt, a, beta, gamma, delta, lambda, mu, t, S, E, I, R, i);
    
        // Step size correction
        [dt, t, S, E, I, R, i] = stepSizeChecker(dt, epsilon, t, S, E, I, R, S1, E1, I1, R1, S2, E2, I2, R2, i);
    }
    
    // Write t, S, I and R to our solution object
    var solution = {
        t: t,
        S: S,
        E: E, 
        I: I,
        R: R
    };
    return solution;
}

/** 
 * Solve the problem using RKF45.
 *
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               A solution object.
 */
function solveProblem(objectOfInputs) {
    // Obtain the parameters of the problem
    var {a, beta, gamma, delta, lambda, mu, t0, tf, S0, E0, I0, R0, dtInitial, epsilon} = objectOfInputs;

    // Solve problem using RKF45 and return result
    var solution = RKF45(dtInitial, epsilon, a, beta, gamma, delta, lambda, mu, t0, tf, S0, E0, I0, R0);
    return solution;
}

/**
 * Tabulates solution data.
 *
 * @param objectOfInputs An object containing the problem parameters.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs) {
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(objectOfInputs);
    var {t, S, E, I, R} = solution;
    // Extract epsilon from objectOfInputs
    var epsilon = objectOfInputs.epsilon;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>S</th>';
    tableContents += '<th>E</th>';
    tableContents += '<th>I</th>';
    tableContents += '<th>R</th>';
    tableContents += "</tr>";
    for (let j = 0; j < S.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + S[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + E[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
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
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generate3DPhasePlot(objectOfInputs) {
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(objectOfInputs);
    var {S, I, R} = solution;

    // Height and width of plot
    adjustPlotHeight("phasePlotXYZ");

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
    var dataXYZ = [plotXYZ];

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
    rmPlot("phasePlotXYZ");
}

/**
 * Generates a XY phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateXYPhasePlot(objectOfInputs) {
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(objectOfInputs);
    var {S, I} = solution;

    // Height and width of plot
    adjustPlotHeight("phasePlotXY");

    // Plot object and data object array
    var plotXY = {
        x: S,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXY = [plotXY];

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
    rmPlot("phasePlotXY");
}

/**
 * Generates a XZ phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateXZPhasePlot(objectOfInputs) {
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(objectOfInputs);
    var {S, R} = solution;
    
    // Height and width of plot
    adjustPlotHeight("phasePlotXZ");
    
    // Plot object and data object array
    var plotXZ = {
        x: S,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXZ = [plotXZ];
    
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
    rmPlot("phasePlotXZ");
}

/**
 * Generates a YZ phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateYZPhasePlot(objectOfInputs) {
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(objectOfInputs);
    var {I, R} = solution;

    // Height and width of plot
    adjustPlotHeight("phasePlotYZ");

    // Plot object and data object array
    var plotYZ = {
        x: I,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataYZ = [plotYZ];

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
    rmPlot("phasePlotYZ");
}

/**
 * Generates a time plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateTimePlot(objectOfInputs) {
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(objectOfInputs);
    var {t, S, E, I, R} = solution;

    // Height and width of plot
    adjustPlotHeight("timePlot");

    // Plot object and data object array
    var plotTS = {
        x: t,
        y: S,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'S'
    };
    var plotTE = {
        x: t,
        y: E,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'E'
    };
    var plotTI = {
        x: t,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'I'
    };
    var plotTR = {
        x: t,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'R'
    };
    var dataTimePlot = [plotTS, plotTE, plotTI, plotTR];

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
    rmPlot("timePlot");
}

/**
 * Generate five plots:
 * - The first is a 3D phase plot of S, I and R.
 * - The second is a 2D phase plot of I against S.
 * - The third is a 2D phase plot of R against S.
 * - The fourth is a 2D phase plot of R against I.
 * - The fifth is a plot of S, I and R.
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generate3DPhasePlot(objectOfInputs);
    generateXYPhasePlot(objectOfInputs);
    generateXZPhasePlot(objectOfInputs);
    generateYZPhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    // Clear HTML and CSS of the plots
    removeTimePlot();
    remove3DPhasePlot();
    removeXYPhasePlot();
    removeXZPhasePlot();
    removeYZPhasePlot();
};