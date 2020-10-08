/**
 * Approximate solution at next step
 * 
 * @param f              A function representing dt times the RHS of the 1st order ODE system we're solving.
 * @param dt             Step size.
 * @param objectOfInputs An object containing all the inputs from the forms.
 * @param t              An array of t values.
 * @param vars           An array of solution values.
 * @param i              Counter variable.
 * @return               Updated [dt, t, vars, i]
 */
function approxRKF45(f, dt, objectOfInputs, t, vars, i) {
    // Initialize variables
    var K = [[]];
    var X1 = [];
    var X2 = [];
    var Rarr = [];
    var sarr = [];
    var {epsilon} = objectOfInputs;

    // K[i][j], i goes from 0 to 5 and represents k1, k2, k3, k4, k5, k6
    // j goes from 0 to the dimensionality of the problem - 1
    K[0] = f(objectOfInputs, t[i], vars[i], dt);
    K[1] = f(objectOfInputs, t[i] + dt/4, arrAdd(vars[i], arrDiv(K[0], 4)), dt);
    K[2] = f(objectOfInputs, t[i] + 3*dt/8, arrAdd(vars[i], arrMult(K[0], 3/32), arrMult(K[1], 9/32)), dt);
    K[3] = f(objectOfInputs, t[i] + 12*dt/13, arrAdd(vars[i], arrMult(K[0], 1932/2197), arrMult(K[1], -7200/2197), arrMult(K[2], 7296/2197)), dt);
    K[4] = f(objectOfInputs, t[i]+dt, arrAdd(vars[i], arrMult(K[0], 439/216), arrMult(K[1], -8), arrMult(K[2], 3680/513), arrMult(K[3], -845/4104)), dt);
    K[5] = f(objectOfInputs, t[i] + dt/2, arrAdd(vars[i], arrMult(K[0], -8/27), arrMult(K[1], 2), arrMult(K[2], -3544/2565), arrMult(K[3], 1859/4104), arrMult(K[4], -11/40)), dt);

    // X1[j] are our fourth order approxs & X2[j] are our 5th order approxs
    for (let j = 0; j < K[0].length; j++) {
        X1[j] = vars[i][j] + 25*K[0][j]/216+1408*K[2][j]/2565+2197*K[3][j]/4104-K[4][j]/5;
        X2[j] = vars[i][j] + 16*K[0][j]/135+6656*K[2][j]/12825+28561*K[3][j]/56430-9*K[4][j]/50+2*K[5][j]/55;
        Rarr[j] = Math.abs(X1[j]-X2[j])/dt;
        sarr[j] = Math.pow(epsilon/(2*Rarr[j]), 0.25);  
    }

    // R is our error estimate
    // s is the factor by which our step size is to be adjusted
    var R = Math.max(...Rarr);
    var s = Math.min(...sarr);
    // If R is less than or equal to epsilon move onto the next step
    if ( R <= epsilon ) {
        t.push(t[i]+dt);
        vars.push(X1);
        i++;
    }
    dt *= s;

    return [dt, t, vars, i];
}

/**
 * Solves the problem using RKF45.
 * 
 * @param f              A function representing the RHS of the 1st order ODE system we're solving.
 * @param objectOfInputs An object containing all the problem parameters.
 * @param vars0          Initial conditions of the problem in a 2D array (e.g. [[x0, y0]]).
 * @return               [t, vars] where vars is an array of solution values.
 */
function RKF45Body(f, objectOfInputs, vars0) {
    // Initialize vars
    var {t0, tf, dtInitial} = objectOfInputs;
    var vars = vars0;
    var t = [t0];
    var dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        dt = Math.min(dt, tf-t[i]);
        [dt, t, vars, i] = approxRKF45(f, dt, objectOfInputs, t, vars, i);
    }

    // Transpose vars
    vars = vars[0].map((_, colIndex) => vars.map(row => row[colIndex]));

    // Return t and vars
    return [t, vars];
}

/** 
 * Solve the problem using RKF45.
 *
 * @param func           RKF45 function to be used.
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               A solution object.
 */
function solveProblem(func, objectOfInputs) {
    // Solve the problem, and time it
    let t0 = performance.now();
    var [t, vars] = func(objectOfInputs);
    let t1 = performance.now();

    // Log execution time to console
    let diff = t1 - t0;
    console.log("Solving the problem took " + diff + " milliseconds.");
    
    // Write t and vars to our solution object
    var solution = {
        t: t,
        vars: vars
    };
    return solution;
}