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
    const { tolType, epsilon } = objectOfInputs;
    const nVars = vars[i].length;

    // Allocate arrays only once
    const K = Array.from({ length: 6 }, () => new Array(nVars));
    const X1 = new Array(nVars);
    const X2 = new Array(nVars);
    const Rarr = new Array(nVars);
    const sarr = new Array(nVars);

    // Precompute the six K vectors without using arrAdd/arrMult
    K[0] = f(objectOfInputs, t[i], vars[i], dt);

    {
        const temp = new Array(nVars);
        for (let j = 0; j < nVars; j++) temp[j] = vars[i][j] + K[0][j] / 4;
        K[1] = f(objectOfInputs, t[i] + dt / 4, temp, dt);
    }

    {
        const temp = new Array(nVars);
        for (let j = 0; j < nVars; j++) temp[j] = vars[i][j] + (3 / 32) * K[0][j] + (9 / 32) * K[1][j];
        K[2] = f(objectOfInputs, t[i] + (3 * dt) / 8, temp, dt);
    }

    {
        const temp = new Array(nVars);
        for (let j = 0; j < nVars; j++) {
            temp[j] = vars[i][j] +
                (1932 / 2197) * K[0][j] -
                (7200 / 2197) * K[1][j] +
                (7296 / 2197) * K[2][j];
        }
        K[3] = f(objectOfInputs, t[i] + (12 * dt) / 13, temp, dt);
    }

    {
        const temp = new Array(nVars);
        for (let j = 0; j < nVars; j++) {
            temp[j] = vars[i][j] +
                (439 / 216) * K[0][j] -
                8 * K[1][j] +
                (3680 / 513) * K[2][j] -
                (845 / 4104) * K[3][j];
        }
        K[4] = f(objectOfInputs, t[i] + dt, temp, dt);
    }

    {
        const temp = new Array(nVars);
        for (let j = 0; j < nVars; j++) {
            temp[j] = vars[i][j] -
                (8 / 27) * K[0][j] +
                2 * K[1][j] -
                (3544 / 2565) * K[2][j] +
                (1859 / 4104) * K[3][j] -
                (11 / 40) * K[4][j];
        }
        K[5] = f(objectOfInputs, t[i] + dt / 2, temp, dt);
    }

    // Compute X1, X2, Rarr, sarr in one loop
    for (let j = 0; j < nVars; j++) {
        const k0 = K[0][j], k2 = K[2][j], k3 = K[3][j], k4 = K[4][j], k5 = K[5][j];

        const x1 = vars[i][j] + (25 / 216) * k0 + (1408 / 2565) * k2 + (2197 / 4104) * k3 - (1 / 5) * k4;
        const x2 = vars[i][j] + (16 / 135) * k0 + (6656 / 12825) * k2 + (28561 / 56430) * k3 - (9 / 50) * k4 + (2 / 55) * k5;

        X1[j] = x1;
        X2[j] = x2;

        const diff = Math.abs(x1 - x2);
        Rarr[j] = tolType && x1 !== 0 ? diff / (dt * Math.abs(x1)) : diff / dt;

        sarr[j] = Math.pow(epsilon / (2 * Rarr[j]), 0.25);
    }

    // Max error estimate
    const R = Math.max(...Rarr);
    const s = Math.min(...sarr);

    // Accept step if within tolerance
    if (R <= epsilon) {
        t.push(t[i] + dt);
        vars.push(X1);
        i++;
    }

    dt *= s;
    return [dt, t, vars, i];
}

/**
 * Test function to determine whether RKF45 should exit
 * @param i              Dummy variable indicating where in loop we are.
 * @param tf             Maximum time value.
 * @param dt             Step size for time.
 * @param t              Time vector for simulation.
 * @param objectOfInputs Object of page inputs.
 * @returns              Boolean.
 */
function testRKF45(i, tf, vars, dt, t, objectOfInputs) {
    test1 = t[i] < tf;
    if (Object.hasOwn(objectOfInputs, 'hMin')) {
        hMin = objectOfInputs.hMin;
    } else {
        hMin = (tf-t[0])/1e8;
    }
    test2 = dt >= hMin;
    if (!test2) {
        var msg = "Exiting RKF45 at t = " + t[i] + " because dt = " + dt + " is <" + hMin + ".";
        msg += "vars = "
        for (let j = 0; j < vars.length; j++) {
            msg += vars[j]
            if (j != vars.length - 1) {
                msg += ", "
            }
        }
        alert(msg);
        console.log(msg);
    }
    test = test1 && test2; 

    
    return test;
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
    var {t0, tf, hInitial} = objectOfInputs;
    var vars = vars0;
    var t = [t0];
    var dt = hInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( testRKF45(i, tf, vars[i], dt, t, objectOfInputs) ) {
        dt = Math.min(dt, tf-t[i]);
        [dt, t, vars, i] = approxRKF45(f, dt, objectOfInputs, t, vars, i);
    }

    // Transpose vars
    vars = vars[0].map((_, colIndex) => vars.map(row => row[colIndex]));

    // Return t and vars
    return [t, vars];
}

/**
 * Solution object constructor
 * 
 * @param t        An array of t values.
 * @param vars     An array of dependent variable values.
 * @return         Nothing. Creates an object containing t and vars, with methods to extract data. 
 */
function SolClass(t, vars) {
    this.t = t;
    this.vars = vars;

    this.extract = function(i) {
        return this.vars[i];
    }

    this.varsLen = function() {
        return (this.vars).length;
    }

    this.tLen = function() {
        return (this.t).length;
    }
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
    var solution = new SolClass(t, vars);
    return solution;
}