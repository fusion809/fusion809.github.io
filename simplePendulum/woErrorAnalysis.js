/**
 * Calculates the 4th and 5th order approximations of theta for the current 
 * step.
 * 
 * @param g             Acceleration due to gravity parameter.
 * @param l             Length of pendulum bob parameter.
 * @param dt            Change in t for current step.
 * @param t             t array.
 * @param theta         theta array.
 * @param thetaDot      thetaDot array.
 * @param i             Iteration counter.
 * @return              Array of 4th and 5th order corrections to theta and 
 * theta dot
 */
function approximatorRKF45(g, l, dt, t, theta, thetaDot, i) {
    // Runge-Kutta-Fehlberg approximations of the change in theta and 
    // thetaDot over the step
    var k1 = dt*f(g, l, t[i], theta[i], thetaDot[i])[0];
    var l1 = dt*f(g, l, t[i], theta[i], thetaDot[i])[1];
    var k2 = dt*f(g, l, t[i]+dt/4, theta[i]+k1/4, thetaDot[i]+l1/4)[0];
    var l2 = dt*f(g, l, t[i]+dt/4, theta[i]+k1/4, thetaDot[i]+l1/4)[1];
    var k3 = dt*f(g, l, t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, thetaDot[i]+3*l1/32+9*l2/32)[0];
    var l3 = dt*f(g, l, t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, thetaDot[i]+3*l1/32+9*l2/32)[1];
    var k4 = dt*f(g, l, t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, thetaDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
    var l4 = dt*f(g, l, t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, thetaDot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
    var k5 = dt*f(g, l, t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, thetaDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
    var l5 = dt*f(g, l, t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, thetaDot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
    var k6 = dt*f(g, l, t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, thetaDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
    var l6 = dt*f(g, l, t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, thetaDot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

    // 4th order approx
    var theta1 = theta[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    var thetaDot1 = thetaDot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    // theta2 and thetaDot2 are our fifth order approximations
    var theta2 = theta[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    var thetaDot2 = thetaDot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

    return [theta1, thetaDot1, theta2, thetaDot2];
}

/**
 * A function that checks whether step size has been adequate and adjusts 
 * it accordingly. If the step size is adequate, it update the t, theta and
 * thetaDot arrays with the current step values.
 * 
 * @param theta1        4th order approximation for theta
 * @param theta2        5th order approximation for theta
 * @param thetaDot1     4th order approximation for theta dot
 * @param thetaDot2     5th order approximation for theta dot
 * @param epsilon       Error tolerance.
 * @param i             Counter variable value.
 * @param dt            dt value.
 * @param t             Array of t values.
 * @param theta         Array of theta values.
 * @param thetaDot      Array of theta dot values.
 * @return              i, dt, and updated t, theta and thetaDot arrays.
 */
function stepSizeChecker(theta1, theta2, thetaDot1, thetaDot2, epsilon, i, dt, t, theta, thetaDot) {
    var RTheta = Math.abs(theta1-theta2)/dt;
    var RThetaDot = Math.abs(thetaDot1-thetaDot2)/dt;
    var sTheta = Math.pow(epsilon/(2*RTheta), 1/4);                
    var sThetaDot = Math.pow(epsilon/(2*RThetaDot), 1/4);
    var R = Math.max(RTheta, RThetaDot);
    var s = Math.min(sTheta, sThetaDot);
    if ( R <= epsilon ) {
        t.push(t[i]+dt);
        theta.push(theta1);
        thetaDot.push(thetaDot1);
        i++;
    }
    dt *= s;

    return [i, dt, t, theta, thetaDot];
}

/**
 * Uses Runge-Kutta-Fehlberg 4/5th order method over the domain of integration.
 * 
 * @param dtInitial     Initial dt value.
 * @param epsilon       Maximum error tolerance.
 * @param g             Acceleration due to gravity.
 * @param l             Length of the pendulum bob.
 * @param t0            Beginning time for the simulation.
 * @param tf            End time of simulation.
 * @param theta0        theta(t0) initial condition.
 * @param thetaDot0     thetaDot(t0) initial condition.
 * @return              Solution object containing t, theta and thetaDot 
 * arrays.
 */
function RKF45(dtInitial, epsilon, g, l, t0, tf, theta0, thetaDot0) {
    // Initialize global variables
    var t = [t0];
    var theta = [theta0];
    var thetaDot = [thetaDot0];
    var i = 0;
    var dt = dtInitial;
    var theta1, theta2, thetaDot1, thetaDot2;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);
    
        // Use approximatorRKF45 to make approximations
        [theta1, thetaDot1, theta2, thetaDot2] = approximatorRKF45(g, l, dt, t, theta, thetaDot, i);
    
        // The following are used to correct the step size
        [i, dt, t, theta, thetaDot] = stepSizeChecker(theta1, theta2, thetaDot1, thetaDot2, epsilon, i, dt, t, theta, thetaDot);
    }

    // Write values to solution object and return it, 
    // a replacement for multiple returns that some languages
    // support.
    var solution = {
        t: t,
        theta: theta,
        thetaDot: thetaDot
    };
    return solution;
}

/** 
 * Solve the problem using RKF45.
 *
 * @param objectOfInputs  An array that contains all the problem parameters.
 * @return                A solution object.
 */
function solveProblemSP(objectOfInputs) {
    // Obtain the parameters of the problem
    var {g, l, t0, tf, theta0, thetaDot0, epsilon, dtInitial} = objectOfInputs;

    // Solve the problem
    var solution = RKF45(dtInitial, epsilon, g, l, t0, tf, theta0, thetaDot0);

    // Write number of steps to table field
    document.getElementById("NRKF45").innerHTML = solution.t.length;

    return solution;
};