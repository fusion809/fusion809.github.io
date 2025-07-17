/** 
* Solve the problem using RKF45
*
* @param objectOfInputs An object containing all the problem parameters.
* @return               [t, vars]
*/
RKF45 = function(objectOfInputs) {
    // Extract initial conditions from object and enter it into RKF45Body
    var {theta10, dtheta10, theta20, dtheta20} = objectOfInputs;
    var vars0 = [[theta10, dtheta10, theta20, dtheta20]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}