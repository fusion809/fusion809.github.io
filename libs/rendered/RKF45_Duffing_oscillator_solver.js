/** 
* Solve the problem using RKF45
*
* @param objectOfInputs An object containing all the problem parameters.
* @return               [t, vars]
*/
RKF45 = function(objectOfInputs) {
    // Extract initial conditions from object and enter it into RKF45Body
    var {"x"0, "dx"0} = objectOfInputs;
    var vars0 = [["x"0, "dx"0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}