/**
 * Add 1d arrays contained with a 2d array
 * 
 * @param arr      A 2d array with each column being a different array to be added.
 * @return         1d array whose elements are the sum of corresponding elements of the 1d subarrays of arr.
 */
function arrAdd(arr) {
    // Initialize relevant variables
    var finalArr = new Array(n1);
    var n1 = arr.length;
    var n2 = arr[0].length;

    // Add elements across each column of arr
    for (let i = 0; i < n2 ; i++) {
        finalArr[i] = arr[0][i];
        for (let j = 1; j < n1; j++) {
            finalArr[i] += arr[j][i];
        }
    }

    // Return the column vector of values
    return finalArr;
}

/**
 * Divide arr by scalar
 * 
 * @param arr      The array whose elements are to be divided.
 * @param scalar   The scalar by which arr's elements are to be divided.
 * @return         Divided array.
 */
function arrDiv(arr, scalar) {
    return arrMult(arr, 1/scalar);
}

/**
 * Array to be multiplied by a scalar
 * 
 * @param arr      Array to be multiplied
 * @param scalar   Scalar by which array is to be multiplied.
 */
function arrMult(arr, scalar) {
    // Initialize relevant variables
    var n = arr.length;
    var finalArr = new Array(n);

    // Multiply each element by scalar
    for (let j = 0; j < n; j++) {
        finalArr[j] = arr[j]*scalar;
    }

    return finalArr;
}