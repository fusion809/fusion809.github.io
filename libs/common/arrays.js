/**
 * Add arrays
 * Previously used an iterative algorithm.
 * This code is from StackOverflow https://stackoverflow.com/a/55533058/1876983.
 * 
 * @param arrs     A list of arrays to be added
 * @return         Array containing elements summed from the input arrays.
 */
function arrAdd(...arrs) {
    const n = arrs.reduce((max, xs) => Math.max(max, xs.length), 0);
    const result = Array.from({ length: n });
    return result.map((_, i) => arrs.map(xs => xs[i] || 0).reduce((sum, x) => sum + x, 0));
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
    arr = arr.map(x => x*scalar);
    return arr;
}