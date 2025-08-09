/**
 * Add arrays
 * 
 * @param arrs   A list of arrays to be added
 * @return       Array containing elements summed from the input arrays.
 */
function arrAdd(...arrs) {
    // Previously used an iterative algorithm (see 
    // https://github.com/fusion809/fusion809.github.io/blob/a8845b185a0a68c283ed456b7a53e93a2ad1aec3/_libs/common/arrays.js).
    // The code below is from StackOverflow:
    // https://stackoverflow.com/a/55533058/1876983.
    const n = arrs.reduce((max, xs) => Math.max(max, xs.length), 0);
    const result = Array.from({ length: n });
    return result.map((_, i) => arrs.map(xs => xs[i] || 0).reduce((sum, x) => sum + x, 0));
}

function arrScAdd(arr, scalar) {
    return arr.map(x => x + scalar);
}
/**
 * Divide arr by scalar
 * 
 * @param arr    The array whose elements are to be divided.
 * @param scalar The scalar by which arr's elements are to be divided.
 * @return       Divided array.
 */
function arrDiv(arr, scalar) {
    return arr.map(x => x/scalar);
}

/**
 * Array to be multiplied by a scalar
 * 
 * @param arr    Array to be multiplied
 * @param scalar Scalar by which array is to be multiplied.
 */
function arrMult(arr, scalar) {
    return arr.map(x => x*scalar);
}

function arrsMult(arr1, arr2) {
    return arr1.map((val, i) => val * arr2[i]);
}

function arrSq(arr) {
    return arrsMult(arr, arr);
}

/**
 * Run a function on two vectors
 * @param f    Function to be run.
 * @param vec1 First vector argument.
 * @param vec2 Second vector argument.
 * @returns    Result from running function on both arguments. 
 */
function func2Vecs(f, vec1, vec2) {
    var val;
    try {
        val = f(...vec1.map((v,i) => f(v, vec2[i])));
    } catch (e) {
        if (e instanceof RangeError) {
            val = vec1.concat(vec2).reduce((a, b) => f(a, b), 0);
        } else {
            throw e; // rethrow if it's not the expected error
        }
    }
    return val;
}

function vecMag(vec) {
    var sum = 0;
    for (let i = 0; i < vec.length; i++) {
        sum += vec[i]^2;
    }
    return Math.sqrt(sum);
}