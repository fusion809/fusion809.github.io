/**
 * Add arrays
 * 
 * @param arrs   A list of arrays to be added
 * @return       Array containing elements summed from the input arrays.
 */
function arrAdd(...arrs) {
    const n = arrs.reduce((max, xs) => Math.max(max, xs.length), 0);
    const result = new Array(n).fill(0);
    for (let a = 0; a < arrs.length; a++) {
        const arr = arrs[a];
        for (let i = 0; i < arr.length; i++) {
            result[i] += arr[i];
        }
    }
    return result;
}

function arrMult(arr, scalar) {
    const result = new Array(arr.length);
    for (let i = 0; i < arr.length; i++) {
        result[i] = arr[i] * scalar;
    }
    return result;
}

function arrDiv(arr, scalar) {
    const result = new Array(arr.length);
    const inv = 1 / scalar; // Faster than dividing each time
    for (let i = 0; i < arr.length; i++) {
        result[i] = arr[i] * inv;
    }
    return result;
}

function arrScAdd(arr, scalar) {
    const result = new Array(arr.length);
    for (let i = 0; i < arr.length; i++) {
        result[i] = arr[i] + scalar;
    }
    return result;
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

function linspace(start, end, n) {
  if (n < 2) {
    return n === 1 ? [start] : [];
  }
  const step = (end - start) / (n - 1);
  return Array.from({ length: n }, (_, i) => start + i * step);
}