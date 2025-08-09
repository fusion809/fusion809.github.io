/**
 * Build spatial grid x_j and spectral wave numbers k for 1D periodic domain.
 * @param {number} N - number of grid points (integer)
 * @param {number} L - domain length
 * @returns {{ xj: Float64Array, k: Float64Array }}
 */
function makeGrid(x0, x1, N) {
  const xj = new Array(N);
  const k = new Array(N);
  const L = x1-x0;
  for (let j = 0; j < N; j++) {
    xj[j] = x0 + j * L / N; // x_j = j * Δx, Δx = L/N

    // frequency in cycles per unit length: fftfreq analogue
    // for j < N/2: positive, else negative wrapped
    const freq = (j < N / 2) ? (j / L) : ((j - N) / L);
    k[j] = 2 * Math.PI * freq; // wave number
  }
  return [xj, k];
}

function toFlatArray(u) {
  if (u == null) throw new Error("input is null/undefined");
  // TypedArray (Float64Array, etc.)
  if (ArrayBuffer.isView(u)) return Array.from(u);
  if (!Array.isArray(u)) throw new Error("expected array or array of [x]");

  return u.map(el => (Array.isArray(el) && el.length === 1 ? el[0] : el));
}
function f(objectOfInputs, t, vars, dt) {
    var {alpha, x0, x1, N} = objectOfInputs;
    var [x, k] = makeGrid(x0, x1, N);
    var u = vars;
    u = toFlatArray(u);
    var dudt = math.dotMultiply(math.ifft(math.dotMultiply(arrSq(k).map(x => -x), math.fft(u))), alpha);
    var du = math.dotMultiply(dudt, dt);
    var du = math.re(du);
    return du;
}

RKF45 = function(objectOfInputs) {
    // Extract initial conditions from object and enter it into RKF45Body
    var {x0, x1, N} = objectOfInputs;
    var [x, k] = makeGrid(x0, x1, N);
    var L = x1-x0;
    var vars0 = [arrScAdd(arrMult(x.map(x => Math.cos(2*Math.PI/L * x)), 10), 10)];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}