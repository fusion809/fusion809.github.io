@def hassim=false;
@def title="Derivation of the equations of motion of the triple pendulum"
@def mintoclevel=1
@def maxtoclevel=2


The following Python code was used to derive the equations of motion

```python
from sympy import symbols, Function, diff, cos, sin, simplify, sqrt, Abs, Eq, solve, latex, Array
from sympy.vector import CoordSys3D
from multiprocessing import Pool, cpu_count
N = CoordSys3D('N');
# Set up symbols
t = symbols('t')
m1r = symbols('m1r'); 
m2r = symbols('m2r'); 
m3r = symbols('m3r'); 
m1b = symbols('m1b'); 
m2b = symbols('m2b'); 
m3b = symbols('m3b'); 
l1 = symbols('11');
l2 = symbols('l2');
l3 = symbols('l3');
g = symbols('g'); 
b1r = symbols('b1r'); 
b2r = symbols('b2r'); 
b3r = symbols('b3r'); 
b1b = symbols('b1b'); 
b2b = symbols('b2b'); 
b3b = symbols('b3b'); 
c1r = symbols("c1r"); 
c2r = symbols("c2r"); 
c3r = symbols("c3r");
c1b = symbols("c1b"); 
c2b = symbols("c2b"); 
c3b = symbols("c3b");

theta1=Function('theta1')(t); 
theta2=Function('theta2')(t); 
theta3=Function('theta3')(t); 
x1b = l1*cos(theta1);
x1r = x1b/2;
y1b = l1*sin(theta1);
y1r = y1b/2;
r1b = x1b*N.i + y1b*N.j;
r1r = x1r*N.i + y1r*N.j;
x2b = x1b + l2*cos(theta2);
y2b = y1b + l2*sin(theta2);
x2r = x1b + l2*cos(theta2)/2;
y2r = y1b + l2*sin(theta2)/2;
r2b = x2b*N.i + y2b*N.j;
r2r = x2r*N.i + y2r*N.j;
x3b = x2b + l3*cos(theta3);
y3b = y2b + l3*sin(theta3);
x3r = x2b + l3*cos(theta3)/2;
y3r = y2b + l3*sin(theta3)/2;
r3b = x3b*N.i + y3b*N.j;
r3r = x3r*N.i + y3r*N.j;

def compute_kinetic(m, l, theta, r):
	T = 0;
	j = 0;
	v = [];
	for i in m:
		v.append(diff(r[j], t))
		T += i*v[j].dot(v[j])/2;
		j += 1
		
	for i in range(len(l)):
		T += m[2*i]*l[i]**2*diff(theta[i], t)**2/24;
	
	return T
	
def compute_potential(g, m, r):
	V = 0;
	for j in range(len(r)):
		y = r[j].dot(N.j)
		V += m[j] * g * y
	
	return V
	
def compute_lagrangian(g, l, m, r, theta):
	L = compute_kinetic(m, l, theta, r) - compute_potential(g, m, r)
	return L

def diss_force(b, c, r):
	v = diff(r, t);
	return -(b+c*sqrt(v.dot(v)))*v
	
def gen_diss_force(b, c, r, coord):
	Q = [0 for i in range(len(coord))]
	for j in range(len(coord)):
		for i in range(len(r)):
			v = diff(r[i], t)
			F = diss_force(b[i], c[i], r[i])
			ehat = diff(r[i], coord[j])
			Q[j] += F.dot(ehat)

	return Q
	
def compute_eq_of_motion(g, m, b, c, r, l, theta):
	t = symbols('t')
	L = compute_lagrangian(g, l, m, r, theta)
	coord = [theta1, theta2, theta3]
	Q = [Function('Qtheta1')(t), Function('Qtheta2')(t), Function('Qtheta3')(t)];
	lhs = [0 for i in range(len(coord))]
	return [Eq(diff(diff(L, diff(coord[i], t)), t) - diff(L, coord[i]), Q[i]) for i in range(len(coord))]
	
m = [m1r, m1b, m2r, m2b, m3r, m3b]
b = [b1r, b1b, b2r, b2b, b3r, b3b]
c = [c1r, c1b, c2r, c2b, c3r, c3b]
r = [r1r, r1b, r2r, r2b, r3r, r3b]
l = [l1, l2, l3];
coord = [theta1, theta2, theta3]
Q = gen_diss_force(b, c, r, coord)
eqns = compute_eq_of_motion(g, m, b, c, r, l, coord)
d2 = [diff(i, (t, 2)) for i in coord];
theta_ddot_syms = symbols('theta1_dd theta2_dd theta3_dd')
#sol = solve(eqn, d2)

# Turn Eq(lhs, rhs) into lhs - rhs == 0
residuals = [eq.lhs - eq.rhs for eq in eqns]

# Substitute second derivatives with dummy symbols (e.g. theta1_dd)
residuals_subs = [res.subs(dict(zip(d2, theta_ddot_syms))) for res in residuals]
from sympy.solvers.solveset import linear_eq_to_matrix

A, b = linear_eq_to_matrix(residuals_subs, theta_ddot_syms)
print(simplify(A))
```