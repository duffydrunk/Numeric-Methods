import math # Import math library for sin and cos functions

def dydt(t, y):
	"""
	f(y,t), y', function defined in the problem statement.
	"""
	return 10*(y-math.sin(t))+math.cos(t)

def rungeKutta4th(t0, y0, t, h):
	"""
	That function applies the Runge-Kutta Order 4 Method
	"""
	n = (int)((t - t0)/h) # Find number of iterations using step size
	y = y0 # Assign initial value to y
	for i in range(1, n + 1): # Start iterations and find k values
		k1 = h * dydt(t0, y)
		k2 = h * dydt(t0 + 0.5 * h, y + 0.5 * k1)
		k3 = h * dydt(t0 + 0.5 * h, y + 0.5 * k2)
		k4 = h * dydt(t0 + h, y + k3)
		y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) # Apply Runge-Kutta Order 4 & Update y_n to y_n+1
		t0 = t0 + h # Update t_n to t_n+1
	return y	

# Define intial conditions, step size and desired value
t0 = 0
y0 = 1
t = 0.2
h = 0.05

rk4 = rungeKutta4th(t0, y0, t, h) # Find y(0.2) value by using Range Kutta Order 4
print ('Using Runge-Kutta Order 4 to find y(0.2):', rungeKutta4th(t0, y0, t, h))

def exact_solution(t):
	"""
	Given exact solution of the f(y,t) in the problem statement.
	"""
	return ((math.e)**(10*t))+math.sin(t)

exact = exact_solution(t)# Find exact solution
print ('Exact solution of find y(0.2):', exact_solution(t))

def dydt_5(t, y):
	"""
	Function that returns fifth derivative of y wrt t
	"""
	return 100000*(y-math.sin(t))+math.cos(t)

def estimated_error(h,zeta,y):
	"""
	Local error for Runge-Kutta Order 4 from lecture notes
	"""
	return ((h**5)/(math.factorial(5)))*dydt_5(zeta,y)

 # Zeta is between 0 and 0.2. Thus estimate Zeta as 0.1 & y(zeta) as 4.3
zeta = 0.1
y = 4.3

print ('Estimated Error:', estimated_error(h, zeta,y)) # Find estimated error

print ('Exact Error:', exact - rk4) # Find exact error
 
 
 
 
