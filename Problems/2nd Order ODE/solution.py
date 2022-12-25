import numpy as np
import matplotlib.pyplot as plt

# f1 = dy/dt
def f1( t, y, z ):
  return z

# f2 = dz/dt
def f2( x, y, z ):
  return (60-y-19.42*(z**2))/77.7


# Function for euler formula
def euler( t0, y0, z0, h, t_final ):
  t = t0
  y = y0
  z = z0
  y_list = [y0]

  n = int((t_final - t)/h) # Find number of iterations using step size

  for i in range(1, n + 1): # Start iterations and find k values
    z = z + h*f2(t,y,z)
    y = y + h*f1(t,y,z)
    t += h
    y_list.append(y)

  # Printing approximation
  print("Approximate euler solution at x = ", t, " is ", "%.6f"% y)
  return y_list


def rungeKutta4th(t0, y0, z0, h, t_final):
  """
  That function applies the Runge-Kutta Order 4 Method to system of two ODE's
  """
  # Create lists to deposit values to tabulate results
  y_list = [y0]
  z_list = [z0]

  k1_list = []
  k2_list = []
  k3_list = []
  k4_list = []

  l1_list = []
  l2_list = []
  l3_list = []
  l4_list = []

  t_list = [t0]

  t = t0
  n = int((t_final - t)/h) # Find number of iterations using step size
  y = y0 # Assign initial value to y values
  z = z0

  for i in range(1, n + 1): # Start iterations and find k values
    k1 = f1(t, y, z)
    l1 = f2(t, y, z)

    k2 = f1(t + 0.5 * h, y + 0.5 *h*k1, z + 0.5 *h*l1)
    l2 = f2(t + 0.5 * h, y + 0.5 *h*k1, z + 0.5 *h*l1)

    k3 = f1(t + 0.5 * h, y + 0.5 *h*k2, z + 0.5 *h*l2)
    l3 = f2(t + 0.5 * h, y + 0.5 *h*k2, z + 0.5 *h*l2)

    k4 = f1(t + 0.5 * h, y + 0.5 *h*k3, z + 0.5 *h*l3)
    l4 = f2(t + 0.5 * h, y + 0.5 *h*k3, z + 0.5 *h*l3)
    
    y = y + (h / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) # Update y_n to y_n+1  
    z = z + (h / 6.0)*(l1 + 2 * l2 + 2 * l3 + l4)

    t += h # Update t_n to t_n+1  

    # Deposit founded values to investigate 
    t_list.append(t)
    y_list.append(y)

  print("Approximate Rk4 solution at x = ", t_list[-1], " is ", "%.6f"% y_list[-1])

  return y_list



# Driver Code
# Initial Values
t0 = 0
y0 = 0
z0 = 0
h = 0.1
 
# Value of x at which we need approximation
t_final = 5
 
euler_list = euler(t0, y0, z0, h, 5)
rk4_list = rungeKutta4th(t0, y0, z0, h, 5)


x_ax = np.linspace(t0,t_final,int((t_final - t0)/h)+1) # Create x axis for plot

fig, ax = plt.subplots()  # Create a figure and an axes.

# ax.plot(x_ax, euler_list,label =" Euler Solution") # Plot euler 
ax.plot(x_ax, rk4_list,label ="y: Current Fluid Level") # Plot rk4 

ax.set_xlabel('Time (s)')  # Add an x-label to the axes.
ax.set_ylabel('Fluid Level')  # Add a y-label to the axes.
ax.set_title("Fluid Level vs Time")  # Add a title to the axes.
ax.legend(loc='center left', bbox_to_anchor=(0.25, -0.25))  # Add a legend.
