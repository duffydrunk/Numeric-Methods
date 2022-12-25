
import math # Import math library for sin and cos functions
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np

def dy1dt(t, y1, y2, y3, c=1, d=5): 
    """
    f(y1,t), y1', function defined in the problem statement.
    """
    return -c * y1 * y2

def dy2dt(t, y1, y2, y3, c=1, d=5): 
    """
    f(y2,t), y2', function defined in the problem statement.
    """
    return (c * y1 * y2) -  (d*y2) 

def dy3dt(t, y1, y2, y3, c=1, d=5): 
    """
    f(y3,t), y3', function defined in the problem statement.
    """
    return d * y2
 
def rungeKutta4th(t_0, y1_0, y2_0, y3_0, t_final, h):
    """
    That function applies the Runge-Kutta Order 4 Method to system of three ODE's
    """
    # Create lists to deposit values to tabulate results
    y1_list = [y1_0]
    y2_list = [y2_0]
    y3_list = [y3_0]

    k1_list = []
    k2_list = []
    k3_list = []
    k4_list = []

    l1_list = []
    l2_list = []
    l3_list = []
    l4_list = []

    m1_list = []
    m2_list = []
    m3_list = []
    m4_list = []

    t_list = [t_0]

    t = t_0
    n = int((t_final - t)/h) # Find number of iterations using step size
    y1 = y1_0 # Assign initial value to y values
    y2 = y2_0
    y3 = y3_0

    for i in range(1, n + 1): # Start iterations and find k values
        k1 = dy1dt(t, y1, y2, y3)
        l1 = dy2dt(t, y1, y2, y3)
        m1 = dy3dt(t, y1, y2, y3)

        k2 = dy1dt(t + 0.5 * h, y1 + 0.5 *h*k1, y2 + 0.5 *h*l1 ,y3 + 0.5 *h*m1)
        l2 = dy2dt(t + 0.5 * h, y1 + 0.5 *h*k1, y2 + 0.5 *h*l1 ,y3 + 0.5 *h*m1)
        m2 = dy3dt(t + 0.5 * h, y1 + 0.5 *h*k1, y2 + 0.5 *h*l1 ,y3 + 0.5 *h*m1)

        k3 = dy1dt(t + 0.5 * h, y1 + 0.5 *h*k2, y2 + 0.5 *h*l2 ,y3 + 0.5 *h*m2)
        l3 = dy2dt(t + 0.5 * h, y1 + 0.5 *h*k2, y2 + 0.5 *h*l2 ,y3 + 0.5 *h*m2)
        m3 = dy3dt(t + 0.5 * h, y1 + 0.5 *h*k2, y2 + 0.5 *h*l2 ,y3 + 0.5 *h*m2)

        k4 = dy1dt(t + 0.5 * h, y1 + 0.5 *h*k3, y2 + 0.5 *h*l3 ,y3 + 0.5 *h*m3)
        l4 = dy2dt(t + 0.5 * h, y1 + 0.5 *h*k3, y2 + 0.5 *h*l3 ,y3 + 0.5 *h*m3)
        m4 = dy3dt(t + 0.5 * h, y1 + 0.5 *h*k3, y2 + 0.5 *h*l3 ,y3 + 0.5 *h*m3)
        
        y1 = y1 + (h / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) # Update y_n to y_n+1  
        y2 = y2 + (h / 6.0)*(l1 + 2 * l2 + 2 * l3 + l4)
        y3 = y3 + (h / 6.0)*(m1 + 2 * m2 + 2 * m3 + m4)

        t += h # Update t_n to t_n+1  

        # Deposit founded values to investigate 
        y1_list.append(y1)
        y2_list.append(y2)
        y3_list.append(y3)

        k1_list.append(k1)
        k2_list.append(k2)
        k3_list.append(k3)
        k4_list.append(k4)

        l1_list.append(l1)
        l2_list.append(l2)
        l3_list.append(l3)
        l4_list.append(l4)

        m1_list.append(m1)
        m2_list.append(m2)
        m3_list.append(m3)
        m4_list.append(m4)

        t_list.append(t)

    return y1_list , y2_list, y3_list, k1_list, k2_list, k3_list, k4_list, \
      l1_list, l2_list, l3_list, l4_list, m1_list, m2_list, m3_list, m4_list, t_list
 
# Define intial conditions, step size and desired value
t_0 = 0
y1_0 = 95
y2_0 = 5
y3_0 = 0
t_final = 1 # 1 months
h = 0.02 # 3 days

y1_list , y2_list, y3_list, k1_list, k2_list, k3_list, k4_list, l1_list, l2_list,\
 l3_list, l4_list, m1_list, m2_list, m3_list, m4_list,\
  t_list = rungeKutta4th(t_0, y1_0, y2_0, y3_0, t_final, h) # Find values by using Range Kutta Order 4

x_ax = np.linspace(t_0,t_final,int((t_final - t_0)/h)+1) # Create x axis for plot

fig, ax = plt.subplots()  # Create a figure and an axes.

ax.plot(x_ax, y1_list,label =" y1: People susceptible to the disease") # Plot y1 
ax.plot(x_ax, y2_list,label =" y2: Infected people that are still in circulation ") # Plot y2 
ax.plot(x_ax, y3_list,label =" y3: Infected people that are removed from the population ") # Plot y3
ax.set_xlabel('Time (month)')  # Add an x-label to the axes.
ax.set_ylabel('Number of People')  # Add a y-label to the axes.
ax.set_title("Number of People vs Time")  # Add a title to the axes.
ax.legend(loc='center left', bbox_to_anchor=(0, -0.35))  # Add a legend.

n = np.arange(x_ax.shape[0]) # n values for table

table = zip(n, t_list, y1_list, y2_list, y3_list) # Create table
header = ["n","t\u2099","y\u2081,\u2099","y\u2082,\u2099","y3\u2083,\u2099"] # Define table header
print(tabulate(table,header,tablefmt="grid" ,floatfmt=(".0f", ".7f",".7f",".7f",".7f"))) # Print table
