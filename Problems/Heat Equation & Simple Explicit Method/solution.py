import numpy as np
from tabulate import tabulate

"""
Formulation of the problem is done with simple explicit method.
"""


def inputs(x_top_limit,x_bot_limit,t_top_limit,t_bot_limit,k,h):
  r = k/(h**2)

  i = int((x_top_limit-x_bot_limit)/h)
  j = int((t_top_limit-t_bot_limit)/k)

  x_values = np.arange(x_bot_limit,x_top_limit+h,h)
  t_values = np.arange(t_bot_limit,t_top_limit+k,k)

  u = np.zeros((t_values.shape[0],x_values.shape[0]))

  return u, t_values, x_values, r

def simple_explicit(u,r,x_values):

  u[0,1:-1] = v_ic_func(x_values[1:-1]) # Apply IC
  u[:,0] = v_bc_func(t_values) # Apply BC
  u[:,10] = v_bc_func(t_values) # Apply BC

  for j in range(u.shape[0]-1):
    
    for i in range(1,u.shape[1]-1):    
      u[j+1,i] = r*u[j,i-1] + (1-2*r)*u[j,i] + r*u[j,i+1]

x_top_limit = 1
x_bot_limit = 0

t_top_limit = 0.5
t_bot_limit = 0

k = 0.0025
h = 0.1

u, t_values, x_values, r = inputs(x_top_limit,x_bot_limit,t_top_limit,t_bot_limit,k,h)

def ic_func(x):
  """
  Inıtial condition funciton
  """
  return np.sin(np.pi*x)

def bc_func(x):
  """
  Boundary condition function
  """
  return 0

v_ic_func = np.vectorize(ic_func)
v_bc_func = np.vectorize(bc_func)

simple_explicit(u,r,x_values)

print("For case 1:")

n = np.arange(t_values.shape[0]) # n values for table
n = np.reshape(n,(t_values.shape[0],1))

table = np.concatenate((n,u),axis=1)

header = ["j","u\u2080,\u2c7c","u\u2081,\u2c7c","u\u2082,\u2c7c","u\u2083,\u2c7c","u\u2084,\u2c7c","u\u2085,\u2c7c", \
          "u\u2086,\u2c7c","u\u2087,\u2c7c","u\u2088,\u2c7c","u\u2089,\u2c7c","u\u2081\u2080,\u2c7c"] # Define table header
print(tabulate(table,header,tablefmt="grid" ,floatfmt=(".0f", ".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f"))) # Print table


print("\nFor case 2:")

k = 0.01
h = 0.1

x_top_limit = 1
x_bot_limit = 0

t_top_limit = 0.5
t_bot_limit = 0


u, t_values, x_values, r = inputs(x_top_limit,x_bot_limit,t_top_limit,t_bot_limit,k,h)

simple_explicit(u,r,x_values)

n = np.arange(t_values.shape[0]) # n values for table
n = np.reshape(n,(t_values.shape[0],1))

table = np.concatenate((n,u),axis=1)

header = ["j","u\u2080,\u2c7c","u\u2081,\u2c7c","u\u2082,\u2c7c","u\u2083,\u2c7c","u\u2084,\u2c7c","u\u2085,\u2c7c", \
          "u\u2086,\u2c7c","u\u2087,\u2c7c","u\u2088,\u2c7c","u\u2089,\u2c7c","u\u2081\u2080,\u2c7c"] # Define table header
print(tabulate(table,header,tablefmt="grid" ,floatfmt=(".0f", ".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f",".7f"))) # Print table
