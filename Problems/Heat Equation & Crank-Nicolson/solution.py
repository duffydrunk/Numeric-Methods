import numpy as np
from tabulate import tabulate

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

def inputs(x_top_limit,x_bot_limit,t_top_limit,t_bot_limit,k,h):
  r = k/(h**2)

  i = int((x_top_limit-x_bot_limit)/h)
  j = int((t_top_limit-t_bot_limit)/k)

  x_values = np.arange(x_bot_limit,x_top_limit+h,h)
  t_values = np.arange(t_bot_limit,t_top_limit+k,k)

  u = np.zeros((t_values.shape[0],x_values.shape[0]))

  return u, t_values, x_values, r

def apply_bc_and_ic(u,x_values,t_values):

  v_ic_func = np.vectorize(ic_func)
  v_bc_func = np.vectorize(bc_func)
  u[0,1:-1] = v_ic_func(x_values[1:-1]) # Apply IC
  u[:,0] = v_bc_func(t_values) # Apply BC
  u[:,-1] = v_bc_func(t_values) # Apply BC

def crank_nicolson(u,r,x_values,t_values):
  N = x_values.shape[0]

  #initialize matrices A, B and b array
  A = np.zeros((N-2,N-2))
  B = np.zeros((N-2,N-2))
  b = np.zeros((N-2))

  for i in range(N-2):
    if i==0:
        A[i,:] = [-2*(1+(1/r)) if j==0 else 1 if j==1 else 0 for j in range(N-2)]
        B[i,:] = [ 2*(1-(1/r)) if j==0 else -1 if j==1 else 0 for j in range(N-2)]
        b[i] = 0 #boundary condition at i=0
    elif i==N-3:
        A[i,:] = [1 if j==N-4 else -2*(1+(1/r)) if j==N-3 else 0 for j in range(N-2)]
        B[i,:] = [-1 if j==N-4 else 2*(1-(1/r)) if j==N-3 else 0 for j in range(N-2)]
        b[i] = 0 #boundary condition at i=N
    else:
        A[i,:] = [1 if j==i-1 or j==i+1 else -2*(1+(1/r)) if j==i else 0 for j in range(N-2)]
        B[i,:] = [-1 if j==i-1 or j==i+1 else 2*(1-(1/r)) if j==i else 0 for j in range(N-2)]

  rhs = B.dot(u[0,1:-1]) + b # Initial RHS 
  counter = 1
  t_values_edited = np.delete(t_values, 0) 
  for t in t_values_edited:
     u[counter,1:-1] = np.linalg.solve(A,rhs) # Solve the system
     rhs = B.dot(u[counter,1:-1]) + b # Update RHS
     counter += 1

# Define parameters
x_top_limit = 1
x_bot_limit = 0

t_top_limit = 0.5
t_bot_limit = 0

k = 0.05
h = 0.2

u, t_values, x_values, r = inputs(x_top_limit,x_bot_limit,t_top_limit,t_bot_limit,k,h) # Transform inputs

apply_bc_and_ic(u,x_values,t_values) # Apply BC and IC to U matrix

crank_nicolson(u,r,x_values,t_values) # Apply Crank-Nicolson and fill the U matrix

n = np.arange(t_values.shape[0]) # n values for table
n = np.reshape(n,(t_values.shape[0],1))

table = np.concatenate((n,u),axis=1) # Form table

header = ["j","u\u2080,\u2c7c","u\u2081,\u2c7c","u\u2082,\u2c7c","u\u2083,\u2c7c","u\u2084,\u2c7c","u\u2085,\u2c7c"] # Define table header
print(tabulate(table,header,tablefmt="grid" ,floatfmt=(".0f", ".7f",".7f",".7f",".7f",".7f",".7f"))) # Print table

print(r)
