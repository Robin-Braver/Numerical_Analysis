
import math

def f1(x):
    return (2*x**3+2)/(3*x**2-1)

def f2(x1,x0):
    up = (x1**3)*x0-x1*(x0**3)+2*x1-2*x0
    down = x1**3-x0**3-x1+x0
    return up/down


x0 = 1.5
x1 = 1.51
for i in range(4):
    x = f2(x1,x0)
    print(x)
    x0 = x1
    x1 = x