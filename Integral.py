"""
以下包含多种积分方法
"""
import math
import numpy as np
"""
复化梯形公式
param: a:积分下限，b：积分上限，func：被积函数，n：分多少份
return:所求积分
"""
def trapezoid(a,b,func,n):
    h = (b-a)/n
    sum = 0
    for i in range(1,n):
        xi = a+i*h
        sum += func(xi)
    sum = sum*2+func(a)+func(b)
    sum *= h/2
    return sum

"""
复化Simpson公式
param: a:积分下限，b：积分上限，func：被积函数，n：分多少份
return:所求积分
"""
def Simpson(a,b,func,n):
    h = (b-a)/n
    sum = 0
    for i in range(1,n):
        xi = a+i*h
        sum += 2*func(xi)
    for i in range(0,n):
        xi = a+(i+0.5)*h
        sum += 4 * func(xi)
    sum = sum+func(a)+func(b)
    sum *= h/6
    return sum

"""
变步长梯形公式
param: a:积分下限，b：积分上限，func：被积函数，error：误差
return:所求积分
"""
def ChangeSteps_trapezoid(a,b,func,error):
    step = 1#初始步长
    T1 = (b-a)/2*(func(b)+func(a))
    while True:
        step *= 2
        h = (b-a)/step
        T2 = 0
        for i in range(1,step):
            xi = a + i * h
            T2 += func(xi)
        T2 = T2 * 2 + func(a) + func(b)
        T2 *= h / 2
        if abs(T2-T1)<=error:
            break
        else:
            T1 = T2
            continue
    return T2

"""
Romberg公式
param: a:积分下限，b：积分上限，func：被积函数，error：误差
return:Romberg公式表格
"""
def Romberg(a,b,func,error):
    #以下四个列表分别记录复化梯形序列、Simpson序列，Cotes序列，Romberg序列
    T = []
    S = []
    C = []
    R = []
    T.append(trapezoid(a,b,func,1))#梯形的第一项
    m = 1
    while True:
        sum_t = trapezoid(a,b,func,2*m)
        T.append(sum_t)
        if m>=1:
            S.append((4*T[-1]-T[-2])/3)
        if m>=2:
            C.append((16*S[-1]-S[-2])/15)
        if m>=3:
            R.append((64*C[-1]-C[-2])/63)
        if m>4 and abs(R[-1]-R[-2])<error:
            break
        m = m *2
    table = []
    for i in [T,S,C,R]:
        i = np.array(i)
        table.append(i)
    table = np.array(table)
    for i in range(m//2+1):
        if i >= 3:
            print(i, i+1 * 2, T[i], S[i-1], C[i-2],R[i-3])
        elif i >= 2:
            print(i, i * 2, T[i],S[i-1], C[i-2])
        elif i >= 1:
            print(i, i * 2, T[i], S[i-1])
        else:
            print(i, i * 2, T[i])
    return table

if __name__ == '__main__':
    def func(x):
        return 4/(1+x**2)
    x = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]
    y = [1, 0.997, 0.989, 0.977, 0.959, 0.946, 0.909, 0.877, 0.841]
    for i in range(1,5):
        print(Simpson(0,1,func,2**i))
    Romberg(0,1,func,0.0001)
