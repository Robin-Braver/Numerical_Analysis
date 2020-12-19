"""
以下包含多种插值法
"""
import numpy as np
from scipy.interpolate import lagrange
""""
lagrange插值法
param: x：x轴坐标，y：y轴坐标
return:差值多项式
"""
def Lagrange(x,y):
    n = len(x)
    p = 0.0
    for i in range(n):
        pt = y[i]
        for j in range(n):
            if i==j:
                continue
            fac = x[i]-x[j]
            pt *= np.poly1d([1.0,-x[j]])/fac
        p += pt
    return p

"""
Newton插值法
param: x：x轴坐标，y：y轴坐标
return:指数向量
"""
def Newton(x,y,_x):
    '''以下步骤求差商'''
    n = len(x)
    p = np.zeros((n, n + 1))  # 差商表
    p[:, 0] = x  # 初始化x列
    p[:, 1] = y  # 初始化y列
    for j in range(2, n+1):  # 依次求每一列
        p[j-1:n,j] = (p[j-1:n,j-1]-p[j-2:n-1,j-1])/(x[j-1:n]-x[:n+1-j])
    q = np.diag(p, k=1)#差商
    print(q)
    xx = [1]
    for i in range(n-1):#
        tmp = _x-x[i]
        xx.append(xx[i]*tmp)
    result = y[0]
    for i in range(1, n):
        a = q[i]
        result += (a*xx[i])
    return result

if __name__ == '__main__':
    x = [0,1,2,4]
    y = [1,9,23,3]
    p = lagrange(x,y)
    print(p(5))
    print(Newton((np.array(x)).T,(np.array(y)).T,5))