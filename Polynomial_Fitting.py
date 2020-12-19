"""
以下包含多种拟合方式
"""
import numpy as np

"""
多项式拟合
param: x：x轴坐标，y：y轴坐标,m:阶次
return:C:系数矩阵,次数由低到高
"""
def PolyFitting(x,y,m):
    n = len(x)
    A = []#计算各阶x的求和
    B = []#计算xy的各阶求和
    for i in range(1,m*2+1):
        x_tmp = 0
        for j in range(0,n):
            x_tmp += x[j]**i
        A.append(x_tmp)
    for i in range(0,m+1):
        xy_tmp = 0
        for j in range(0,n):
            xy_tmp += y[j]*x[j]**i
        B.append(xy_tmp)
    #接下来构造系数矩阵即可。
    X = []
    A.insert(0,n)#将个数插进去，方便构造矩阵
    for i in range(m+1):
        row = []
        for j in range(i,i+m+1):
            row.append(A[j])
        X.append(row)
    X = np.mat(X)
    B = np.mat(B).T
    C = X.I*B
    C = C.getA().tolist()
    p = []
    for i in range(len(C)):
        p.append(C[i][0])
    p.reverse()
    return p




if __name__ == '__main__':
    x = [0.75, 0.86, 0.96, 1.08, 1.12, 1.26, 1.35, 1.51, 1.55, 1.60, 1.63, 1.67, 1.71, 1.78, 1.85]
    y = np.array([10, 12, 15, 17, 20, 27, 35, 41, 48, 50, 51, 54, 59, 66, 75])
    p = np.polyfit(x,y,3)
    p = np.poly1d(p)
    print(p)
    PolyFitting(x,y,3)