"""
以下包含多种线性方程组解法
"""

import numpy as np
""""
Guass选主元解方程组
param: A:系数矩阵，B：常数向量
return:求得的方程组解x向量
"""
def Guass(A,B):
    augA = np.concatenate((A,B),axis=1)#拼接系数矩阵和向量
    rows = A.shape[0]
    for k in range(rows):#消元
        prow = np.argmax(np.abs(augA[k:,k]))+k#找到最大行
        pivot = augA[prow,k]#最大值
        if prow != k:
            augA[[k,prow],:] = augA[[prow,k],:]
        for i in range(k+1,rows):
            mi = augA[i,k]/pivot
            augA[i,:] = augA[i,:] - augA[k,:]*mi
    x = np.zeros(rows)
    k = rows-1
    x[k] = augA[k,-1]/augA[k,k]
    for k in range(rows-2,-1,-1):
        tx = x[k+1:]
        ta = augA[k,k+1:-1].flatten()
        x[k] = (augA[k,-1]-np.sum(tx*ta))/augA[k,k]
    return x

def GaussianElimination(A,B):
    N = len(A)
    for i in range(1,N):
        for j in range(i,N):
            # 计算消元因子delta
            delta = A[j][i-1]/A[i-1][i-1]
            # 从第i-1行开始消元
            for k in range(i-1,N):
                # 对A进行消元
                A[j][k] = A[j][k] - A[i-1][k]*delta
            # 对B进行消元
            B[j] = B[j]-B[i-1]*delta
    # 进行回代，直接将方程的解保留在B中
    B[N-1] = B[N-1]/A[N-1][N-1]
    for i in range(N-2,-1,-1):
        for j in range(N-1,i,-1):
            B[i] = B[i]- A[i][j]*B[j]
        B[i] = B[i]/A[i][i]
    # 返回所有解的列表
    return B


""""
LU分解解方程组
param: A:系数矩阵，B：常数向量
return:L,U矩阵求得的方程组解x向量
"""
def LU(A,B):
    L = np.eye(len(A)) #生成对角矩阵
    U = np.zeros(np.shape(A))
    for r in range(1,len(A)):#计算U第一列和L第一行
       U[0,r-1]  = A[0,r-1]
       L[r,0] = A[r,0]/A[0,0]
    U[0,-1] = A[0,-1]
    for r in range(1,len(A)):
        for i in range(r,len(A)):#求U
            sum = 0
            for k in range(0,r):
                sum += L[r,k]*U[k,i]
            U[r,i] = A[r,i]-sum
        for i in range(r+1,len(A)):#求L
            sum = 0
            for k in range(0, r):
                sum += L[i, k] * U[k, r]
            L[i, r] = (A[i, r] - sum)/U[r,r]

    '''以上过程求解L,U矩阵'''
    L = np.mat(L)
    U = np.mat(U)
    print(L)
    print(U)
    y = np.matmul(L.I,B)
    print(y)
    x = np.matmul(U.I,y.T)
    x = x.tolist()
    x = [i[0] for i in x]
    result = []
    for i in x:
        if i>10:
            result.append(round(i,4))
        else:
            result.append(round(i,5))
    return result

""""
Jacobi迭代解方程组
param: A:系数矩阵，B：常数向量,X:初始迭代点,max_iteration:最大迭代次数,err:误差
return:求得的方程组解x向量,迭代次数，误差
"""
def Jacobi(A,B,X,max_iteration,err):
    n = len(A)
    x0 = X
    x1 = np.zeros(n)
    for time in range(max_iteration):
        for i in range(n):
            sum = 0
            for j in range(n):
                if i != j:
                    sum += x0[j]*A[i][j]
            x1[i] = (B[i]-sum)/A[i][i]#更新每个x
        error = max(abs(x1-x0))
        if error<err:#小于指定误差
            break
        else:
            x0 = x1.copy()#更新继续迭代,不可直接复制，否则为一个变量
        #print(x1)
    result = []
    for i in x1:
        if i > 10:
            result.append(round(i, 4))
        else:
            result.append(round(i, 5))
    return result

""""
Gauss-Seidel迭代解方程组
param: A:系数矩阵，B：常数向量,X:初始迭代点,max_iteration:最大迭代次数,err:误差
return:求得的方程组解x向量,迭代次数，误差
"""
def Gauss_Seidel(A,B,X,max_iteration,err):
    n = len(A)
    x0 = X
    x1 = np.zeros(n)
    for time in range(max_iteration):
        temp_x = x0.copy()
        for i in range(n):
            sum = 0
            for j in range(n):
                if i != j:
                    sum += x0[j]*A[i][j]
            x1[i] = (B[i]-sum)/A[i][i]#更新每个x
            x0[i] = x1[i].copy()#立即更新x0
        error = max(abs(x1-temp_x))
        if error<err:#小于指定误差
            break
        else:
            x0 = x1.copy()#更新继续迭代,不可直接复制，否则为一个变量
        #print(x1)
    result = []
    for i in x1:
        if i > 10:
            result.append(round(i, 4))
        else:
            result.append(round(i, 5))
    return result


if __name__ == '__main__':
    # A = np.array([[-5,2,-1],[1,0,3],[3,1,6]], dtype='float')
    # A1 = np.array([[-5,2,-1],[3,1,6],[1,0,3]], dtype='float')
    # B = np.array([-1,7,18], dtype='float')
    # B1 = np.array([-1, 18, 7], dtype='float')
    #x0 = np.array([0.9,2.9,1.9],dtype='float')
    A = np.array([[10, 3, 1], [2, -10, 3], [1,3, 10]], dtype='float')
    B = np.array([14,-5,14], dtype='float')
    A1 = np.array([[10, 3, 1], [1, 3, 10], [2, -10, 3]], dtype='float')
    B1 = np.array([14,14,-5], dtype='float')
    x0 = np.array([0,0, 0], dtype='float')
    # print(LU(A,B))
    print(Guass(A1,B1.reshape(3,1)))
    print(Jacobi(A1,B1,x0,100,0.0001))
    print(Gauss_Seidel(A1,B1,x0,100,0.0001))