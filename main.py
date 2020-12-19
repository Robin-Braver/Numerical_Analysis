import math
import matplotlib.pyplot as plt
import  Polynomial_Fitting
import Integral
from Linear_equations import *
from Interpolation_method import *
from Nonlinear_equation import *


'''Q1:运用不动点迭代法求解 根号5的值，精确到小数点后 9 位   2.236067978'''
def Q1():
    def func1(x):
        return x-(x**2-5)/4
    result1 = fixed_Iteration(func1,2.5,0.000000001,10000)
"""Q2 """
def Q2():
    def func(x):
        return x**3-3*x
    def dfunc(x):
        return 3*x**2-3
    for i in range(700000,800000):
        i = i/1000000
        result = newton_Iteration(func,dfunc,i,0.001,0.001,100)
        if result[-1] != 0:
            index1 = i
            break
    for i in range(-800000,-700000):
        i = i/1000000
        result = newton_Iteration(func,dfunc,i,0.001,0.001,100)
        if result[-1] == 0:
            index2 = i
            break
    x1 = -2
    x2 = -0.8
    x3 = 0.5
    x4 = 0.8
    x5 = 2
    r1 = newton_Iteration(func,dfunc,x1,0.001,0.001,100)[-1]
    r2 = newton_Iteration(func,dfunc,x2,0.001,0.001,100)[-2]
    r3 = newton_Iteration(func,dfunc,x3,0.001,0.001,100)[-3]
    r4 = newton_Iteration(func,dfunc,x4,0.001,0.001,100)[-4]
    r5 = newton_Iteration(func,dfunc,x5,0.001,0.001,100)[-5]
    print("临界值为：",index1)
    print("取-2、-0.8、-0.5、0.8、2五个值，其分别属于:",r1,r2,r3,r4,r5)


"""Q3"""
def Q3():
    A = [[4,-1,0,0,0],[-1,4,-1,0,0],[0,-1,4,-1,0],[0,0,-1,4,-1],[0,0,0,-1,4]]
    B = [100,0,0,0,200]
    A = np.array(A)
    B = np.array(B)
    x = np.array([1,1,1,1,1],dtype='float32')
    Lu = LU(A,B)
    Ja = Jacobi(A,B,x,100,0.000001)
    Ga = Gauss_Seidel(A,B,x,100,0.000001)
    print("追赶法:",Lu)
    print("Jacobi:",Ja)
    print("Guass-Seidel:",Ga)

'''Q4'''
#(1)
def Q4_1():
    x = [1.00,1.02,1.04,1.06]
    y = []
    def func2(x):
        return math.exp(x)*(3*x-math.exp(x))
    for i in x:
        y.append(func2(i))
    p = Lagrange(x,y)
    print("插值多项式为：")
    print(p)
    print("1.03的函数值为：",round(p(1.03),4))
#(2)
def Q4_2():
    def f(x):
        return math.exp(x)*(3*x-math.exp(x))
    def df(x):
        return math.exp(x)*(3*x+3-2*math.exp(x))
    def cof_matirx(x,y):
        result = []
        result1=[1,x,x**2,x**3]
        result.append(result1)
        result2=[1,y,y**2,y**3]
        result.append(result2)
        result3=[0,1,2*x,3*x*x]
        result.append(result3)
        result4=[0,1,2*y,3*y*y]
        result.append(result4)
        result=np.array(result)
        return result
    y=[f(1),f(1.05),df(1),df(1.05)]
    y=np.array(y)
    z=GaussianElimination(cof_matirx(1,1.05),y)
    z = z.tolist()
    z.reverse()
    p = np.poly1d(z)
    def Her(z,x):
        return z[0]+z[1]*x+z[2]*x*x+z[3]*x*x*x
    print("插值多项式：",p)
    print("1.03近似值为：",round(p(1.03),4))

"""Q5"""
def Q5():
    def func(x):
        return 1/(1+x**2)
    def generateDate(n,a,b):
        result = []
        h = (b-a)/n
        for i in range(n):
            result.append(a+i*h)
        result.append(b)
        return result
    x1 = generateDate(2,-5,5)
    x2 = generateDate(10,-5,5)
    print(x1,x2)
    y1 = [func(i) for i in x1]
    y2 = [func(i) for i in x2]
    p1 = Lagrange(x1,y1)
    p2 = Lagrange(x2,y2)
    print("真值，5个，10个对应结果：",func(4.98),p1(4.98),p2(4.98))
    print("绝对误差分别为：",abs(p1(4.98)-func(4.98)),abs(p2(4.98)-func(4.98)))

    x = generateDate(100,-5,5)
    yt = [func(i) for i in x]
    py1 = [p1(i) for i in x]
    py2 = [p2(i) for i in x]
    plt.plot(x, yt, color='r', label='truth')
    plt.plot(x, py1, color='g', label='less')
    plt.plot(x, py2, '-', color='k', label='more')  # label每个plot指定一个字符串标签
    plt.legend(loc='best')
    plt.show()


'''Q6'''
def Q6():
    x = [i for i in range(25)]
    y = [14,13,13,13,13,14,15,17,19,21,22,24,27,30,31,30,28,26,24,23,21,19,17,16,14]
    y_ = [math.log(i) for i in y]
    c2 = Polynomial_Fitting.PolyFitting(x,y,2)
    c3 = Polynomial_Fitting.PolyFitting(x,y,3)
    c4 = Polynomial_Fitting.PolyFitting(x,y,4)
    c_ = Polynomial_Fitting.PolyFitting(x,y_,2)
    p2 = np.poly1d(c2)
    p3 = np.poly1d(c3)
    p4 = np.poly1d(c4)
    p_ = np.poly1d(c_)
    print("多项式分别为：")
    print(p2)
    print(p3)
    print(p4)
    print(p_)
    x = np.linspace(0,25,25)
    y2 = [p2(i) for i in x]
    y3 = [p3(i) for i in x]
    y4 = [p4(i) for i in x]
    y_ = [math.exp(p_(i)) for i in x]
    plt.plot(x,y,color='r', label='truth')
    plt.plot(x, y2, '-',color='y', label='2')  # label每个plot指定一个字符串标签
    plt.plot(x, y3, '-.', color='b', label='3')
    plt.plot(x, y4, color='g', label='4')
    plt.plot(x,y_, color='k', label='special')
    plt.legend(loc='best')
    plt.show()
    err2 = 0
    err4 = 0
    err3 = 0
    err_ = 0
    for i in range(len(x)):
        err2 += abs(p2(x[i])-y[i])
        err3 += abs(p3(x[i])-y[i])
        err4 += abs(p4(x[i]) - y[i])
        err_ += abs(math.exp(p_(x[i])) - y[i])
    print("2次、3次、4次指数绝对误差之和分别为：",err2,err3,err4,err_)

'''Q7'''
def Q7():
    def func1(x):
        if x ==0:
            return 1
        else:
            return math.sin(x)/x
    def func(x):
        return math.exp(-(x**2))
    result = Integral.trapezoid(0, 1, func1, 10)
    print("复化梯形公式求得的结果为：",result)
    S10 = Integral.Simpson(0, 1, func, 10)
    S50 = Integral.Simpson(0, 1, func, 50)
    S100 = Integral.Simpson(0, 1, func, 100)
    print("n=10,50,100的结果分别为：",S10,S50,S100)

Q6()


