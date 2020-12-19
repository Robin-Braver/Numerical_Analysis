"""
以下包含多种非线性方程解法
"""

""""
以下为二分法
param: a,b代表起始区间,func：求解函数，error：误差
return:计算的值
"""
def dichotomy(a,b,func,error):
    if(a>=b):
        print("a,b选择有误")
        return None
    if(func(a)*func(b)>0):
        print("选择区间无解")
        return None
    while(b-a>=error):
        mid = (a + b) / 2
        if func(mid)*func(a)<0:
            b = mid
        else:
            a = mid
    return mid

""""
以下为试值法
param: a,b代表起始区间,func：求解函数，error：误差
return:计算的值
"""
def shizhi(a,b,func,error):
    max = 10 #最大迭代次数
    if (a >= b):
        print("a,b选择有误")
        return None
    if (func(a) * func(b) > 0):
        print("选择区间无解")
        return None
    time = 0
    while True:
        time += 1
        c = b-(func(b)*(b-a)/(func(b)-func(a)))
        if func(c) == 0:#c就是解
            return c
        elif func(b)*func(c)>0:#解在ac间
            b = c
        else:#解在bc间
            a = c
        if(b-a)<error or time==10:
            break
    return c

""""
不动点迭代法
param: func：迭代函数,p0：迭代开始点,tolerance：容忍误差,max_iteration:迭代次数
return:各个迭代点及误差
"""
def fixed_Iteration(func,p0,tolerance,max_iteration):
    count = 1
    for i in range(max_iteration):
        p_tmp = func(p0)
        print("第{0}次计算的结果:{1}".format(count, round(p_tmp,9)))
        error = abs(p_tmp-p0)
        reerror = error/p_tmp
        if error<tolerance or reerror<tolerance:
            break
        p0 = p_tmp
        count+=1

""""
牛顿迭代法
param: func：迭代函数,dfun:迭代函数倒数,p0：迭代开始点,tolerance1：p0容忍误差,tolerance2：y容忍误差,max_iteration:迭代次数
return:最终迭代点，误差，迭代次数，迭代点函数值
"""
def newton_Iteration(func,dfunc,p0,tolerance1,tolerance2,max_iteration):
    p = [p0]
    for i in range(max_iteration):
        p_tmp = p[i]-func(p[i])/dfunc(p[i])
        p.append(p_tmp)
        # error1 = abs(p[i+1]-p[i])
        # error2 = 2*error1/(abs(p[i+1])+tolerance1)
        # if(error1<tolerance1 or error2<tolerance1 or abs(func(p[i+1])<tolerance2)):
        #     break
    return p

""""
割线法
param: func：迭代函数,p0：迭代开始点,p1:跌点开始点,tolerance1：p0容忍误差,tolerance2：y容忍误差,max_iteration:迭代次数
return:最终迭代点，误差，迭代次数，迭代点函数值
"""
def gexian_Iteration(func,p0,p1,tolerance1,tolerance2,max_iteration):
    p = [p0,p1]
    for i in range(max_iteration):
        p_tmp = p[i+1]-func(p[i+1])/(func(p[i+1]-p[i]))*(p[i+1]-p[i])
        p.append(p_tmp)
        error1 = abs(p[i+2]-p[i+1])
        error2 = 2*error1/(abs(p[i+2])+tolerance1)
        if(error1<tolerance1 or error2<tolerance1 or abs(func(p[i+2])<tolerance2)):
            break
    return p[-1],error1,i,func(p[-1])



if __name__ == '__main__':
    def func(x):
        return x**3-3*x-2
    def dfunc(x):
        return 3*(x**2)-3
    print(gexian_Iteration(func,2.1,0.001,0.001,10))
