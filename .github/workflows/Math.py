
PI = 3.14159265358979323846
TOL = 0.000000000000001
BERNOULLI = [1]
import math as math
import timeit
def genBernoulli(n):
    for k in range(0,n):
        if len(BERNOULLI)> k:
            BERNOULLI[k] = bernoulliNum(k)
        else:
            BERNOULLI.append(bernoulliNum(k))

def abs(x):
    if x<0:
        return -x
    else:
        return x
    
def fac(x):
    r = 1
    t = 2
    while t<=x:
        r *= t
        t += 1
    return r

binCoef = lambda n, k :  fac(n)/(fac(n-k)*fac(k))

def exp(x):
    t = 1.0
    r = 1.0
    n = 1
    while abs(t)>TOL:
        t *= x/n
        n += 1
        r += t
    #print("Exp " + str(x) + ": " + str(n))
    return r

def ln(x):
    if x <= 0:
        return
    t = 1.0
    r = 1.0
    n = 1
    #print("---------ln(" + str(x) + ")----------")
    while abs(t)>TOL:
        t *= (2*n-1)*(x-1)*(x-1)/((x+1)*(x+1)*(2*n+1))
        r += t#/(2*n+1)
        n += 1
        #print(str(n) + ": " + str(t))
    r *= 2*(x-1)/(x+1)
    #print("Ln " + str(x) + ": " + str(n))
    return r


E = exp(1)

def pow(x, n):
    if x<0 and n%1 != 0:
        return
    if x == 0:
        if n<=0:
            return
        return 0
    return (-1 if x<0 and n%2 == 1 else 1)*exp(n*ln(abs(x)))

def root(x, n):
    if x<0:
        return
    return exp(ln(x)/n)

def sin(x):
    x = x%(2*PI)
    t = x
    r = x
    n = 3
    while abs(t)>TOL:
        t *= x*x/(n*(n-1))
        if n%4 == 3:
            r -= t
        else:
            r += t
        n += 2
    #print(n)
    return r

def cos(x):
    x = x%(2*PI)
    t = 1.0
    r = 1.0 
    n = 2
    while abs(t)>TOL:
        t *= x*x/(n*(n-1))
        if n%4 == 2:
            r -= t
        else:
            r += t
        n += 2
    #print(n)
    return r

def bernoulliNum(n):
    if len(BERNOULLI)>n:
        return BERNOULLI[n]
    if n%2 == 1 and  n!= 1:
        return 0
    B = 1
    k = 0
    while k < n:
        B -= binCoef(n, k)*bernoulliNum(k)/(n-k+1)
        k += 1
    return B

def eulerNum(n):
    if n == 0:
        return 1
    if n%2 == 1:
        return 0
    Sk = 0
    Tk = 1
    k = 1
    while k<= n:
        Sl = 0
        l = 0
        while l <= 2*k:
            Sl += (-1 if l%2==1 else 1)*binCoef(2*k, l)*pow(k-l, n)
            l += 1
        Tk *= -0.5
        Sk += Tk*Sl
        print(Sk)
        k += 1
    return Sk
        
        

def tan(x):
    neg = -1 if x%PI>PI/2 else 1
    x = x%(2*PI)
    if x > PI:
        x = x - (2*PI)
    x = x%(PI/2)
    #print(x)
    r = x
    t = x/2
    n = 2
    while abs(t)>TOL:
        #print(r)
        #print(bernoulliNum(2*n))
        t *= (-4)*x*x/(2*n*(2*n-1))
        r += bernoulliNum(2*n)*t*(1-pow(4,n))
        n += 1
    return r*neg
    #return sin(x)/cos(x)

def sec(x):
    neg = -1 if PI/2< x%(2*PI) < 1.5*PI else 1
    x = x%(PI/2)
    r = 1+(x*x/2)
    t = -1*x*x/2
    n = 2
    while abs(t)>TOL:
        t *= -1*x*x/(2*n*(2*n-1))
        r += t*eulerNum(2*n)
        n += 1
    return r * neg

def csc(x):
    neg = -1 if x%(2*PI) > PI else 1
    x = x%PI
    r = 1/x
    t = -2/x
    n = 1
    while abs(t)>TOL:
        t *= -1*x*x/(2*n*(2*n-1))
        r += bernoulliNum(2*n)*t*(pow(2, 2*n-1)-1)
        n += 1
    return r

def cot(x):
    neg = -1 if x%(2*PI) > PI else 1
    x = x%PI
    t = 1/x
    r = 1/x
    n = 1
    while abs(t)>TOL:
        t *= -1*x*x*4/(2*n*(2*n-1))
        r += bernoulliNum(2*n)*t
        n += 1
    return r

def arcsin(x):
    if abs(x)>1:
        return
    elif abs(x) == 1:
        return x*PI/2
    r = x
    t = x
    n = 1
    while abs(t) > TOL:
        t *= (2*n-1)*x*x*(2*n*(2*n-1))/((4*n*n)*(2*n+1))
        r+= t
        n+=1
    return r

def arccos(x):
    if abs(x)>1:
        return
    return PI/2-arcsin(x)

def arctan(x):
    if abs(x)>1:
        return
    if abs(x) == 1:
        return x*PI/4
    r = x
    t = x
    n = 1
    while abs(t) > TOL:
        t *= -(2*n-1)*x*x/(2*n+1)
        r += t
        n += 1
    return r

genBernoulli(100)
def testTrig():
    pt = timeit.repeat(lambda: (sin(n/100) for n in range(1,628)),repeat=100,number=10000)
    mt = timeit.repeat(lambda: (math.sin(n/100) for n in range(1,628)),repeat=100,number=10000)
    print("Min. runtime of sin over 100 trials")
    print(min(pt))
    print("Min. runtime of math.sin over 100 trials:")
    print(min(mt))
    print("Avg. runtime of sin over 100 trials")
    print(sum(pt)/100)
    print("Avg. runtime of math.sin over 100 trials:")
    print(sum(mt)/100)

    diff = []
    for n in range(1, 628):
        diff.append(abs(sin(n/100)-math.sin(n/100))/abs(math.sin(n/100)))
    mean = sum(diff)/626
    print("Mean % error of sin: " + str(100*mean))
    for n in range(0, 627):
        diff[n] = math.pow(diff[n]-mean, 2)
    print("St. Dev of % error: " + str(math.sqrt(sum(diff)/6.26)))
    print("------------------------------------------")
    pt=timeit.repeat(lambda: (cos(n/100) for n in range(1,628)),repeat=100,number=10000)
    mt=timeit.repeat(lambda: (math.cos(n/100) for n in range(1,1001)),repeat=100,number=10000)
    print("Min. runtime of cos over 100 trials")
    print(min(pt))
    print("Min. runtime of math.cos over 100 trials:")
    print(min(mt))
    print("Avg. runtime of cos over 100 trials")
    print(sum(pt)/100)
    print("Avg. runtime of math.cos over 100 trials:")
    print(sum(mt)/100)

    diff = []
    for n in range(1, 628):
        diff.append(abs(cos(n/100)-math.cos(n/100))/abs(math.cos(n/100)))
    mean = sum(diff)/626
    print("Mean % error of cos: " + str(100*mean))
    for n in range(0, 627):
        diff[n] = math.pow(diff[n]-mean, 2)
    print("St. Dev of % error: " + str(math.sqrt(sum(diff)/6.26)))
    print("------------------------------------------")
    pt =timeit.repeat(lambda: (tan(n/100) for n in range(1,628)),repeat=100,number=10000)
    mt =timeit.repeat(lambda: (math.tan(n/100) for n in range(1,628)),repeat=100,number=10000)
    print("Min. runtime of tan over 100 trials")
    print(min(pt))
    print("Min. runtime of math.tan over 100 trials:")
    print(min(mt))
    print("Avg. runtime of tan over 100 trials")
    print(sum(pt)/100)
    print("Avg. runtime of math.tan over 100 trials:")
    print(sum(mt)/100)

    diff = []
    for n in range(1, 628):
        diff.append(abs(tan(n/100)-math.tan(n/100))/abs(math.tan(n/100)))
    mean = sum(diff)/626
    print("Mean % error of tan: " + str(100*mean))
    for n in range(0, 627):
        diff[n] = math.pow(diff[n]-mean, 2)
    print("St. Dev of % error: " + str(math.sqrt(sum(diff)/6.26)))
    print("------------------------------------------")
    pt=timeit.repeat(lambda: (pow(2,n/100) for n in range(1,1001)),repeat=100,number=10000)
    mt=timeit.repeat(lambda: (math.pow(2,n/100) for n in range(1,1001)),repeat=100,number=10000)
    print("Min. runtime of pow(2,n) over 100 trials")
    print(min(pt))
    print("Min. runtime of math.pow(2,n) over 100 trials:")
    print(min(mt))
    print("Avg. runtime of pow(2,n) over 100 trials")
    print(sum(pt)/100)
    print("Avg. runtime of math.pow(2,n) over 100 trials:")
    print(sum(mt)/100)
    diff = []
    for n in range(1, 1001):
        diff.append(abs(pow(2,n/100)-math.pow(2,n/100))/abs(math.pow(2,n/100)))
    mean = sum(diff)/1000
    print("Mean % error of pow: " + str(100*mean))
    for n in range(0, 1000):
        diff[n] = math.pow(diff[n]-mean, 2)
    print("St. Dev of % error: " + str(math.sqrt(sum(diff)/10.00)))

#print(tan(0.01))
#print(math.tan(0.01))
#for n in range(0,100):
#    print(str(n) + ": " + str(BERNOULLI[n]))
testTrig()
