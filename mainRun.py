import scipy.special as spc
import numpy as np
import scipy as sp
import plotly.offline as py
import numbers
import matplotlib.pyplot as plt


parameter_array = {
    'R': 4,
    'l': 0.5,
    'k': 0.011,
    'alf': 0.005,
    'c': 1.6,
    'betta': 0.008,
    'P': 40,
    'a': 0.8,
    'T': 180}

array_norm_N = {'0.1':7,
                '0.01':18,
                '0.001':52,
                '0.0001':137,
                '0.00001':347,
                '0.000001':877,
                '0.0000001':2202}

dataChanged=[0]
def updateParameret_array(dan):
    parameter_array['R'] = dan['R']
    parameter_array['l'] = dan['l']
    parameter_array['k'] = dan['k']
    parameter_array['alf'] = dan['alf']
    parameter_array['c'] = dan['c']
    parameter_array['betta'] = dan['betta']
    parameter_array['P'] = dan['P']
    parameter_array['a'] = dan['a']
    parameter_array['T'] = dan['T']

# надо как нибудь оптимизировать
bessel0 = [0, ]
temp = spc.jn_zeros(1, 100000)

for t in temp:
    bessel0.append(t)
##
def Bn(n):
    j0 = spc.j0(bessel0[n])
    j1 = spc.j1(bessel0[n] / 5)
    numerator = 2 * parameter_array['betta'] * parameter_array['P']*j1
    denominator = 5*parameter_array['c']*np.pi * (parameter_array['a'])**2 *bessel0[n]*(j0)**2
    return numerator/denominator

def An(n,t):
    numerator = Bn(n)*parameter_array['c']*parameter_array['l']*(parameter_array['R'])**2
    denominator = (2*parameter_array['alf'] * (parameter_array['R'])**2 )+ (parameter_array['k']* \
                  bessel0[n]**2 * parameter_array['l'])

    x =(((-2*parameter_array['alf'] * parameter_array['R']**2) - (parameter_array['k']* \
                  bessel0[n]**2 * parameter_array['l']))*t)/( parameter_array['c']*parameter_array['l']*parameter_array['R']**2)

    all = numerator/denominator
    return all * (1- np.exp(x))

def analitic_eps(rank):
    #rank = 10**(-rank)
    N = ((2*parameter_array['betta']*parameter_array['R']*parameter_array['P']*10**(1/2))/ (15*(parameter_array['a'])**2 *parameter_array['k']*
                                                                      (np.pi)**3 * rank))**(2/3)
    return round(N)

flag=[]
NMAX=[0]


def experimental_eps(r, t, eps, date=None):
    if date != None:
        updateParameret_array(date)

    numerator = (parameter_array['l'] * parameter_array['betta'] * parameter_array['P'])
    denominator = (50 * sp.pi * parameter_array['alf'] * parameter_array['a'] ** 2)
    all = numerator / denominator
    x = ((-2 * parameter_array['alf'] * t) / (parameter_array['c'] * parameter_array['l']))
    u = all * (1 - np.exp(x))

    N = analitic_eps(eps)

    for i in np.arange(1, N+1, 1):
        u += An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])
        flag.append(u)

    i=N-2
    analit_eps = flag[N-1]
    while (abs ( analit_eps - flag[i])<=eps)and(i>0):
        i=i-1
    if ((i+2)>NMAX[0]):
        NMAX[0] = i+2
    flag.clear()
    return NMAX[0]

def u(r,t, eps,colslag,date=None):
    if date != None:
        updateParameret_array(date)

    numerator = (parameter_array['l']*parameter_array['betta']*parameter_array['P'] )
    denominator = (50*sp.pi*parameter_array['alf']*parameter_array['a']**2)
    all = numerator/denominator
    x = ((-2*parameter_array['alf']*t)/(parameter_array['c']*parameter_array['l']))
    u = all * (1- np.exp(x))
    if (colslag==0):
        if dataChanged[0]==1:

            N = analitic_eps(eps)
            for i in np.arange(1, N, 1):
                u += An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])

        else:
            N = array_norm_N[str(eps)]

            for i in np.arange(1,N,1):
                u +=An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])

    else:
        for i in np.arange(1,colslag,1):
            u +=An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])

    return u

def getparams():
    return parameter_array


def drawUIT(da, eps,colSlag):
    data = [dict(
    visible=False,
    line=dict(color='#ff0000', width=2),
    name='t = ' + str(step),
    x=np.arange(0, da['R']+0.1, 0.1),
    y=[u(t, step, eps,colSlag, da) for t in np.arange(0, da['R']+0.1, 0.1)]) for step in
    np.arange(da['T']+1)]
    steps = []
    for i in range(len(data)):
        step = dict(
            method='restyle',
            args=['visible', [False] * len(data)])
        step['args'][1][i] = True # Toggle i'th trace to "visible"
        steps.append(step)
    sliders = [dict(
    active=10,
    currentvalue={"prefix": "Time: "},
    pad={"t": 50},
    steps=steps
    )]
    layout = dict(sliders=sliders)
    fig = dict(data=data, layout=layout)
    py.plot(fig, filename='Temperature in time.html')

def drawUIR(da, eps,colSlag):
    data = [dict(
    visible=False,
    line=dict(color='#ff0000', width=2),
    name='r = ' + str(step),
    x=np.arange(da['T']),
    y=[u(step, t, eps,colSlag,da) for t in np.arange(da['T']+1)]) for step in np.arange(0, da['R']+0.1, 0.1)]

    steps = []

    for i in range(len(data)):
        step = dict(
            method='restyle',
            args=['visible', [False] * len(data)])
        step['args'][1][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)
    sliders = [dict(
    active=10,
    currentvalue={"prefix": "Radius: "},
    pad={"t": 50},
    steps=steps
    )]
    layout = dict(sliders=sliders)
    fig = dict(data=data, layout=layout)
    py.plot(fig, filename='Temperature in radius.html')

# Для построения графика погрешности
Xnew = []
Ynew=[]

def plot_eps(x,y):

    fig=  plt.figure("График погрешности ")
    ax = fig.add_subplot(111)
    plt.grid(True)
    ax.set_yticks([2202,1900,1700,1500,1250,1000,877,650,500,347,137,7])
    plt.plot(x, y, marker='o',markersize = 2.5,markeredgecolor = '#0000FF')
    plt.plot(x, y,linewidth = 1,color = '#00BFFF')
    plt.plot(Xnew,Ynew,color='#FF0000')
    plt.xlabel("accuracy")
    plt.ylabel("number of members")

    ax.set_xscale('log')



    plt.show()

# X = [0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1]
# Y=[2202,877,347,137,52,18,7]
#
# yt = 0.0000001
# st = -8
# i =0
# while (yt<=10**-1):
#     Ynew.append(experimental_eps(0,180,yt,parameter_array))
#     NMAX[0]=0
#     Xnew.append(yt)
#     yt= yt + 9*10**st
#     i=i+1
#     if(i==9):
#         st=st+1
#         i=0
#         yt=10**(st+1)
#
#
# plot_eps(X,Y)

exit = 0
while(exit==0):
    da = getparams()
    print("Текущие параметры")
    print("R = {0} l = {1} k = {2}".format(da['R'], da['l'], da['k']))
    print("alfa = {0} c = {1} betta = {2}".format(da['alf'], da['c'], da['betta']))
    print("P = {0}, a = {1}, T={2}".format(da['P'], da['a'], da['T']))
    yes = input("Введите +, если хотите изменить параметры  иначе - ")

    if yes.upper() == "+":
        dataChanged[0]=1
        while True:
            try:
                R = float(input("Введите R "))
                if((R<=0) or R>=10):
                    print("Вы должны ввести положительное число больше, 0 и меньше 10, попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести положительное число больше, 0 и меньше 10, попробуйте снова.")

        while True:
            try:
                l = float(input("Введите l "))
                if  (l<=0 or l>=5):
                    print("Вы должны ввести положительное число, больше 0 и меньше 5 , попробуйте снова.")
                else:
                     break
            except ValueError:
                print("Вы должны ввести положительное число, больше 0 и меньше 5 , попробуйте снова.")

        while True:
            try:
                k = float(input("Введите k "))
                if  ((k<=0) or (k>1)):
                    print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")

        while True:
            try:
                alfa = float(input("Введите alfa "))
                if((alfa<=0) or (alfa>1)):
                    print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")

        while True:
            try:
                c = float(input("Введите c "))
                if ((c<=0) or (c>10)):
                    print("Вы должны ввести положительное число, больше 0 и меньше 10 , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести положительное число, больше 0 и меньше 10 , попробуйте снова.")

        while True:
            try:
                betta = float(input("Введите betta "))
                if ((betta<=0) or (betta>1)):
                    print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")

        while True:
            try:
                P = float(input("Введите P "))
                if(P<=0):
                    print("Вы должны ввести положительное число, попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести положительное число, попробуйте снова.")

        while True:
            try:
                T = int(input("Введите T "))
                if(T <= 0):
                     print("Вы должны ввести целое положительное число , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести целое положительное число , попробуйте снова.")

        while True:
            try:
                a = float(input("Введите a "))
                if(a <= 0):
                    print("Вы должны ввести целое положительное число , попробуйте снова.")
                else:
                     break
            except ValueError:
                print("Вы должны ввести целое положительное число , попробуйте снова.")


        da['R'] = float(R)
        da['l'] = float(l)
        da['k'] = float(k)
        da['alf'] = float(alfa)
        da['c'] = float(c)
        da['betta'] = float(betta)
        da['P'] = float(P)
        da['a'] = float(a)
        da['T'] = float(T)
        updateParameret_array(da)
        da = getparams()
        print("Текущие параметры")
        print("R = {0} l = {1} k = {2}".format(da['R'], da['l'], da['k']))
        print("alfa = {0} c = {1} betta = {2}".format(da['alf'], da['c'], da['betta']))
        print("P = {0}, a = {1}, T={2}".format(da['P'], da['a'], da['T']))

    scoreOrSeries = (input("Хотите вычислить ряд ( введите - 1) или погрешность (введите 2) ? "))

    if scoreOrSeries == "1":
        yes_no = (input(
            "Хотите посчитать ряд используя точность( введите - 1) или количество слагаемых для посчета(введите 2) ? "))

        if yes_no == "1":

            while True:
                try:
                    eps = float(input("Введите точность (пример 0.001 или 10e-4)"))
                    if ((eps<0) or (eps>=0.9)):
                        print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")

            print("Подождите, пожалуйста, выполняются вычисления...")
            drawUIR(da, eps, 0)
            drawUIT(da, eps, 0)

        if yes_no == "2":

            while True:
                try:
                    n = int(input("Введите количество слагаемых N = "))
                    if ((n<=3) or (n>10000)):
                        print("Вы должны ввести целое число > 3 и < 10 000, попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести целое число > 3 и < 10 000, попробуйте снова.")

            print("Подождите, пожалуйста, выполняются вычисления...")
            drawUIR(da, 0, n)
            drawUIT(da, 0, n)

    if scoreOrSeries == "2":
        error = (input(
            "Рассчитать кол -во слагаемых для достижени конкретной точности (введите 1 ) или для всего ряда (введите 2 (это займет много времени ))"))
        if error == "1":
            while True:
                try:
                    ep = float(input("Введите точность (пример 0.001 или 10e-4)"))
                    if ((ep<0) or (ep>=0.9)):
                        print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")

            print("Подождите, пожалуйста, выполняются вычисления...")
            for t in np.arange(1, da['T'] + 1):
                for r in np.arange(0, da['R'] + 0.1, 0.1):
                    experimental_eps(r, t, ep, da)
            print("NMAX избыточое для eps = 10^-", ep, "равно ", analitic_eps(ep))
            print("NMAX достаточное для eps = 10^-", ep, "равно ", NMAX)
            NMAX[0] = 0
            flag.clear()

        if error == "2":
            print("Подождите, пожалуйста, выполняются вычисления...")
            Y = []
            EPS = [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1] # степени
            for ep in EPS:
                for t in np.arange(1, da['T'] + 1):
                    for r in np.arange(0, da['R'] + 0.1, 0.1):
                        experimental_eps(r, t, ep, da)
                print("N избыточое для eps = ", ep, "равно ", analitic_eps(ep))
                print("N достаточное для eps = ", ep, "равно ", NMAX)
                Y.append(NMAX[0])
                NMAX[0] = 0
                flag.clear()

            plot_eps(EPS, Y)
    exit = int(input("Для выхода нажмите - 1 , чтобы продолжить - 0"))
    dataChanged[0]=0
    print("-------------------------------------------------------------------------------------")
