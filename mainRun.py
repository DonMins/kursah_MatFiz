import scipy.special as spc
import numpy as np
import  scipy as sp
import plotly.offline as py



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
temp = spc.jn_zeros(1, 10000)
print(temp)
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
    numerator = Bn(n)*parameter_array['c']*parameter_array['l']*parameter_array['R']**2
    denominator = 2*parameter_array['alf'] * parameter_array['R']**2 + parameter_array['k']* \
                  bessel0[n]**2 * parameter_array['l']

    x =(-2*parameter_array['alf'] * parameter_array['R']**2 - parameter_array['k']* \
                  bessel0[n]**2 * parameter_array['l'] *t) /( parameter_array['c']*parameter_array['l']*parameter_array['R']**2)

    all = numerator/denominator
    return  all * (1- np.exp(x))

def u(r,t, eco, date=None):
    if date != None:
        updateParameret_array(date)

    u = (parameter_array['l']*parameter_array['betta']*parameter_array['P'] *\
         ( 1-sp.exp(-2*parameter_array['alf']*t/(parameter_array['c']*parameter_array['l']))))/\
         (25*sp.pi*2*parameter_array['alf']*parameter_array['a']**2)
    for i in np.arange(1,120,1):
        mn = bessel0[i]
        u +=An(i, t) * spc.j0(mn * r / parameter_array['R'])
    return u

def getparams():
    return parameter_array


def drawUIT(da, eps):
    data = [dict(
    visible=False,
    line=dict(color='#ff0000', width=3),
    name='t = ' + str(step),
    x=np.arange(0, 5, 0.1),
    y=[u(t, step, eps, da) for t in np.arange(0, 5, 0.1)]) for step in
    range(100)]
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

def drawUIR(da, eps):
    data = [dict(
    visible=False,
    line=dict(color='#ff0000', width=3),
    name='r = ' + str(step),
    x=np.arange(100),
    y=[u(step, t, eps, da) for t in np.arange(100)]) for step in np.arange(0, 5, 0.1)]

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


da = getparams()
print("Текущие параметры")
print("R = {0} l = {1} k = {2}".format(da['R'], da['l'], da['k']))
print("alfa = {0} c = {1} betta = {2}".format(da['alf'], da['c'],da['betta']))
print("P = {0}, a = {1}, T={2}".format(da['P'], da['a'],da['T']))
yes = input("Введите Yes, если хотите изменить параметры ")
if yes.upper() == "YES":
    R = input("Введите R ")
    l = input("Введите l ")
    k = input("Введите k ")
    alfa = input("Введите alfa ")
    c = input("Введите c ")
    betta = input("Введите betta ")
    P = input("Введите P ")
    T = input("Введите T ")
    da['R'] = float(R)
    da['l'] = float(l)
    da['k'] = float(k)
    da['alf'] = float(alfa)
    da['c'] = float(c)
    da['betta'] = float(betta)
    da['P'] = float(P)
    da['a'] = float(R) / 5
    da['T'] = float(T)

eps = int(input("Введите точность (целоечисло -количество знаков после запятой)"))
print("Подождите, пожалуйста, выполняются вычисления...")
drawUIR(da, eps)
drawUIT(da, eps)
