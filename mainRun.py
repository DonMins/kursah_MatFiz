import scipy.special as spc
import numpy as np
import  scipy as sp


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
# надо как нибудь оптимизировать
bessel0 = [0, ]
temp = spc.jn_zeros(1, 100)
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

    x = denominator *t / parameter_array['c']*parameter_array['l']*parameter_array['R']**2

    all = numerator/denominator
    return  all * (1- np.exp(-x))

print(An(1,1))
