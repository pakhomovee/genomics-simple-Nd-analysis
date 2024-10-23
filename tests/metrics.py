import with_mutations
import numpy as np
import math
from scipy import integrate
Ne = 5000 #same for now
mu = 10 ** -8
t_ND = 1500 # neanderthal date
LEN = 100000
def Pois(k, t):
    L = mu * t * LEN
    return np.power(L, k) / math.factorial(k) * np.exp(-L)

def estimate1(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2) * Pois(k2, t2 - t_ND) * Pois(k3, 2 * t1 - t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/1000 *np.exp(-(t2-2000)/1000) * 1/11000 * np.exp(-(t1-28000)/11000)
    def f1(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/11000 * np.exp(-(t2-28000)/11000 * 3) * 1/11000 * np.exp(-(t1-t2)/11000)
    I1 = integrate.dblquad(f, 2000, 28000, 28000, np.inf)
    I2 = integrate.dblquad(f1, 28000, np.inf, lambda x : x, np.inf)
    return I1[0] + I2[0] * np.exp(-26000/1000)

def estimate2(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, 2 * t1 - t2 - t_ND) * Pois(k2, t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/11000 * np.exp(-(t2-28000)/11000 * 3) * 1/11000 * np.exp(-(t1-t2)/11000)
    I1 = integrate.dblquad(f, 28000, np.inf, lambda x : x, np.inf)
    return I1[0]

def estimate3(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2 - t_ND) * Pois(k2, 2 * t1 - t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/11000 * np.exp(-(t2-28000)/11000 * 3) * 1/11000 * np.exp(-(t1-t2)/11000)
    I1 = integrate.dblquad(f, 28000, np.inf, lambda x : x, np.inf)
    return I1[0]

def estimateno1(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2) * Pois(k2, t2 - t_ND) * Pois(k3, 2 * t1 - t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/1000 *np.exp(-(t2-2000)/1000) * 1/11000 * np.exp(-(t1-28000)/11000)
    def f1(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/11000 * np.exp(-(t2-28000)/11000 * 3) * 1/11000 * np.exp(-(t1-t2)/11000)
    I1 = integrate.dblquad(f, 2000, 28000, 28000, np.inf)
    I2 = integrate.dblquad(f1, 28000, np.inf, lambda x : x, np.inf)
    return I1[0] + I2[0] * np.exp(-26000/1000)

def estimateno2(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, 2 * t1 - t2 - t_ND) * Pois(k2, t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/11000 * np.exp(-(t2-28000)/11000 * 3) * 1/11000 * np.exp(-(t1-t2)/11000)
    I1 = integrate.dblquad(f, 28000, np.inf, lambda x : x, np.inf)
    return I1[0]

def estimateno3(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2 - t_ND) * Pois(k2, 2 * t1 - t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/11000 * np.exp(-(t2-28000)/11000 * 3) * 1/11000 * np.exp(-(t1-t2)/11000)
    I1 = integrate.dblquad(f, 28000, np.inf, lambda x : x, np.inf)
    return I1[0]

def cumulative_admix(k1, k2, k3):
    return estimate1(k1, k2, k3) * (1 - 2/3 * np.exp(-26000/10000)) + estimate2(k1, k2, k3) * np.exp(-26000/10000) / 3 + estimate3(k1, k2, k3) * np.exp(-26000/10000) / 3

def cumulative_no_admix(k1, k2, k3):
    return estimateno1(k1, k2, k3) * (1 - 2/3 * np.exp(-26000/10000)) + estimateno2(k1, k2, k3) * np.exp(-26000/10000) / 3 + estimateno3(k1, k2, k3) * np.exp(-26000/10000) / 3

testid = 1
for test in with_mutations.generate_tests(30):
    k1, k2, k3 = with_mutations.get_mutations(test.ts)
    t1, t2 = with_mutations.get_times(test.ts)
    print('=' * 10, f'Test #{testid}', '=' * 10)
    print(f"#k1: {k1}, #k2: {k2}, #k3: {k3}, #admix_score: {cumulative_admix(k1, k2, k3)}")
    print(f"#k1: {k1}, #k2: {k2}, #k3: {k3}, #no_admix_score: {cumulative_no_admix(k1, k2, k3)}")
    print('=' * 20)
    testid += 1