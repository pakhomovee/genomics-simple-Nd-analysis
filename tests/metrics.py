import with_mutations
import numpy as np
import math
from scipy import integrate
import params
Params = params.PARAMS()
def Pois(k, t):
    L = Params.mu * t * Params.seq_len
    return np.power(L, k) / math.factorial(k) * np.exp(-L)

def estimate1(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2) * Pois(k2, t2 - Params.t_ND) * Pois(k3, 2 * t1 - t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_ND *np.exp(-(t2-Params.t_admix)/Params.N_ND) * 1/Params.N_ANC * np.exp(-(t1-Params.t_all)/Params.N_ANC)
    def f1(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_ANC * np.exp(-(t2-Params.t_all)/Params.N_ANC * 3) * 1/Params.N_ANC * np.exp(-(t1-t2)/Params.N_ANC)
    I1 = integrate.dblquad(f, Params.t_admix, Params.t_all, Params.t_all, np.inf)
    I2 = integrate.dblquad(f1, Params.t_all, np.inf, lambda x : x, np.inf)
    return I1[0] + I2[0] * np.exp(-26000/Params.N_ND)

def estimate2(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, 2 * t1 - t2 - Params.t_ND) * Pois(k2, t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_ANC * np.exp(-(t2-Params.t_all)/Params.N_ANC * 3) * 1/Params.N_ANC * np.exp(-(t1-t2)/Params.N_ANC)
    I1 = integrate.dblquad(f, Params.t_all, np.inf, lambda x : x, np.inf)
    return I1[0]

def estimate3(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2 - Params.t_ND) * Pois(k2, 2 * t1 - t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_ANC * np.exp(-(t2-Params.t_all)/Params.N_ANC * 3) * 1/Params.N_ANC * np.exp(-(t1-t2)/Params.N_ANC)
    I1 = integrate.dblquad(f, Params.t_all, np.inf, lambda x : x, np.inf)
    return I1[0]

def estimateno1(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, 2 * t1 - t2 - Params.t_ND) * Pois(k2, t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_NND * np.exp(-t2/Params.N_NND) * 1/Params.N_ANC * np.exp(-(t1-Params.t_all)/Params.N_ANC)
    def f1(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_ANC * np.exp(-(t2-Params.t_all)/Params.N_ANC * 3) * 1/Params.N_ANC * np.exp(-(t1-t2)/Params.N_ANC)
    I1 = integrate.dblquad(f, 0, Params.t_all, Params.t_all, np.inf)
    I2 = integrate.dblquad(f1, Params.t_all, np.inf, lambda x : x, np.inf)
    return I1[0] + I2[0] * np.exp(-24/10)

def estimateno2(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2 - Params.t_ND) * Pois(k2, t2) * Pois(k3, 2 * t1 - t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1/Params.N_ANC * np.exp(-(t2-Params.t_all)/Params.N_ANC * 3) * 1/Params.N_ANC * np.exp(-(t1-t2)/Params.N_ANC)
    I1 = integrate.dblquad(f, Params.t_all, np.inf, lambda x : x, np.inf)
    return I1[0]

def estimateno3(k1, k2, k3):
    def prob(k1, k2, k3, t1, t2):  # returns P{k1, k2, k3 | t1, t2}
        return Pois(k1, t2 - Params.t_ND) * Pois(k2, 2 * t1 - t2) * Pois(k3, t2)
    def f(t1, t2):
        return prob(k1, k2, k3, t1, t2) * 1 / Params.N_ANC * np.exp(-(t2 - Params.t_all) / Params.N_ANC * 3) * 1 / Params.N_ANC * np.exp(-(t1 - t2) / Params.N_ANC)
    I1 = integrate.dblquad(f, Params.t_all, np.inf, lambda x: x, np.inf)
    return I1[0]

def cumulative_admix(k1, k2, k3):
    return estimate1(k1, k2, k3) * (1 - 2/3 * np.exp(-26000/Params.N_NND)) + estimate2(k1, k2, k3) * np.exp(-26000/Params.N_NND) / 3 + estimate3(k1, k2, k3) * np.exp(-26000/Params.N_NND) / 3

def cumulative_no_admix(k1, k2, k3):
    return estimateno1(k1, k2, k3) * (1 - 2/3 * np.exp(-24/10)) + estimateno2(k1, k2, k3) * np.exp(-24/10) / 3 + estimateno3(k1, k2, k3) * np.exp(-24/10) / 3

stats = [[0, 0], [0, 0]]
testid = 1
processed = 0
total_admix = 0
for test in with_mutations.generate_tests(Params.TESTCOUNT, False):
    k1, k2, k3 = with_mutations.get_mutations(test.ts)
    t1, t2 = with_mutations.get_times(test.ts)
    total_admix += test.has_nd_ancestry
    if testid % 100 == 0:
        print('=' * 10, f'Test #{testid}', '=' * 10)
    if testid % 1000 == 0:
        print(stats[0][0] / Params.TESTCOUNT, stats[0][1] / Params.TESTCOUNT)
        print(stats[1][0] / Params.TESTCOUNT, stats[1][1] / Params.TESTCOUNT)
    try:
        x = cumulative_admix(k1, k2, k3)
        y = cumulative_no_admix(k1, k2, k3)
        #print(f"#k1: {k1}, #k2: {k2}, #k3: {k3}, #admix_score: {x}")
        #print(f"#k1: {k1}, #k2: {k2}, #k3: {k3}, #no_admix_score: {y}")
        #print("PREDICTED PROBABILITY:", prob)
        #print('=' * 20)
        #print(f"predicted {x} vs {y}, answer {test.has_nd_ancestry}")
        prob = x / (x + y)
        stats[test.has_nd_ancestry][prob >= Params.threshold] += 1
        processed += 1
    except Exception as e:
        pass
    testid += 1
print(stats[0][0] / Params.TESTCOUNT, stats[0][1] / Params.TESTCOUNT)
print(stats[1][0] / Params.TESTCOUNT, stats[1][1] / Params.TESTCOUNT)
print(f"processed {processed} tests")
print(f"total admixed {total_admix} tests")