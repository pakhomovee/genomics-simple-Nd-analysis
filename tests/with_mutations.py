'''
verified, correct
generate_demography(seed, do_draw, cnt)
> generates a coalescent tree of ancestry, population structure:
         ANC
    /          \
   NND           ND
  /   \           |
AFR    EUR_PURE   | 0.03 admixture rate
          \       /
             EUR
with fixed seed and cnt samples of AFR, EUR and ND. (whatever that means)
if do_draw is True, saves the tree to tree.svg
'''
import time
from random import randint as rnd
import msprime
import numpy as np
from IPython.display import SVG, display
import params
Params = params.PARAMS()
def generate_demography(seed, do_draw=False, cnt=1):
    # Define demography structure
    demography = msprime.Demography()
    demography.add_population(name="ANC", initial_size=Params.N_ANC)
    demography.add_population(name="EUR_PURE", initial_size=Params.N_EUR)
    demography.add_population(name="EUR", initial_size=Params.N_EUR)
    demography.add_population(name="AFR", initial_size=Params.N_AFR)
    demography.add_population(name="NND", initial_size=Params.N_NND)
    demography.add_population(name="ND", initial_size=Params.N_ND)
    demography.add_population_split(time=Params.t_all, derived=["ND", "NND"], ancestral="ANC")
    demography.add_population_split(time=Params.t_AFR_EUR, derived=["EUR_PURE", "AFR"], ancestral="NND")
    demography.add_admixture(Params.t_admix, derived="EUR", ancestral=["EUR_PURE", "ND"], proportions=[0.97, 0.03])
    demography.sort_events()
    samples = [
        msprime.SampleSet(cnt, population="AFR", time=0),
        msprime.SampleSet(cnt, population="EUR", time=0),
        msprime.SampleSet(cnt, population="ND", time=Params.t_ND)
    ]
    ts = msprime.sim_ancestry(samples=samples, demography=demography, random_seed=seed, ploidy=1, sequence_length=Params.seq_len)
    ts = msprime.sim_mutations(ts, rate=Params.mu)
    if do_draw: # Draw
        styles = []
        # Create a style for each population, programmatically (or just type the string by hand)
        for colour, p in zip(['red', 'green', 'blue', 'purple', 'black', 'orange'], ts.populations()):
            # target the symbols only (class "sym")
            s = f".node.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
            styles.append(s)
            #print(f'"{s}" applies to nodes from population {p.metadata["name"]} (id {p.id})')
        css_string = " ".join(styles)
        nd_labels = {}  # An array of labels for the nodes
        for n in ts.nodes():
            # Set sample node labels from metadata. Here we use the population name, but you might want
            # to use the *individual* name instead, if the individuals in your tree sequence have names
            if n.is_sample():
                nd_labels[n.id] = ts.population(n.population).metadata["name"]
                #print(ts.population(n.population))
        with open('../tree.svg', 'w') as outfile: # Save
            outfile.write(SVG(ts.first().draw_svg(y_axis=True, size=(1000,400), time_scale='rank', y_label=' ', style=css_string, node_labels=nd_labels)).data)
    return ts

def HAS_ND_ANCESTRY(ts):
    nd_sample = -1
    eur_sample = -1
    afr_sample = -1
    for n in ts.nodes():
        if n.is_sample():
            if ts.population(n.population).metadata["name"] == "ND":
                nd_sample = n.id
            elif ts.population(n.population).metadata["name"] == "EUR":
                eur_sample = n.id
            else:
                afr_sample = n.id
    #print("EUR-ND:", ts.nodes()[ts.first().mrca(eur_sample, nd_sample)].time)
    #print("EUR-AFR:", ts.nodes()[ts.first().mrca(eur_sample, afr_sample)].time)
    if (ts.nodes()[ts.first().mrca(eur_sample, nd_sample)]).time < ts.nodes()[ts.first().mrca(eur_sample, afr_sample)].time:
        return True
    return False

config_probs = [
    0.97 * (1 - 2/3 * np.exp(-2.4)),
    0.97 * np.exp(-2.4) / 3,
    0.97 * np.exp(-2.4) / 3,
    0.03 * (1 - 2/3 * np.exp(-26)),
    0.03 * np.exp(-26) / 3,
    0.03 * np.exp(-26) / 3
]

def get_mutations(ts): # this uses order!!!! AFR, EUR, ND
    k1, k2, k3 = 0, 0, 0
    for v in ts.variants():
        v = v.genotypes
        if v[0] == v[1] and v[0] != v[2]:
            k1 += 1
        elif v[0] == v[2] and v[0] != v[1]:
            k2 += 1
        elif v[1] == v[2] and v[0] != v[1]:
            k3 += 1
    return (k1, k2, k3)

def get_times(ts):
    nd_sample = -1
    eur_sample = -1
    afr_sample = -1
    for n in ts.nodes():
        if n.is_sample():
            if ts.population(n.population).metadata["name"] == "ND":
                nd_sample = n.id
            elif ts.population(n.population).metadata["name"] == "EUR":
                eur_sample = n.id
            else:
                afr_sample = n.id
    # print("EUR-ND:", ts.nodes()[ts.first().mrca(eur_sample, nd_sample)].time)
    # print("EUR-AFR:", ts.nodes()[ts.first().mrca(eur_sample, afr_sample)].time)
    t2 = (ts.nodes()[ts.first().mrca(eur_sample, nd_sample)]).time
    t1 = ts.nodes()[ts.first().mrca(eur_sample, afr_sample)].time
    return t1, t2

def calc_mutations_probs(ts, i):
    k1, k2, k3 = get_mutations(ts)
    prob = 0
    
    return prob
def estimate(ts):
    res = 0
    for i in range(3, 6):
        #  res += P(mutation|config)P(config)
        res += config_probs[i] * calc_mutations_probs(ts, i)
    return res

class Test:
    def __init__(self, ts):
        self.ts = ts
        self.has_nd_ancestry = HAS_ND_ANCESTRY(ts)

def generate_tests(cnt, require_admix = True):
    tests = []
    i = 0
    while len(tests) < cnt:
        tests.append(Test(generate_demography(i + 1)))
        i += 1
        if not HAS_ND_ANCESTRY(tests[-1].ts) and require_admix:
            tests.pop()
    print("GENERATED:")
    for test in tests:
        print(f"{test.has_nd_ancestry}: {estimate(test.ts)}, {get_mutations(test.ts)}")
    print('=' * 10)
    return tests

generate_demography(42, True, 1)