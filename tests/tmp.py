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
import tskit
from IPython.display import SVG, display
import params
from params import Params


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
    #print(demography)
    samples = [
        msprime.SampleSet(cnt, population="AFR", time=0),
        msprime.SampleSet(cnt, population="EUR", time=0),
        msprime.SampleSet(cnt, population="ND", time=Params.t_ND)
    ]
    ts = msprime.sim_ancestry(samples=samples, demography=demography, random_seed=seed, ploidy=1,
                              sequence_length=Params.seq_len)
    ts = msprime.sim_mutations(ts, rate=Params.mu)
    if do_draw:  # Draw
        styles = []
        # Create a style for each population, programmatically (or just type the string by hand)
        for colour, p in zip(['red', 'green', 'blue', 'purple', 'black', 'orange'], ts.populations()):
            # target the symbols only (class "sym")
            s = f".node.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
            styles.append(s)
            # print(f'"{s}" applies to nodes from population {p.metadata["name"]} (id {p.id})')
        css_string = " ".join(styles)
        nd_labels = {}  # An array of labels for the nodes
        for n in ts.nodes():
            # Set sample node labels from metadata. Here we use the population name, but you might want
            # to use the *individual* name instead, if the individuals in your tree sequence have names
            if n.is_sample():
                nd_labels[n.id] = ts.population(n.population).metadata["name"]
                # print(ts.population(n.population))
        with open('../tree.svg', 'w') as outfile:  # Save
            outfile.write(
                SVG(ts.first().draw_svg(y_axis=True, size=(1000, 400), time_scale='rank', y_label=' ', style=css_string,
                                        node_labels=nd_labels)).data)

    return ts


def connected(m):
    for i in range(len(m) - 1):
        if m[i][1] == m[i + 1][0]:
            return True
    return False


def remove_one(m):
    mas = m
    while connected(mas) == True:
        for i in range(len(mas) - 1):
            if mas[i][1] == mas[i + 1][0]:
                mas[i][1] = mas[i + 1][1]
                mas.pop(i + 1)
                break
    return mas


def get_migrating_tracts_ind(ts, pop_name, ind, T_anc):
    pop = -1
    for i in ts.populations():
        if i.metadata['name'] == pop_name:
            pop = i.id

    mig = ts.tables.migrations
    migration_int = []

    for tree in ts.trees():  # перебираем все деревья. Как известно, каждому дереву отвечает участок днк
        anc_node = ind  # chose observable node
        while tree.time(tree.parent(
                anc_node)) <= T_anc:  # идем в прошлое до вершины anc_node по предкам нашего мексиканца, пока не наткнемся на миграцию
            anc_node = tree.parent(anc_node)
        migs = np.where(mig.node == anc_node)[0]  # выбирем все строки, соответствующие заданному узлу

        # идем по таблице миграций с anc_node и проверяем, чтобы миграции попадали в тот самый участок днк
        for i in migs:

            stroka = mig[i]
            if stroka.time == T_anc and stroka.dest == pop and tree.interval.left >= stroka.left and tree.interval.right <= stroka.right:
                migration_int.append([tree.interval.left, tree.interval.right])

    migration_int2 = []
    for i in range(len(migration_int)):
        if migration_int[i][0] != migration_int[i][1]:
            migration_int2.append(migration_int[i])
    migration_int = migration_int2

    mi = remove_one(migration_int)
    mi.sort()

    return mi


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
    # print("EUR-ND:", ts.nodes()[ts.first().mrca(eur_sample, nd_sample)].time)
    # print("EUR-AFR:", ts.nodes()[ts.first().mrca(eur_sample, afr_sample)].time)
    print(ts.nodes()[ts.first().mrca(eur_sample, nd_sample)].time, end= ' ')
    if (ts.nodes()[ts.first().mrca(eur_sample, nd_sample)].time < ts.nodes()[
        ts.first().mrca(eur_sample, afr_sample)].time):
        return True
    return False


config_probs = [
    0.97 * (1 - 2 / 3 * np.exp(-2.4)),
    0.97 * np.exp(-2.4) / 3,
    0.97 * np.exp(-2.4) / 3,
    0.03 * (1 - 2 / 3 * np.exp(-26)),
    0.03 * np.exp(-26) / 3,
    0.03 * np.exp(-26) / 3
]


def get_mutations(ts):  # this uses order!!!! AFR, EUR, ND
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
    eur_sample = 1
    afr_sample = -1
    n = ts.nodes()[eur_sample]
    ok = False
    while ts.first().parent(n.id) != tskit.NULL:
        n = ts.nodes()[ts.first().parent(n.id)]
        if n.population == 5:
            ok = True
    return ok


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


def generate_tests(cnt, require_admix=0, seed=0):
    tests = []
    i = 0
    while len(tests) < cnt:
        tests.append(Test(generate_demography(seed + i + 1, True)))
        i += 1
        ts = tests[-1].ts
        #print(ts.tables.populations)
        if not get_times(ts):
            tests.pop()
    print(get_migrating_tracts_ind(ts, "ND", 1, Params.t_admix), tests[-1].has_nd_ancestry, get_times(ts))
    return tests


generate_tests(1, 1, 1)
