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
from IPython.display import SVG, display
def generate_demography(seed, do_draw=False, cnt=1):
    # Define demography structure
    demography = msprime.Demography()
    demography.add_population(name="ANC", initial_size=11_000)
    demography.add_population(name="EUR_PURE", initial_size=5_000)
    demography.add_population(name="EUR", initial_size=5_000)
    demography.add_population(name="AFR", initial_size=5_000)
    demography.add_population(name="NND", initial_size=10_000)
    demography.add_population(name="ND", initial_size=1_000)
    demography.add_population_split(time=28_000, derived=["ND", "NND"], ancestral="ANC")
    demography.add_population_split(time=4_000, derived=["EUR_PURE", "AFR"], ancestral="NND")
    demography.add_admixture(2_000, derived="EUR", ancestral=["EUR_PURE", "ND"], proportions=[0.97, 0.03])
    demography.sort_events()
    ts = msprime.sim_ancestry(samples={"AFR": cnt, "EUR": cnt, "ND": cnt}, demography=demography, random_seed=seed, ploidy=1)
    if do_draw: # Draw
        styles = []
        # Create a style for each population, programmatically (or just type the string by hand)
        for colour, p in zip(['red', 'green', 'blue', 'purple', 'black', 'orange'], ts.populations()):
            # target the symbols only (class "sym")
            s = f".node.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
            styles.append(s)
            print(f'"{s}" applies to nodes from population {p.metadata["name"]} (id {p.id})')
        css_string = " ".join(styles)
        nd_labels = {}  # An array of labels for the nodes
        for n in ts.nodes():
            # Set sample node labels from metadata. Here we use the population name, but you might want
            # to use the *individual* name instead, if the individuals in your tree sequence have names
            if n.is_sample():
                nd_labels[n.id] = ts.population(n.population).metadata["name"]
                print(ts.population(n.population))
        with open('../tree.svg', 'w') as outfile: # Save
            outfile.write(SVG(ts.first().draw_svg(y_axis=True, size=(1000,400), time_scale='rank', y_label=' ', style=css_string, node_labels=nd_labels)).data)
    return ts

'''
kind of verified
def NO_ADMIX_CONF1(ts):
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
    if (ts.nodes()[ts.first().mrca(eur_sample, afr_sample)]).time < ts.nodes()[ts.first().mrca(eur_sample, nd_sample)].time:
        return True
    return False

ITERS = 100000
cnt = 0
for seed in range(ITERS):
    ts = generate_demography(seed + 1, False, 1)
    if NO_ADMIX_CONF1(ts):
        cnt += 1
print(cnt / ITERS)
'''

'''
*** check that admixture rate is 0.03 indeed ***
*** check passed ***
total = 0
drifted = 0
for seed in range(100):
    ts = generate_demography(seed + 1, False, 1)
    nd_sample = -1
    eur_sample = -1
    for n in ts.nodes():
        if n.is_sample():
            if ts.population(n.population).metadata["name"] == "ND":
                nd_sample = n.id
    for n in ts.nodes():
        if n.is_sample():
            if ts.population(n.population).metadata["name"] == "EUR":
                total += 1
                if ts.population((ts.nodes()[ts.first().mrca(nd_sample, n.id)]).population).metadata["name"] == "ND":
                    drifted += 1
print(drifted / total)
'''