from matplotlib import pyplot as plt
import math
from metrics import  *
stats = [[0, 0], [0, 0]]
testid = 1
processed = 0
total_admix = 0
all = []
pts = []
pts1 = []
d = {}
c = {}
d1 = {}
c1 = {}
for ts in range(Params.TESTCOUNT):
    test = with_mutations.generate_tests(1, 1, ts)[0]
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
        prob = x / (x + y)
        if math.isnan(prob):
            continue
        stats[test.has_nd_ancestry][prob >= Params.threshold] += 1
        processed += 1
        if test.has_nd_ancestry:
            k = int(t2) // 1000
            if k not in d:
                d[k] = 0
                c[k] = 0
            d[k] += prob
            c[k] += 1
            all.append((t2, prob))
        else:
            k = int(t2) // 1000
            if k not in d1:
                d1[k] = 0
                c1[k] = 0
            d1[k] += prob
            c1[k] += 1
    except Exception as e:
        pass
    testid += 1
print(stats[0][0] / processed, stats[0][1] / processed)
print(stats[1][0] / processed, stats[1][1] / processed)
print(f"processed {processed} tests")
print(f"total admixed {total_admix} tests")
for key in d.keys():
    pts.append((key, d[key] / c[key]))
for key in d1.keys():
    pts1.append((key, d1[key] / c1[key]))
pts.sort()
pts1.sort()
plt.plot([i[0] for i in pts], [i[1] for i in pts])
plt.plot([i[0] for i in pts1], [i[1] for i in pts1])
print(pts)
print("========")
print(pts1)
plt.savefig('graph.jpeg')
plt.clf()
plt.scatter([i[0] for i in all], [i[1] for i in all])
plt.savefig('scat.jpeg')
'''plt.clf()
plt.style.use('_mpl-gallery')
plt.hist(ptsx, bins=100, linewidth=0.5, edgecolor="white")
plt.savefig('hist.jpeg')'''
# threshold = 0.2 => 50/50
# threshold = 0.1
'''
ignoring t2 > 40k
0.9193548387096774 0.03735144312393888
0.004244482173174873 0.03904923599320883



'''