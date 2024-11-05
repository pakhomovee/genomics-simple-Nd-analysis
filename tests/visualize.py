from metrics import  *
stats = [[0, 0], [0, 0]]
testid = 1
processed = 0
total_admix = 0
for test in with_mutations.generate_tests(Params.TESTCOUNT, 0):
    k1, k2, k3 = with_mutations.get_mutations(test.ts)
    t1, t2 = with_mutations.is_ok(test.ts)
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
