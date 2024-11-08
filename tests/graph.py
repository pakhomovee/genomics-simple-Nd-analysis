from matplotlib import pyplot as plt
import math
from multiprocessing import Process, Queue
from metrics import *


def get_data(start, end, step, output_queue):
    testid = 1
    processed = 0
    total_admix = 0
    stats = [[0, 0], [0, 0]]
    pts = {}
    pts1 = {}

    for ts in range(start, end, step):
        test = with_mutations.generate_tests(1, 0, 500 * ts + 1)[0]
        k1, k2, k3 = with_mutations.get_mutations(test.ts)
        t1, t2 = with_mutations.get_times(test.ts)
        total_admix += test.has_nd_ancestry

        if testid % 100 == 0:
            print('=' * 10, f'Test #{testid}', '=' * 10)
        if testid % 1000 == 0:
            print(stats[0][0] / Params.TESTCOUNT, stats[0][1] / Params.TESTCOUNT)
            print(stats[1][0] / Params.TESTCOUNT, stats[1][1] / Params.TESTCOUNT)

        try:
            x = max(cumulative_admix(k1, k2, k3), 0)
            y = max(cumulative_no_admix(k1, k2, k3), 0)
            prob = x / (x + y)
            if math.isnan(prob):
                continue
            stats[test.has_nd_ancestry][prob >= Params.threshold] += 1
            processed += 1
            if test.has_nd_ancestry:
                k = int(t2) // 1000
                if k not in pts:
                    pts[k] = [0, 0]
                pts[k][0] += prob
                pts[k][1] += 1
            else:
                k = int(t2) // 1000
                if k not in pts1:
                    pts1[k] = [0, 0]
                pts1[k][0] += prob
                pts1[k][1] += 1
        except Exception as e:
            pass

        testid += 1

    output_queue.put((pts, pts1, stats, processed))  # Put results in queue

if __name__ == '__main__':
    # Set up multiprocessing
    proc = []
    output_queue = Queue()

    for i in range(Params.threads):
        p = Process(target=get_data, args=(i, Params.TESTCOUNT, Params.threads, output_queue))
        p.start()
        proc.append(p)

    # Collect results from each process
    results = []
    for p in proc:
        p.join()
        results.append(output_queue.get())  # Get data from queue

    # Process the collected results (example: accumulate results)
    all_pts = []
    all_pts1 = []
    total_stats = [[0, 0], [0, 0]]
    total_processed = 0

    for pts, pts1, stats, processed in results:
        all_pts.append(pts)
        all_pts1.append(pts1)
        total_stats[0][0] += stats[0][0]
        total_stats[0][1] += stats[0][1]
        total_stats[1][0] += stats[1][0]
        total_stats[1][1] += stats[1][1]
        total_processed += processed

    print("Final results:")
    print("Total processed:", total_processed)
    #print("Points with ancestry:", all_pts)
    #print("Points without ancestry:", all_pts1)
    #print("Statistics:", total_stats)
    print(total_stats[0][0] / total_processed, total_stats[0][1] / total_processed)
    print(total_stats[1][0] / total_processed, total_stats[1][1] / total_processed)
    pts = []
    pts1 = []
    d = {}
    d1 = {}
    c = {}
    c1 = {}
    for arr in all_pts:
        for x in arr.keys():
            if x not in d:
                d[x] = 0
                c[x] = 0
            d[x] += arr[x][0]
            c[x] += arr[x][1]
    for arr in all_pts1:
        for x in arr.keys():
            if x not in d1:
                d1[x] = 0
                c1[x] = 0
            d1[x] += arr[x][0]
            c1[x] += arr[x][1]
    plt.ticklabel_format(useOffset=False)
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