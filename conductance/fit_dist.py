# Usage:
# Ex: python fit_dist.py dist/merged_SST1_N_20.csv dist/SST1-PR_distribution-1x49.csv dist/SST1-THR_distribution-1x41.csv dist/SST1-CN_vs_PR-15x49.csv dist/SST1-PR_vs_TH-49x49x41.csv 3 1 fits/SST1.fit.csv
#
import sys
import string
import numpy as np
import random
from conductance import readCSV, writeCSV

bin_size = 3

DO_PRINT = False

def read_PR_file(filename):
    rs = []
    dist = []
    with open(filename, 'r') as f:
        line = f.readline()
        line = line.split(',')
        line = [l.strip() for l in line]
        for i in range(len(line)):
            r = i + 1; # radius in microns
            rs.append(r * 1e-6)
            dist.append(float(line[i]))
        assert abs(1 - sum(dist)) < 1e-4
    return rs, dist

def read_THR_file(filename):
    return read_PR_file(filename)

def read_CN_vs_PR_file(filename):
    cn_pr = dict()
    with open(filename, 'r') as f:
        cn = 0
        for line in f:
            line = line.split(',')
            line = [l.strip() for l in line]
            cn_pr[cn] = [float(l) for l in line]
            assert abs(1 - sum(cn_pr[cn])) < 1e-4
            cn += 1
    return cn_pr

def read_PR_vs_TH_file(filename):
    pr_th = dict()
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(',')
            line = [l.strip() for l in line]

            prs = (int(line[0]) / bin_size, int(line[1]) / bin_size)
            line = line[2:]

            #prs = (float(prs[0]) * bin_size * 1e-6, float(prs[1]) * bin_size * 1e-6)
            if prs in pr_th:
                for i in range(len(pr_th[prs])):
                    pr_th[prs][i] += float(line[i])
            else:
                pr_th[prs] = [float(l) for l in line]

    # normalize them
    for prs, _ in pr_th.iteritems():
        s = sum(pr_th[prs])
        if s != 0:
            pr_th[prs] = [x / s for x in pr_th[prs]]

        s = sum(pr_th[prs])
        assert abs(1 - s) < 1e-4 or abs(s) < 1e-4

    return pr_th



def solve(infile, PR_file, THR_file, CN_vs_PR_file, PR_vs_TH_file, bin_size, n_fits, outfile):
    for fit in range(n_fits):
        print '\n\n\n ----------------------- FIT ', fit, ' ----------------\n\n\n'

        # Read inputs
        #
        pores, throats, left, right, ucs, n, permeability = readCSV(infile)

        pr, pr_dist = read_PR_file(PR_file)
        assert len(pr) == len(pr_dist)

        thr, thr_dist = read_THR_file(THR_file)
        assert len(thr) == len(thr_dist)

        cn_pr = read_CN_vs_PR_file(CN_vs_PR_file)
        assert len(cn_pr) == 15
        assert len(cn_pr[0]) == len(pr)

        pr_th = read_PR_vs_TH_file(PR_vs_TH_file)
        assert len(pr_th) == (len(pr) / bin_size + (bin_size != 1))**2
        assert len(pr_th[(1, 1)]) == len(thr)

        # Assign pore radii according to CN_vs_PR
        #
        cns = [0] * len(pores)
        for throat in throats:
            cns[throat[0]] += 1
            cns[throat[1]] += 1

        for i in range(len(pores)):
            dist = cn_pr[cns[i]] # pore radius distribution
            new_r = -1 # pore radius
            while new_r == -1:
                x = random.random()
                for j in range(len(pr)):
                    if x < dist[j]:
                        new_r = pr[j]
                        break
                    x -= dist[j]
            pores[i] = (pores[i][0], pores[i][1], pores[i][2], new_r)

        # sanity check -- should be same as CN_vs_PR
        #
        if DO_PRINT:
            cn_pr_sanity = dict()
            for cn in range(15):
                cn_pr_sanity[cn] = [0] * len(pr)
            for i in range(len(pores)):
                cn_pr_sanity[cns[i]][int(pores[i][3] * 1e6) - 1] += 1
            for cn in range(15):
                s = sum(cn_pr_sanity[cn])
                cn_pr_sanity[cn] = [float(x) / s for x in cn_pr_sanity[cn]]
                print cn, ': ', ','.join(["%.3f" % x for x in cn_pr_sanity[cn]])

        # Assign throat radii
        #
        generics = 0
        for i in range(len(throats)):
            p1 = throats[i][0]
            p2 = throats[i][1]
            pr1 = pores[p1][3]
            pr2 = pores[p2][3]
            key = (int(pr1 * 1e6) / bin_size, int(pr2 * 1e6) / bin_size)

            dist = pr_th[key] # throat radius distribution
            if abs(sum(dist)) < 1e-4:
                generics += 1
                dist = thr_dist

            new_r = -1 # throat radius
            x = random.random()
            for j in range(len(pr)):
                if x < dist[j]:
                    new_r = pr[j]
                    break
                x -= dist[j]
            assert new_r != -1

            throats[i] = (throats[i][0], throats[i][1], new_r)

        print 'Used THR for %d out of %d throats' % (generics, len(throats))

        #
        # get stats for plotting
        #

        # Plot #1: Average pore radius vs coordination number
        #
        avg_pr_vs_cn_plot = []
        for cn in range(15):
            radii = []
            for i in range(len(pores)):
                if cns[i] == cn:
                    radii.append(pores[i][3])
            avg_pr_vs_cn_plot.append(np.mean(radii))
        print '\n\nfigure;'
        print 'x = [', ' '.join(str(x) for x in range(15)), '];'
        print 'y = [', ' '.join(str(x) for x in avg_pr_vs_cn_plot), '];'
        print 'scatter(x, y);'
        print "xlabel('Coordination number');"
        print "ylabel('Avg pore radius');"

        # Plot #2: Average throat radius vs pore radius
        #
        avg_thr_vs_pr_plot = dict()
        for i in range(len(throats)):
            p1 = throats[i][0]
            p2 = throats[i][1]
            pr1 = pores[p1][3]
            pr2 = pores[p2][3]

            if pr1 not in avg_thr_vs_pr_plot:
                avg_thr_vs_pr_plot[pr1] = []
            if pr2 not in avg_thr_vs_pr_plot:
                avg_thr_vs_pr_plot[pr2] = []

            avg_thr_vs_pr_plot[pr1].append(throats[i][2])
            avg_thr_vs_pr_plot[pr2].append(throats[i][2])
        print '\n\nfigure;'
        print 'x = [', ' '.join(str(x) for x, _ in avg_thr_vs_pr_plot.iteritems()), '];'
        print 'y = [', ' '.join(str(np.mean(x)) for _, x in avg_thr_vs_pr_plot.iteritems()), '];'
        print 'scatter(x, y);'
        print "xlabel('Pore radius');"
        print "ylabel('Average throat radius');"

        if n_fits == 1:
            fit_outfile = outfile
        else:
            fit_outfile = outfile + "." + str(fit)
        writeCSV(pores, throats, left, right, ucs, n, permeability, fit_outfile)


if __name__ == '__main__':
    infile = sys.argv[1]
    PR_file = sys.argv[2]
    THR_file = sys.argv[3]
    CN_vs_PR_file = sys.argv[4]
    PR_vs_TH_file = sys.argv[5]
    bin_size = int(sys.argv[6])
    n_fits = int(sys.argv[7])
    outfile = sys.argv[8]

    solve(infile, PR_file, THR_file, CN_vs_PR_file, PR_vs_TH_file, bin_size, n_fits, outfile)
