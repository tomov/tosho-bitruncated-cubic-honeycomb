# Usage:
# Ex: python fit_dist.py dist/merged_SST1_N_20.csv dist/SST1-PR_distribution-1x49.csv dist/SST1-THR_distribution-1x41.csv dist/SST1-CN_vs_PR-15x49.csv dist/SST1-PR_vs_TH-49x49x41.csv 1 fits/SST1.fit.csv
#
import sys
import string
import numpy as np
from conductance import readCSV

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

            prs = (float(line[0]) * 1e-6, float(line[1]) * 1e-6)
            line = line[2:]
            pr_th[prs] = [float(l) for l in line]

            s = sum(pr_th[prs]) # TODO why are they not normalized?
            if s != 0:
                pr_th[prs] = [x / s for x in pr_th[prs]]

            s = sum(pr_th[prs])
            assert abs(1 - s) < 1e-4 or abs(s) < 1e-4
    return pr_th


if __name__ == '__main__':
    infile = sys.argv[1]
    PR_file = sys.argv[2]
    THR_file = sys.argv[3]
    CN_vs_PR_file = sys.argv[4]
    PR_vs_TH_file = sys.argv[5]
    n_fits = int(sys.argv[6])
    outfile = sys.argv[7]

    pores, throats, left, right, ucs, n, permeability = readCSV(infile)

    pr, pr_dist = read_PR_file(PR_file)

    thr, thr_dist = read_THR_file(THR_file)

    cn_pr = read_CN_vs_PR_file(CN_vs_PR_file)
    assert len(cn_pr) == 15
    assert len(cn_pr[0]) == len(pr)

    pr_th = read_PR_vs_TH_file(PR_vs_TH_file)
    assert len(pr_th) == len(pr)**2
    assert len(pr_th[(1e-6, 1e-6)]) == len(thr)
