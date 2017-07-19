# Usage for single file:
# python fit_dist.py [topology] [PR] [THR] [CN_vs_PR] [PR_vs_TH] [bin size] [# of fits] [output file]
#    bin size is for grouping the pore sizes in PR_vs_TH
#    if # of fits > 1, each output is named [output file].X, where X = the index of the fit
#
# Ex: python fit_dist.py top/SST1-topology-N_20.csv dist/SST1-PR_distribution-1x49.csv dist/SST1-THR_distribution-1x41.csv dist/SST1-CN_vs_PR-15x49.csv dist/SST1-PR_vs_TH-49x49x41.csv 3 1 fits/SST1.fit.csv
#
# Ex: python fit_dist.py SST1-topologies dist/SST1-PR_distribution-1x49.csv dist/SST1-THR_distribution-1x41.csv dist/SST1-CN_vs_PR-15x49.csv dist/SST1-PR_vs_TH-49x49x41.csv 3 1 SST1-fits
# Ex: python fit_dist.py SST1-topologies dist/SST1-PR_distribution-1x49.csv dist/SST1-THR_distribution-1x41.csv dist/SST1-CN_vs_PR-15x49.csv dist/SST1-PR_vs_TH-49x49x41.csv 1 1 SST1-fits

# Usage for directories:
# python fit_dist.py [input dir] [bin size] [# of fits] [output dir]
#    files in the input dir must have the following format:
#    [material name]-topology-....csv
#    [material name]-PR_distribution-....csv
#    [material name]-TH_distribution-....csv
#    [material name]-CN_vs_PR-....csv
#    [material name]-PR_vs_TH-....csv
#    where the topology is generated by run_tosho.py and merge_outs.py
#
# Ex: python fit_dist.py dist 3 1 fits
# Ex: python fit_dist.py dist 3 1 fits > wtf.m
# Ex: python fit_dist.py dist 3 5 fits
#
import sys
import os
import string
import numpy as np
import random
from conductance import readCSV, writeCSV

bin_size = 3

DO_PRINT = False 

one_pr_per_th = False


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
    pr_pr_th = dict()
    pr_th = dict()
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(',')
            line = [l.strip() for l in line]

            prs = (int(line[0]) / bin_size, int(line[1]) / bin_size)
            line = line[2:]

            #prs = (float(prs[0]) * bin_size * 1e-6, float(prs[1]) * bin_size * 1e-6)
            if prs in pr_pr_th:
                for i in range(len(pr_pr_th[prs])):
                    pr_pr_th[prs][i] += float(line[i])
            else:
                pr_pr_th[prs] = [float(l) for l in line]

            if prs[0] in pr_th:
                for i in range(len(pr_th[prs[0]])):
                    pr_th[prs[0]][i] += float(line[i])
            else:
                pr_th[prs[0]] = [float(l) for l in line]


    # normalize them
    for prs, _ in pr_pr_th.iteritems():
        s = sum(pr_pr_th[prs])
        if s != 0:
            pr_pr_th[prs] = [x / s for x in pr_pr_th[prs]]

        s = sum(pr_pr_th[prs])
        assert abs(1 - s) < 1e-4 or abs(s) < 1e-4

    for pr, _ in pr_th.iteritems():
        s = sum(pr_th[pr])
        if s != 0:
            pr_th[pr] = [x / s for x in pr_th[pr]]

        s = sum(pr_th[pr])
        assert abs(1 - s) < 1e-4 or abs(s) < 1e-4

    return pr_pr_th, pr_th


def sample(dist):
    assert abs(sum(dist) - 1) < 1e-4
    while True:
        x = random.random()
        for j in range(len(dist)):
            if x < dist[j]:
                return j
            x -= dist[j]
    assert False


def draw(hist): # pass by reference! modifies hist in place
    for j in range(len(hist)):
        if hist[j] > 0:
            hist[j] -= 1
            return j
    return None 

def draw2(hist1, hist2): # pass by reference! modifies hist in place
    assert len(hist1) == len(hist2)
    for j in range(len(hist1)):
        if hist1[j] > 0 and hist2[j] > 0:
            hist1[j] -= 1
            hist2[j] -= 1
            return j
    return None 

def solve(infile, PR_file, THR_file, CN_vs_PR_file, PR_vs_TH_file, bin_size, n_fits, outfile):
    print '\n\n ======================= solving for ', infile, ' ================\n\n'

    for fit in range(n_fits):
        print '\n\n -------------- FIT ', fit, ' -------------\n\n'

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

        pr_pr_th, pr_th = read_PR_vs_TH_file(PR_vs_TH_file)
        assert len(pr_pr_th) == (len(pr) / bin_size + (bin_size != 1))**2
        assert len(pr_pr_th[(1, 1)]) == len(thr)

        # Assign pore radii according to CN_vs_PR
        #
        cns = [0] * len(pores)
        for throat in throats:
            cns[throat[0]] += 1
            cns[throat[1]] += 1

        # build histogram and draw from there
        #
        cn_pr_hist = []
        for i in range(len(cn_pr)):
            cnt = sum([1 for x in cns if x == i])
            cn_pr_hist.append([int(x * cnt) for x in cn_pr[i]])
            if DO_PRINT:
                print 'for cn ', i, ' (', cnt, '): ', cn_pr[i]

        # assign pore radii in random order
        #
        idxs = range(len(pores))
        random.shuffle(idxs)

        randoms = 0
        for i in idxs:
            cn = cns[i]        

            _ = draw(cn_pr_hist[cn]) # crucial to pass by reference!
            if _ is not None:
                new_r = pr[_]
            else: # out of pore radii to give away
                randoms += 1
                dist = cn_pr[cn]
                new_r = pr[sample(dist)]

            pores[i] = (pores[i][0], pores[i][1], pores[i][2], new_r)

        print 'Assigned', randoms, 'out of', len(pores), 'pore radii randomly'

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

        # Build histogram for throat radii
        #
        thr_hist = [int(x * len(throats)) for x in thr_dist]

        # Assign throat radii according to PR_vs_TH
        #
        if not one_pr_per_th:

            # get scaling factors (throat counts for pairs of pore radii)
            #
            cnts = dict()
            for i in range(len(throats)):
                p1 = throats[i][0]
                p2 = throats[i][1]
                pr1 = pores[p1][3]
                pr2 = pores[p2][3]
                key = (int(pr1 * 1e6) / bin_size, int(pr2 * 1e6) / bin_size)
    
                if key not in cnts:
                    cnts[key] = 0
                cnts[key] += 1
    
            # build histogram and draw from there
            #
            pr_pr_th_hist = dict()
            for key, hist in pr_pr_th.iteritems():
                cnt = cnts[key] if key in cnts else 0
                pr_pr_th_hist[key] = [int(x * cnt) for x in hist]
                if DO_PRINT:
                    print 'for throat pr-pr ', key, ' (', cnt, '): ', pr_pr_th_hist[key]
    
            # assign throat radii in random order
            #
            idxs = range(len(throats))
            random.shuffle(idxs)
    
            generics = 0
            single_prs = 0
            randoms = 0
            generic_hists = 0
            for i in idxs:
                p1 = throats[i][0]
                p2 = throats[i][1]
                pr1 = pores[p1][3]
                pr2 = pores[p2][3]
                key = (int(pr1 * 1e6) / bin_size, int(pr2 * 1e6) / bin_size)
    
                #_ = draw(pr_pr_th_hist[key]) # crucial to pass by reference!
                _ = draw2(pr_pr_th_hist[key], thr_hist) # crucial to pass by reference!
                #_ = draw(thr_hist) # sanity
                if _ is not None:
                    new_r = thr[_]
                else: # out of throat radii to give away
                    _ = draw(thr_hist)
                    if _ is not None:
                        new_r = thr[_]
                        generic_hists += 1
                    else:
                        randoms += 1
                        dist = pr_pr_th[key] # throat radius distribution
                        if abs(sum(dist)) < 1e-4:
                            dist = pr_th[int((pr1 + pr2) / 2 * 1e6) / bin_size]
                            if abs(sum(dist)) < 1e-4:
                                dist = thr_dist
                                generics += 1
                            else:
                                single_prs += 1
        
                        assert len(dist) == len(thr)
                        new_r = thr[sample(dist)] # throat radius
    
                throats[i] = (throats[i][0], throats[i][1], new_r)
    
            print 'Out of %d throat radii, assigned %d from generic histogram and %d randomly; out of the random ones, used single-PR distirbution for %d and prior distribution for %d' % (len(throats), generic_hists, randoms, single_prs, generics)


        else:

            # get scaling factors (throat counts for pairs of pore radii)
            #
            cnts = dict()
            for i in range(len(throats)):
                p1 = throats[i][0]
                p2 = throats[i][1]
                pr1 = pores[p1][3]
                pr2 = pores[p2][3]
                key = int((pr1 + pr2) / 2 * 1e6) / bin_size

                if key not in cnts:
                    cnts[key] = 0
                cnts[key] += 1

            # build histogram and draw from there
            #
            pr_th_hist = dict()
            for key, hist in pr_th.iteritems():
                cnt = cnts[key] if key in cnts else 0
                pr_th_hist[key] = [int(x * cnt) for x in hist]
                if DO_PRINT:
                    print 'for throat pr ', key, ' (', cnt, '): ', pr_th_hist[key]

            # assign throat radii in random order
            #
            idxs = range(len(throats))
            random.shuffle(idxs)

            generics = 0
            randoms = 0
            generic_hists = 0
            for i in idxs:
                p1 = throats[i][0]
                p2 = throats[i][1]
                pr1 = pores[p1][3]
                pr2 = pores[p2][3]
                key = int((pr1 + pr2) / 2 * 1e6) / bin_size

                #_ = draw(pr_th_hist[key]) # crucial to pass by reference!
                _ = draw2(pr_th_hist[key], thr_hist) # crucial to pass by reference!
                #_ = draw(thr_hist) # sanity
                if _ is not None:
                    new_r = thr[_]
                else: # out of throat radii to give away
                    _ = draw(thr_hist)
                    if _ is not None:
                        new_r = thr[_]
                        generic_hists += 1
                    else:
                        randoms += 1
                        dist = pr_th[key] # throat radius distribution
                        if abs(sum(dist)) < 1e-4:
                            dist = thr_dist
                            generics += 1

                        assert len(dist) == len(thr)
                        new_r = thr[sample(dist)] # throat radius

                throats[i] = (throats[i][0], throats[i][1], new_r)

            print 'Out of %d throat radii, assigned %d from generic histogram and %d randomly; out of these, used prior distribution for %d' % (len(throats), generic_hists, randoms, generics)


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
            fit_outfile = outfile[:-4] + "." + str(fit) + ".csv"
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

    if infile.lower().endswith('.csv'):
        solve(infile, PR_file, THR_file, CN_vs_PR_file, PR_vs_TH_file, bin_size, n_fits, outfile)

    else:
        indir = infile
        outdir = outfile

        all_files = []
        for (path, dirs, files) in os.walk(indir):
            print '\n============ EXPLORING DIRECTORY', path, '=================\n'
            for filename in files:
                if filename.lower().endswith('.csv'):
                    all_files.append((path, filename))
                    print 'FILE: ', path, filename

        for path, filename in all_files:
            infile = os.path.join(path, filename)
            outfile = os.path.join(outdir, filename[:-4] + '.fit.csv')
            solve(infile, PR_file, THR_file, CN_vs_PR_file, PR_vs_TH_file, bin_size, n_fits, outfile)
