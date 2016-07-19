1) Find a network with a given size and target CN distribution:

    python run_tosho.py input.txt solutions/ 1

2) Merge the coords...csv and neigh...csv files into format acceptable for conductance.py:

    python merge_outs.py solutions/coords_N_10_cost.....csv solutions/neigh_N_10_cost....csv conductance/datadir/blabla.csv

3) Find max conductance, num of critical pores / throats, etc

    cd conductance
    python conductance.py datadir/ output.csv
