#include "critical_features.h"

void solve(const char* infile, const char* outfile)
{
    std::vector<int> left, right;
    Network net = readCSV(infile, left, right);

    printf("solving %s\n", infile);
    printf("N = %d\n", net.n);
    printf("UCS = %lf\n", net.ucs);
    printf("perm = %e\n", net.permeability);
    printf("# pores = %d\n", (int)net.pores.size());
    printf("# throats = %d\n", (int)net.throats.size());
    printf("# pores on left boundary = %d\n", (int)left.size());
    printf("# pores on right boundary = %d\n", (int)right.size());

    // Get boundary pores, if none were provided
    //
    if (left.size() == 0 || right.size() == 0)
    {
        getBoundaries(net, direction, left, right);
        printf("Computing boundaries: # left = %d, # right = %d\n", (int)left.size(), (int)right.size());
    }

    // Remove critical pores until reaching a given target histogram
    //
    getBoundaries(net, 0 /*direction*/, left, right);
    Network new_net = fitCriticalPoreCNs(net, left, right, {2, 3, 4, 5, 9, 9, 10, 12, 14, 14});
    
    getBoundaries(net, 1 /*direction*/, left, right);
    new_net = fitCriticalPoreCNs(new_net, left, right, {3, 3, 5, 5, 8, 8, 8, 9, 9, 10, 14, 14, 14});

    getBoundaries(net, 2 /*direction*/, left, right);
    new_net = fitCriticalPoreCNs(new_net, left, right, {1, 2, 3, 3, 5, 5, 6, 8});
}


int main(int argc, char* argv[])
{
    srand(time(NULL));

    const char *infile = argv[1];
    const char *outfile = argv[1];

    int direction = 0;

    if (streq(infile + strlen(infile) - 4, ".csv"))
    {
        solve(infile, outfile);
    }
    else
    {
        const char* dirname = infile;
    }


    return 0;
}
