// Same as conductance.py but faster
// usage: ./conductance [input file/dir] [output file] [direction]
// example #1: ./conductance.py input.csv output.csv 0
// example #2: ./conductance.py datadir output.csv 1
//
// Unlike conductance.py, this does NOT do the maximum conductance path (max_g = 0)
//
#include "critical_features.h"

void solve(const char* infile, const char* outfile, int direction)
{
    printf("\n\n------------------- solving %s ------------------\n\n", infile);

    std::vector<int> left, right;
    Network net = readCSV(infile, left, right);

    printf("N = %d\n", net.n);
    printf("UCS = %lf\n", net.ucs);
    printf("perm = %e\n", net.permeability);
    printf("# pores = %d\n", (int)net.pores.size());
    printf("# throats = %d\n", (int)net.throats.size());
    printf("direction = %d\n", direction);
    printf("# pores on left boundary = %d\n", (int)left.size());
    printf("# pores on right boundary = %d\n", (int)right.size());

    // Get boundary pores, if none were provided
    //
    if (left.size() == 0 || right.size() == 0)
    {
        getBoundaries(net, direction, left, right);
        printf("Computing boundaries: # left = %d, # right = %d\n", (int)left.size(), (int)right.size());
    }

    int flow;
    // Critical throats
    //
    std::vector<Edge> criticalThroats;
    flow = maxFlow(net, left, right, false /*doubleVertices*/, criticalThroats);
    std::cout<<"Critical throats (count = "<<flow<<") = [";
    for (auto it : criticalThroats)
    {
        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
    }
    std::cout<<"]\n\n";

    // Critical pores
    //
    std::vector<Edge> criticalPores;
    flow = maxFlow(net, left, right, true /*doubleVertices*/, criticalPores);
    std::cout<<"Critical pores (count = "<<flow<<") = [";
    for (auto it : criticalPores)
    {
        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
    }
    std::cout<<"]\n\n";

    // Append to output file
    //
    FILE *f = fopen(outfile, "a");
    fprintf(f, "%s,", infile);
    fprintf(f, "%e,", net.permeability);
    fprintf(f, "0,0, ,"); // no max conductance
    fprintf(f, "%d,", (int)criticalThroats.size());
    for (auto edge : criticalThroats)
    {
        fprintf(f, "%d-%d ", edge.first + 1, edge.second + 1);
    }
    fprintf(f, ",%d,", (int)criticalPores.size());
    for (auto edge : criticalPores)
    {
        assert(edge.first < edge.second);
        fprintf(f, "%d ", edge.first + 1);
    }
    fprintf(f, "\n");
    fclose(f);
}


int main(int argc, char* argv[])
{
    srand(time(NULL));

    const char *infile = argv[1];
    const char *outfile = argv[2];
    int direction = atoi(argv[3]);

    if (streq(infile + strlen(infile) - 4, ".csv"))
    {
        solve(infile, outfile, direction);
    }
    else
    {
        const char* dirname = infile;
        std::vector<std::string> infiles = ls_csv(dirname);

        for (auto infile : infiles)
        {
            solve(infile.c_str(), outfile, direction);
        }
    }

    return 0;
}
