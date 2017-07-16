// Finds 
//
// Ex: ./reduce merged/merged_N_20.csv cp_cns.txt reduced/merged_N_20.reduced.csv
// Ex: ./reduce merged cp_cns.txt reduced 
// Ex: ./reduce reduced/merged_N_20.reduced.csv cp_cns.txt /dev/null
//
#include "critical_features.h"

std::vector<std::vector<int>> readCN(const char* infile)
{
    FILE* f = fopen(infile, "r");
    char line[3][100000];

    std::vector<std::vector<int>> cp_cns(3);
    for (int direction = 0; direction < 3; direction++)
    {
        fgets(line[direction], sizeof(line[direction]), f);

        int offset;
        int cn;
        const char* s = line[direction];
        while (sscanf(s, " %d%n", &cn, &offset) == 1)
        {
            cp_cns[direction].push_back(cn);
            s += offset;
        }

        printf("Target CNs in X direction: ");
        for (auto cn : cp_cns[direction])
        {
            printf("%d, ", cn);
        }
        printf("\n");
    }

    return cp_cns;
}


void solve(const char* infile, const char* cnfile, const char* outfile)
{
    printf("\n\n------------------- solving %s ------------------\n\n", infile);

    std::vector<int> left, right;
    Network net = readCSV(infile, left, right);

    printf("N = %d\n", net.n);
    printf("UCS = %lf\n", net.ucs);
    printf("perm = %e\n", net.permeability);
    printf("# pores = %d\n", (int)net.pores.size());
    printf("# throats = %d\n", (int)net.throats.size());
    printf("# pores on left boundary = %d\n", (int)left.size());
    printf("# pores on right boundary = %d\n", (int)right.size());

    std::vector<std::vector<int>> cp_cns = readCN(cnfile);

    // Remove critical pores until reaching a given target histogram
    //
    Network new_net = net;
    for (int direction = 0; direction < 3; direction++)
    {
        getBoundaries(net, direction, left, right);
        new_net = fitCriticalPoreCNs(new_net, left, right, cp_cns[direction]);
    }

    new_net.save(outfile);
}


int main(int argc, char* argv[])
{
    srand(time(NULL));

    const char *infile = argv[1];
    const char *cnfile = argv[2];
    const char *outfile = argv[3];

    if (streq(infile + strlen(infile) - 4, ".csv"))
    {
        solve(infile, cnfile, outfile);
    }
    else
    {
        const char* indir = infile;
        const char* outdir = outfile;

        std::vector<std::string> infiles = ls_csv(indir);
        for (auto infile : infiles)
        {
            std::size_t slash = infile.find_last_of("/\\");
            std::string infilename;
            if (slash == std::string::npos)
            {
                infilename = infile;
            }
            else
            {
                infilename = infile.substr(slash + 1);
            }

            std::string outfile = (std::string)(outdir) + "/" + infilename.substr(0, infilename.length() - 4) + ".reduced.csv";

            solve(infile.c_str(), cnfile, outfile.c_str());
        }
    }

    return 0;
}
