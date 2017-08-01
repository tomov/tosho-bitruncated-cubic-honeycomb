// Keep removing critical pores until we match a target critical pore coordination number histogram & critical throat count
//
// Usage: ./reduce [input file / directory] [coordination # file] [output file / directory] [fraction of pores to remove on each iteration] [how many terminal iterations] [coordScale]
// 
// Ex: ./reduce merged/merged_N_20.csv cp_cns.txt reduced/merged_N_20.reduced.csv 0.5 2 1000
// Ex: ./reduce merged cp_cns.txt reduced 0.5 2 1000
// Ex: ./reduce reduced/merged_N_20.reduced.csv cp_cns.txt /dev/null 0.5 2 1000
//
#include "critical_features.h"

std::vector<std::vector<int>> readCN(const char* infile, std::vector<int> &ct /*out*/, std::vector<int> &cp /*out*/)
{
    // Format: Each line = # of critical throats [space] # of critical pores [space] coordination numbers of critical pores, separated by spaces
    // One line per direction; directions are X, Y, Z
    //
    FILE* f = fopen(infile, "r");
    char line[3][100000];

    std::vector<std::vector<int>> cp_cns(3);
    ct.clear();
    cp.clear();
    for (int direction = 0; direction < 3; direction++)
    {
        fgets(line[direction], sizeof(line[direction]), f);

        int offset;
        int _ct, _cp, cn;
        const char* s = line[direction];

        sscanf(s, " %d%d%n", &_ct, &_cp, &offset); 
        ct.push_back(_ct);
        cp.push_back(_cp);
        s += offset;

        while (sscanf(s, " %d%n", &cn, &offset) == 1)
        {
            cp_cns[direction].push_back(cn);
            s += offset;
        }
        assert(_cp == cp_cns[direction].size());

        printf("Target CNs in direction %d: ", direction);
        for (auto cn : cp_cns[direction])
        {
            printf("%d, ", cn);
        }
        printf("\n");
    }

    return cp_cns;
}


void solve(const char* infile, const char* cnfile, const char* outfile, double frac, int terminal, double coordScale)
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

    std::vector<int> ct, cp;
    std::vector<std::vector<int>> cp_cns = readCN(cnfile, ct /*out*/, cp /*out*/);

    // Remove critical pores until reaching a given target histogram
    //
    Network new_net = net;
    std::vector<std::vector<Edge>> criticalThroats(3), criticalPores(3);
    for (int round = 0, remaining = 1000000000; remaining > 0; round++, remaining--)
    {
        printf("\n\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ROUND %d >>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n\n", round);
        for (int direction = 0; direction < 3; direction++)
        {
            printf("\n\n>>> direction = %d\n", direction);
            getBoundaries(net, direction, left, right, coordScale);
            new_net = fitCriticalPoreCNs(new_net, left, right, cp_cns[direction], frac, criticalThroats[direction] /*out*/, criticalPores[direction] /*out*/);

            // Once we reach the # of critical features in either direction, give it one more iteration
            //
            if (criticalThroats[direction].size() <= ct[direction] || criticalPores[direction].size() <= cp[direction])
            {
                if (remaining > terminal)
                {
                    remaining = terminal;
                    printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!! TERMINAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                }
            }
            if (frac == 0) // special case -- instead of doing an infinite loop, just exit (this is for sanity checking stuff)
            {
                remaining = 1;
            }
        }
    }

    new_net.writeCSV(outfile);
}


int main(int argc, char* argv[])
{
    srand(time(NULL));

    const char *infile = argv[1];
    const char *cnfile = argv[2];
    const char *outfile = argv[3];
    double frac = atof(argv[4]);
    int terminal = atoi(argv[5]);
    double coordScale = atof(argv[6]);

    if (streq(infile + strlen(infile) - 4, ".csv"))
    {
        solve(infile, cnfile, outfile, frac, terminal, coordScale);
    }
    else
    {
        const char* indir = infile;
        const char* outdir = outfile;

        std::vector<std::string> infiles = ls_csv(indir);
        for (auto infile : infiles)
        {
            printf("main loop at file %s\n", infile.c_str());
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
            printf("main loop truncated filename to %s\n", infilename.c_str());

            std::string outfile = (std::string)(outdir) + "/" + infilename.substr(0, infilename.length() - 4) + ".reduced.csv";
            printf("main loop outfile %s\n", outfile.c_str());

            solve(infile.c_str(), cnfile, outfile.c_str(), frac, terminal, coordScale);
        }
    }

    return 0;
}
