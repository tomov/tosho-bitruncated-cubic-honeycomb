// Same as conductance.py but faster
//
#include <boost/config.hpp>
#include <iostream>
#include <cstdio>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>

const bool DEBUG = false;
const int MAX_CN = 14;

void replace(std::string& line, char what, char with)
{
    for (int i = 0; i < line.size(); i++)
    {
        if (line[i] == what)
        {
            line[i] = with;
        }
    }
}

bool streq(const char *str1, const char *str2)
{
    return strcmp(str1, str2) == 0;
}

struct Pore
{
    double x, y, z, r;

    Pore() {}
    Pore(double _x, double _y, double _z, double _r)
        : x(_x), y(_y), z(_z), r(_r) {}

    double operator[](int direction) const
    {
        switch (direction)
        {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                return 1.0/0;
        }
    }
};

struct Throat
{
    int p1, p2;
    double r;

    Throat() {}
    Throat(int _p1, int _p2, double _r)
        : p1(_p1), p2(_p2), r(_r) {}
};

struct Network
{
    std::vector<Pore> pores;
    std::vector<Throat> throats;
    double ucs;
    int n;
    double permeability;
};

Network readCSV(const char* infile, std::vector<int> &left /*out*/, std::vector<int> &right)
{
    Network net;
    left.clear();
    right.clear();

    // Row format:
    //  0 1 2  3     4     5   6  7  8  9  10  11  12  13       14  15       16  17
    // [X,Y,Z,Pore1,Pore2,ThR,Pr,LB,RB,UCS,N,Perm, Qx, QXtotal, Qy, QYtotal, Qz, QZtotal] 
    //

    FILE* f = fopen(infile, "r");
    char l[10000];
    char buf[10000];
    char* line[100];
    for (int i = 0; i < sizeof(line) / sizeof(line[0]); i++)
    {
        line[i] = &buf[i * 100];
    }

    // Read each l from the input file
    //
    while (fgets(l, sizeof(l), f))
    {
        char *s = l, *e;
        size_t len = strlen(l);
        int n = 0;

        // Parse the l
        //
        while (s < l + len)
        {
            e = strchr(s, ',');
            if (e == nullptr)
            {
                e = l + len;
            }
            e[0] = '\0';
            sscanf(s, "%s", line[n]);
            s = e + 1;
            n++;
        }

        // Convert each input to the right data structure
        //
        //for (int i = 0; i < n; i++) printf("%s,", line[i]);
        //
        if (!streq(line[0], "") && !streq(line[0], "0"))
        {
            Pore pore(atof(line[0]) * 1e-6, atof(line[1]) * 1e-6, atof(line[2]) * 1e-6, atof(line[6]) * 1e-6);
            net.pores.push_back(pore);
            if (DEBUG) printf("%lf %lf %lf: %lf,", pore.x, pore.y, pore.z, pore.r);
        }
        if (!streq(line[3], "") && !streq(line[3], "0"))
        {
            Throat throat(atoi(line[3]) - 1, atoi(line[4]) - 1, atof(line[5]) * 1e-6);
            net.throats.push_back(throat);
            if (DEBUG) printf("%d %d: %lf,", throat.p1, throat.p2, throat.r);
        }
        if (!streq(line[7], "") && !streq(line[7], "0"))
        {
            int idx = atoi(line[7]) - 1;
            left.push_back(idx);
            if (DEBUG) printf("%d,", idx);
        }
        if (!streq(line[8], "") && !streq(line[8], "0"))
        {
            int idx = atoi(line[8]) - 1;
            right.push_back(idx);
            if (DEBUG) printf("%d,", idx);
        }
        if (!streq(line[9], "") && !streq(line[9], "0"))
        {
            net.ucs = atof(line[9]);
            if (DEBUG) printf("%lf,", net.ucs);
        }
        if (!streq(line[10], "") && !streq(line[10], "0"))
        {
            net.n = atoi(line[10]);
            if (DEBUG) printf("%d,", net.n);
        }
        if (!streq(line[11], "") && !streq(line[11], "0"))
        {
            net.permeability = atof(line[11]);
            if (DEBUG) printf("%lf,", net.permeability);
        }

        if (DEBUG) printf("\n");
    }
    fclose(f);

    return net;
}


void getBoundaries(const Network &net, int direction, std::vector<int> &starting /*out*/, std::vector<int> &ending /*out*/)
{
    starting.clear();
    ending.clear();
    std::vector<double> left_b(net.pores.size()), right_b(net.pores.size());
    double min_b = net.pores[0][direction], max_b = net.pores[0][direction];
    for (int i = 0; i < net.pores.size(); i++)
    {
        min_b = std::min(min_b, net.pores[i][direction]);
        max_b = std::max(max_b, net.pores[i][direction]);
        left_b[i] = net.pores[i][direction] - net.pores[i].r;
        right_b[i] = net.pores[i][direction] + net.pores[i].r;
    }
    for (int i = 0; i < net.pores.size(); i++)
    {
        if (left_b[i] <= min_b)
        {
            starting.push_back(i);
        }
        if (right_b[i] >= max_b)
        {
            ending.push_back(i);
        }
    }
}

typedef std::pair<int, int> Edge;

using namespace boost;

template <class Graph, class ResCapMap>
filtered_graph<Graph, is_residual_edge<ResCapMap> >
residual_graph(Graph& g, ResCapMap residual_capacity) {
  return filtered_graph<Graph, is_residual_edge<ResCapMap> >
	(g, is_residual_edge<ResCapMap>(residual_capacity));
}

int maxFlow(const Network &net, const std::vector<int> &left, const std::vector<int> &right, bool doubleVertices, std::vector<Edge> &critical /*out*/)
{

    // Transform into graph and double the edges
    //
    std::vector<Pore> verts = net.pores;
    std::vector<Edge> edges;
    for (auto throat : net.throats)
    {
        edges.push_back(Edge(throat.p1, throat.p2));
        edges.push_back(Edge(throat.p2, throat.p1));
    }

    // Add source & sink
    //
    verts.push_back(Pore(0, 0, 0, 0));
    verts.push_back(Pore(0, 0, 0, 0));
    int source = verts.size() - 2;
    int sink = verts.size() - 1;
    printf("source = %d, sink = %d\n", source, sink);

    for (int i = 0; i < left.size(); i++)
    {
        edges.push_back(Edge(source, left[i]));
    }
    for (int i = 0; i < right.size(); i++)
    {
        edges.push_back(Edge(right[i], sink));
    }

    // Add rev edges
    //
    int edges_n = edges.size();
    std::vector<int> cap, /*flow, */rev;
    for (int i = 0; i < edges_n; i++)
    {
        cap.push_back(1);
        //flow.push_back(0);
        rev.push_back(i + edges_n);
    }
    for (int i = 0; i < edges_n; i++)
    {
        int u = edges[i].first;
        int v = edges[i].second;
        edges.push_back(Edge(v, u));
        //flow.push_back(0);
        cap.push_back(0);
        rev.push_back(i);
    }

    // Make capacities to/from source/sink infinite
    //
    for (int i = 0; i < edges.size(); i++)
    {
        int u = edges[i].first;
        int v = edges[i].second;
        if (u == sink || u == source || v == sink || v == source)
        {
            cap[i] = 1000000000;
        }
    }

    // Double vertices
    //
    int verts_n = verts.size();
    if (doubleVertices)
    {
        for (int i = 0; i < verts_n; i++)
        {
            verts.push_back(verts[i]);
        }
        
        // Make existing edges u_out -> v_in
        for (int i = 0; i < edges_n; i++)
        {
            Edge e = edges[i];
            Edge new_edge = Edge(e.first + verts_n, e.second);
            edges[i] = new_edge;

            // Reverse edge
            e = edges[i + edges_n];
            assert(e.second + verts_n == new_edge.first);
            assert(e.first == new_edge.second);
            new_edge = Edge(e.first, e.second + verts_n);
            edges[i + edges_n] = new_edge;
        }

        source += verts_n;

        // Add in-out edges
        for (int u = 0; u < verts_n; u++)
        {
            int v = u + verts_n;

            edges.push_back(Edge(u, v));
            //flow.push_back(0);
            cap.push_back(1);
            rev.push_back(edges.size());

            // Add reverse edge
            edges.push_back(Edge(v, u));
            //flow.push_back(0);
            cap.push_back(0);
            rev.push_back(edges.size() - 2);
        }
    }

    // Build adjacency list
    //
    std::vector<std::vector<int> > adj(verts.size());
    for (int i = 0; i < edges.size(); i++)
    {
        int u = edges[i].first;
        adj[u].push_back(i);
    }

    //
    //
    // Convert to Boost Graph structure 
    //
    //

    typedef adjacency_list_traits<vecS, vecS, directedS> Traits;
    typedef adjacency_list<vecS, vecS, directedS, 
        property<vertex_name_t, std::string,
            property<vertex_index_t, long,
                  property<vertex_color_t, boost::default_color_type,
                          property<vertex_distance_t, long,
                                    property<vertex_predecessor_t, Traits::edge_descriptor > > > > >,
        property<edge_capacity_t, int,
            property<edge_residual_capacity_t, int,
                property<edge_reverse_t, Traits::edge_descriptor> > >
    > Graph;
    
    Graph g;
    int flow;
    
    property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
    property_map<Graph, edge_reverse_t>::type reverse = get(edge_reverse, g);
    property_map<Graph, edge_residual_capacity_t>::type residual_capacity = get(edge_residual_capacity, g);

    // Add verts
    std::vector<Traits::vertex_descriptor> vertex_descriptors(verts.size());
    for (int i = 0; i < verts.size(); i++)
    {
        vertex_descriptors[i] = add_vertex(g);
    }
    Traits::vertex_descriptor s = vertex_descriptors[source];
    Traits::vertex_descriptor t = vertex_descriptors[sink];
   
    // Add edges
    std::vector<Traits::edge_descriptor> edge_descriptors(edges.size());
    for (int i = 0; i < edges.size(); i++)
    {
        Traits::vertex_descriptor u = vertex_descriptors[edges[i].first];
        Traits::vertex_descriptor v = vertex_descriptors[edges[i].second];
        edge_descriptors[i] = add_edge(u, v, g).first;
    }
    
    // Add capacities & reverse edges
    for (int i = 0; i < edges.size(); i++)
    {
        //put(edge_capacity, g, edge_descriptors[i], cap[i]);
        capacity[edge_descriptors[i]] = cap[i];
        residual_capacity[edge_descriptors[i]] = cap[i];
        reverse[edge_descriptors[i]] = edge_descriptors[rev[i]];
    }

    graph_traits<Graph>::vertex_iterator u_iter, u_end;
    graph_traits<Graph>::out_edge_iterator ei, e_end;
    if (DEBUG)
    {
        for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
          for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
            if (capacity[*ei] > 0)
              std::cout << "f " << *u_iter << " " << target(*ei, g) << " " 
                        << capacity[*ei] << ", " << residual_capacity[*ei] << std::endl;
    }

    // Run max flow
    auto start = std::chrono::steady_clock::now();

    //flow = push_relabel_max_flow(g, s, t); // NOTE: push-relabel doesn't work well with the residual capacities; cannot run BFS on residual network afterwards
    //flow = edmonds_karp_max_flow(g, s, t); 
    flow = boykov_kolmogorov_max_flow(g, s, t);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout<<"Elapsed time is :  "<< std::chrono::duration_cast<std::chrono::microseconds>(diff).count() / 1000000.0<<" s "<<std::endl;

    // Print output
    std::cout << "The total flow: " << flow << std::endl;

    if (DEBUG)
    {
        std::cout << "c flow values:" << std::endl;
        for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
          for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
            if (capacity[*ei] > 0)
              std::cout << "f " << *u_iter << " " << target(*ei, g) << " " 
                        << (capacity[*ei] - residual_capacity[*ei]) << std::endl;
    }

    // Find the critical edges (throats / pores)
    //
	std::vector<default_color_type> color(num_vertices(g));
	std::vector<Traits::edge_descriptor> pred(num_vertices(g));

    boost::queue<Traits::vertex_descriptor> Q;
    breadth_first_search(detail::residual_graph(g, residual_capacity), s, Q,
         make_bfs_visitor(record_edge_predecessors(&pred[0], on_tree_edge())),
         &color[0]);

    critical.clear();
    for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
    {
        if (DEBUG) std::cout<<*u_iter<<": color "<<color[*u_iter]<<", pred "<<pred[*u_iter]<<std::endl;
        for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
        {
            int u = *u_iter;
            int v = target(*ei, g);
            int cf = residual_capacity[*ei];
            if (find(critical.begin(), critical.end(), Edge(u, v)) != critical.end() || find(critical.begin(), critical.end(), Edge(v, u)) != critical.end())
            {
                continue;
            }
            if (doubleVertices && v != u + verts_n)
            {
                continue;
            }
            if (cf == 0 && color[u] != color[v])
            {
                critical.push_back(Edge(u, v));
                if (DEBUG) std::cout<<"             critical edge "<<u<<" "<<v<<std::endl;
            }
        }
    }
    assert(critical.size() == flow); // NOTE: fails w/ push-relabel

    return flow;
}



Network fitCriticalPoreCNs(const Network &net, const std::vector<int> &left, const std::vector<int> &right, const std::vector<int> &cp_cns)
{
    // find current critical pores
    //
    std::vector<Edge> criticalPores;
    int flow = maxFlow(net, left, right, true /*doubleVertices*/, criticalPores);
    std::cout<<"Critical pores (count = "<<flow<<") = [";
    for (auto it : criticalPores)
    {
        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
    }
    std::cout<<"]\n\n";

    // Define target CP coordination #s
    //
    std::vector<int> target(MAX_CN + 1); // target CP coordination # histogram
    for (auto cn : cp_cns)
    {
        target[cn]++;
    }

    // Get all coordination #s
    //
    std::vector<int> cns(net.pores.size());
    for (auto throat : net.throats)
    {
        cns[throat.p1]++;
        cns[throat.p2]++;
    }

    // Find CP coordination #s
    //
    std::vector<int> hist(MAX_CN + 1); // CP coordination # histogram
    std::vector<std::vector<int> > cps(MAX_CN + 1); // the actual CPs
    for (auto it : criticalPores)
    {
        int cp = it.first;
        hist[cns[cp]]++;
        cps[cns[cp]].push_back(cp);
    }

    std::cout<<"CP CN hist:            ";
    for (int cn = 0; cn <= MAX_CN; cn++)
    {
        std::cout<<std::setw(3)<<hist[cn]<<" ";
    }
    std::cout<<"\n";
    std::cout<<"target CP CN hist:     ";
    for (int cn = 0; cn <= MAX_CN; cn++)
    {
        std::cout<<std::setw(3)<<target[cn]<<" ";
    }
    std::cout<<"\n";


    // Modify target histogram so it's a subset of the actual histogram
    //
    srand(time(NULL));
    for (int cn = 0; cn <= MAX_CN; cn++)
    {
        if (DEBUG) std::cout<<cn<<": "<<hist[cn]<<" vs. "<<target[cn]<<"\n";
        // As long as there are more CPs with the given CN in the target distribution,
        // keep reassigning them to different CNs randomly (closer CNs are preferred)
        //
        while (hist[cn] < target[cn])
        {
            bool done = false;
            // Must make sure to reassign to a new CN that is feasible
            //
            do
            {
                std::vector<double> cn_weights(MAX_CN + 1);
                int sum = 0;
                if (DEBUG) std::cout<<" weights for "<<cn<<":\n";
                for (int i = 0; i <= MAX_CN; i++)
                {
                    cn_weights[i] = (i == cn) ? 0 : 1.0 / ((i - cn) * (i - cn));
                    sum += cn_weights[i];
                    if (DEBUG) std::cout<<"      "<<i<<" -> "<<cn_weights[i]<<"\n";
                }

                double r = sum * ((double) rand() / (RAND_MAX));
                int new_cn = -1;
                for (int i = 0; i <= MAX_CN; i++)
                {
                    if (r < cn_weights[i])
                    {
                        new_cn = i;
                        break;
                    }
                    r -= cn_weights[i];
                }
                assert(new_cn != -1);
                if (DEBUG) std::cout<<" new cn = "<<new_cn<<"\n";
                if (hist[new_cn] >= target[new_cn] + 1)
                {
                    target[cn]--;
                    target[new_cn]++;
                    done = true;
                    if (DEBUG) std::cout<<"                WORKS!\n";
                }
            } while (!done);
        }
    }

    std::cout<<"new target CP CN hist: ";
    for (int cn = 0; cn <= MAX_CN; cn++)
    {
        std::cout<<std::setw(3)<<target[cn]<<" ";
    }
    std::cout<<"\n";

    // Mark extra CPs for removal (randomly)
    //
    std::vector<bool> removed(net.pores.size());
    for (int cn = 0; cn <= MAX_CN; cn++)
    {
        int to_remove = hist[cn] - target[cn];
        assert(to_remove >= 0);

        std::vector<int> temp = cps[cn];
        assert(temp.size() == hist[cn]);
        std::random_shuffle(temp.begin(), temp.end());

        for (int i = 0; i < to_remove; i++)
        {
            int p = temp[i];
            removed[p] = true;
            if (DEBUG) std::cout<<" remove "<<p<<" (cn = "<<cns[p]<<")\n";
        }
    }

    // Create new network with the removed CPs disconnected
    //
    std::vector<Throat> new_throats;
    for (auto throat : net.throats)
    {
        if (!removed[throat.p1] && !removed[throat.p2])
        {
            new_throats.push_back(throat);
        }
    }

    Network new_net = net;
    new_net.throats = new_throats;
    std::cout<<"new net has "<<new_net.throats.size()<<" throats\n";

    // Find critical pores in new net
    //
    std::vector<Edge> newCriticalPores;
    flow = maxFlow(new_net, left, right, true /*doubleVertices*/, newCriticalPores);
    std::cout<<"New critical pores (count = "<<flow<<") = [";
    for (auto it : newCriticalPores)
    {
        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
    }
    std::cout<<"]\n\n";

    // Find new histogram
    std::vector<int> new_hist(MAX_CN + 1); // CP coordination # histogram
    for (auto it : newCriticalPores)
    {
        int cp = it.first;
        new_hist[cns[cp]]++;
    }

    std::cout<<"new CP CN hist:        ";
    for (int cn = 0; cn <= MAX_CN; cn++)
    {
        std::cout<<std::setw(3)<<new_hist[cn]<<" ";
    }
    std::cout<<"\n";

    return new_net;
}




int main(int argc, char* argv[])
{
    using namespace boost;

    int direction = 0;
    const char *infile = argv[1];

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

    int flow;
    // Critical throats
    //
//    std::vector<Edge> criticalThroats;
//    flow = maxFlow(net, left, right, false /*doubleVertices*/, criticalThroats);
//    std::cout<<"Critical throats (count = "<<ct<<") = [";
//    for (auto it : criticalThroats)
//    {
//        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
//    }
//    std::cout<<"]\n\n";
//
//    // Critical pores
//    //
//    std::vector<Edge> criticalPores;
//    flow = maxFlow(net, left, right, true /*doubleVertices*/, criticalPores);
//    std::cout<<"Critical pores (count = "<<flow<<") = [";
//    for (auto it : criticalPores)
//    {
//        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
//    }
//    std::cout<<"]\n\n";

    // Remove critical pores until reaching a given target histogram
    //
    Network new_net = fitCriticalPoreCNs(net, left, right, {2, 3, 4, 5, 9, 9, 10, 12, 14, 14});

    return 0;
}
