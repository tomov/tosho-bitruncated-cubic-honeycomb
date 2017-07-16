#include <boost/config.hpp>
#include <iostream>
#include <cstdio>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <chrono>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>

const bool DEBUG = false;

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
    std::vector<int> edges; // corresponding edges

    Throat() {}
    Throat(int _p1, int _p2, double _r)
        : p1(_p1), p2(_p2), r(_r) {}
};

std::vector<Pore> pores;
std::vector<Throat> throats;
std::vector<int> left, right;
double ucs;
int N;
double permeability;

int readCSV(const char* infile)
{
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
            pores.push_back(pore);
            if (DEBUG) printf("%lf %lf %lf: %lf,", pore.x, pore.y, pore.z, pore.r);
        }
        if (!streq(line[3], "") && !streq(line[3], "0"))
        {
            Throat throat(atoi(line[3]) - 1, atoi(line[4]) - 1, atof(line[5]) * 1e-6);
            throats.push_back(throat);
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
            ucs = atof(line[9]);
            if (DEBUG) printf("%lf,", ucs);
        }
        if (!streq(line[10], "") && !streq(line[10], "0"))
        {
            N = atoi(line[10]);
            if (DEBUG) printf("%d,", N);
        }
        if (!streq(line[11], "") && !streq(line[11], "0"))
        {
            permeability = atof(line[11]);
            if (DEBUG) printf("%lf,", permeability);
        }

        if (DEBUG) printf("\n");
    }
    fclose(f);

    return 0;
}


void getBoundaries(const std::vector<Pore> &pores, int direction, std::vector<int> &starting, std::vector<int> &ending)
{
    starting.clear();
    ending.clear();
    std::vector<double> left_b(pores.size()), right_b(pores.size());
    double min_b = pores[0][direction], max_b = pores[0][direction];
    for (int i = 0; i < pores.size(); i++)
    {
        min_b = std::min(min_b, pores[i][direction]);
        max_b = std::max(max_b, pores[i][direction]);
        left_b[i] = pores[i][direction] - pores[i].r;
        right_b[i] = pores[i][direction] + pores[i].r;
    }
    for (int i = 0; i < pores.size(); i++)
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

int maxFlow(std::vector<int> left, std::vector<int> right, bool doubleVertices, std::vector<Edge> &critical)
{

    // Transform into graph and double the edges
    //
    std::vector<Pore> verts = pores;
    std::vector<Edge> edges;
    std::vector<Throat*> edgeThroats; // pointers to the corresponding throats
    for (int i = 0; i < throats.size(); i++)
    {
        edges.push_back(Edge(throats[i].p1, throats[i].p2));
        edges.push_back(Edge(throats[i].p2, throats[i].p1));

        edgeThroats.push_back(&throats[i]);
        edgeThroats.push_back(&throats[i]);

        throats[i].edges.clear();
        throats[i].edges.push_back(edges.size() - 2);
        throats[i].edges.push_back(edges.size() - 1);
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
        edgeThroats.push_back(nullptr);
    }
    for (int i = 0; i < right.size(); i++)
    {
        edges.push_back(Edge(right[i], sink));
        edgeThroats.push_back(nullptr);
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

    // Double verts TODO
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





int main(int argc, char* argv[])
{
    using namespace boost;

    int direction = 0;
    const char *infile = argv[1];

    readCSV(infile);
    printf("solving %s\n", infile);
    printf("N = %d\n", N);
    printf("UCS = %lf\n", ucs);
    printf("perm = %e\n", permeability);
    printf("# pores = %d\n", (int)pores.size());
    printf("# throats = %d\n", (int)throats.size());
    printf("# pores on left boundary = %d\n", (int)left.size());
    printf("# pores on right boundary = %d\n", (int)right.size());

    // Get boundary pores
    //
    if (left.size() == 0 || right.size() == 0)
    {
        getBoundaries(pores, direction, left, right);
        printf("Computing boundaries: # left = %d, # right = %d\n", (int)left.size(), (int)right.size());
    }

    // Critical throats
    //
    std::vector<Edge> criticalThroats;
    int ct = maxFlow(left, right, false /*doubleVertices*/, criticalThroats);
    std::cout<<"Critical throats (count = "<<ct<<") = [";
    for (auto it : criticalThroats)
    {
        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
    }
    std::cout<<"]\n\n";

    // Critical pores
    //
    std::vector<Edge> criticalPores;
    int cp = maxFlow(left, right, true /*doubleVertices*/, criticalPores);
    std::cout<<"Critical pores (count = "<<cp<<") = [";
    for (auto it : criticalPores)
    {
        std::cout<<"("<<it.first<<", "<<it.second<<"), ";
    }
    std::cout<<"]\n\n";

    return 0;
}
