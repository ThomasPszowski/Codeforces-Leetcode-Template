#pragma once
#include "declarations.h"

struct Edge {
    int to;
    int weight;
};

using wGraph = vector<vector<Edge>>; // Weighted graph
using uGraph = vector<vector<int>>;  // Unweighted graph


// Overload the << operator for Edge
ostream& operator<<(ostream& os, const Edge& edge) {
    os << edge.to << " (" << edge.weight << ")";
    return os;
}

// Overload the << operator for wGraph (Weighted Graph)
ostream& operator<<(ostream& os, const wGraph& graph) {
    os << endl;
    for (int i = 0; i < graph.size(); ++i) {
        os << "Node " << i << ": ";
        for (const auto& edge : graph[i]) {
            os << edge << " "; // Print each edge
        }
        os << endl;
    }
    return os;
}

// Overload the << operator for uGraph (Unweighted Graph)
ostream& operator<<(ostream& os, const uGraph& graph) {
    os << endl;
    for (int i = 0; i < graph.size(); ++i) {
        os << "Node " << i << ": ";
        for (int neighbor : graph[i]) {
            os << neighbor << " "; // Print each neighbor
        }
        os << endl;
    }
    return os;
}

// Function to convert adjacency matrix to wGraph
wGraph adjMatrixToWGraph(const vector<vector<int>>& adjMatrix) {
    int n = adjMatrix.size();
    wGraph graph(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (adjMatrix[i][j] != 0) {  // Edge exists if weight is non-zero
                graph[i].push_back({ j, adjMatrix[i][j] });
            }
        }
    }

    return graph;
}

// Function to convert wGraph to adjacency matrix
vector<vector<int>> wGraphToAdjMatrix(const wGraph& graph) {
    int n = graph.size();
    vector<vector<int>> adjMatrix(n, vector<int>(n, 0));  // Initialize matrix with zeros

    for (int i = 0; i < n; ++i) {
        for (const auto& edge : graph[i]) {
            adjMatrix[i][edge.to] = edge.weight;  // Set the weight in the matrix
        }
    }

    return adjMatrix;
}

// Convert an adjacency matrix to an unweighted graph (uGraph)
uGraph adjMatrixToUGraph(const vector<vector<int>>& adj_matrix) {
    int n = adj_matrix.size();
    uGraph graph(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (adj_matrix[i][j]) {
                graph[i].push_back(j);
            }
        }
    }

    return graph;
}

// Convert an unweighted graph (uGraph) to an adjacency matrix
vector<vector<int>> uGraphToAdjMatrix(const uGraph& graph) {
    int n = graph.size();
    vector<vector<int>> adj_matrix(n, vector<int>(n, 0));

    for (int i = 0; i < n; ++i) {
        for (int neighbor : graph[i]) {
            adj_matrix[i][neighbor] = 1;  // Set to 1 for an edge
        }
    }

    return adj_matrix;
}

// Dijkstra's Algorithm for shortest paths in weighted graphs
vector<int> dijkstra(const wGraph& graph, int start) {
    int n = graph.size();
    vector<int> distances(n, numeric_limits<int>::max());
    distances[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({ 0, start });

    while (!pq.empty()) {
        int dist = pq.top().first;
        int node = pq.top().second;
        pq.pop();

        if (dist > distances[node]) continue;

        for (const auto& edge : graph[node]) {
            int next_node = edge.to;
            int next_dist = dist + edge.weight;
            if (next_dist < distances[next_node]) {
                distances[next_node] = next_dist;
                pq.push({ next_dist, next_node });
            }
        }
    }
    return distances;
}

// DFS for unweighted graphs, marking visited nodes
void dfs(int node, const uGraph& graph, vector<bool>& visited) {
    visited[node] = true;
    for (int neighbor : graph[node]) {
        if (!visited[neighbor]) {
            dfs(neighbor, graph, visited);
        }
    }
}

// BFS for unweighted graphs to find shortest path in terms of edges
vector<int> bfs(int start, const uGraph& graph) {
    int n = graph.size();
    vector<int> distance(n, -1);  // -1 for unvisited
    queue<int> q;
    q.push(start);
    distance[start] = 0;

    while (!q.empty()) {
        int node = q.front();
        q.pop();

        for (int neighbor : graph[node]) {
            if (distance[neighbor] == -1) {
                distance[neighbor] = distance[node] + 1;
                q.push(neighbor);
            }
        }
    }
    return distance;
}

// Connected Components Count for undirected unweighted graphs
int find_connected_components(const uGraph& graph) {
    int n = graph.size();
    vector<bool> visited(n, false);
    int components = 0;

    for (int node = 0; node < n; ++node) {
        if (!visited[node]) {
            dfs(node, graph, visited);
            components++;
        }
    }
    return components;
}

// Cycle Detection for undirected unweighted graphs
bool has_cycle_undirected(int node, int parent, const uGraph& graph, vector<bool>& visited) {
    visited[node] = true;
    for (int neighbor : graph[node]) {
        if (!visited[neighbor]) {
            if (has_cycle_undirected(neighbor, node, graph, visited)) return true;
        }
        else if (neighbor != parent) {
            return true;
        }
    }
    return false;
}

bool contains_cycle_undirected(const uGraph& graph) {
    int n = graph.size();
    vector<bool> visited(n, false);
    for (int node = 0; node < n; ++node) {
        if (!visited[node]) {
            if (has_cycle_undirected(node, -1, graph, visited)) return true;
        }
    }
    return false;
}

// Cycle Detection in Directed Graphs (unweighted) with color-marking
bool has_cycle_directed(int node, const uGraph& graph, vector<int>& visited) {
    visited[node] = 1;
    for (int neighbor : graph[node]) {
        if (visited[neighbor] == 1) return true;
        if (visited[neighbor] == 0 && has_cycle_directed(neighbor, graph, visited)) return true;
    }
    visited[node] = 2;
    return false;
}

bool contains_cycle_directed(const uGraph& graph) {
    int n = graph.size();
    vector<int> visited(n, 0);
    for (int node = 0; node < n; ++node) {
        if (visited[node] == 0) {
            if (has_cycle_directed(node, graph, visited)) return true;
        }
    }
    return false;
}

// Bipartiteness Check for undirected unweighted graphs
// Returns true if the graph is bipartite, false otherwise
// A bipartite graph is a graph whose vertices can be divided into two disjoint sets
// U and V such that every edge connects a vertex in U to one in V
bool is_bipartite(const uGraph& graph) {
    int n = graph.size();
    vector<int> color(n, -1);
    queue<int> q;

    for (int start = 0; start < n; ++start) {
        if (color[start] == -1) {
            q.push(start);
            color[start] = 0;

            while (!q.empty()) {
                int node = q.front();
                q.pop();

                for (int neighbor : graph[node]) {
                    if (color[neighbor] == -1) {
                        color[neighbor] = 1 - color[node];
                        q.push(neighbor);
                    }
                    else if (color[neighbor] == color[node]) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// Prim's Algorithm for MST in weighted graphs
int prims_mst_cost(const wGraph& graph) {
    int n = graph.size();
    vector<int> min_edge(n, INT_MAX);
    vector<bool> in_mst(n, false);
    min_edge[0] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({ 0, 0 });
    int mst_cost = 0;

    while (!pq.empty()) {
        int weight = pq.top().first;
        int node = pq.top().second;
        pq.pop();

        if (in_mst[node]) continue;
        in_mst[node] = true;
        mst_cost += weight;

        for (const auto& edge : graph[node]) {
            int neighbor = edge.to;
            int edge_weight = edge.weight;
            if (!in_mst[neighbor] && edge_weight < min_edge[neighbor]) {
                min_edge[neighbor] = edge_weight;
                pq.push({ edge_weight, neighbor });
            }
        }
    }
    return mst_cost;
}


// Function to return both MST cost and the MST edges
pair<int, vector<pair<int, int>>> prims_mst(const wGraph& graph) {
    int n = graph.size();
    vector<int> min_edge(n, INT_MAX);   // To track the minimum edge weight to a node
    vector<bool> in_mst(n, false);      // To track nodes that are included in the MST
    vector<int> parent(n, -1);          // To track the parent node of each node in the MST
    min_edge[0] = 0;                    // Start from node 0
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({ 0, 0 });                    // {weight, node}

    int mst_cost = 0;
    vector<pair<int, int>> mst_edges;   // To store the edges of the MST

    while (!pq.empty()) {
        int weight = pq.top().first;
        int node = pq.top().second;
        pq.pop();

        if (in_mst[node]) continue;      // Skip if node is already in MST

        in_mst[node] = true;
        mst_cost += weight;

        // If the parent is not -1, the edge (parent[node], node) is part of the MST
        if (parent[node] != -1) {
            mst_edges.push_back({ parent[node], node });
        }

        for (const auto& edge : graph[node]) {
            int neighbor = edge.to;
            int edge_weight = edge.weight;

            // If the neighbor is not in MST and we find a smaller weight edge
            if (!in_mst[neighbor] && edge_weight < min_edge[neighbor]) {
                min_edge[neighbor] = edge_weight;
                parent[neighbor] = node;  // Set the parent of the neighbor
                pq.push({ edge_weight, neighbor });
            }
        }
    }

    return { mst_cost, mst_edges };  // Return the total cost and the MST edges
}

// Floyd-Warshall Algorithm for shortest paths in weighted graphs
// Returns a 2D vector of shortest distances between all pairs of nodes
vector<vector<int>> floyd_warshall(const wGraph& graph) {
    int n = graph.size();
    vector<vector<int>> dist(n, vector<int>(n, INT_MAX));

    for (int i = 0; i < n; ++i) {
        dist[i][i] = 0;
        for (const auto& edge : graph[i]) {
            dist[i][edge.to] = edge.weight;
        }
    }

    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist[i][k] != INT_MAX && dist[k][j] != INT_MAX) {
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }
    return dist;
}


// Topological Sort for directed unweighted graphs (DAG)
// Returns an empty vector if a cycle is detected
// Time complexity: O(V + E)
// Returns a topological ordering of the nodes
vector<int> topological_sort(const uGraph& graph) {
    int n = graph.size();
    vector<int> in_degree(n, 0);
    vector<int> result;
    queue<int> q;

    for (const auto& edges : graph) {
        for (int neighbor : edges) {
            in_degree[neighbor]++;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (in_degree[i] == 0) q.push(i);
    }

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        result.push_back(node);

        for (int neighbor : graph[node]) {
            in_degree[neighbor]--;
            if (in_degree[neighbor] == 0) q.push(neighbor);
        }
    }

    if (result.size() != n) result.clear(); // Cycle detected
    return result;
}
