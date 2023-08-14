#ifndef COMPONENT_FINDER_H
#define COMPONENT_FINDER_H

#include <array>
#include <vector>
#include <map>
#include <deque>

typedef unsigned int uint;

struct Face {
    Face() = default;
    Face(uint a, uint b, uint c) : data({ a,b,c }) {}
    Face(const std::array<uint, 3>& d) : data(d) {}
    Face(const Face& f) : data(f.data) {}
    std::array<uint, 3> data;
    bool visited = false;
};

std::map<uint, Face> query_faces_with_index(unsigned int i, std::map<uint, Face> &faces) {

    std::map<uint, Face> queried_faces;

    for (std::pair<const uint, Face>& f : faces) {

        // If face has already been visited, ignore
        if (f.second.visited)
            continue;

        for (const auto &fi : f.second.data) {
            // Found the index
            if (fi == i) {
                f.second.visited = true;
                queried_faces.insert(f);
                break;
            }
        }
    }

    return queried_faces;

}

std::map<uint, Face> query_faces_intersection(const Face &f0, std::map<uint, Face>& faces) {

    std::map<uint, Face> queried_faces;

    for (std::pair<const uint, Face>& f : faces) {

        // If face has already been visited, ignore
        if (f.second.visited)
            continue;

        bool found = false;

        for (const auto& fi : f.second.data) {
            for (const auto& fj : f0.data) {
                // Found an intersection
                if (fi == fj) {
                    f.second.visited = true;
                    queried_faces.insert(f);
                    found = true;
                    break;
                }
            }

            if (found)
                break;
        }
    }

    return queried_faces;

}

// unvisited_faces is a copy of all faces where the key specified is the index found in all_faces
// and visited faces will be ignored
// component should initially be empty and will add faces with respective indices
void closed_face_list(unsigned int f0, std::map<uint, Face> &unvisited_faces, std::map<uint, Face>& component) {
    Face f = unvisited_faces[f0];
    // Search for all faces with each index in face

    // Search for faces that intersect face f
    std::map<uint, Face> intersecting_faces = query_faces_intersection(f, unvisited_faces);

    // Add these faces to the component
    for (const auto& face : intersecting_faces) {
        component.insert(face);
    }

    // Recall closed_face_list on each face found
    for (const auto& face : intersecting_faces) {
        closed_face_list(face.first, unvisited_faces, component);
    }
    
}

std::vector<std::map<uint, Face>> find_all_components(std::map<uint, Face>& unvisited_faces, uint minComponentSize = 0) {

    std::vector<std::map<uint, Face>> components;
    if (unvisited_faces.empty())
        return components;

    do {
        components.push_back(std::map<uint, Face>{});

        // Choose any face to start, first for convenience
        uint f0 = unvisited_faces.begin()->first;

        closed_face_list(f0, unvisited_faces, components[components.size() - 1]);

        // Find all faces that were visited
        std::vector<uint> visited;
        uint componentSize = components.size();
        uint numVisited = components[componentSize - 1].size();
        visited.reserve(numVisited);
        for (auto& f : components[componentSize - 1])
            visited.push_back(f.first);

        // Now remove entries that have been visited
        for (auto& i : visited) {
            std::map<uint, Face>::iterator iter = unvisited_faces.find(i);
            if (iter != unvisited_faces.end())
                unvisited_faces.erase(iter);
        }

        // Remove the component if it is too small
        if (components[componentSize - 1].size() < minComponentSize)
            components.pop_back();


    } while (unvisited_faces.size());


    // Sort components by largest to smallest
    auto largest = [](const std::map<uint, Face>& c1, const std::map<uint, Face>& c2) {
        return (c1.size() > c2.size());
    };

    std::sort(components.begin(), components.end(), largest);

    return components;
}


class Graph {
    uint V; // No. of vertices

    // Pointer to an array containing adjacency lists
    std::list<uint>* adj;

    // A function used by DFS
    void DFSUtil(uint v, bool visited[], std::vector<uint>& vertices);
    void DFSIterative(uint v, bool visited[], std::vector<uint>& vertices);

public:
    Graph(uint V); // Constructor
    ~Graph();
    void addEdge(uint v, uint w);
    std::list<uint>* adjacencyList() { return adj; }

    std::vector<std::vector<uint>> connectedComponents(uint minComponentSize);
};

// Method to print connected components in an
// undirected graph
std::vector<std::vector<uint>> Graph::connectedComponents(uint minComponentSize)
{
    // Mark all the vertices as not visited
    bool* visited = new bool[V];
    for (uint v = 0; v < V; v++)
        visited[v] = false;

    std::vector<std::vector<uint>> components;

    for (uint v = 0; v < V; v++) {
        if (visited[v] == false) {
            // This gives all components from vertex v
            std::vector<uint> verts;
            DFSIterative(v, visited, verts);
            // Store as a component if bigger than desired size
            if (minComponentSize <= verts.size())
                components.push_back(verts);
        }
    }
    delete[] visited;
    return components;
}

void Graph::DFSUtil(uint v, bool visited[], std::vector<uint> &vertices)
{
    // Mark the current node as visited and print it
    visited[v] = true;
    vertices.push_back(v);

    // Recur for all the vertices
    // adjacent to this vertex not already visited
    std::list<uint>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFSUtil(*i, visited, vertices);
}

void Graph::DFSIterative(uint src, bool visited[], std::vector<uint>& vertices) {

    // src : source node
    std::deque<uint> stack;
    stack.push_back(src);

    while (stack.size()) {
        src = stack.front();
        stack.pop_front();
        if (!visited[src]) {
            visited[src] = true;
            vertices.push_back(src);
            std::list<uint>::iterator i;
            for (i = adj[src].begin(); i != adj[src].end(); ++i)
                if (!visited[*i])
                    stack.push_front(*i);
        }
    }

}

Graph::Graph(uint V)
{
    this->V = V;
    adj = new std::list<uint>[V];
}

Graph::~Graph() { delete[] adj; }

// method to add an undirected edge
void Graph::addEdge(uint v, uint w)
{
    // Check if edge already exists
    for (const auto& e : adj[v])
        if (e == w)
            return;

    adj[v].push_back(w);
    adj[w].push_back(v);
}

class GraphVector {
    uint V; // No. of vertices

    // Pointer to the positions in the graph (in 3D)
    const double* points = 0;

    // threshold to accept point
    double threshold;

    // Pointer to an array containing adjacency lists
    std::list<uint>* adj;

    // A function used by DFS
    void DFSUtil(uint v, bool visited[], std::vector<uint>& vertices);

public:
    GraphVector(uint V, const double * points, double thresh = 0.5); // Constructor
    ~GraphVector();
    void addEdge(uint v, uint w);
    std::list<uint>* adjacencyList() { return adj; }

    std::vector<std::vector<uint>> connectedComponents(uint minComponentSize);
};

// Method to print connected components in an
// undirected graph
std::vector<std::vector<uint>> GraphVector::connectedComponents(uint minComponentSize)
{
    // Mark all the vertices as not visited
    bool* visited = new bool[V];
    for (uint v = 0; v < V; v++)
        visited[v] = false;

    std::vector<std::vector<uint>> components;

    for (uint v = 0; v < V; v++) {
        if (visited[v] == false) {
            // This gives all components from vertex v
            std::vector<uint> verts;
            DFSUtil(v, visited, verts);
            // Store as a component if bigger than desired size
            if (minComponentSize <= verts.size())
                components.push_back(verts);
        }
    }
    delete[] visited;
    return components;
}

void GraphVector::DFSUtil(uint v, bool visited[], std::vector<uint>& vertices)
{
    // Mark the current node as visited and print it
    visited[v] = true;
    vertices.push_back(v);

    // Recur for all the vertices
    // adjacent to this vertex
    std::list<uint>::iterator i;
    uint vp = vertices.size() - 2;
    double tx0 = points[3 * v] - points[3 * vp];
    double ty0 = points[3 * v + 1] - points[3 * vp + 1];
    double tz0 = points[3 * v + 2] - points[3 * vp + 2];
    for (i = adj[v].begin(); i != adj[v].end(); ++i) {

        uint n = *i;

        if (!visited[n]) {
            // Compute the current tangent and the previous tangent
            if (vertices.size() > 2) {

                double tx = points[3 * n] - points[3 * v];
                double ty = points[3 * n + 1] - points[3 * v + 1];
                double tz = points[3 * n + 2] - points[3 * v + 2];

                // Compute the dot product
                double dprod = tx * tx0 + ty * ty0 + tz * tz0;


                if (dprod < threshold)
                    continue;
            }

            DFSUtil(n, visited, vertices);

        }
    }
}


GraphVector::GraphVector(uint V, const double* pts, double thresh)
{
    this->V = V;
    this->points = pts;
    adj = new std::list<uint>[V];
}

GraphVector::~GraphVector() { delete[] adj; }

// method to add an undirected edge
void GraphVector::addEdge(uint v, uint w)
{
    // Check if edge already exists
    for (const auto& e : adj[v])
        if (e == w)
            return;

    adj[v].push_back(w);
    adj[w].push_back(v);
}


std::vector<std::vector<uint>> find_all_components_graph(std::map<uint, Face>& unvisited_faces, uint nVerts, uint minComponentSize = 0) {

    std::vector<std::vector<uint>> components;
    if (unvisited_faces.empty())
        return components;


    Graph G(nVerts);

    for (auto& f : unvisited_faces) {
        // Each face contributes 3 edges for a triangle
        for (int d = 0; d < 3; d++) {
            uint dp1 = (d + 1) % 3;
            G.addEdge(f.second.data[d], f.second.data[dp1]);
        }
    }

    components = G.connectedComponents(minComponentSize);

    // Sort components by largest to smallest
    auto largest = [](const std::vector<uint>& c1, const std::vector<uint>& c2) {
        return (c1.size() > c2.size());
    };

    std::sort(components.begin(), components.end(), largest);

    // For now just return a single component
    return components;
}
#endif