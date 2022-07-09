#ifndef LC_GRAPH_H
#define LC_GRAPH_H

#include <vector>
#include <list>
#include <memory>

/*
	Directed graph data structure
*/

namespace LC { namespace Math {

    // A class to represent a graph object and detect distinct closed cycles (also computing their length)
    class Graph {
    public:
        struct Edge {
            unsigned int src, dest;
            float weight;
        };
        struct Pair {
            Pair() = default;
            Pair(const unsigned int& pid, const float& w) : id(pid), weight(w) {}
            unsigned int id;
            float weight;
        };
        struct Generator {
            Generator() = default;
            Generator(const unsigned int& i, const float& cLen) : id(i), cycleLength(cLen) {}
            unsigned int id;
            float cycleLength;
        };

        // a vector of vectors of Pairs to represent an adjacency list
        std::vector<std::vector<Pair>> adjList;

        // Graph Constructor
        /*
            vector of edges contains all edges between nodes with given weight
            n is the total number of nodes (even if there are no edges for the given node)
        */
        Graph(const std::vector<Edge>& edges, int n = -1) {
            // Empty graph
            if (n == 0)
                return;
            // Auto detect the number of nodes
            else if (n == -1) {
                n = DetectNodeCount(edges);
            }


            // resize the vector to hold `n` elements of type vector<Edge>
            adjList.resize(n);

            // add edges to the directed graph
            for (const auto& edge : edges) {
                Insert(edge);
            }
        }

        Graph(unsigned int n = 0) {
            // Empty graph
            if (n == 0)
                return;

            // resize the vector to hold `n` elements of type vector<Edge>
            adjList.resize(n);
        }

        void Insert(unsigned int src, unsigned int dst, float w) {
            adjList[src].emplace_back(dst, w);
        }

        void Insert(const Edge& e) {
            int src = e.src;
            int dest = e.dest;
            float weight = e.weight;

            // insert at the end
            adjList[src].emplace_back(dest, weight);
        }

        unsigned int DetectNodeCount(const std::vector<Edge>& edges) {
            unsigned int n = 0;
            for (const auto& e : edges) {
                if (n < e.dest) n = e.dest;
                if (n < e.src) n = e.src;
            }
            return n;
        }

        std::vector<Generator> detectDisjointCycles() {
            // Mark all the vertices as not visited and not part of recursion
            // stack
            unsigned int V = adjList.size();
            std::unique_ptr<bool[]> visited(new bool[V]);
            std::unique_ptr<bool[]>  recStack(new bool[V]);
            std::vector<Generator> generators;

            // Reset visited and stack
            for (int i = 0; i < V; i++) {
                visited[i] = false;
                recStack[i] = false;
            }

            // Call the recursive helper function to detect cycle in different
            // DFS trees
            for (int i = 0; i < V; i++) {

                float length = 0.f;
                if (!visited[i] && isCyclicUtil(i, visited.get(), recStack.get(), length)) {
                    generators.emplace_back(i, length);
                }
            }

            return generators;
        }

    private:
        bool isCyclicUtil(unsigned int v, bool visited[], bool* recStack, float& length) {
            if (visited[v] == false) {
                // Mark the current node as visited and part of recursion stack
                visited[v] = true;
                recStack[v] = true;

                // Recur for all the vertices adjacent to this vertex
                std::vector<Pair>::iterator i;
                for (i = adjList[v].begin(); i != adjList[v].end(); ++i) {
                    if (!visited[i->id] && isCyclicUtil(i->id, visited, recStack, length)) {
                        length += i->weight;
                        return true;
                    }
                    else if (recStack[i->id]) {
                        length += i->weight;
                        return true;
                    }

                }

            }
            recStack[v] = false;  // remove the vertex from recursion stack
            return false;
        }

    };
	
}}

#endif