#ifndef COMPONENT_FINDER_H
#define COMPONENT_FINDER_H

#include <array>
#include <vector>
#include <map>

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

// unvisited_faces is a copy of all faces where the key specified is the index found in all_faces
// and visited faces will be ignored
// component should initially be empty and will add faces with respective indices
void closed_face_list(unsigned int f0, std::map<uint, Face> &unvisited_faces, std::map<uint, Face>& component) {
    Face f = unvisited_faces[f0];
    // Search for all faces with each index in face
    for (const auto& i : f.data) {
        std::map<uint, Face> intersecting_faces = query_faces_with_index(i, unvisited_faces);
        // Add these faces to query_list and recall closed_face_list on each face found
        for (const auto& f : intersecting_faces) {
            component.insert(f);
        }
        for (const auto& f : intersecting_faces) {
            closed_face_list(f.first, unvisited_faces, component);
        }
    }
}

std::vector<std::map<uint, Face>> find_all_components(std::map<uint, Face>& unvisited_faces, uint minComponentSize = 0) {

    std::vector<std::map<uint, Face>> components;
    if (unvisited_faces.empty())
        return components;

    do {
        components.push_back(std::map<uint, Face>{});

        // Choose component
        uint f0 = unvisited_faces.begin()->first;

        closed_face_list(f0, unvisited_faces, components[components.size() - 1]);
        // Erase all entries that were visited
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
    auto compare = [](const std::map<uint, Face>& c1, const std::map<uint, Face>& c2) {
        return (c1.size() > c2.size());
    };

    std::sort(components.begin(), components.end(), compare);

    return components;
}

#endif