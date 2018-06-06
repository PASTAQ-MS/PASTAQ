#ifndef TAPP_GRID_WEIGHTEDMESH_HPP
#define TAPP_GRID_WEIGHTEDMESH_HPP

#include <vector>

#include "tapp-grid.hpp"
#include "tapp-mesh.hpp"

class WeightedMesh : public Mesh {
   protected:
    std::vector<double> mWeights;
    std::vector<unsigned int> mCounts;

   public:
    WeightedMesh(Grid::Dimensions dimensions = {}, Grid::Bounds bounds = {});

    // This method will set the value stored at the given index to be proportional
    // to the existing value on that position and the given weight, incrementing
    // in turn the counter on the same position.
    bool set(unsigned int i, unsigned int j, double value, double weight);

    // Prints the mesh to std::out for debugging purposes.
    void printAll();
};

#endif /* TAPP_GRID_WEIGHTEDMESH_HPP */
