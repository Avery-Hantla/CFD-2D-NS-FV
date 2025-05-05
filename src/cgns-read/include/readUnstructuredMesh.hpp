#pragma once

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  #if defined(COMPILEDLL)
    #define MESHREADERLIB_EXPORT __declspec(dllexport)
  #elif defined(COMPILELIB)
    #define MESHREADERLIB_EXPORT
  #else
    #define MESHREADERLIB_EXPORT __declspec(dllimport)
  #endif
#else
  #define MESHREADERLIB_EXPORT
#endif

#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <tuple>

#include <cgnslib.h>

#include "../include/readMeshBase.hpp"
#include "../include/types.hpp"

class MESHREADERLIB_EXPORT ReadUnstructuredMesh : public ReadMeshBase {
public:
  using CoordinateType = typename std::vector<std::vector<std::vector<double>>>;
  using IndexType = std::vector<cgsize_t>;
  using InterfaceConnectivityType = typename std::vector<std::vector<InterfaceConnectivity<IndexType>>>;
  using BoundaryConditionInformationType = typename std::vector<std::vector<BoundaryConditionInformation<IndexType>>>;
  
  using InternalElementConnectivityType = typename std::vector<std::vector<IndexType>>;
  using InterfaceElementConnectivityType = typename std::vector<std::vector<std::vector<IndexType>>>;
  using BoundaryElementConnectivityType = typename std::vector<std::vector<std::vector<IndexType>>>;
  using StartLocationType = typename std::vector<IndexType>;
  using StartLocationForElementConnectivityType = typename std::tuple<std::vector<std::vector<cgsize_t>>,
    std::vector<std::vector<cgsize_t>>, std::vector<std::vector<cgsize_t>>>;

public:
  ReadUnstructuredMesh(std::filesystem::path cgnsFilePath);
  void readMesh() override final;

  const CoordinateType& getCoordinates() const { return _coordinates; };
  const InternalElementConnectivityType& getInternalCells() const { return _internalElementConnectivity; };
  const InterfaceConnectivityType getInterfaceConnectivity() const { return _interfaceConnectivity; };
  const BoundaryConditionInformationType& getBoundaryConditions() const { return _boundaryConditions; };

protected:
  void readCoorinates() override final;
  void readInterfaceConnectivity() override final;
  void readBoundaries() override final;

private:
  void readCellConnectivity();
  StartLocationForElementConnectivityType getStartLocationForElementConnectivity() const;

private:
  CoordinateType _coordinates;
  InterfaceConnectivityType _interfaceConnectivity;
  BoundaryConditionInformationType _boundaryConditions;

  InternalElementConnectivityType _internalElementConnectivity;
  InterfaceElementConnectivityType _interfaceElementConnectivity;
  BoundaryElementConnectivityType _boundaryElementConnectivity;
};