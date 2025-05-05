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
#include <array>
#include <unordered_map>
#include <string>
#include <cmath>

#include <cgnslib.h>

#include "../include/readMeshBase.hpp"
#include "../include/types.hpp"

class MESHREADERLIB_EXPORT ReadStructuredMesh : public ReadMeshBase {
public:
  using ZoneType = typename std::vector<std::vector<cgsize_t>>;
  using CoordinateType = typename std::vector<std::vector<std::vector<std::vector<double>>>>;
  using IndexType = typename std::array<unsigned, 2>;
  using InterfaceConnectivityType = typename std::vector<InterfaceConnectivity<IndexType>>;
  using BoundaryConditionInformationType = typename std::vector<std::vector<BoundaryConditionInformation<IndexType>>>;

public:
  ReadStructuredMesh(std::filesystem::path cgnsFilePath);
  void readMesh() override final;

  const CoordinateType& getCoordinates() const { return _coordinates; }
  const InterfaceConnectivityType& getInterfaceConnectivity() const { return _interfaceConnectivity; }
  const BoundaryConditionInformationType& getBoundaryConditions() const { return _boundaryConditions; }

protected:
  void readCoorinates() override final;
  void readInterfaceConnectivity() override final;
  void readBoundaries() override final;

private:
  CoordinateType _coordinates;
  InterfaceConnectivityType _interfaceConnectivity;
  BoundaryConditionInformationType _boundaryConditions;
};