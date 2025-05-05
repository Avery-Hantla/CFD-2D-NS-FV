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

#include <iostream>
#include <filesystem>
#include <cassert>
#include <vector>
#include <array>
#include <unordered_map>

#include <cgnslib.h>

#include "../include/types.hpp"

template<typename IndexType>
struct MESHREADERLIB_EXPORT InterfaceConnectivity {
  std::array<int, 2> zones;
  std::vector<IndexType> ownerIndex;
  std::vector<IndexType> neighbourIndex;
};

template<typename IndexType>
struct MESHREADERLIB_EXPORT BoundaryConditionInformation {
  std::string boundaryName;
  unsigned boundaryType;
  std::vector<IndexType> indices;
};

class MESHREADERLIB_EXPORT ReadMeshBase {
public:
  using ZoneType = typename std::vector<std::vector<cgsize_t>>;

public:
  ReadMeshBase(std::filesystem::path cgnsFilePath);
  ~ReadMeshBase();

// interface / API
public:
  virtual void readMesh() = 0;

// internal implementation hidden from the user
protected:
  void readBases();
  void readZones(ZoneType_t zoneType);
  void createBoundaryNameLookupTable();

// internal implementations that must be implemented by derived classes
protected:
  virtual void readCoorinates() = 0;
  virtual void readInterfaceConnectivity() = 0;
  virtual void readBoundaries() = 0;
  unsigned getBoundaryConditionID(BCType_t bcType);

protected:
  std::filesystem::path _cgnsFilePath;
  int _numberOfZones = 0;
  int _numberOfBoundaries = 0;
  int _fileIndex = 0;

  ZoneType _zoneSizesMax, _zoneSizesMin;
  std::unordered_map<std::string, cgsize_t> _zoneNamesMap;
  std::unordered_map<std::string, BCType_t> _boundaryConditionMap;
};