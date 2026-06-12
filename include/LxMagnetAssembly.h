//
/// \brief LxMagnetAssembly base class (DD4hep version)
//
/// Concrete magnet model classes derive from this and implement:
///   GetAssembly()    — returns the shared assembly, building it once and
///                      caching it in TGeoManager by (modelName + "_assembly")
///   GetFieldVolume() — static; returns a fresh field volume with a unique
///                      name for each placement
///   ConstructSupport() — places support table + pedestal into a mother volume
//

#ifndef LxMagnetAssembly_dd4hep_h
#define LxMagnetAssembly_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxMagnetAssembly
{
public:
  LxMagnetAssembly(const std::string& modelName) : fModelName(modelName) {};
  virtual ~LxMagnetAssembly() {};

  // Returns the shared assembly for this magnet model.
  // Builds it on first call and registers it in TGeoManager;
  // on subsequent calls retrieves it from TGeoManager directly.
  virtual dd4hep::Assembly GetAssembly(dd4hep::Detector& description) = 0;

  // Returns a fresh field volume with the given unique name.
  // Must be static so it can be called without an instance when the
  // assembly is retrieved from TGeoManager on subsequent placements.
  // Concrete classes implement a static version; this non-static
  // virtual is provided for the case where the caller does have an instance.
  virtual dd4hep::Volume GetFieldVolume(dd4hep::Detector&  description,
                                        const std::string& lvname) = 0;

  // Places support structure into motherVol at pos.
  virtual void ConstructSupport(dd4hep::Detector&       description,
                                dd4hep::Volume&         motherVol,
                                const dd4hep::Position& pos,
                                const std::string&      mname,
                                bool                    rotate = false) = 0;

  const std::string& ModelName() const { return fModelName; }

protected:
  std::string fModelName;
};

#endif
