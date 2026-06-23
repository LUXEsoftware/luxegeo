//
/// \brief LxTargetChamber class (DD4hep version)
//
/// Architecture mirrors LxMagnetAssembly:
///   - GetChamberVolume()   — builds the chamber body once, caches in TGeoManager,
///                            returns cached on subsequent calls.
///   - GetTargetVolume()    — always returns a fresh vacuum box (the "field volume"
///                            equivalent) with a unique name for each placement.
///   - ConstructSupport()   — places table + pedestal into a given mother volume.
//

#ifndef LxTargetChamber_dd4hep_h
#define LxTargetChamber_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxTargetChamber
{
public:
  LxTargetChamber() {};
  ~LxTargetChamber() {};

  // Returns shared chamber body — builds once, caches in TGeoManager.
  dd4hep::Volume GetChamberVolume(dd4hep::Detector& description);

  // Returns a fresh target volume (vacuum box) with a unique name per placement.
  // Static so it can be called without an instance.
  static dd4hep::Volume GetTargetVolume(dd4hep::Detector&  description,
                                        const std::string& volname);

  // Places table + pedestal support into motherVol at (TargetChamberXPos, 0, zpos).
  static void ConstructSupport(dd4hep::Detector&  description,
                               dd4hep::Volume&    motherVol,
                               const std::string& tcname,
                               double             zpos);

protected:
  dd4hep::Volume BuildChamberVolume(dd4hep::Detector& description);
  dd4hep::Assembly ConstructTargetSupportAssembly(dd4hep::Detector& description);
};

#endif
