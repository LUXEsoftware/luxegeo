//
/// \brief LxAux class (DD4hep version)
//

#ifndef LxAux_dd4hep_h
#define LxAux_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxAux
{
public:
  LxAux() {};
  ~LxAux() {};

  // Build a table assembly: top plate + legs.
  // ytop < 0 falls back to OPPPBasePlateY constant; rleg < 0 falls back to OPPPBasePlateR constant.
  static dd4hep::Assembly BuildTable(dd4hep::Detector&   description,
                                     const std::string&  tname,
                                     double xs, double ys, double zs,
                                     int nleg,
                                     double ytop = -1.0,
                                     double rleg = -1.0);

  // Build a rectangular concrete pedestal volume (not an assembly — placed directly by caller).
  static dd4hep::Volume BuildPedestal(dd4hep::Detector&  description,
                                      const std::string& pname,
                                      double xs, double ys, double zs);

  // Merge all placed volumes from srcAssembly (offset by translation+rotation)
  // into dstAssembly. Mirrors LxAux::AddAssmblyVolumes behaviour.
  static void AddAssemblyVolumes(dd4hep::Assembly&   dstAssembly,
                                 dd4hep::Assembly&   srcAssembly,
                                 const dd4hep::Position&    translation,
                                 const dd4hep::RotationZYX& rotation = dd4hep::RotationZYX(0,0,0));
};

#endif
