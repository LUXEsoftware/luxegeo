//
/// \brief LxECal class (DD4hep version)
//

#ifndef LxECal_dd4hep_h
#define LxECal_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxECal
{
public:
  LxECal() : fECalLayerZ(0.) {};
  virtual ~LxECal() {};

  void Construct(dd4hep::Detector&          description,
                 dd4hep::DetElement&        sdet,
                 xml_h&                     e,
                 dd4hep::SensitiveDetector& sd);

protected:
  dd4hep::Assembly  ConstructSupportAssembly (dd4hep::Detector& description);
  dd4hep::Assembly  ConstructCasingAssembly  (dd4hep::Detector& description, double ecalz);
  dd4hep::Volume    ConstructPCB             (dd4hep::Detector& description,
                                              double ecalx, double ecalz);
  void              ConstructShielding       (dd4hep::Detector& description,
                                              dd4hep::Volume&   motherVol);
  void              ConstructDumpShielding   (dd4hep::Detector& description,
                                              dd4hep::Volume&   motherVol);

  double fECalLayerZ;   // half-layer pitch * 2 + extra — set during Construct
};

#endif
