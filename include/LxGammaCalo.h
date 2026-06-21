//
/// \brief LxDetectorGammaCalo8 class (DD4hep version)
//

#ifndef LxDetectorGammaCalo8_dd4hep_h
#define LxDetectorGammaCalo8_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxDetectorGammaCalo8
{
public:
  LxDetectorGammaCalo8() {};
  virtual ~LxDetectorGammaCalo8() {};

  void Construct(dd4hep::Detector&   description,
                 dd4hep::DetElement& sdet,
                 xml_h&              e,
                 dd4hep::SensitiveDetector& sd);

protected:
  void ConstructGammaBeamDump(dd4hep::Detector& description,
                              dd4hep::Volume&   motherVol);

  void ConstructGammaBeamDumpMagnetized(dd4hep::Detector& description,
                                        dd4hep::Volume&   motherVol);

  void ConstructGammaBeamDumpShielding(dd4hep::Detector& description,
                                       dd4hep::Volume&   motherVol);

  void ConstructGammaBeamDumpShieldingTight(dd4hep::Detector& description,
                                            dd4hep::Volume&   motherVol);
};

#endif
