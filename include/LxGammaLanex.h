//
/// \brief LxDetectorComptonFluka class (DD4hep version)
//

#ifndef LxGammaLanex_h
#define LxGammaLanex_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxDetectorComptonFluka
{
public:
  LxDetectorComptonFluka() {};
  virtual ~LxDetectorComptonFluka() {};

  void Construct(dd4hep::Detector&          description,
                 dd4hep::DetElement&        sdet,
                 xml_h&                     e,
                 dd4hep::SensitiveDetector& sd);

protected:
  dd4hep::Assembly ConstructSupportAssembly(dd4hep::Detector& description,
                                            double&           sphight);

  void ConstructComptShielding(dd4hep::Detector& description,
                               dd4hep::Volume&   motherVol);
};

#endif
