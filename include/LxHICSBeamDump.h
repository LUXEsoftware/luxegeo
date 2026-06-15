//
/// \brief LxHICSBeamDump class (DD4hep version)
//

#ifndef LxHICSBeamDump_dd4hep_h
#define LxHICSBeamDump_dd4hep_h 1

#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxHICSBeamDump
{
public:
  LxHICSBeamDump() : fHICSDumpAngle(0.0) {};
  virtual ~LxHICSBeamDump() {};

  void Construct(dd4hep::Detector&   description,
                 dd4hep::DetElement& sdet,
                 xml_h&              e);

protected:
  void ConstructElectronShielding(dd4hep::Detector& description,
                                  dd4hep::Volume&   motherVol);

  void ConstructNeutronAbsorber(dd4hep::Detector& description,
                                dd4hep::Volume&   motherVol);

  dd4hep::Assembly ConstructSupportAssembly(dd4hep::Detector& description,
                                            double&           supporthight);

  dd4hep::Assembly ConstructHICSElDetSupportAssembly(dd4hep::Detector& description,
                                                     dd4hep::Volume&   motherVol,
                                                     double&           sphight);

private:
  double fHICSDumpAngle;
};

#endif
