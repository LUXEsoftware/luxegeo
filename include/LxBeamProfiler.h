//
/// \brief LxBeamProfiler class (DD4hep version)
//

#ifndef LxBeamProfiler_h
#define LxBeamProfiler_h 1

#include <string>
#include <vector>
#include "DD4hep/DetFactoryHelper.h"

class LxBeamProfiler
{
public:
  LxBeamProfiler() {};
  virtual ~LxBeamProfiler() {};

  void Construct(dd4hep::Detector&          description,
                 dd4hep::DetElement&        sdet,
                 xml_h&                     e,
                 dd4hep::SensitiveDetector& sd);

protected:
  void CreateMaterial(dd4hep::Detector& description);

  // Returns the PCB+sensor assembly and the sensor PlacedVolume (for DetElement cloning)
  dd4hep::Assembly ConstructPCBAssembly(dd4hep::Detector&          description,
                                        dd4hep::SensitiveDetector& sd,
                                        dd4hep::PlacedVolume&      pvSensor);

  dd4hep::Assembly ConstructPCBSupportAssembly(dd4hep::Detector& description);
  dd4hep::Assembly ConstructMotorsAssembly(dd4hep::Detector& description);
};

#endif
