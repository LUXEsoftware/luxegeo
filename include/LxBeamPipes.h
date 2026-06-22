//
/// \brief LxBeamPipes class (DD4hep version)
//

#ifndef LxBeamPipes_dd4hep_h
#define LxBeamPipes_dd4hep_h 1

#include <map>
#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxBeamPipes
{
public:
  LxBeamPipes(dd4hep::Detector& description, dd4hep::Volume& motherVol);
  virtual ~LxBeamPipes() {};

  // Build and place the named beam pipe section into the mother volume.
  // Checks TGeoManager first — skips construction if already placed.
  void Construct(const std::string& pipeName);

private:
  void ConstructBeamPipeToDump();
  void ConstructBeamPipeToIP();
  void ConstructGammaVacuumChamber();
  void ConstructBeamPipeTM();
  void ConstructBeamPipeInc();

  typedef void (LxBeamPipes::*ProcessFT)();
  std::map<std::string, ProcessFT> fFunctionMap;

  dd4hep::Detector& fDescription;
  dd4hep::Volume&   fMotherVol;
};

#endif
