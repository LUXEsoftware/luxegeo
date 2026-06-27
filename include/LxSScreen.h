//
/// \brief LxSScreen class (DD4hep version)
//
/// Architecture mirrors LxCherenkov:
///   - Single plugin DECLARE_DETELEMENT(LxSScreen,...) called once per
///     <detector type="LxSScreen"> XML entry.
///   - First call builds shared volumes (scintArmLogical, camera assembly)
///     and stores the LxSScreen pointer via description.addExtension<LxSScreen>.
///   - Each call places the appropriate mother volume based on the detector
///     name attribute: "LxHICSScint" → scintArmMotherLogical,
///                     "LxBremsScint" → BremScintArmMotherLogical.
///   - A function pointer map dispatches to the correct placement method.
//

#ifndef LxSScreen_dd4hep_h
#define LxSScreen_dd4hep_h 1

#include <functional>
#include <map>
#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxSScreen
{
public:
  LxSScreen() : fBuilt(false) {};
  virtual ~LxSScreen() {};

  // Build all shared volumes — called once on first plugin invocation.
  void Build(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd);

  // Place the detector named by detName into motherVol.
  // Returns a DetElement for this placement.
  dd4hep::DetElement Place(dd4hep::Detector&   description,
                           dd4hep::Volume&     motherVol,
                           const std::string&  detName,
                           int                 detID);
  void BuildSupportAssembly(dd4hep::Detector& description, dd4hep::Volume& motherVol);

private:
  // Placement methods — one per detector type
  dd4hep::DetElement PlaceHICSScint (dd4hep::Detector& description,
                                     dd4hep::Volume&   motherVol, int detID);
  dd4hep::DetElement PlaceBremsScint(dd4hep::Detector& description,
                                     dd4hep::Volume&   motherVol, int detID);

  // Build sub-assemblies
  void BuildScintArm    (dd4hep::Detector& description, dd4hep::SensitiveDetector& sd);
  void BuildCamera      (dd4hep::Detector& description);

  bool fBuilt;

  // Shared volumes built once
  dd4hep::Volume fScintArmMotherVol;
  dd4hep::Volume fBremScintArmMotherVol;
  dd4hep::Volume fScintArmVol;           // shared between both mothers
  dd4hep::Volume fCameraMotherVol;       // ScintCameraMotherBoxLogical

  // Sensor DetElement template — cloned per placement
  dd4hep::DetElement fPhosphorTemplate;

  // Dispatch map: detector name → placement method pointer
  using PlaceFn = std::function<dd4hep::DetElement(dd4hep::Detector&,
                                                    dd4hep::Volume&, int)>;
  std::map<std::string, PlaceFn> fPlaceMap;
};

#endif
