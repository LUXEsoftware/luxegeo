//
/// \brief LxTrackerOPPP class (DD4hep version)
//
/// Architecture mirrors LxCherenkov/LxSScreen:
///   - Single plugin DECLARE_DETELEMENT(LxTrackerOPPP,...) called once per
///     <detector type="LxTrackerOPPP"> XML entry.
///   - First call builds all shared volumes and stores the instance via
///     description.addExtension<LxTrackerOPPP>.
///   - Each call places OPPPTrackerL or OPPPTrackerR based on detector name.
///   - Dispatch map routes detector name to placement method.
//

#ifndef LxTrackerOPPP_dd4hep_h
#define LxTrackerOPPP_dd4hep_h 1

#include <functional>
#include <map>
#include <string>
#include "DD4hep/DetFactoryHelper.h"

class LxTrackerOPPP
{
public:
  LxTrackerOPPP();
  virtual ~LxTrackerOPPP() {};

  // Build all shared volumes — called once on first plugin invocation.
  void Build(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd);

  // Place the detector named by detName into motherVol.
  dd4hep::DetElement Place(dd4hep::Detector&  description,
                           dd4hep::Volume&    motherVol,
                           const std::string& detName,
                           int                detID);

private:
  // Placement methods
  dd4hep::DetElement PlaceTrackerL(dd4hep::Detector& description,
                                   dd4hep::Volume&   motherVol, int detID);
  dd4hep::DetElement PlaceTrackerR(dd4hep::Detector& description,
                                   dd4hep::Volume&   motherVol, int detID);

  // Sub-assembly builders
  dd4hep::Assembly  ConstructStaveAssembly         (dd4hep::Detector& description);
  dd4hep::Assembly  ConstructStaveHolderAssembly   (dd4hep::Detector& description);
  dd4hep::Volume    ConstructSensorFlex            (dd4hep::Detector& description,
                                                    dd4hep::SensitiveDetector& sd);
  dd4hep::Assembly  ConstructSupportAssembly       (dd4hep::Detector& description, double& sphight);
  dd4hep::Assembly  ConstructServiceSupportAssembly(dd4hep::Detector& description, double& width);
  dd4hep::Assembly  ConstructCoolingPipeConnectorAssembly(dd4hep::Detector& description);
  dd4hep::Assembly  ConstructInnerServiceLinesAssembly   (dd4hep::Detector& description);
  dd4hep::Assembly  ConstructOuterServiceLinesAssembly   (dd4hep::Detector& description);
  dd4hep::Assembly  ConstructInnerServiceLinesTermAssembly(dd4hep::Detector& description);
  dd4hep::Assembly  ConstructFlatCableTerminator1Assembly(dd4hep::Detector& description);
  dd4hep::Assembly  ConstructFlatCableTerminator2Assembly(dd4hep::Detector& description);

  void ConstructSideServiceLines(dd4hep::Detector& description, dd4hep::Volume& motherVol);
  void ConstructElectronicsRack (dd4hep::Detector& description, dd4hep::Volume& motherVol);

  void ConstructSideServiceLinesOutAssemblies(dd4hep::Detector& description,
                                              dd4hep::Assembly& sideServiceAssembly,
                                              dd4hep::Assembly& sideFlatCableTermAssemblyL,
                                              dd4hep::Assembly& sideFlatCableTermAssemblyR);

  void ConstructSideServiceLinesInAssemblies(dd4hep::Detector& description,
                                             dd4hep::Assembly& sideServiceAssembly,
                                             dd4hep::Assembly& sideFlatCableTermAssemblyL,
                                             dd4hep::Assembly& sideFlatCableTermAssemblyR);

  void FillSensorDE(dd4hep::DetElement &staveDE);

  bool fBuilt;

  // Shared volumes built once
  dd4hep::Volume   fTrackerLVol;
  dd4hep::Volume   fTrackerRVol;
  dd4hep::Volume   fStaveFlexVol;       // logicStaveSensorContainer
  std::vector<dd4hep::PlacedVolume> fSensorPVs;
  dd4hep::DetElement fTrackerLDE;
  dd4hep::DetElement fTrackerRDE;

  // Cached assemblies (built once, reused)
  dd4hep::Assembly fCoolingPipeConnectorAssembly;
  dd4hep::Assembly fCableTerm1Assembly;
  dd4hep::Assembly fCableTerm2Assembly;
  bool fCoolingPipeConnBuilt;
  bool fCableTerm1Built;
  bool fCableTerm2Built;

  // Z positions stored for service line continuity
  double fOutServFlatCblAssmblyZpos;
  double fInServCoolPipeAssmblyZpos;
  double fInServFlatCbl0AssmblyZpos;
  double fInServFlatCbl1AssmblyZpos;
  double fInServFlatCbl2AssmblyZpos;

  // Tracker geometry dimensions (stored for reuse in placement methods)
  double fdxdetc, fdydetc, fdzdetc;
  double fdetxpos, fdetypos, fdetzpos;

  // Dispatch map
  using PlaceFn = std::function<dd4hep::DetElement(dd4hep::Detector&,
                                                    dd4hep::Volume&, int)>;
  std::map<std::string, PlaceFn> fPlaceMap;
};

#endif
