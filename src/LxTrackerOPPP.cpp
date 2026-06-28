//
/// \brief Implementation of the LxTrackerOPPP class (DD4hep version)
//

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoPara.h"
#include "TGeoArb8.h"   // TGeoTrd2 for Trapezoid

#include "LxAux.h"
#include "LxTrackerOPPP.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
LxTrackerOPPP::LxTrackerOPPP()
  : fBuilt(false),
    fCoolingPipeConnBuilt(false), fCableTerm1Built(false), fCableTerm2Built(false),
    fOutServFlatCblAssmblyZpos(0.), fInServCoolPipeAssmblyZpos(0.),
    fInServFlatCbl0AssmblyZpos(0.), fInServFlatCbl1AssmblyZpos(0.),
    fInServFlatCbl2AssmblyZpos(0.),
    fdxdetc(0.), fdydetc(0.), fdzdetc(0.),
    fdetxpos(0.), fdetypos(0.), fdetzpos(0.)
{}


void LxTrackerOPPP::FillSensorDE(dd4hep::DetElement &staveDE)
{
  for (int idet = 0; idet < fSensorPVs.size(); ++idet) {
    dd4hep::DetElement sensorDE(staveDE, "sensor_" + std::to_string(idet), idet);
    sensorDE.setPlacement(fSensorPVs[idet]);
  }
}


// ---------------------------------------------------------------------------
void LxTrackerOPPP::Build(dd4hep::Detector& description, dd4hep::SensitiveDetector& sd)
{
  if (fBuilt) return;
  fBuilt = true;

  Material opppContainerMaterial = description.material(
      description.constant<std::string>("EnvironmentMaterial"));

  double OPPPDetTopPlateX  = description.constant<double>("OPPPDetTopPlateX");
  double OPPPStaveLiftY    = description.constant<double>("OPPPStaveLiftY");
  double OPPPColdPlateY    = description.constant<double>("OPPPColdPlateY");
  double OPPPDetTopPlateZ  = description.constant<double>("OPPPDetTopPlateZ");

  fdxdetc = OPPPDetTopPlateX;
  fdydetc = OPPPStaveLiftY + OPPPColdPlateY;
  fdzdetc = OPPPDetTopPlateZ;

  Box solidOPPPTrackerL("solidOPPPTrackerL", fdxdetc/2., fdydetc/2., fdzdetc/2.);
  Box solidOPPPTrackerR("solidOPPPTrackerR", fdxdetc/2., fdydetc/2., fdzdetc/2.);
  fTrackerLVol = Volume("logicOPPPTrackerL", solidOPPPTrackerL, opppContainerMaterial);
  fTrackerRVol = Volume("logicOPPPTrackerR", solidOPPPTrackerR, opppContainerMaterial);

  // Build sensor flex (contains sensitive volumes)
  fStaveFlexVol = ConstructSensorFlex(description, sd);

  // Build sub-assemblies
  Assembly staveAssembly     = ConstructStaveAssembly(description);
  Assembly holderAssembly    = ConstructStaveHolderAssembly(description);
  Assembly innerServAssembly = ConstructInnerServiceLinesAssembly(description);
  Assembly outerServAssembly = ConstructOuterServiceLinesAssembly(description);
  Assembly innerTermAssembly = ConstructInnerServiceLinesTermAssembly(description);

  int OPPPTrackerNLayers      = description.constant<int>("OPPPTrackerNLayers");
  double OPPPDetOffsetZ       = description.constant<double>("OPPPDetOffsetZ");
  double OPPPColdPlateEpoxyThikness       = description.constant<double>("OPPPColdPlateEpoxyThikness");
  double OPPPColdPlateCarbonFiberThikness = description.constant<double>("OPPPColdPlateCarbonFiberThikness");
  double OPPPStaveSideZpos    = description.constant<double>("OPPPStaveSideZpos");
  double OPPPStaveX           = description.constant<double>("OPPPStaveX");
  double OPPPTrussZ           = description.constant<double>("OPPPTrussZ");
  double OPPPTrackerActiveX   = description.constant<double>("OPPPTrackerActiveX");
  double OPPPStaveInOutX      = description.constant<double>("OPPPStaveInOutX");
  double OPPPStaveInOutZ      = description.constant<double>("OPPPStaveInOutZ");
  double OPPPTrackerInterLayerZ = description.constant<double>("OPPPTrackerInterLayerZ");
  double OPPPSensorZ          = description.constant<double>("OPPPSensorZ");
  double OPPPSesorFPCLayer1Z  = description.constant<double>("OPPPSesorFPCLayer1Z");
  double OPPPSesorFPCLayer2Z  = description.constant<double>("OPPPSesorFPCLayer2Z");
  double OPPPSensorY          = description.constant<double>("OPPPSensorY");
  double OPPPDetYLowPos       = description.constant<double>("OPPPDetYLowPos");
  double OPPPServiceCablTermBase1X = description.constant<double>("OPPPServiceCablTermBase1X");
  bool   OverlapTest          = (description.constant<int>("OverlapTest") != 0);

  double opppdetmodz = OPPPDetOffsetZ;
  double staveextraz = OPPPColdPlateEpoxyThikness + OPPPColdPlateCarbonFiberThikness + OPPPStaveSideZpos;
  double sasxpos  = (OPPPStaveX - fdxdetc)/2.;
  double sasypos  = (OPPPColdPlateY - fdydetc)/2. + OPPPStaveLiftY;
  double saszpos  = (OPPPTrussZ - fdzdetc)/2. + staveextraz + opppdetmodz;
  double osasxpos = sasxpos + OPPPTrackerActiveX - OPPPStaveInOutX;

  Position trmin (sasxpos,  sasypos, saszpos + OPPPStaveInOutZ);
  Position trmout(osasxpos, sasypos, saszpos);

  double sensstavez = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  double flexxpos   = (OPPPStaveX - fdxdetc)/2.;
  double flexoxpos  = flexxpos + OPPPTrackerActiveX - OPPPStaveInOutX;
  double flexzpos   = -fdzdetc/2. + OPPPStaveSideZpos + opppdetmodz - sensstavez/2.;

  Position trflex   (flexxpos,  (fdydetc - OPPPSensorY)/2., flexzpos + OPPPStaveInOutZ);
  Position trflexout(flexoxpos, (fdydetc - OPPPSensorY)/2., flexzpos);

  fTrackerLDE = DetElement("detector_0", 0);
  fTrackerRDE = DetElement("detector_1", 1);

  // G4RotationMatrix(G4ThreeVector(0,0,1), pi): axis-angle around Z → RotationZYX(pi, 0, 0)
  RotationZYX rotPi(M_PI, 0., 0.);

  int nStaveFlexR = 2 * OPPPTrackerNLayers;  // stave copy counter for R tracker

  for (int il = 0; il < OPPPTrackerNLayers; ++il) {
    // Left tracker — assemblies placed into regular Volume: use placeVolume directly
    fTrackerLVol.placeVolume(staveAssembly,     trmin);
    fTrackerLVol.placeVolume(holderAssembly,    trmin);
    fTrackerLVol.placeVolume(staveAssembly,     trmout);
    fTrackerLVol.placeVolume(holderAssembly,    trmout);
    fTrackerLVol.placeVolume(innerServAssembly, trmin);
    fTrackerLVol.placeVolume(outerServAssembly, trmout);
    fTrackerLVol.placeVolume(innerTermAssembly, trmin);

    // Right tracker
    Position trminr (-trmin.X(),  trmin.Y(),  trmin.Z());
    Position trmoutr(-trmout.X(), trmout.Y(), trmout.Z());
    fTrackerRVol.placeVolume(staveAssembly,
        Transform3D(rotPi, trminr));
    fTrackerRVol.placeVolume(holderAssembly,    trminr);
    fTrackerRVol.placeVolume(staveAssembly,
        Transform3D(rotPi, trmoutr));
    fTrackerRVol.placeVolume(holderAssembly,    trmoutr);
    fTrackerRVol.placeVolume(innerServAssembly,
        Transform3D(rotPi, trminr));
    fTrackerRVol.placeVolume(outerServAssembly,
        Transform3D(rotPi, trmoutr));
    Position trminrtcterm(trmin.X() - fdxdetc + OPPPServiceCablTermBase1X,
                          trmin.Y(), trmin.Z());
    fTrackerRVol.placeVolume(innerTermAssembly, trminrtcterm);

    // Stave flex placements into L tracker (regular volume: direct placeVolume)
    PlacedVolume pvFlexL0 = fTrackerLVol.placeVolume(fStaveFlexVol, trflex);
    pvFlexL0.addPhysVolID("stave", 2*il);
    DetElement staveDE_L0(fTrackerLDE, "stave_" + std::to_string(2*il), 2*il);
    staveDE_L0.setPlacement(pvFlexL0);
    FillSensorDE(staveDE_L0);
    if (OverlapTest) pvFlexL0.ptr()->CheckOverlaps();

    PlacedVolume pvFlexL1 = fTrackerLVol.placeVolume(fStaveFlexVol, trflexout);
    pvFlexL1.addPhysVolID("stave", 2*il+1);
    DetElement staveDE_L1(fTrackerLDE, "stave_" + std::to_string(2*il+1), 2*il+1);
    staveDE_L1.setPlacement(pvFlexL1);
    FillSensorDE(staveDE_L1);
    if (OverlapTest) pvFlexL1.ptr()->CheckOverlaps();

    // Stave flex placements into R tracker
    Position trflexr   (-trflex.X(),    trflex.Y(),    trflex.Z());
    Position trflexoutr(-trflexout.X(), trflexout.Y(), trflexout.Z());
    PlacedVolume pvFlexR0 = fTrackerRVol.placeVolume(fStaveFlexVol,
        Transform3D(rotPi, trflexr));
    pvFlexR0.addPhysVolID("stave", nStaveFlexR + 2*il);
    DetElement staveDE_R0(fTrackerRDE, "stave_" + std::to_string(2*il), 2*il);
    staveDE_R0.setPlacement(pvFlexR0);
    FillSensorDE(staveDE_R0);
    if (OverlapTest) pvFlexR0.ptr()->CheckOverlaps();

    PlacedVolume pvFlexR1 = fTrackerRVol.placeVolume(fStaveFlexVol,
        Transform3D(rotPi, trflexoutr));
    pvFlexR1.addPhysVolID("stave", nStaveFlexR + 2*il + 1);
    DetElement staveDE_R1(fTrackerRDE, "stave_" + std::to_string(2*il+1), 2*il+1);
    staveDE_R1.setPlacement(pvFlexR1);
    FillSensorDE(staveDE_R1);
    if (OverlapTest) pvFlexR1.ptr()->CheckOverlaps();

    // Advance Z for next layer
    trmin    = Position(trmin.X(),    trmin.Y(),    trmin.Z()    + OPPPTrackerInterLayerZ);
    trmout   = Position(trmout.X(),   trmout.Y(),   trmout.Z()   + OPPPTrackerInterLayerZ);
    trflex   = Position(trflex.X(),   trflex.Y(),   trflex.Z()   + OPPPTrackerInterLayerZ);
    trflexout= Position(trflexout.X(),trflexout.Y(),trflexout.Z()+ OPPPTrackerInterLayerZ);
  }

  // Compute and store world placement positions for use in PlaceTrackerL/R
  double IPMagnetZpos      = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength = description.constant<double>("FlashMFieldLength");
  double OPPPDetZtoMagnet  = description.constant<double>("OPPPDetZtoMagnet");
  double OPPPDetXPos       = description.constant<double>("OPPPDetXPos");

  fdetxpos = fdxdetc/2. + OPPPDetXPos;
  fdetypos = (OPPPColdPlateY - fdydetc)/2.;
  fdetzpos = IPMagnetZpos + FlashMFieldLength/2. + OPPPDetZtoMagnet
             + fdzdetc/2. - opppdetmodz;

  // Register dispatch map
  fPlaceMap["OPPPTrackerL"] = [this](dd4hep::Detector& d, dd4hep::Volume& v, int id)
                              { return PlaceTrackerL(d, v, id); };
  fPlaceMap["OPPPTrackerR"] = [this](dd4hep::Detector& d, dd4hep::Volume& v, int id)
                              { return PlaceTrackerR(d, v, id); };

  printout(INFO, "LxTrackerOPPP::Build", "Built OPPP Tracker volumes.");
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxTrackerOPPP::PlaceTrackerL(dd4hep::Detector& description,
                                                  dd4hep::Volume&   motherVol,
                                                  int               detID)
{
  double OPPPColdPlateY    = description.constant<double>("OPPPColdPlateY");
  double OPPPDetTopPlateY  = description.constant<double>("OPPPDetTopPlateY");
  double OPPPDetBottomSupportZ = description.constant<double>("OPPPDetBottomSupportZ");
  double OPPPBasePlateX    = description.constant<double>("OPPPBasePlateX");
  double OPPPBasePlateZ    = description.constant<double>("OPPPBasePlateZ");
  double HICSBasePlateY    = description.constant<double>("HICSBasePlateY");
  double FloorSurfaceYpos  = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest       = (description.constant<int>("OverlapTest") != 0);

  PlacedVolume pvL = motherVol.placeVolume(fTrackerLVol,
      Position(fdetxpos, fdetypos, fdetzpos));
  pvL.addPhysVolID("detector", detID);
  if (OverlapTest) pvL.ptr()->CheckOverlaps();

  fTrackerLDE.setPlacement(pvL);

  // Support assembly
  double supporthight;
  Assembly supportAssembly = ConstructSupportAssembly(description, supporthight);
  double supypos = (OPPPColdPlateY - OPPPDetTopPlateY)/2. - fdydetc;
  motherVol.placeVolume(supportAssembly, Position(fdetxpos, supypos, fdetzpos));

  // Service support assembly
  double servicesupportwidth;
  Assembly serviceSupportAssembly = ConstructServiceSupportAssembly(description, servicesupportwidth);
  motherVol.placeVolume(serviceSupportAssembly,
      Position(fdxdetc + description.constant<double>("OPPPDetXPos") + servicesupportwidth/2.,
               supypos, fdetzpos));

  // Hexapod L
  double hexhightl = 0.;
  Assembly hexAssemblyL = LxAux::BuildHexapod(description, "TrackerL", hexhightl);
  double hexzpos = fdetzpos + (OPPPDetBottomSupportZ - fdzdetc)/2.;
  motherVol.placeVolume(hexAssemblyL,
      Position(fdetxpos,
               supypos + OPPPDetTopPlateY/2. - supporthight - hexhightl/2.,
               hexzpos));

  supporthight += hexhightl;

  // Table and pedestal
  double baseypos = supypos - supporthight + (OPPPColdPlateY - description.constant<double>("OPPPBasePlateY"))/2.;
  double basezpos = fdetzpos + (OPPPDetBottomSupportZ - fdzdetc)/2.;
  double ypestal  = 0.5*dd4hep::m;
  double ylevel   = baseypos + HICSBasePlateY/2.;
  double tblhight = ylevel - ypestal - FloorSurfaceYpos;

  Assembly tablesupport = LxAux::BuildTable(description, "TrackerEcalTable",
                                            OPPPBasePlateX, tblhight, OPPPBasePlateZ, 4);
  motherVol.placeVolume(tablesupport, Position(0., ylevel, basezpos));

  Volume pedestal = LxAux::BuildPedestal(description, "TrackerEcal",
                                         OPPPBasePlateX, ypestal, OPPPBasePlateZ);
  PlacedVolume pvPed = motherVol.placeVolume(pedestal,
      Position(0., 0.5*ypestal + FloorSurfaceYpos, basezpos));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();

  ConstructElectronicsRack(description, motherVol);

  return fTrackerLDE;
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxTrackerOPPP::PlaceTrackerR(dd4hep::Detector& description,
                                                  dd4hep::Volume&   motherVol,
                                                  int               detID)
{
  double OPPPColdPlateY    = description.constant<double>("OPPPColdPlateY");
  double OPPPDetTopPlateY  = description.constant<double>("OPPPDetTopPlateY");
  double OPPPDetBottomSupportZ = description.constant<double>("OPPPDetBottomSupportZ");
  double OPPPDetYLowPos    = description.constant<double>("OPPPDetYLowPos");
  bool   OverlapTest       = (description.constant<int>("OverlapTest") != 0);

  PlacedVolume pvR = motherVol.placeVolume(fTrackerRVol,
      Position(-fdetxpos, fdetypos - OPPPDetYLowPos, fdetzpos));
  pvR.addPhysVolID("detector", detID);
  if (OverlapTest) pvR.ptr()->CheckOverlaps();

  fTrackerRDE.setPlacement(pvR);

  // Support assembly
  double supporthight;
  Assembly supportAssembly = ConstructSupportAssembly(description, supporthight);
  double supypos = (OPPPColdPlateY - OPPPDetTopPlateY)/2. - fdydetc;
  motherVol.placeVolume(supportAssembly,
      Position(-fdetxpos, supypos - OPPPDetYLowPos, fdetzpos));

  // Service support assembly
  double servicesupportwidth;
  Assembly serviceSupportAssembly = ConstructServiceSupportAssembly(description, servicesupportwidth);
  // G4RotationMatrix(G4ThreeVector(0,0,1), pi) → RotationZYX(pi, 0, 0)
  motherVol.placeVolume(serviceSupportAssembly,
      Transform3D(RotationZYX(M_PI, 0., 0.),
                  Position(-(fdxdetc + description.constant<double>("OPPPDetXPos") + servicesupportwidth/2.),
                           supypos - OPPPDetYLowPos, fdetzpos)));

  // Hexapod R
  double hexhightr = OPPPDetYLowPos;
  Assembly hexAssemblyR = LxAux::BuildHexapod(description, "TrackerR", hexhightr);
  double hexhightl = 0.;  // from L build — approximate as 0 for R
  double hexzpos = fdetzpos + (OPPPDetBottomSupportZ - fdzdetc)/2.;
  motherVol.placeVolume(hexAssemblyR,
      Position(-fdetxpos,
               supypos + OPPPDetTopPlateY/2. - supporthight - hexhightl/2.
               + (hexhightl - hexhightr)/2. - OPPPDetYLowPos,
               hexzpos));

  ConstructSideServiceLines(description, motherVol);

  return fTrackerRDE;
}


// ---------------------------------------------------------------------------
dd4hep::DetElement LxTrackerOPPP::Place(dd4hep::Detector&  description,
                                         dd4hep::Volume&    motherVol,
                                         const std::string& detName,
                                         int                detID)
{
  auto it = fPlaceMap.find(detName);
  if (it == fPlaceMap.end())
    dd4hep::except("LxTrackerOPPP::Place", "Unknown detector name '%s'", detName.c_str());
  return it->second(description, motherVol, detID);
}


// ---------------------------------------------------------------------------
dd4hep::Volume LxTrackerOPPP::ConstructSensorFlex(dd4hep::Detector&          description,
                                                   dd4hep::SensitiveDetector& sd)
{
  Material sensorContainerMat = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material fpc1Mat            = description.material(description.constant<std::string>("OPPPSesorFPCMaterial1"));
  Material fpc2Mat            = description.material(description.constant<std::string>("OPPPSesorFPCMaterial2"));
  Material sensorMat          = description.material(description.constant<std::string>("OPPPTrackerMaterial"));

  double OPPPNStaveSensors    = description.constant<double>("OPPPNStaveSensors");
  double OPPPSensorX          = description.constant<double>("OPPPSensorX");
  double OPPPSensorY          = description.constant<double>("OPPPSensorY");
  double OPPPSensorZ          = description.constant<double>("OPPPSensorZ");
  double OPPPSensorPixelX     = description.constant<double>("OPPPSensorPixelX");
  double OPPPSensorPixelY     = description.constant<double>("OPPPSensorPixelY");
  double OPPPSensorPixelZ     = description.constant<double>("OPPPSensorPixelZ");
  double OPPPSensorNCellX     = description.constant<double>("OPPPSensorNCellX");
  double OPPPSensorNCellY     = description.constant<double>("OPPPSensorNCellY");
  double OPPPStaveSensorsGapX = description.constant<double>("OPPPStaveSensorsGapX");
  double OPPPSesorFPCLayer1Z  = description.constant<double>("OPPPSesorFPCLayer1Z");
  double OPPPSesorFPCLayer2Z  = description.constant<double>("OPPPSesorFPCLayer2Z");
  bool   OverlapTest          = (description.constant<int>("OverlapTest") != 0);

  int nSensors = (int)OPPPNStaveSensors;
  double sensstavex = nSensors * (OPPPStaveSensorsGapX + OPPPSensorX) - OPPPStaveSensorsGapX;
  double sensstavez = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  double sensx = OPPPSensorNCellX * OPPPSensorPixelX;
  double sensy = OPPPSensorNCellY * OPPPSensorPixelY;

  Box    solidStaveSensorContainer(sensstavex/2., OPPPSensorY/2., sensstavez/2.);
  Volume logicStaveSensorContainer("logicStaveSensorContainer",
                                    solidStaveSensorContainer, sensorContainerMat);

  // FPC layers
  Box    solidFPCLayer1(sensstavex/2., OPPPSensorY/2., OPPPSesorFPCLayer1Z/2.);
  Volume logicFPCLayer1("logicFPCLayer1", solidFPCLayer1, fpc1Mat);
  Box    solidFPCLayer2(sensstavex/2., OPPPSensorY/2., OPPPSesorFPCLayer2Z/2.);
  Volume logicFPCLayer2("logicFPCLayer2", solidFPCLayer2, fpc2Mat);

  double lzpos = (sensstavez - OPPPSesorFPCLayer1Z)/2.;
  logicStaveSensorContainer.placeVolume(logicFPCLayer1,
      Position(0., 0., lzpos)).addPhysVolID("FPCLayer1", 0);
  lzpos -= (OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z)/2.;
  logicStaveSensorContainer.placeVolume(logicFPCLayer2,
      Position(0., 0., lzpos)).addPhysVolID("FPCLayer2", 0);

  // Sensor (non-sensitive outer shell + sensitive inner)
  Box solidOPPPSensor1(OPPPSensorX/2., OPPPSensorY/2., OPPPSensorZ/2.);
  Box solidOPPPSensorSensitive(sensx/2., sensy/2., OPPPSensorPixelZ/2.);

  SubtractionSolid solidOPPPSensor("solidOPPPSensor", solidOPPPSensor1, solidOPPPSensorSensitive,
      Position(0., (OPPPSensorY - sensy)/2., (OPPPSensorPixelZ - OPPPSensorZ)/2.));

  Volume logicOPPPSensor    ("logicOPPPSensor",    solidOPPPSensor,          sensorMat);
  Volume logicOPPPSensitive ("logicOPPPSensitive", solidOPPPSensorSensitive, sensorMat);
  logicOPPPSensitive.setSensitiveDetector(sd);

  // Place sensors into container
  double detxpos = (OPPPSensorX - sensstavex)/2.;
  lzpos -= (OPPPSesorFPCLayer2Z + OPPPSensorZ)/2.;

  for (int idet = 0; idet < nSensors; ++idet) {
//cld     logicStaveSensorContainer.placeVolume(logicOPPPSensor,
//cld         Position(detxpos, 0., lzpos)).addPhysVolID("OPPPSensor", idet);
    logicStaveSensorContainer.placeVolume(logicOPPPSensor,
        Position(detxpos, 0., lzpos));
    PlacedVolume pvSens = logicStaveSensorContainer.placeVolume(logicOPPPSensitive,
        Position(detxpos,
                 (OPPPSensorY - sensy)/2.,
                 lzpos - (OPPPSensorZ - OPPPSensorPixelZ)/2.));
    pvSens.addPhysVolID("sensor", idet);
    if (OverlapTest) pvSens.ptr()->CheckOverlaps();
    fSensorPVs.push_back(pvSens);
    detxpos += OPPPSensorX + OPPPStaveSensorsGapX;
  }

  return logicStaveSensorContainer;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructStaveAssembly(dd4hep::Detector& description)
{
  Material opppColdPlateMat   = description.material(description.constant<std::string>("OPPPColdPlateMaterial"));
  Material opppSpaceFrameMat  = description.material(description.constant<std::string>("OPPPFrameMaterial"));
  Material opppCoolPipeMat    = description.material(description.constant<std::string>("OPPPCoolingPipeMaterial"));
  Material opppCoolingMat     = description.material(description.constant<std::string>("OPPPCoolingMaterial"));
  Material opppStaveSideMat   = description.material(description.constant<std::string>("OPPPStaveSideMaterial"));

  int    NTruss    = (int)description.constant<double>("OPPPNTruss");
  double dxbot     = description.constant<double>("OPPPTrussX");
  double dybot     = description.constant<double>("OPPPColdPlateY");
  double dxtop     = description.constant<double>("OPPPTrussTopXY");
  double dytop     = description.constant<double>("OPPPTrussTopXY");
  double dz        = description.constant<double>("OPPPTrussZ");
  double dtop      = description.constant<double>("OPPPTrussTopThickness");
  double dside     = 0.5*dd4hep::mm;
  double coolpiper = description.constant<double>("OPPPCoolingPipeR");
  double coolwaterr= coolpiper - description.constant<double>("OPPPCoolingPipeWall");
  double dxcplate  = description.constant<double>("OPPPColdPlateX");
  double dzcplate  = description.constant<double>("OPPPColdPlateEpoxyThikness")
                     + description.constant<double>("OPPPColdPlateCarbonFiberThikness");
  double OPPPStaveX       = description.constant<double>("OPPPStaveX");
  double OPPPStaveSideZ   = description.constant<double>("OPPPStaveSideZ");
  double OPPPStaveSideZpos= description.constant<double>("OPPPStaveSideZpos");
  bool   OverlapTest      = (description.constant<int>("OverlapTest") != 0);

  // Truss (Trapezoid)
  Trapezoid solidTruss1(dxbot/2., dxtop/2., dybot/2., dytop/2., dz/2.);
  // Cuts: solidTrussCut1 and solidTrussCut2 with Z offset -dtop
  Trapezoid solidTrussCut1(dxbot, dxbot, dybot/2., dytop/2., dz/2.);
  Trapezoid solidTrussCut2((dxbot-dside)/2., (dxtop-dside)/2., dybot, dybot, dz/2.);

  SubtractionSolid solidTrussS1("solidTrussS1", solidTruss1, solidTrussCut1,
      Position(0., 0., -dtop));
  SubtractionSolid solidTruss("solidTruss", solidTrussS1, solidTrussCut2,
      Position(0., 0., -dtop));
  Volume logicTruss("logicTruss", solidTruss, opppSpaceFrameMat);

  Assembly staveAssembly("staveAssembly");

  for (int ii = 0; ii < NTruss; ++ii)
    staveAssembly.placeVolume(logicTruss, Position(dxbot*(ii+0.5-0.5*NTruss), 0., 0.));

  Box    solidTrussBar(NTruss*dxbot/2., dytop/2., dtop/2.);
  Volume logicTrussBar("logicTrussBar", solidTrussBar, opppSpaceFrameMat);
  staveAssembly.placeVolume(logicTrussBar, Position(0., 0., (dz-dtop)/2.-dtop));

  // Cooling pipes
  double coolpipeypos = 3.*dd4hep::mm;
  Tube   solidTrkCoolPipe (0., coolpiper,  NTruss*dxbot/2., 0., 2.*M_PI);
  Tube   solidTrkCoolWater(0., coolwaterr, NTruss*dxbot/2., 0., 2.*M_PI);
  Volume logicTrkCoolPipe ("logicTrkCoolPipe",  solidTrkCoolPipe,  opppCoolPipeMat);
  Volume logicTrkCoolWater("logicTrkCoolWater", solidTrkCoolWater, opppCoolingMat);
  logicTrkCoolPipe.placeVolume(logicTrkCoolWater, Position(0.,0.,0.)).addPhysVolID("CoolWater",0);

  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2): axis-angle around Y → RotationZYX(0, pi/2, 0)
  RotationZYX rotY90(0., M_PI/2., 0.);
  staveAssembly.placeVolume(logicTrkCoolPipe,
      Transform3D(rotY90, Position(0.,  coolpipeypos, -dz/2.+coolpiper)));
  staveAssembly.placeVolume(logicTrkCoolPipe,
      Transform3D(rotY90, Position(0., -coolpipeypos, -dz/2.+coolpiper)));

  // Cold plate
  Box    solidColdPlate(dxcplate/2., dybot/2., dzcplate/2.);
  Volume logicColdPlate("logicColdPlate", solidColdPlate, opppColdPlateMat);
  staveAssembly.placeVolume(logicColdPlate, Position(0., 0., -(dzcplate+dz)/2.));

  // End caps
  double dxendcap = 0.5*(dxcplate - NTruss*dxbot);
  Trapezoid solidTrussEndCap(dxendcap/2., dxendcap/2., dybot/2., dytop/2., dz/2.);
  Volume logicTrussEndCap("logicTrussEndCap", solidTrussEndCap, opppStaveSideMat);
  staveAssembly.placeVolume(logicTrussEndCap, Position( 0.5*(dxendcap+NTruss*dxbot), 0., 0.));
  staveAssembly.placeVolume(logicTrussEndCap, Position(-0.5*(dxendcap+NTruss*dxbot), 0., 0.));

  // Stave sides
  double stavex         = OPPPStaveX;
  double stavesidey     = 10.*dd4hep::mm;
  double stavesidez     = OPPPStaveSideZ;
  double stavesidezpos  = OPPPStaveSideZpos;
  double stavesidepincutr    = 0.6*dd4hep::mm;
  double stavesidescrewcutr  = 1.0*dd4hep::mm;
  double stavesidepincutx    = 9.0*dd4hep::mm;
  double stavesidescrewcutx  = 5.0*dd4hep::mm;
  double dxstaveside = 0.5*(stavex - dxcplate);

  Box  solidStaveSide1(dxstaveside/2., stavesidey/2., stavesidez/2.);
  Tube solidStvSidePinCut  (0., stavesidepincutr,   stavesidez, 0., 2.*M_PI);
  Tube solidStvSideScrewCut(0., stavesidescrewcutr, stavesidez, 0., 2.*M_PI);

  SubtractionSolid solidStaveSide2("solidStaveSide2", solidStaveSide1, solidStvSidePinCut,
      Position(stavesidepincutx - dxstaveside/2., 0., 0.));
  SubtractionSolid solidStaveSideIn("solidStaveSideIn", solidStaveSide2, solidStvSideScrewCut,
      Position(stavesidescrewcutx - dxstaveside/2., 0., 0.));

  Volume logicStaveSideIn ("logicStaveSideIn",  solidStaveSideIn,  opppStaveSideMat);
  Volume logicStaveSideOut("logicStaveSideOut", solidStaveSide2,   opppStaveSideMat);

  double ecdz = 0.5*(stavesidez - dz) - dzcplate - stavesidezpos;
  // G4RotationMatrix(G4ThreeVector(0,0,1), pi) → RotationZYX(pi, 0, 0)
  RotationZYX rotZ180(M_PI, 0., 0.);
  staveAssembly.placeVolume(logicStaveSideIn,
      Position(-0.5*(dxcplate+dxstaveside), 0., ecdz));
  staveAssembly.placeVolume(logicStaveSideOut,
      Transform3D(rotZ180, Position(0.5*(dxcplate+dxstaveside), 0., ecdz)));

  return staveAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructSupportAssembly(dd4hep::Detector& description,
                                                          double& sphight)
{
  Material opppDetSupportMat = description.material(
      description.constant<std::string>("OPPPDetSupportMaterial"));

  double OPPPDetBottomSupportX = description.constant<double>("OPPPDetBottomSupportX");
  double OPPPDetBottomSupportY = description.constant<double>("OPPPDetBottomSupportY");
  double OPPPDetBottomSupportZ = description.constant<double>("OPPPDetBottomSupportZ");
  double OPPPDetTopPlateX      = description.constant<double>("OPPPDetTopPlateX");
  double OPPPDetTopPlateY      = description.constant<double>("OPPPDetTopPlateY");
  double OPPPDetTopPlateZ      = description.constant<double>("OPPPDetTopPlateZ");
  double OPPPDetSupportBallH   = description.constant<double>("OPPPDetSupportBallH");

  Box    solidOPPPBottomSupport(OPPPDetBottomSupportX/2., OPPPDetBottomSupportY/2., OPPPDetBottomSupportZ/2.);
  Volume logicOPPPBottomSupport("logicOPPPBottomSupport", solidOPPPBottomSupport, opppDetSupportMat);
  Box    solidOPPPTopPlate(OPPPDetTopPlateX/2., OPPPDetTopPlateY/2., OPPPDetTopPlateZ/2.);
  Volume logicOPPPTopPlate("logicOPPPTopPlate", solidOPPPTopPlate, opppDetSupportMat);
  Sphere solidOPPPDetSupportBall(0., OPPPDetSupportBallH/2., 0., M_PI, 0., 2.*M_PI);
  Volume logicOPPPDetSupportBall("logicOPPPDetSupportBall", solidOPPPDetSupportBall, opppDetSupportMat);

  Assembly oppSupportAssembly("OPPPSupportAssembly");
  oppSupportAssembly.placeVolume(logicOPPPTopPlate, Position(0., 0., 0.));

  double by = -(OPPPDetBottomSupportY + OPPPDetTopPlateY)/2. - OPPPDetSupportBallH;
  double bz = (OPPPDetBottomSupportZ - OPPPDetTopPlateZ)/2.;
  oppSupportAssembly.placeVolume(logicOPPPBottomSupport, Position(0., by, bz));

  double dbz = 0.5*OPPPDetTopPlateZ - 1.05*OPPPDetSupportBallH;
  double dbx = 0.5*OPPPDetTopPlateX - 1.05*OPPPDetSupportBallH;
  double dby = by + (OPPPDetBottomSupportY + OPPPDetSupportBallH)/2.;
  oppSupportAssembly.placeVolume(logicOPPPDetSupportBall, Position( dbx, dby,  dbz));
  oppSupportAssembly.placeVolume(logicOPPPDetSupportBall, Position(-dbx, dby,  dbz));
  oppSupportAssembly.placeVolume(logicOPPPDetSupportBall, Position( 0.,  dby, -dbz));

  sphight = OPPPDetBottomSupportY + OPPPDetTopPlateY + OPPPDetSupportBallH;
  return oppSupportAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructStaveHolderAssembly(dd4hep::Detector& description)
{
  Material mat = description.material(description.constant<std::string>("OPPPStaveHolderMaterial"));

  double BaseX    = description.constant<double>("OPPPAngleHolderBaseX");
  double BaseY    = description.constant<double>("OPPPAngleHolderBaseY");
  double HolderZ  = description.constant<double>("OPPPAngleHolderZ");
  double HolderH  = description.constant<double>("OPPPAngleHolderH");
  double BaseT    = description.constant<double>("OPPPAngleHolderBaseT");
  double LCutX    = description.constant<double>("OPPPAngleHolderLCutX");
  double LCutY    = description.constant<double>("OPPPAngleHolderLCutY");
  double RCutR    = description.constant<double>("OPPPAngleHolderRCutR");
  double HolderX  = description.constant<double>("OPPPAngleHolderX");
  double HolderY  = description.constant<double>("OPPPAngleHolderY");
  double HolderT  = description.constant<double>("OPPPAngleHolderT");
  double TopCutXY = description.constant<double>("OPPPAngleHolderTopCutXY");
  double ColdPlateY   = description.constant<double>("OPPPColdPlateY");
  double ColdPlateX   = description.constant<double>("OPPPColdPlateX");
  double ColdPlateEpoxyThikness        = description.constant<double>("OPPPColdPlateEpoxyThikness");
  double ColdPlateCarbonFiberThikness  = description.constant<double>("OPPPColdPlateCarbonFiberThikness");
  double StaveX       = description.constant<double>("OPPPStaveX");
  double StaveSideZ   = description.constant<double>("OPPPStaveSideZ");
  double StaveSideZpos= description.constant<double>("OPPPStaveSideZpos");
  double TrussZ       = description.constant<double>("OPPPTrussZ");
  double StaveLiftY   = description.constant<double>("OPPPStaveLiftY");

  // -----------------------------------------------------------------------
  // Base body: box with two cuts (one rectangular, one L+round)
  // -----------------------------------------------------------------------
  Box solidAngleHold1(BaseX/2., BaseY/2., HolderZ/2.);
  Box solidAngleHoldCut(BaseX, BaseY/2., HolderZ/2.);

  SubtractionSolid solidAngleHold2("solidAngleHold2", solidAngleHold1, solidAngleHoldCut,
      Position(0., HolderH, BaseT));

  double trlcuty = 0.5*(LCutY - BaseY) + HolderH + RCutR;
  double trlcutz = -0.5*(HolderZ - BaseT);

  Box solidAngleHoldLCut(LCutX/2., LCutY/2., BaseT);
  SubtractionSolid solidAngleHold3("solidAngleHold3", solidAngleHold2, solidAngleHoldLCut,
      Position(0., trlcuty, trlcutz));

  Tube solidAngleHoldRCut(0., RCutR, BaseT, 0., 2.*M_PI);
  SubtractionSolid solidAngleHold("solidAngleHold4", solidAngleHold3, solidAngleHoldRCut,
      Position(0., trlcuty - LCutY/2., trlcutz));

  // -----------------------------------------------------------------------
  // Tip with filet — full chain matching original active code
  // -----------------------------------------------------------------------
  double atipy      = HolderY - BaseY;
  double trtipfcutr = BaseT - HolderT;
  double trtipfcutz = HolderT;
  double trtiplcuty = trlcuty - BaseY;
  double topcutx    = std::sqrt(2.) * TopCutXY;

  Box solidAngleHoldTip0(HolderX/2., atipy/2., BaseT/2.);
  Box solidAngleHoldTipFCut(HolderX, atipy/2., BaseT/2.);

  // Cut filet corner box
  SubtractionSolid solidAngleHoldTip01("solidAngleHoldTip01",
      solidAngleHoldTip0, solidAngleHoldTipFCut,
      Position(0., trtipfcutr, trtipfcutz));

  // Round filet
  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2) → RotationZYX(0, pi/2, 0)
  Tube solidAngleHoldRFCut(0., trtipfcutr, HolderX, 0., 2.*M_PI);
  SubtractionSolid solidAngleHoldTip1("solidAngleHoldTip1",
      solidAngleHoldTip01, solidAngleHoldRFCut,
      Transform3D(RotationZYX(0., M_PI/2., 0.),
                  Position(0., trtipfcutr - atipy/2., BaseT/2.)));

  // L-cut and round cut (same as base body)
  SubtractionSolid solidAngleHoldTip2("solidAngleHoldTip2",
      solidAngleHoldTip1, solidAngleHoldLCut,
      Position(0., trtiplcuty, 0.));

  SubtractionSolid solidAngleHoldTip3("solidAngleHoldTip3",
      solidAngleHoldTip2, solidAngleHoldRCut,
      Position(0., trtiplcuty + LCutY/2., 0.));

  // Two diagonal corner cuts
  // G4RotationMatrix(G4ThreeVector(0,0,1), pi/4) → RotationZYX(pi/4, 0, 0)
  // G4RotationMatrix(G4ThreeVector(0,0,-1), pi/4) → RotationZYX(-pi/4, 0, 0)
  Box solidAngleHoldTopCut(topcutx/2., topcutx, BaseT);

  SubtractionSolid solidAngleHoldTip4("solidAngleHoldTip4",
      solidAngleHoldTip3, solidAngleHoldTopCut,
      Transform3D(RotationZYX(M_PI/4., 0., 0.),
                  Position(HolderX/2., atipy/2., 0.)));

  SubtractionSolid solidAngleHoldTip("solidAngleHoldTip",
      solidAngleHoldTip4, solidAngleHoldTopCut,
      Transform3D(RotationZYX(-M_PI/4., 0., 0.),
                  Position(-HolderX/2., atipy/2., 0.)));

  Volume logicAngleHold   ("logicAngleHold",    solidAngleHold,    mat);
  Volume logicAngleHoldTip("logicAngleHoldTip", solidAngleHoldTip, mat);

  // -----------------------------------------------------------------------
  // Assembly placement
  // -----------------------------------------------------------------------
  double dzcplate    = ColdPlateEpoxyThikness + ColdPlateCarbonFiberThikness;
  double dxstaveside = 0.5*(StaveX - ColdPlateX);
  double ahxpos      = 0.5*(ColdPlateX + dxstaveside);
  double ahypos      = (BaseY - ColdPlateY)/2. - StaveLiftY;
  double ahzpos      = 0.5*(HolderZ - TrussZ) - dzcplate - StaveSideZpos + StaveSideZ;
  double ahtipy      = ahypos + 0.5*(BaseY + atipy);
  double ahtipz      = ahzpos + 0.5*(BaseT - HolderZ);

  Assembly holderAssembly("holderAssembly");
  holderAssembly.placeVolume(logicAngleHold,    Position( ahxpos, ahypos, ahzpos));
  holderAssembly.placeVolume(logicAngleHold,    Position(-ahxpos, ahypos, ahzpos));
  holderAssembly.placeVolume(logicAngleHoldTip, Position( ahxpos, ahtipy, ahtipz));
  holderAssembly.placeVolume(logicAngleHoldTip, Position(-ahxpos, ahtipy, ahtipz));

  return holderAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructServiceSupportAssembly(dd4hep::Detector& description,
                                                                 double& width)
{
  Material opppServiceSupMat = description.material(description.constant<std::string>("OPPPDetSupportMaterial"));
  Material opppSSTMat        = description.material(description.constant<std::string>("OPPPServiceSupportTubeMaterial"));

  double OPPPServiceSupportX          = description.constant<double>("OPPPServiceSupportX");
  double OPPPServiceSupportY          = description.constant<double>("OPPPServiceSupportY");
  double OPPPServiceSupportZ          = description.constant<double>("OPPPServiceSupportZ");
  double OPPPServiceSupportTubeRin    = description.constant<double>("OPPPServiceSupportTubeRin");
  double OPPPServiceSupportTubeThickness = description.constant<double>("OPPPServiceSupportTubeThickness");
  double OPPPServiceSupportTubeX     = description.constant<double>("OPPPServiceSupportTubeX");

  Box    solidServiceSupport(OPPPServiceSupportX/2., OPPPServiceSupportY/2., OPPPServiceSupportZ/2.);
  Volume logicOPPPServiceSupport("logicOPPPServiceSupport", solidServiceSupport, opppServiceSupMat);

  double rout = OPPPServiceSupportTubeRin + OPPPServiceSupportTubeThickness;
  Tube   solidOPPPServiceSupportTube(OPPPServiceSupportTubeRin, rout,
                                     OPPPServiceSupportTubeX/2., 0., 2.*M_PI);
  Volume logicOPPPServiceSupportTube("logicOPPPServiceSupportTube", solidOPPPServiceSupportTube, opppSSTMat);

  Assembly serviceSupportAssembly("ServiceSupportAssembly");
  serviceSupportAssembly.placeVolume(logicOPPPServiceSupport,
      Position(OPPPServiceSupportTubeX/2., 0., 0.));

  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2) → RotationZYX(0, pi/2, 0)
  RotationZYX rotY90(0., M_PI/2., 0.);
  serviceSupportAssembly.placeVolume(logicOPPPServiceSupportTube,
      Transform3D(rotY90, Position(-OPPPServiceSupportX/2., 0.,  OPPPServiceSupportZ/4.)));
  serviceSupportAssembly.placeVolume(logicOPPPServiceSupportTube,
      Transform3D(rotY90, Position(-OPPPServiceSupportX/2., 0., -OPPPServiceSupportZ/4.)));

  width = OPPPServiceSupportX + OPPPServiceSupportTubeX;
  return serviceSupportAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructCoolingPipeConnectorAssembly(dd4hep::Detector& description)
{
  if (fCoolingPipeConnBuilt) return fCoolingPipeConnectorAssembly;
  fCoolingPipeConnBuilt = true;

  Material coolPipeCnctrMat = description.material(description.constant<std::string>("OPPPServiceCoolPipeCnctrMaterial"));
  Material opppCoolantMat   = description.material(description.constant<std::string>("OPPPCoolingMaterial"));
  Material coolPipeMat      = description.material(description.constant<std::string>("OPPPServiceCoolPipeMaterial"));

  double OPPPServiceCoolPipeCnctrRin  = description.constant<double>("OPPPServiceCoolPipeCnctrRin");
  double OPPPServiceCoolPipeCnctrRout = description.constant<double>("OPPPServiceCoolPipeCnctrRout");
  double OPPPServiceCoolPipeCnctrL    = description.constant<double>("OPPPServiceCoolPipeCnctrL");
  double OPPPServiceCoolPipeCnctrSRin = description.constant<double>("OPPPServiceCoolPipeCnctrSRin");
  double OPPPServiceCoolPipeCnctrSRout= description.constant<double>("OPPPServiceCoolPipeCnctrSRout");
  double OPPPServiceCoolPipeCnctrSL   = description.constant<double>("OPPPServiceCoolPipeCnctrSL");
  double OPPPServiceCoolPipeRout      = description.constant<double>("OPPPServiceCoolPipeRout");

  Tube solidInServCoolPipeCnctr   (OPPPServiceCoolPipeCnctrRin,   OPPPServiceCoolPipeCnctrRout,
                                   OPPPServiceCoolPipeCnctrL/2.,  0., 2.*M_PI);
  Tube solidInServCoolPipeCnctrW  (0., OPPPServiceCoolPipeCnctrRin,  OPPPServiceCoolPipeCnctrL/2.,  0., 2.*M_PI);
  Tube solidInServCoolPipeCnctrS  (OPPPServiceCoolPipeCnctrSRin,  OPPPServiceCoolPipeCnctrSRout,
                                   OPPPServiceCoolPipeCnctrSL/2., 0., 2.*M_PI);
  Tube solidInServCoolPipeCnctrSW (0., OPPPServiceCoolPipeCnctrSRin, OPPPServiceCoolPipeCnctrSL/2., 0., 2.*M_PI);
  Tube solidInServCoolPipeCnctrSP (OPPPServiceCoolPipeCnctrSRout, OPPPServiceCoolPipeRout,
                                   OPPPServiceCoolPipeCnctrSL/2., 0., 2.*M_PI);

  Volume logicInServCoolPipeCnctr  ("logicInServCoolPipeCnctr",   solidInServCoolPipeCnctr,   coolPipeCnctrMat);
  Volume logicInServCoolPipeCnctrW ("logicInServCoolPipeCnctrW",  solidInServCoolPipeCnctrW,  opppCoolantMat);
  Volume logicInServCoolPipeCnctrS ("logicInServCoolPipeCnctrS",  solidInServCoolPipeCnctrS,  coolPipeCnctrMat);
  Volume logicInServCoolPipeCnctrSW("logicInServCoolPipeCnctrSW", solidInServCoolPipeCnctrSW, opppCoolantMat);
  Volume logicInServCoolPipeCnctrSP("logicInServCoolPipeCnctrSP", solidInServCoolPipeCnctrSP, coolPipeMat);

  fCoolingPipeConnectorAssembly = Assembly("CoolingPipeConnectorAssembly");
  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2) → RotationZYX(0, pi/2, 0)
  RotationZYX rotY90(0., M_PI/2., 0.);
  double cx = OPPPServiceCoolPipeCnctrL/2.;
  fCoolingPipeConnectorAssembly.placeVolume(logicInServCoolPipeCnctr,
      Transform3D(rotY90, Position(cx, 0., 0.)));
  fCoolingPipeConnectorAssembly.placeVolume(logicInServCoolPipeCnctrW,
      Transform3D(rotY90, Position(cx, 0., 0.)));
  cx += (OPPPServiceCoolPipeCnctrL + OPPPServiceCoolPipeCnctrSL)/2.;
  fCoolingPipeConnectorAssembly.placeVolume(logicInServCoolPipeCnctrS,
      Transform3D(rotY90, Position(cx, 0., 0.)));
  fCoolingPipeConnectorAssembly.placeVolume(logicInServCoolPipeCnctrSW,
      Transform3D(rotY90, Position(cx, 0., 0.)));
  fCoolingPipeConnectorAssembly.placeVolume(logicInServCoolPipeCnctrSP,
      Transform3D(rotY90, Position(cx, 0., 0.)));

  return fCoolingPipeConnectorAssembly;
}


////////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------
// Helper: ecdz used throughout service line functions
// ---------------------------------------------------------------------------
static double ecdz(dd4hep::Detector& d)
{
  return 0.5*(d.constant<double>("OPPPStaveSideZ") - d.constant<double>("OPPPTrussZ"))
         - d.constant<double>("OPPPColdPlateEpoxyThikness")
         - d.constant<double>("OPPPColdPlateCarbonFiberThikness")
         - d.constant<double>("OPPPStaveSideZpos");
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructFlatCableTerminator1Assembly(dd4hep::Detector& description)
{
  if (fCableTerm1Built) return fCableTerm1Assembly;
  fCableTerm1Built = true;

  Material mat = description.material(description.constant<std::string>("OPPPServiceCablTermMaterial"));

  double ctbasey   = description.constant<double>("OPPPStaveLiftY");
  double Base1X    = description.constant<double>("OPPPServiceCablTermBase1X");
  double Base1Z    = description.constant<double>("OPPPServiceCablTermBase1Z");
  double Mid1X     = description.constant<double>("OPPPServiceCablTermMid1X");
  double Mid1Y     = description.constant<double>("OPPPServiceCablTermMid1Y");
  double Mid1Z     = description.constant<double>("OPPPServiceCablTermMid1Z");

  Box    solidBase1(Base1X/2., ctbasey/2., Base1Z/2.);
  Volume logicBase1("logicOPPPCableTermBase1", solidBase1, mat);
  Box    solidMid1 (Mid1X/2.,  Mid1Y/2.,  Mid1Z/2.);
  Volume logicMid1 ("logicOPPPCableTermMid1",  solidMid1,  mat);

  fCableTerm1Assembly = Assembly("CableTerm1Assembly");
  fCableTerm1Assembly.placeVolume(logicBase1,
      Position(Base1X/2., ctbasey/2., Mid1Z - Base1Z/2.));
  fCableTerm1Assembly.placeVolume(logicMid1,
      Position(Base1X/2., ctbasey/2. + 0.5*(ctbasey + Mid1Y), Mid1Z/2.));

  return fCableTerm1Assembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructFlatCableTerminator2Assembly(dd4hep::Detector& description)
{
  if (fCableTerm2Built) return fCableTerm2Assembly;
  fCableTerm2Built = true;

  Material mat = description.material(description.constant<std::string>("OPPPServiceCablTermMaterial"));

  double ctbasey   = description.constant<double>("OPPPStaveLiftY");
  double Base2X    = description.constant<double>("OPPPServiceCablTermBase2X");
  double Base2Z    = description.constant<double>("OPPPServiceCablTermBase2Z");
  double Mid2X     = description.constant<double>("OPPPServiceCablTermMid2X");
  double Mid2Y     = description.constant<double>("OPPPServiceCablTermMid2Y");
  double Mid2Z     = description.constant<double>("OPPPServiceCablTermMid2Z");

  Box    solidBase2(Base2X/2., ctbasey/2., Base2Z/2.);
  Volume logicBase2("logicOPPPCableTermBase2", solidBase2, mat);
  Box    solidMid2 (Mid2X/2.,  Mid2Y/2.,  Mid2Z/2.);
  Volume logicMid2 ("logicOPPPCableTermMid2",  solidMid2,  mat);

  fCableTerm2Assembly = Assembly("CableTerm2Assembly");
  fCableTerm2Assembly.placeVolume(logicBase2,
      Position(Base2X/2., ctbasey/2., 0.));
  fCableTerm2Assembly.placeVolume(logicMid2,
      Position(Base2X/2., ctbasey/2. + 0.5*(ctbasey + Mid2Y), 0.));

  return fCableTerm2Assembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructInnerServiceLinesAssembly(dd4hep::Detector& description)
{
  Material coolantMat  = description.material(description.constant<std::string>("OPPPCoolingMaterial"));
  Material coolPipeMat = description.material(description.constant<std::string>("OPPPServiceCoolPipeMaterial"));
  Material cableMat    = description.material(description.constant<std::string>("OPPPServiceCablMaterial"));

  double CoolPipeRout      = description.constant<double>("OPPPServiceCoolPipeRout");
  double CoolPipeThick     = description.constant<double>("OPPPServiceCoolPipeThick");
  double CnctrL            = description.constant<double>("OPPPServiceCoolPipeCnctrL");
  double CnctrSL           = description.constant<double>("OPPPServiceCoolPipeCnctrSL");
  double StaveX            = description.constant<double>("OPPPStaveX");
  double DetTopPlateX      = description.constant<double>("OPPPDetTopPlateX");
  double BentR             = description.constant<double>("OPPPInnerServiceCoolPipeBentR");
  double BentDZ            = description.constant<double>("OPPPInnerServiceCoolPipeDZ");
  double CoolPipeYdist     = description.constant<double>("OPPPServiceCoolPipeYdist");
  double StaveSideZpos     = description.constant<double>("OPPPStaveSideZpos");
  double StaveSideZ        = description.constant<double>("OPPPStaveSideZ");
  // flat cable sec0
  double Sec0X   = description.constant<double>("OPPPInnerServiceCablSec0X");
  double Sec0Y   = description.constant<double>("OPPPInnerServiceCablSec0Y");
  double Sec0Z   = description.constant<double>("OPPPInnerServiceCablSec0Z");
  double Sec0Th  = description.constant<double>("OPPPInnerServiceCablSec0Thick");
  double Sec1X   = description.constant<double>("OPPPInnerServiceCablSec1X");
  double Sec1Y   = description.constant<double>("OPPPInnerServiceCablSec1Y");
  double Sec1Th  = description.constant<double>("OPPPInnerServiceCablSec1Thick");
  double Sec20Y  = description.constant<double>("OPPPInnerServiceCablSec20Y");
  double Sec20Th = description.constant<double>("OPPPInnerServiceCablSec20Thick");
  double Sec21X  = description.constant<double>("OPPPInnerServiceCablSec21X");
  double Sec21Y  = description.constant<double>("OPPPInnerServiceCablSec21Y");
  double Sec21Z  = description.constant<double>("OPPPInnerServiceCablSec21Z");
  double Sec21Th = description.constant<double>("OPPPInnerServiceCablSec21Thick");
  double Sec22X  = description.constant<double>("OPPPInnerServiceCablSec22X");
  double Sec22Y  = description.constant<double>("OPPPInnerServiceCablSec22Y");
  double Sec22Z  = description.constant<double>("OPPPInnerServiceCablSec22Z");
  double Sec22Th = description.constant<double>("OPPPInnerServiceCablSec22Thick");
  double Sec31Y  = description.constant<double>("OPPPInnerServiceCablSec31Y");
  double Sec31Th = description.constant<double>("OPPPInnerServiceCablSec31Thick");
  double Sec32Y  = description.constant<double>("OPPPInnerServiceCablSec32Y");
  double Sec32Th = description.constant<double>("OPPPInnerServiceCablSec32Thick");
  double NStaveSensors   = description.constant<double>("OPPPNStaveSensors");
  double SensorX         = description.constant<double>("OPPPSensorX");
  double StaveSensorsGapX= description.constant<double>("OPPPStaveSensorsGapX");

  double coolprin = CoolPipeRout - CoolPipeThick;
  double inStaveCoolPipeX = DetTopPlateX - StaveX - CnctrL - CnctrSL;
  double pbphi = std::acos(1. - BentDZ / (2.*BentR));
  double bentx = 2.*BentR * std::sin(pbphi);
  double L0 = 0.2*(inStaveCoolPipeX - bentx);
  double L1 = inStaveCoolPipeX - bentx - L0;

  Tube   s0 (coolprin, CoolPipeRout, L0/2., 0., 2.*M_PI);
  Tube   s0W(0., coolprin,           L0/2., 0., 2.*M_PI);
  Tube   s1 (coolprin, CoolPipeRout, L1/2., 0., 2.*M_PI);
  Tube   s1W(0., coolprin,           L1/2., 0., 2.*M_PI);
//cld   Torus  sR (coolprin, CoolPipeRout, BentR, 0., pbphi);
//cld   Torus  sRW(0., coolprin,           BentR, 0., pbphi);
  Torus  sR (BentR, coolprin, CoolPipeRout, 0., pbphi);
  Torus  sRW(BentR, 0., coolprin,           0., pbphi);

  Volume pipe0 ("logicInServCoolPipe0",  s0,  coolPipeMat);
  Volume pipe0W("logicInServCoolPipe0W", s0W, coolantMat);
  Volume pipe1 ("logicInServCoolPipe1",  s1,  coolPipeMat);
  Volume pipe1W("logicInServCoolPipe1W", s1W, coolantMat);
  Volume pipeR ("logicInServCoolPipeR",  sR,  coolPipeMat);
  Volume pipeRW("logicInServCoolPipeRW", sRW, coolantMat);

  // Cooling pipe sub-assembly
  Assembly coolingPipeAssembly("InnerServiceCoolingPipeAssembly");
  Assembly connAssembly = ConstructCoolingPipeConnectorAssembly(description);
  LxAux::AddAssemblyVolumes(coolingPipeAssembly, connAssembly, Position(StaveX/2., 0., 0.));

  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2) → RotationZYX(0, pi/2, 0)
  RotationZYX rotY90(0., M_PI/2., 0.);
  double cx = StaveX/2. + CnctrL + CnctrSL/2. + (CnctrSL + L0)/2.;
  coolingPipeAssembly.placeVolume(pipe0,  Transform3D(rotY90, Position(cx, 0., 0.)));
  coolingPipeAssembly.placeVolume(pipe0W, Transform3D(rotY90, Position(cx, 0., 0.)));
  cx += bentx + (L0 + L1)/2.;
  coolingPipeAssembly.placeVolume(pipe1,  Transform3D(rotY90, Position(cx, 0., BentDZ)));
  coolingPipeAssembly.placeVolume(pipe1W, Transform3D(rotY90, Position(cx, 0., BentDZ)));
  fInServCoolPipeAssmblyZpos = BentDZ;

  double trpiperx = StaveX/2. + CnctrL + CnctrSL + L0;
  // rotateX(pi/2) then rotateY(pi/2) → RotationZYX(0, pi/2, pi/2)
//cld   RotationZYX rotpiper(0., M_PI/2., M_PI/2.);
  RotationZYX rotpiper = RotationZYX(0.0, M_PI/2.0, 0.0) * RotationZYX(0.0, 0.0, M_PI/2.0);
  coolingPipeAssembly.placeVolume(pipeR,
      Transform3D(rotpiper, Position(trpiperx, 0., BentR)));
  coolingPipeAssembly.placeVolume(pipeRW,
      Transform3D(rotpiper, Position(trpiperx, 0., BentR)));
  // rotpiper then rotateY(pi) → RotationZYX(0, 3pi/2, pi/2)
//cld   RotationZYX rotpiperflip(0., 3.*M_PI/2., M_PI/2.);
  RotationZYX rotpiperflip = RotationZYX(0.0, M_PI, 0.0) * rotpiper;
  coolingPipeAssembly.placeVolume(pipeR,
      Transform3D(rotpiperflip, Position(trpiperx + bentx, 0., BentDZ - BentR)));
  coolingPipeAssembly.placeVolume(pipeRW,
      Transform3D(rotpiperflip, Position(trpiperx + bentx, 0., BentDZ - BentR)));

  // Flat cables
  double sensstavex = (int)NStaveSensors * (StaveSensorsGapX + SensorX) - StaveSensorsGapX;
  double flatcabledx = DetTopPlateX - 0.5*(StaveX + sensstavex);
  double sec20dx = flatcabledx - Sec0X - Sec1X;
  double sec31dx = sec20dx - Sec21X;
  double sec32dx = sec20dx - Sec22X;

  // Sec0 Para: G4Para(dx, dy, dz, 0, theta, 0)
  double thet0 = std::atan2(Sec0Z, Sec0X);
  double s0dx  = Sec0Th / std::cos(thet0);
  TGeoPara* tp0 = new TGeoPara("solidInServFlatCblSec0",
      s0dx/2./dd4hep::cm, Sec0Y/2./dd4hep::cm, Sec0X/2./dd4hep::cm,
      0., thet0*180./M_PI, 0.);
  Volume vSec0("logicInServFlatCblSec0", Solid(tp0), cableMat);

  Box    bSec1(Sec1X/2., Sec1Y/2., Sec1Th/2.);
  Volume vSec1("logicInServFlatCblSec1", bSec1, cableMat);

  Box    bSec20(sec20dx/2., Sec20Y/2., Sec20Th/2.);
  Volume vSec20("logicInServFlatCblSec20", bSec20, cableMat);

  // Sec21 Para: phi = M_PI (180 degrees)
  double thet21 = std::atan2(Sec21Z, Sec21X);
  double s21dx  = Sec21Th / std::cos(thet21);
  TGeoPara* tp21 = new TGeoPara("solidInServFlatCblSec21",
      s21dx/2./dd4hep::cm, Sec21Y/2./dd4hep::cm, Sec21X/2./dd4hep::cm,
      0., thet21*180./M_PI, 180.);
  Volume vSec21("logicInServFlatCblSec21", Solid(tp21), cableMat);

  // Sec22 Para: phi = M_PI
  double thet22 = std::atan2(Sec22Z, Sec22X);
  double s22dx  = Sec22Th / std::cos(thet22);
  TGeoPara* tp22 = new TGeoPara("solidInServFlatCblSec22",
      s22dx/2./dd4hep::cm, Sec22Y/2./dd4hep::cm, Sec22X/2./dd4hep::cm,
      0., thet22*180./M_PI, 180.);
  Volume vSec22("logicInServFlatCblSec22", Solid(tp22), cableMat);

  Box    bSec31(sec31dx/2., Sec31Y/2., Sec31Th/2.);
  Volume vSec31("logicInServFlatCblSec31", bSec31, cableMat);
  Box    bSec32(sec32dx/2., Sec32Y/2., Sec32Th/2.);
  Volume vSec32("logicInServFlatCblSec32", bSec32, cableMat);

  Assembly cablesAssembly("InnerServiceCablesAssembly");
  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2) → RotationZYX(0, pi/2, 0)
  RotationZYX rotcbl(0., M_PI/2., 0.);

  double x0 = 0.5*(sensstavex + Sec0X);
  double z0 = StaveSideZpos - StaveSideZ/2. - 0.5*(Sec0Z + s0dx);
  cablesAssembly.placeVolume(vSec0,  Transform3D(rotcbl, Position(x0, 0., z0)));

  double x1 = x0 + 0.5*(Sec0X + Sec1X);
  double z1 = z0 - 0.5*(Sec0Z + s0dx - Sec1Th);
  cablesAssembly.placeVolume(vSec1, Position(x1, 0., z1));

  double x20 = x1 + 0.5*(Sec1X + sec20dx);
  double z20 = z1 + 0.5*(Sec20Th - Sec1Th);
  cablesAssembly.placeVolume(vSec20, Position(x20, 0., z20));
  fInServFlatCbl0AssmblyZpos = z20 - Sec20Th/2.;

  double x21 = x1 + 0.5*(Sec1X + Sec21X);
  double z21 = z20 + 0.5*(Sec21Z + s21dx + Sec20Th);
  cablesAssembly.placeVolume(vSec21, Transform3D(rotcbl, Position(x21, 0., z21)));

  double x22 = x1 + 0.5*(Sec1X + Sec22X);
  double z22 = z21 + 0.5*(Sec22Z + s22dx + s21dx - Sec21Z);
  cablesAssembly.placeVolume(vSec22, Transform3D(rotcbl, Position(x22, 0., z22)));

  double x31 = x21 + 0.5*(Sec21X + sec31dx);
  double z31 = z21 + 0.5*(Sec21Z + s21dx - Sec31Th);
  cablesAssembly.placeVolume(vSec31, Position(x31, 0., z31));
  fInServFlatCbl1AssmblyZpos = z31 - Sec31Th/2.;

  double x32 = x22 + 0.5*(Sec22X + sec32dx);
  double z32 = z22 + 0.5*(Sec22Z + s22dx - Sec32Th);
  cablesAssembly.placeVolume(vSec32, Position(x32, 0., z32));
  fInServFlatCbl2AssmblyZpos = z32 - Sec32Th/2.;

  // Final inner service lines assembly
  Assembly innerServiceLinesAssembly("InnerServiceLinesAssembly");
  double ez = ecdz(description);

  LxAux::AddAssemblyVolumes(innerServiceLinesAssembly, coolingPipeAssembly,
                            Position(0.,  CoolPipeYdist/2., ez));
  LxAux::AddAssemblyVolumes(innerServiceLinesAssembly, coolingPipeAssembly,
                            Position(0., -CoolPipeYdist/2., ez));
  LxAux::AddAssemblyVolumes(innerServiceLinesAssembly, cablesAssembly,
                            Position(0., 0., ez));

  return innerServiceLinesAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructInnerServiceLinesTermAssembly(dd4hep::Detector& description)
{
  Assembly cableTerm1Assembly = ConstructFlatCableTerminator1Assembly(description);

  double DetTopPlateX = description.constant<double>("OPPPDetTopPlateX");
  double StaveX       = description.constant<double>("OPPPStaveX");
  double Base1X       = description.constant<double>("OPPPServiceCablTermBase1X");
  double StaveLiftY   = description.constant<double>("OPPPStaveLiftY");
  double ColdPlateY   = description.constant<double>("OPPPColdPlateY");
  double Sec20Th      = description.constant<double>("OPPPInnerServiceCablSec20Thick");

  double ez = ecdz(description);
  double trterm1x = DetTopPlateX - StaveX/2. - Base1X;
  double trterm1y = -StaveLiftY - ColdPlateY/2.;
  double trterm1z = fInServFlatCbl0AssmblyZpos + Sec20Th + ez;

  Assembly innerFlatCableTermAssembly("InnerServiceLinesTermAssembly");
  LxAux::AddAssemblyVolumes(innerFlatCableTermAssembly, cableTerm1Assembly,
                            Position(trterm1x, trterm1y, trterm1z));
  return innerFlatCableTermAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxTrackerOPPP::ConstructOuterServiceLinesAssembly(dd4hep::Detector& description)
{
  Material coolantMat  = description.material(description.constant<std::string>("OPPPCoolingMaterial"));
  Material coolPipeMat = description.material(description.constant<std::string>("OPPPServiceCoolPipeMaterial"));
  Material cableMat    = description.material(description.constant<std::string>("OPPPServiceCablMaterial"));

  double CoolPipeRout  = description.constant<double>("OPPPServiceCoolPipeRout");
  double CoolPipeThick = description.constant<double>("OPPPServiceCoolPipeThick");
  double CnctrL        = description.constant<double>("OPPPServiceCoolPipeCnctrL");
  double CnctrSL       = description.constant<double>("OPPPServiceCoolPipeCnctrSL");
  double StaveX        = description.constant<double>("OPPPStaveX");
  double DetTopPlateX  = description.constant<double>("OPPPDetTopPlateX");
  double TrackerActiveX= description.constant<double>("OPPPTrackerActiveX");
  double StaveInOutX   = description.constant<double>("OPPPStaveInOutX");
  double CoolPipeYdist = description.constant<double>("OPPPServiceCoolPipeYdist");
  double StaveSideZpos = description.constant<double>("OPPPStaveSideZpos");
  double StaveSideZ    = description.constant<double>("OPPPStaveSideZ");
  double Sec0X  = description.constant<double>("OPPPOuterServiceCablSec0X");
  double Sec0Y  = description.constant<double>("OPPPOuterServiceCablSec0Y");
  double Sec0Z  = description.constant<double>("OPPPOuterServiceCablSec0Z");
  double Sec0Th = description.constant<double>("OPPPOuterServiceCablSec0Thick");
  double Sec1Y  = description.constant<double>("OPPPOuterServiceCablSec1Y");
  double Sec1Th = description.constant<double>("OPPPOuterServiceCablSec1Thick");

  double coolprin = CoolPipeRout - CoolPipeThick;
  double outPipeX = DetTopPlateX - StaveX - TrackerActiveX + StaveInOutX - CnctrL - CnctrSL;
  double flatcabledx = DetTopPlateX - 0.5*(StaveX + TrackerActiveX) - TrackerActiveX + StaveInOutX;
  double sec1dx = flatcabledx - Sec0X;

  Tube   sP(coolprin, CoolPipeRout, outPipeX/2., 0., 2.*M_PI);
  Tube   sPW(0., coolprin,          outPipeX/2., 0., 2.*M_PI);
  Volume logicOutPipe ("logicOutServCoolPipe0",  sP,  coolPipeMat);
  Volume logicOutPipeW("logicOutServCoolPipe0W", sPW, coolantMat);

  Assembly connAssembly = ConstructCoolingPipeConnectorAssembly(description);
  Assembly outerCoolingPipeAssembly("OuterServiceCoolingPipeAssembly");
  LxAux::AddAssemblyVolumes(outerCoolingPipeAssembly, connAssembly, Position(StaveX/2., 0., 0.));

  RotationZYX rotY90(0., M_PI/2., 0.);
  double cx = StaveX/2. + CnctrL + CnctrSL/2. + (CnctrSL + outPipeX)/2.;
  outerCoolingPipeAssembly.placeVolume(logicOutPipe,  Transform3D(rotY90, Position(cx, 0., 0.)));
  outerCoolingPipeAssembly.placeVolume(logicOutPipeW, Transform3D(rotY90, Position(cx, 0., 0.)));

  // Flat cables
  double thet0 = std::atan2(Sec0Z, Sec0X);
  double s0dx  = Sec0Th / std::cos(thet0);
  TGeoPara* tp0 = new TGeoPara("solidOutServFlatCblSec0",
      s0dx/2./dd4hep::cm, Sec0Y/2./dd4hep::cm, Sec0X/2./dd4hep::cm,
      0., thet0*180./M_PI, 0.);
  Volume vSec0("logicOutServFlatCblSec0", Solid(tp0), cableMat);

  Box    bSec1(sec1dx/2., Sec1Y/2., Sec1Th/2.);
  Volume vSec1("logicOutServFlatCblSec1", bSec1, cableMat);

  Assembly outerCablesAssembly("OuterServiceCablesAssembly");
  RotationZYX rotcbl(0., M_PI/2., 0.);

  double x0 = 0.5*(TrackerActiveX + Sec0X);
  double z0 = StaveSideZpos - StaveSideZ/2. - 0.5*(Sec0Z + s0dx);
  outerCablesAssembly.placeVolume(vSec0, Transform3D(rotcbl, Position(x0, 0., z0)));

  double x1 = x0 + 0.5*(Sec0X + sec1dx);
  double z1 = z0 - 0.5*(Sec0Z + s0dx - Sec1Th);
  outerCablesAssembly.placeVolume(vSec1, Position(x1, 0., z1));
  fOutServFlatCblAssmblyZpos = z1 - Sec1Th/2.;

  // Final outer service lines assembly
  Assembly outerServiceLinesAssembly("OuterServiceLinesAssembly");
  double ez = ecdz(description);

  LxAux::AddAssemblyVolumes(outerServiceLinesAssembly, outerCoolingPipeAssembly,
                            Position(0.,  CoolPipeYdist/2., ez));
  LxAux::AddAssemblyVolumes(outerServiceLinesAssembly, outerCoolingPipeAssembly,
                            Position(0., -CoolPipeYdist/2., ez));
  LxAux::AddAssemblyVolumes(outerServiceLinesAssembly, outerCablesAssembly,
                            Position(0., 0., ez));

  return outerServiceLinesAssembly;
}


// ---------------------------------------------------------------------------
void LxTrackerOPPP::ConstructSideServiceLinesOutAssemblies(dd4hep::Detector& description,
    dd4hep::Assembly& sideServiceAssembly,
    dd4hep::Assembly& sideFlatCableTermAssemblyL,
    dd4hep::Assembly& sideFlatCableTermAssemblyR)
{
  Material coolantMat  = description.material(description.constant<std::string>("OPPPCoolingMaterial"));
  Material coolPipeMat = description.material(description.constant<std::string>("OPPPServiceCoolPipeMaterial"));
  Material cableMat    = description.material(description.constant<std::string>("OPPPServiceCablMaterial"));

  double CoolPipeRout  = description.constant<double>("OPPPServiceCoolPipeRout");
  double CoolPipeThick = description.constant<double>("OPPPServiceCoolPipeThick");
  double ServiceSupportX  = description.constant<double>("OPPPServiceSupportX");
  double ServiceSupportTubeX = description.constant<double>("OPPPServiceSupportTubeX");
  double CoolPipeYdist = description.constant<double>("OPPPServiceCoolPipeYdist");
  double StaveLiftY    = description.constant<double>("OPPPStaveLiftY");
  double ColdPlateY    = description.constant<double>("OPPPColdPlateY");
  double Base1X        = description.constant<double>("OPPPServiceCablTermBase1X");
  double Base2X        = description.constant<double>("OPPPServiceCablTermBase2X");
  double Sec0X  = description.constant<double>("OPPPSideOutServiceCablSec0X");
  double Sec0Y  = description.constant<double>("OPPPSideOutServiceCablSec0Y");
  double Sec0Z  = description.constant<double>("OPPPSideOutServiceCablSec0Z");
  double Sec0Th = description.constant<double>("OPPPSideOutServiceCablSec0Thick");
  double Sec1X  = description.constant<double>("OPPPSideOutServiceCablSec1X");
  double Sec1Y  = description.constant<double>("OPPPSideOutServiceCablSec1Y");
  double Sec1Th = description.constant<double>("OPPPSideOutServiceCablSec1Thick");
  double Sec20X = description.constant<double>("OPPPSideOutServiceCablSec20X");
  double Sec20Y = description.constant<double>("OPPPSideOutServiceCablSec20Y");
  double Sec20Th= description.constant<double>("OPPPSideOutServiceCablSec20Thick");
  double Sec21X = description.constant<double>("OPPPSideOutServiceCablSec21X");
  double Sec21Y = description.constant<double>("OPPPSideOutServiceCablSec21Y");
  double Sec21Z = description.constant<double>("OPPPSideOutServiceCablSec21Z");
  double Sec21Th= description.constant<double>("OPPPSideOutServiceCablSec21Thick");
  double Sec22X = description.constant<double>("OPPPSideOutServiceCablSec22X");
  double Sec22Y = description.constant<double>("OPPPSideOutServiceCablSec22Y");
  double Sec22Z = description.constant<double>("OPPPSideOutServiceCablSec22Z");
  double Sec22Th= description.constant<double>("OPPPSideOutServiceCablSec22Thick");
  double Sec31X = description.constant<double>("OPPPSideOutServiceCablSec31X");
  double Sec31Y = description.constant<double>("OPPPSideOutServiceCablSec31Y");
  double Sec31Th= description.constant<double>("OPPPSideOutServiceCablSec31Thick");
  double Sec32X = description.constant<double>("OPPPSideOutServiceCablSec32X");
  double Sec32Y = description.constant<double>("OPPPSideOutServiceCablSec32Y");
  double Sec32Th= description.constant<double>("OPPPSideOutServiceCablSec32Thick");

  double coolprin     = CoolPipeRout - CoolPipeThick;
  double sideL        = ServiceSupportX + ServiceSupportTubeX;

  Tube   sPipe (0., CoolPipeRout, sideL/2., 0., 2.*M_PI);
  Tube   sPipeW(0., coolprin,     sideL/2., 0., 2.*M_PI);
  Volume logicSideServCoolPipe ("logicSideServCoolPipe",  sPipe,  coolPipeMat);
  Volume logicSideServCoolPipeW("logicSideServCoolPipeW", sPipeW, coolantMat);
  logicSideServCoolPipe.placeVolume(logicSideServCoolPipeW, Position(0.,0.,0.))
      .addPhysVolID("SideServCoolPipeCoolant", 0);

  // Flat cables
  double thet0  = std::atan2(Sec0Z, Sec0X);
  double s0dx   = Sec0Th / std::cos(thet0);
  TGeoPara* tp0 = new TGeoPara("solidSdOutServFlatCblSec0",
      s0dx/2./dd4hep::cm, Sec0Y/2./dd4hep::cm, Sec0X/2./dd4hep::cm,
      0., thet0*180./M_PI, 0.);
  Volume vSec0("logicSdOutServFlatCblSec0", Solid(tp0), cableMat);

  Box bS1(Sec1X/2., Sec1Y/2., Sec1Th/2.);   Volume vSec1 ("logicSdOutServFlatCblSec1",  bS1,  cableMat);
  Box bS20(Sec20X/2.,Sec20Y/2.,Sec20Th/2.); Volume vSec20("logicSdOutServFlatCblSec20", bS20, cableMat);

  double thet21 = std::atan2(Sec21Z, Sec21X);
  double s21dx  = Sec21Th / std::cos(thet21);
  TGeoPara* tp21 = new TGeoPara("solidSdOutServFlatCblSec21",
      s21dx/2./dd4hep::cm, Sec21Y/2./dd4hep::cm, Sec21X/2./dd4hep::cm,
      0., thet21*180./M_PI, 0.);
  Volume vSec21("logicSdOutServFlatCblSec21", Solid(tp21), cableMat);

  double thet22 = std::atan2(Sec22Z, Sec22X);
  double s22dx  = Sec22Th / std::cos(thet22);
  TGeoPara* tp22 = new TGeoPara("solidSdOutServFlatCblSec22",
      s22dx/2./dd4hep::cm, Sec22Y/2./dd4hep::cm, Sec22X/2./dd4hep::cm,
      0., thet22*180./M_PI, 0.);
  Volume vSec22("logicSdOutServFlatCblSec22", Solid(tp22), cableMat);

  Box bS31(Sec31X/2.,Sec31Y/2.,Sec31Th/2.); Volume vSec31("logicSdOutServFlatCblSec31", bS31, cableMat);
  Box bS32(Sec32X/2.,Sec32Y/2.,Sec32Th/2.); Volume vSec32("logicSdOutServFlatCblSec32", bS32, cableMat);

  Assembly cablesAssembly("OuterSideServiceCablesAssembly");
  RotationZYX rotcbl(0., M_PI/2., 0.);

  double x0 = 0.5*Sec0X;
  double z0 = fOutServFlatCblAssmblyZpos + 0.5*(s0dx - Sec0Z);
  cablesAssembly.placeVolume(vSec0, Transform3D(rotcbl, Position(x0, 0., z0)));

  double x1  = x0 + 0.5*(Sec0X + Sec1X);
  double z1  = z0 - 0.5*(Sec0Z + s0dx - Sec1Th);
  cablesAssembly.placeVolume(vSec1, Position(x1, 0., z1));

  double x20 = x1 + 0.5*(Sec1X + Sec20X);
  double z20 = z1 - 0.5*(Sec20Th - Sec1Th);
  cablesAssembly.placeVolume(vSec20, Position(x20, 0., z20));

  double x21 = x1 + 0.5*(Sec1X + Sec21X);
  double z21 = z20 - 0.5*(Sec21Z + s21dx + Sec20Th);
  cablesAssembly.placeVolume(vSec21, Transform3D(rotcbl, Position(x21, 0., z21)));

  double x22 = x1 + 0.5*(Sec1X + Sec22X);
  double z22 = z21 - 0.5*(Sec22Z + s22dx + s21dx - Sec21Z);
  cablesAssembly.placeVolume(vSec22, Transform3D(rotcbl, Position(x22, 0., z22)));

  double x31 = x21 + 0.5*(Sec21X + Sec31X);
  double z31 = z21 - 0.5*(Sec21Z + s21dx - Sec31Th);
  cablesAssembly.placeVolume(vSec31, Position(x31, 0., z31));

  double x32 = x21 + 0.5*(Sec22X + Sec32X);
  double z32 = z22 - 0.5*(Sec22Z + s22dx - Sec32Th);
  cablesAssembly.placeVolume(vSec32, Position(x32, 0., z32));

  // Assemble service components
  sideServiceAssembly = Assembly("SideOutServiceAssembly");
  double ez = ecdz(description);
  RotationZYX rotY90(0., M_PI/2., 0.);
  sideServiceAssembly.placeVolume(logicSideServCoolPipe,
      Transform3D(rotY90, Position(sideL/2.,  CoolPipeYdist/2., ez)));
  sideServiceAssembly.placeVolume(logicSideServCoolPipe,
      Transform3D(rotY90, Position(sideL/2., -CoolPipeYdist/2., ez)));
  LxAux::AddAssemblyVolumes(sideServiceAssembly, cablesAssembly, Position(0., 0., ez));

  // Terminators
  sideFlatCableTermAssemblyL = Assembly("SideOutFlatCableTermAssemblyL");
  sideFlatCableTermAssemblyR = Assembly("SideOutFlatCableTermAssemblyR");

  Assembly cTerm1 = ConstructFlatCableTerminator1Assembly(description);
  double trm1y = -StaveLiftY - ColdPlateY/2.;
  double trm1z = z32 + Sec32Th/2. + ez;
  double trm1x = x32 + Sec32X/2. - Base1X;
  LxAux::AddAssemblyVolumes(sideFlatCableTermAssemblyL, cTerm1, Position(trm1x, trm1y, trm1z));
  LxAux::AddAssemblyVolumes(sideFlatCableTermAssemblyR, cTerm1,
                            Position(-trm1x - Base1X, trm1y, trm1z));

  Assembly cTerm2 = ConstructFlatCableTerminator2Assembly(description);
  double trm2x = x31 + Sec31X/2.;
  double trm2z = (z31 + z20)/2. + ez;
  LxAux::AddAssemblyVolumes(sideFlatCableTermAssemblyL, cTerm2, Position(trm2x, trm1y, trm2z));
  LxAux::AddAssemblyVolumes(sideFlatCableTermAssemblyR, cTerm2,
                            Position(-trm2x - Base2X, trm1y, trm2z));
}


// ---------------------------------------------------------------------------
void LxTrackerOPPP::ConstructSideServiceLinesInAssemblies(dd4hep::Detector& description,
    dd4hep::Assembly& sideServiceAssembly,
    dd4hep::Assembly& sideFlatCableTermAssemblyL,
    dd4hep::Assembly& sideFlatCableTermAssemblyR)
{
  Material coolantMat  = description.material(description.constant<std::string>("OPPPCoolingMaterial"));
  Material coolPipeMat = description.material(description.constant<std::string>("OPPPServiceCoolPipeMaterial"));
  Material cableMat    = description.material(description.constant<std::string>("OPPPServiceCablMaterial"));

  double CoolPipeRout  = description.constant<double>("OPPPServiceCoolPipeRout");
  double CoolPipeThick = description.constant<double>("OPPPServiceCoolPipeThick");
  double ServiceSupportX   = description.constant<double>("OPPPServiceSupportX");
  double ServiceSupportTubeX = description.constant<double>("OPPPServiceSupportTubeX");
  double CoolPipeYdist = description.constant<double>("OPPPServiceCoolPipeYdist");
  double StaveLiftY    = description.constant<double>("OPPPStaveLiftY");
  double ColdPlateY    = description.constant<double>("OPPPColdPlateY");
  double Base2X        = description.constant<double>("OPPPServiceCablTermBase2X");
  double Sec01X  = description.constant<double>("OPPPSideInServiceCablSec01X");
  double Sec01Y  = description.constant<double>("OPPPSideInServiceCablSec01Y");
  double Sec01Th = description.constant<double>("OPPPSideInServiceCablSec01Thick");
  double Sec02X  = description.constant<double>("OPPPSideInServiceCablSec02X");
  double Sec02Y  = description.constant<double>("OPPPSideInServiceCablSec02Y");
  double Sec02Th = description.constant<double>("OPPPSideInServiceCablSec02Thick");

  double coolprin = CoolPipeRout - CoolPipeThick;
  double sideL    = ServiceSupportX + ServiceSupportTubeX;

  Tube   sPipe (0., CoolPipeRout, sideL/2., 0., 2.*M_PI);
  Tube   sPipeW(0., coolprin,     sideL/2., 0., 2.*M_PI);
  Volume logicSideInServCoolPipe ("logicSideInServCoolPipe",  sPipe,  coolPipeMat);
  Volume logicSideInServCoolPipeW("logicSideInServCoolPipeW", sPipeW, coolantMat);
  logicSideInServCoolPipe.placeVolume(logicSideInServCoolPipeW, Position(0.,0.,0.))
      .addPhysVolID("SideInServCoolPipeCoolant", 0);

  Box    bS01(Sec01X/2., Sec01Y/2., Sec01Th/2.);
  Volume vSec01("logicSdInServFlatCblSec01", bS01, cableMat);
  Box    bS02(Sec02X/2., Sec02Y/2., Sec02Th/2.);
  Volume vSec02("logicSdInServFlatCblSec02", bS02, cableMat);

  Assembly cablesAssembly("InnerSideServiceCablesAssembly");
  double x01 = 0.5*Sec01X;
  double z01 = fInServFlatCbl1AssmblyZpos + 0.5*Sec01Th;
  cablesAssembly.placeVolume(vSec01, Position(x01, 0., z01));

  double x02 = 0.5*Sec02X;
  double z02 = fInServFlatCbl2AssmblyZpos + 0.5*Sec02Th;
  cablesAssembly.placeVolume(vSec02, Position(x02, 0., z02));

  // Assemble service components
  sideServiceAssembly = Assembly("SideInServiceAssembly");
  double ez = ecdz(description);
  RotationZYX rotY90(0., M_PI/2., 0.);
  sideServiceAssembly.placeVolume(logicSideInServCoolPipe,
      Transform3D(rotY90, Position(sideL/2.,  CoolPipeYdist/2., fInServCoolPipeAssmblyZpos + ez)));
  sideServiceAssembly.placeVolume(logicSideInServCoolPipe,
      Transform3D(rotY90, Position(sideL/2., -CoolPipeYdist/2., fInServCoolPipeAssmblyZpos + ez)));
  LxAux::AddAssemblyVolumes(sideServiceAssembly, cablesAssembly, Position(0., 0., ez));

  // Terminators
  sideFlatCableTermAssemblyL = Assembly("SideInFlatCableTermAssemblyL");
  sideFlatCableTermAssemblyR = Assembly("SideInFlatCableTermAssemblyR");

  Assembly cTerm2 = ConstructFlatCableTerminator2Assembly(description);
  double trm2x = x01 + Sec01X/2.;
  double trm2y = -StaveLiftY - ColdPlateY/2.;
  double trm2z = (z01 + z02)/2. + ez;
  LxAux::AddAssemblyVolumes(sideFlatCableTermAssemblyL, cTerm2, Position(trm2x, trm2y, trm2z));
  LxAux::AddAssemblyVolumes(sideFlatCableTermAssemblyR, cTerm2,
                            Position(-trm2x - Base2X, trm2y, trm2z));
}


// ---------------------------------------------------------------------------
void LxTrackerOPPP::ConstructSideServiceLines(dd4hep::Detector& description,
                                               dd4hep::Volume&   motherVol)
{
  Material envMat = description.material(description.constant<std::string>("EnvironmentMaterial"));

  double ServiceSupportX   = description.constant<double>("OPPPServiceSupportX");
  double ServiceSupportTubeX = description.constant<double>("OPPPServiceSupportTubeX");
  double ServiceSupportZ   = description.constant<double>("OPPPServiceSupportZ");
  double StaveLiftY        = description.constant<double>("OPPPStaveLiftY");
  double ColdPlateY        = description.constant<double>("OPPPColdPlateY");
  double DetTopPlateX      = description.constant<double>("OPPPDetTopPlateX");
  double DetTopPlateZ      = description.constant<double>("OPPPDetTopPlateZ");
  double VacChambertoOPPPDetZGap = description.constant<double>("VacChambertoOPPPDetZGap");
  double OPPPDetXPos       = description.constant<double>("OPPPDetXPos");
  double OPPPDetYLowPos    = description.constant<double>("OPPPDetYLowPos");
  double TrackerNLayers    = description.constant<double>("OPPPTrackerNLayers");
  double TrackerInterLayerZ= description.constant<double>("OPPPTrackerInterLayerZ");
  double StaveInOutZ       = description.constant<double>("OPPPStaveInOutZ");
  double OPPPDetOffsetZ    = description.constant<double>("OPPPDetOffsetZ");
  bool   OverlapTest       = (description.constant<int>("OverlapTest") != 0);

  double sservconx = ServiceSupportX + ServiceSupportTubeX;
  double sservcony = StaveLiftY + ColdPlateY;
  double ssccutx   = 0.2*sservconx;
  double ssccutz   = 0.5*(ServiceSupportZ - DetTopPlateZ - VacChambertoOPPPDetZGap);

  Box  solidContainer(sservconx/2., sservcony/2., ServiceSupportZ/2.);
  Box  solidCut(ssccutx, sservcony, ssccutz);

  SubtractionSolid solidContL("solidOPPPServSideContainerL",
      solidContainer, solidCut, Position(-sservconx/2., 0., -ServiceSupportZ/2.));
  SubtractionSolid solidContR("solidOPPPServSideContainerR",
      solidContainer, solidCut, Position( sservconx/2., 0., -ServiceSupportZ/2.));

  Volume logicContL("logicOPPPServSideContainerL", solidContL, envMat);
  Volume logicContR("logicOPPPServSideContainerR", solidContR, envMat);

  // Build assemblies
  Assembly sidesOutServAssembly(""), sidesOutServCableTermAssemblyL(""), sidesOutServCableTermAssemblyR("");
  ConstructSideServiceLinesOutAssemblies(description, sidesOutServAssembly,
                                         sidesOutServCableTermAssemblyL,
                                         sidesOutServCableTermAssemblyR);

  Assembly sidesInServAssembly(""), sidesInServCableTermAssemblyL(""), sidesInServCableTermAssemblyR("");
  ConstructSideServiceLinesInAssemblies(description, sidesInServAssembly,
                                        sidesInServCableTermAssemblyL,
                                        sidesInServCableTermAssemblyR);

  // Placement loop into container volumes
  double staveextraz = description.constant<double>("OPPPColdPlateEpoxyThikness")
                       + description.constant<double>("OPPPColdPlateCarbonFiberThikness")
                       + description.constant<double>("OPPPStaveSideZpos");
  double sasxpos = -sservconx/2.;
  double sasypos =  StaveLiftY/2.;
  double saszpos = (description.constant<double>("OPPPTrussZ") - DetTopPlateZ)/2.
                   + staveextraz + OPPPDetOffsetZ;

  Position trmin (sasxpos, sasypos, saszpos + StaveInOutZ);
  Position trmout(sasxpos, sasypos, saszpos);
  RotationZYX rotPi(M_PI, 0., 0.);

  for (int il = 0; il < (int)TrackerNLayers; ++il) {
    // Into L container (regular volume — direct placeVolume)
    logicContL.placeVolume(sidesInServAssembly,         trmin);
    logicContL.placeVolume(sidesInServCableTermAssemblyL,trmin);
    logicContL.placeVolume(sidesOutServAssembly,         trmout);
    logicContL.placeVolume(sidesOutServCableTermAssemblyL,trmout);

    // Into R container
    Position trminr (-trmin.X(),  trmin.Y(),  trmin.Z());
    Position trmoutr(-trmout.X(), trmout.Y(), trmout.Z());
    logicContR.placeVolume(sidesInServAssembly,
        Transform3D(rotPi, trminr));
    logicContR.placeVolume(sidesInServCableTermAssemblyR,  trminr);
    logicContR.placeVolume(sidesOutServAssembly,
        Transform3D(rotPi, trmoutr));
    logicContR.placeVolume(sidesOutServCableTermAssemblyR, trmoutr);

    trmin  = Position(trmin.X(),  trmin.Y(),  trmin.Z()  + TrackerInterLayerZ);
    trmout = Position(trmout.X(), trmout.Y(), trmout.Z() + TrackerInterLayerZ);
  }

  // Place containers into world
  double detxpos = DetTopPlateX/2. + OPPPDetXPos;
  double detypos = -StaveLiftY/2.;
  double scposx  = detxpos + 0.5*(DetTopPlateX + sservconx);

  PlacedVolume pvL = motherVol.placeVolume(logicContL,
      Position(scposx, detypos, fdetzpos));
  pvL.addPhysVolID("OPPPServSideContainerL", 0);
  if (OverlapTest) pvL.ptr()->CheckOverlaps();

  PlacedVolume pvR = motherVol.placeVolume(logicContR,
      Position(-scposx, detypos - OPPPDetYLowPos, fdetzpos));
  pvR.addPhysVolID("OPPPServSideContainerR", 0);
  if (OverlapTest) pvR.ptr()->CheckOverlaps();
}

////////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------
void LxTrackerOPPP::ConstructElectronicsRack(dd4hep::Detector& description,
                                              dd4hep::Volume&   motherVol)
{
  Material rackContainerMat = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material rackMat          = description.material(description.constant<std::string>("TrackerElectronicsRackMaterial"));
  Material pcbMat           = description.material("Tracker_FR4");

  double TrackerElectronicsRackX     = description.constant<double>("TrackerElectronicsRackX");
  double TrackerElectronicsRackY     = description.constant<double>("TrackerElectronicsRackY");
  double TrackerElectronicsRackZ     = description.constant<double>("TrackerElectronicsRackZ");
  double TrackerElectronicsRackThick = description.constant<double>("TrackerElectronicsRackThick");
  double TrackerElectronicsPCBThick  = description.constant<double>("TrackerElectronicsPCBThick");
  double TrackerElectronicsRackXpos  = description.constant<double>("TrackerElectronicsRackXpos");
  double OPPPBasePlateX              = description.constant<double>("OPPPBasePlateX");
  double OPPPDetBottomSupportZ       = description.constant<double>("OPPPDetBottomSupportZ");
  double FloorSurfaceYpos            = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest                 = (description.constant<int>("OverlapTest") != 0);

  Box solidRackContainer(TrackerElectronicsRackX/2., TrackerElectronicsRackY/2., TrackerElectronicsRackZ/2.);
  Volume logicRackContainer("logicTrackerElectronicsRackContainer", solidRackContainer, rackContainerMat);

  // Rack walls: outer box minus inner box
  std::vector<double> racksize{TrackerElectronicsRackX, TrackerElectronicsRackY, TrackerElectronicsRackZ};
  Box solidRack1(racksize[0]/2., racksize[1]/2., racksize[2]/2.);
  std::for_each(racksize.begin(), racksize.end(), [&](double& x){ x -= 2.*TrackerElectronicsRackThick; });
  Box solidRackCut(racksize[0]/2., racksize[1]/2., racksize[2]/2.);
  SubtractionSolid solidRack("solidTrackerElectronicsRack", solidRack1, solidRackCut,
      Position(0., 0., 0.));
  Volume logicRack("logicTrackerElectronicsRack", solidRack, rackMat);

  // PCBs
  Box    solidXYPCB(racksize[0]/2., racksize[1]/2., TrackerElectronicsPCBThick/2.);
  Volume logicXYPCB("logicTrackerElectronicsXYPCB", solidXYPCB, pcbMat);

  double yzpcbz = 0.5*(racksize[2] - TrackerElectronicsPCBThick);
  Box    solidYZPCB(TrackerElectronicsPCBThick/2., racksize[1]/2., yzpcbz/2.);
  Volume logicYZPCB("logicTrackerElectronicsYZPCB", solidYZPCB, pcbMat);

  double xzpcbx = 0.5*(racksize[0] - TrackerElectronicsPCBThick);
  Box    solidXZPCB(xzpcbx/2., TrackerElectronicsPCBThick/2., yzpcbz/2.);
  Volume logicXZPCB("logicTrackerElectronicsXZPCB", solidXZPCB, pcbMat);

  logicRackContainer.placeVolume(logicRack, Position(0.,0.,0.)).addPhysVolID("RackWalls", 0);
  logicRackContainer.placeVolume(logicXYPCB,Position(0.,0.,0.)).addPhysVolID("XYPCB", 0);

  double yzpcbzpos = 0.5*(yzpcbz + TrackerElectronicsPCBThick);
  logicRackContainer.placeVolume(logicYZPCB, Position(0.,0., yzpcbzpos)).addPhysVolID("YZPCB", 0);
  logicRackContainer.placeVolume(logicYZPCB, Position(0.,0.,-yzpcbzpos)).addPhysVolID("YZPCB", 1);

  double xzpcbxpos = 0.5*(xzpcbx + TrackerElectronicsPCBThick);
  std::vector<double> xsign{1., 1., -1., -1.};
  std::vector<double> zsign{1.,-1.,  1., -1.};
  for (int ii = 0; ii < 4; ++ii)
    logicRackContainer.placeVolume(logicXZPCB,
        Position(xzpcbxpos*xsign[ii], 0., yzpcbzpos*zsign[ii])).addPhysVolID("XZPCB", ii);

  double rackxpos = OPPPBasePlateX/2. - TrackerElectronicsRackX/2. - TrackerElectronicsRackXpos;
  double basezpos = fdetzpos + (OPPPDetBottomSupportZ - fdzdetc)/2. - TrackerElectronicsRackZ;

  PlacedVolume pvRack = motherVol.placeVolume(logicRackContainer,
      Position(rackxpos, 0.5*TrackerElectronicsRackY + FloorSurfaceYpos, basezpos));
  pvRack.addPhysVolID("TrackerElectronicsRack", 0);
  if (OverlapTest) pvRack.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------
static Ref_t create_LxTrackerOPPP(dd4hep::Detector& description,
                                   xml_h             e,
                                   dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  sd.setType("tracker");
  // Retrieve or create shared LxTrackerOPPP instance
  LxTrackerOPPP* trk = nullptr;
  try {
    trk = description.extension<LxTrackerOPPP>();
    printout(INFO, "LxTrackerOPPP", "Reusing existing instance for '%s'.", detName.c_str());
  } catch (...) {
    printout(INFO, "LxTrackerOPPP", "Creating new instance for '%s'.", detName.c_str());
    trk = new LxTrackerOPPP();
    trk->Build(description, sd);
    description.addExtension<LxTrackerOPPP>(trk);
  }

  // Envelope
  Volume       fLogicWorld = description.worldVolume();
  Assembly     envelope(detName + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.attr<int>(xml_tag_t("system_id")));

  dd4hep::DetElement sdet(detName, detID);
  sdet.setPlacement(envPV);

  // Dispatch to placement method
  dd4hep::DetElement detDE = trk->Place(description, envelope, detName, detID);

  // Clone sensor DetElements for all stave placements of this tracker
  // (stave DetElements are attached inside PlaceTrackerL/R)
  sdet.add(detDE);

  return sdet;
}

DECLARE_DETELEMENT(LxTrackerOPPP, create_LxTrackerOPPP)
