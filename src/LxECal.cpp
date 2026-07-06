//
/// \brief Implementation of the LxECal class (DD4hep version)
//

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxAux.h"
#include "LxECal.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxECal::Construct(dd4hep::Detector&          description,
                       dd4hep::DetElement&        sdet,
                       xml_h&                     e,
                       dd4hep::SensitiveDetector& sd)
{
  Volume fLogicWorld = description.worldVolume();

  xml_comp_t x_det(e);
  Assembly   envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  sd.setType("calorimeter");

  // -----------------------------------------------------------------------
  // Materials
  // -----------------------------------------------------------------------
  Material ecalContainerMaterial = description.material(
      description.constant<std::string>("EnvironmentMaterial"));
  Material Air         = ecalContainerMaterial;
  Material Silicon     = description.material("G4_Si");
//   Material Aluminium   = description.material("G4_Al");
  Material Aluminium   = description.material("Aluminium");

//   Material PLASTIC_SC  = description.material("G4_PLASTIC_SC_VINYLTOLUENE");
  Material Wabsorber_PL  = description.material("Wabsorber_PL");
  Material Wabsorber_MGS = description.material("Wabsorber_MGS");
  Material FanoutMatB  = description.material("BackFanoutMaterial");
  Material FanoutMatF  = description.material("FrontFanoutMaterial");
  Material C_fiber     = description.material("ECalCarbonFiber");

  // -----------------------------------------------------------------------
  // Global constants (from lxs->...)
  // -----------------------------------------------------------------------
  double base_airx   = description.constant<double>("ECalX");
  double base_airy   = description.constant<double>("ECalY");
  int    n_layers    = (int)description.constant<double>("ECalNLayers");
  double ECalXpos    = description.constant<double>("ECalXpos");
  double ECalCasingGapZ    = description.constant<double>("ECalCasingGapZ");
  double OPPPDetTopPlateX  = description.constant<double>("OPPPDetTopPlateX");
  double OPPPDetTopPlateZ  = description.constant<double>("OPPPDetTopPlateZ");
  double OPPPDetXPos       = description.constant<double>("OPPPDetXPos");
  double OPPPTrackerECalZ  = description.constant<double>("OPPPTrackerECalZ");
  double OPPPECalPCBY      = description.constant<double>("OPPPECalPCBY");
  double ECalTopPlateZ     = description.constant<double>("ECalTopPlateZ");
  double OPPPDetTopPlateY  = description.constant<double>("OPPPDetTopPlateY");
  double IPMagnetZpos      = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength = description.constant<double>("FlashMFieldLength");
  double OPPPDetZtoMagnet  = description.constant<double>("OPPPDetZtoMagnet");
  double ECalSensorPixelX  = description.constant<double>("ECalSensorPixelX");
  double ECalSensorPixelY  = description.constant<double>("ECalSensorPixelY");
  double ECalSensorNCellX  = description.constant<double>("ECalSensorNCellX");
  double ECalSensorNCellY  = description.constant<double>("ECalSensorNCellY");
  bool   OverlapTest       = (description.constant<int>("OverlapTest") != 0);

  // TBeam scenario string — from XML string constant
  std::string TBeam_senrio = description.constant<std::string>("TBeam_senrio");

  // -----------------------------------------------------------------------
  // Local constants — direct numerical assignments moved to XML <define>
  // -----------------------------------------------------------------------
  double Lcal_absorber_pitch      = description.constant<double>("Lcal_absorber_pitch");
  double Lcal_tungsten_thickness  = description.constant<double>("Lcal_tungsten_thickness");
  double Lcal_silicon_thickness   = description.constant<double>("Lcal_silicon_thickness");
  double Lcal_epoxy_heightF       = description.constant<double>("Lcal_epoxy_heightF");
  double Lcal_kapton_heightF      = description.constant<double>("Lcal_kapton_heightF");
  double Lcal_copper_heightF      = description.constant<double>("Lcal_copper_heightF");
  double Lcal_epoxy_heightB       = description.constant<double>("Lcal_epoxy_heightB");
  double Lcal_kapton_heightB      = description.constant<double>("Lcal_kapton_heightB");
  double Lcal_copper_heightB      = description.constant<double>("Lcal_copper_heightB");
  double Lcal_pad_metal_thickness = description.constant<double>("Lcal_pad_metal_thickness");
  double CFhz                     = description.constant<double>("Lcal_CF_halfz");
  double DUTextrahz               = description.constant<double>("Lcal_DUT_extra_halfz");

  // Auxiliary derived quantities — remain in C++
  double hAbsorberPitchDZ   = Lcal_absorber_pitch / 2.;
  double hTungstenDZ        = Lcal_tungsten_thickness / 2.;
  double hSiliconDZ         = Lcal_silicon_thickness / 2.;
  double Lcal_fanoutF_thickness = Lcal_epoxy_heightF + Lcal_kapton_heightF + Lcal_copper_heightF;
  double hFanoutFrontDZ     = Lcal_fanoutF_thickness / 2.;
  double Lcal_fanoutB_thickness = Lcal_epoxy_heightB + Lcal_kapton_heightB + Lcal_copper_heightB;
  double hFanoutBackDZ      = Lcal_fanoutB_thickness / 2.;
  double hMetalDZ           = Lcal_pad_metal_thickness / 2.;

  double airhz      = hAbsorberPitchDZ + hTungstenDZ;
  double airhz1mm   = hAbsorberPitchDZ;
  double airhz_PSC  = 10.*dd4hep::mm;
  double airhz_TR_PSC = 3.75*dd4hep::mm;

  fECalLayerZ = 2.*airhz + DUTextrahz;
  double base_airz = n_layers * fECalLayerZ;

  double sensor_x = ECalSensorPixelX * ECalSensorNCellX;
  double sensor_y = ECalSensorPixelY * ECalSensorNCellY;

  // -----------------------------------------------------------------------
  // Container
  // -----------------------------------------------------------------------
  Box    solidECalContainer(base_airx/2., base_airy/2., base_airz/2.);
  Volume logicECalContainer("logicECalContainer", solidECalContainer, ecalContainerMaterial);

  // -----------------------------------------------------------------------
  // Base unit volumes (various layer types)
  // -----------------------------------------------------------------------
  Box solidBaseUnit      ("solidBaseUnit",      base_airx/2., base_airy/2., airhz);
  Box solidBaseUnitFor1mm("solidBaseUnitFor1mm",base_airx/2., base_airy/2., airhz1mm);
//   Box solidBase_PSC      ("solidBase_PSC",      base_airx/2., base_airy/2., airhz_PSC);
//   Box solidBase_TR_PSC   ("solidBase_TR_PSC",   base_airx/2., base_airy/2., airhz_TR_PSC);

  Volume logicBaseUnit         ("logicBaseUnit",         solidBaseUnit,       Air);
  Volume logicBaseUnit_T       ("logicBaseUnit_T",       solidBaseUnit,       Air);
  Volume logicBaseUnit_S       ("logicBaseUnit_S",       solidBaseUnit,       Air);
  Volume logicBaseUnit_PL      ("logicBaseUnit_PL",      solidBaseUnit,       Air);
  Volume logicBaseUnitONLYAB   ("logicBaseUnitONLYAB",   solidBaseUnit,       Air);
  Volume logicBaseUnitONLYAB_PL("logicBaseUnitONLYAB_PL",solidBaseUnit,       Air);
  Volume logicBaseUnitFor1mm   ("logicBaseUnitFor1mm",   solidBaseUnitFor1mm, Air);
//   Volume logicBaseUnit_PSC     ("logicBaseUnit_PSC",     solidBase_PSC,       PLASTIC_SC);
//   Volume logicBaseUnit_TR_PSC  ("logicBaseUnit_TR_PSC",  solidBase_TR_PSC,    PLASTIC_SC);

  // -----------------------------------------------------------------------
  // Absorber placements
  // -----------------------------------------------------------------------
  double zposAbs0     = -airhz + hTungstenDZ;
  double zposAbs0_MSG = -airhz + hTungstenDZ*1.02;

  Box    solidAbs0    ("solidAbs0",     base_airx/2., base_airy/2., hTungstenDZ);
  Box    solidAbs0_MGS("solidAbs0_MGS", base_airx/2., base_airy/2., 1.02*hTungstenDZ);
  Volume logicAbs0_PL ("logicAbs0_PL",  solidAbs0,     Wabsorber_PL);
  Volume logicAbs0_MGS("logicAbs0_MGS", solidAbs0_MGS, Wabsorber_MGS);

  logicBaseUnit.placeVolume(logicAbs0_MGS,
      Position(0.,0.,zposAbs0_MSG)).addPhysVolID("Absorber0_MGS",0);
  logicBaseUnit_PL.placeVolume(logicAbs0_PL,
      Position(0.,0.,zposAbs0)).addPhysVolID("Absorber0_PL",0);
  logicBaseUnitONLYAB.placeVolume(logicAbs0_MGS,
      Position(0.,0.,zposAbs0_MSG)).addPhysVolID("Absorber0_MGS_AB",0);
  logicBaseUnitONLYAB_PL.placeVolume(logicAbs0_PL,
      Position(0.,0.,zposAbs0)).addPhysVolID("Absorber0_PL_AB",0);

  // -----------------------------------------------------------------------
  // Carbon fiber support with sensor stack
  // -----------------------------------------------------------------------
  Box    solidCF("solidCF", base_airx/2., base_airy/2., CFhz);
  Volume logicCF("logicCF", solidCF, C_fiber);

  Box    solidECalSensor  ("solidECalSensor",   sensor_x/2., sensor_y/2., hSiliconDZ);
  Box    solidMetal       ("solidMetal",        sensor_x/2., sensor_y/2., hMetalDZ);
  Box    solidFanoutFrnt  ("solidFanoutFrnt",   sensor_x/2., sensor_y/2., hFanoutFrontDZ);
  Box    solidFanoutBack  ("solidFanoutBack",   sensor_x/2., sensor_y/2., hFanoutBackDZ);

  Volume logicECalSensor  ("logicECalSensor",   solidECalSensor,  Silicon);
  Volume logicMetalV      ("logicMetalV",       solidMetal,       Aluminium);
  Volume logicFanoutFrnt  ("logicFanoutFront",  solidFanoutFrnt,  FanoutMatF);
  Volume logicFanoutBack  ("logicFanoutBack2",  solidFanoutBack,  FanoutMatB);

  // Mark sensor as sensitive
  logicECalSensor.setSensitiveDetector(sd);

  // Place sensor stack into CF volume (once)
  double SensorAtCF = -CFhz + hFanoutFrontDZ;
  logicCF.placeVolume(logicFanoutFrnt, Position(0.,0.,SensorAtCF)).addPhysVolID("FanOut0",0);
  SensorAtCF += hFanoutFrontDZ + hMetalDZ;
  logicCF.placeVolume(logicMetalV,     Position(0.,0.,SensorAtCF)).addPhysVolID("PadMetal0",0);
  SensorAtCF += hMetalDZ + hSiliconDZ;
  PlacedVolume pvSensorInCF = logicCF.placeVolume(logicECalSensor, Position(0.,0.,SensorAtCF));
  // No addPhysVolID here — layer ID assigned at base unit placement level
  SensorAtCF += hSiliconDZ + hMetalDZ;
  logicCF.placeVolume(logicMetalV,     Position(0.,0.,SensorAtCF)).addPhysVolID("PadMetal1",0);
  SensorAtCF += hMetalDZ + hFanoutBackDZ;
  logicCF.placeVolume(logicFanoutBack, Position(0.,0.,SensorAtCF)).addPhysVolID("FanOut1",0);

  // Sensor DetElement template — cloned per layer placement
  DetElement sensorTemplate("sensor", 0);
  sensorTemplate.setPlacement(pvSensorInCF);

  // CF into base units
  logicBaseUnit.placeVolume(logicCF,
      Position(0.,0., airhz - CFhz)); //.addPhysVolID("CF0",0);
  logicBaseUnit_PL.placeVolume(logicCF,
      Position(0.,0., airhz - CFhz)); //.addPhysVolID("CF0",0);
  logicBaseUnitFor1mm.placeVolume(logicCF,
      Position(0.,0., airhz1mm - CFhz)); //.addPhysVolID("CF0",0);
  logicBaseUnit_S.placeVolume(logicCF,
      Position(0.,0., airhz - CFhz)); //.addPhysVolID("CF0",0);

  // -----------------------------------------------------------------------
  // DetElement for the container (system level)
  // -----------------------------------------------------------------------
  DetElement containerDE("ECalContainer", 0);

  // -----------------------------------------------------------------------
  // Layer placement loop driven by TBeam_senrio string
  // -----------------------------------------------------------------------
  double zpos_PSC    = 1130.*dd4hep::mm;
  double zpos_TR_PSC = 20.*dd4hep::mm;
  double ypos_stag[30] = {};  // zero-initialised
  double zyposLC = 0.;
  double zposLC  = -base_airz/2.;
  int iplacelayer   = 0;
  int iplaceElements= 0;
  int Is_SC = 0, Is_TSC = 0, Is_stag = 0;

  std::istringstream iss(TBeam_senrio);
  while (iss && iplacelayer < n_layers) {
    std::string Osub;
    iss >> Osub;
    if (Osub.empty()) break;

    // Parse optional air gap suffix  "TOKEN:N"
    std::string sub;
    int i_air = 0;
    auto colon = Osub.find(':');
    if (colon != std::string::npos) {
      i_air = std::atoi(Osub.c_str() + colon + 1) * (int)dd4hep::mm;
      sub = Osub.substr(0, colon);
    } else {
      sub = Osub;
    }

    // Helper lambda — place a base unit into container and assign DetElement
    auto placeSensorLayer = [&](Volume& vol, const std::string& prefix, double zoffset) {
      std::string pname = prefix + std::to_string(iplacelayer);
      PlacedVolume pv = logicECalContainer.placeVolume(vol,
          Position(0., zyposLC + ypos_stag[iplacelayer], zposLC + zoffset));
      pv.addPhysVolID("layer", iplacelayer);
      if (OverlapTest) pv.ptr()->CheckOverlaps();

      DetElement layerDE(containerDE, "layer_" + std::to_string(iplacelayer), iplacelayer);
      layerDE.setPlacement(pv);
      DetElement sensorDE = sensorTemplate.clone("sensor", 0);
      layerDE.add(sensorDE);
    };

    auto placeAbsLayer = [&](Volume& vol, const std::string& prefix, double zoffset) {
      std::string pname = prefix + std::to_string(iplaceElements);
      PlacedVolume pv = logicECalContainer.placeVolume(vol,
          Position(0., zyposLC, zposLC + zoffset));
      pv.addPhysVolID("absorber", iplaceElements);
      if (OverlapTest) pv.ptr()->CheckOverlaps();
    };

    if (sub == "S") {
      placeSensorLayer(logicBaseUnit_S,  "DUTAS", airhz);
      zposLC += fECalLayerZ; ++iplacelayer;
    } else if (sub == "T") {
      placeSensorLayer(logicBaseUnit_T,  "DUTAS", airhz);
      zposLC += fECalLayerZ; ++iplacelayer;
    } else if (sub == "AS") {
      placeSensorLayer(logicBaseUnit,    "DUTAS", airhz);
      zposLC += fECalLayerZ; ++iplacelayer;
    } else if (sub == "AS_PL") {
      placeSensorLayer(logicBaseUnit_PL, "DUTASP", airhz);
      zposLC += fECalLayerZ; ++iplacelayer;
    } else if (sub == "A_PL") {
      placeAbsLayer(logicBaseUnitONLYAB_PL, "ABP", airhz);
      zposLC += fECalLayerZ;
    } else if (sub == "A") {
      placeAbsLayer(logicBaseUnitONLYAB, "AB", airhz);
      zposLC += fECalLayerZ;
    } else if (sub == "JSM") {
      placeSensorLayer(logicBaseUnitFor1mm, "DUTJSM", airhz1mm);
      zposLC += 2.*airhz1mm + DUTextrahz + i_air; ++iplacelayer;
    } else if (sub == "SC") {
      if (Is_SC == 0) {
//         PlacedVolume pv = logicECalContainer.placeVolume(logicBaseUnit_PSC,
//             Position(0., zyposLC, zpos_PSC));
//         pv.addPhysVolID("SC", 1);
//         if (OverlapTest) pv.ptr()->CheckOverlaps();
        Is_SC = 1;
      }
    } else if (sub == "TSC") {
      if (Is_TSC == 0) {
//         PlacedVolume pv = logicECalContainer.placeVolume(logicBaseUnit_TR_PSC,
//             Position(0., zyposLC, zpos_TR_PSC));
//         pv.addPhysVolID("TSC", 1);
//         if (OverlapTest) pv.ptr()->CheckOverlaps();
        Is_TSC = 1;
      }
    } else if (sub == "stag") {
      if (Is_stag == 0) {
        double TB_stag[8] = {0., 0., 0.2*dd4hep::mm, -0.7*dd4hep::mm,
                             1.5*dd4hep::mm, -1.0*dd4hep::mm, 0., 0.};
        for (int istag = 0; istag < 8; ++istag)
          ypos_stag[istag] = TB_stag[istag];
        Is_stag = 1;
      }
    } else {
      break;
    }
    ++iplaceElements;
  }

  // -----------------------------------------------------------------------
  // Compute world placement positions
  // -----------------------------------------------------------------------
  double shift_x      = base_airx/2. + ECalXpos;
  double trackerzpos  = IPMagnetZpos + FlashMFieldLength/2. + OPPPDetZtoMagnet
                        + OPPPDetTopPlateZ/2.;
  double shift_z      = trackerzpos + OPPPDetTopPlateZ/2. + OPPPTrackerECalZ + base_airz/2.;

  // Place container into world
  PlacedVolume pvCont = envelope.placeVolume(logicECalContainer,
      Position(shift_x, 0., shift_z));
  pvCont.addPhysVolID("system", x_det.id());
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();
  containerDE.setPlacement(pvCont);
  sdet.add(containerDE);

  // -----------------------------------------------------------------------
  // PCB container
  // -----------------------------------------------------------------------
  Volume logicECalPCBContainer = ConstructPCB(description, base_airx, base_airz);
  PlacedVolume pvPCB = envelope.placeVolume(logicECalPCBContainer,
      Position(shift_x, 0.5*(base_airy + OPPPECalPCBY), shift_z));
  pvPCB.addPhysVolID("ECalPCB", 0);
  if (OverlapTest) pvPCB.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Casing assembly
  // -----------------------------------------------------------------------
  Assembly casingAssembly = ConstructCasingAssembly(description, base_airz);
  double casingx = shift_x;
  double casingz = shift_z + ECalCasingGapZ/2.;
  LxAux::AddAssemblyVolumes(envelope, casingAssembly, Position(casingx, 0., casingz));

  // -----------------------------------------------------------------------
  // Support assembly
  // -----------------------------------------------------------------------
  Assembly supportAssembly = ConstructSupportAssembly(description);
  double detxpos  = OPPPDetTopPlateX/2. + OPPPDetXPos;
  double supypos  = -(description.constant<double>("ECalY") + OPPPDetTopPlateY)/2.;
  double detzpos  = shift_z - base_airz/2. + ECalTopPlateZ/2.;
  LxAux::AddAssemblyVolumes(envelope, supportAssembly, Position(detxpos, supypos, detzpos));

  // -----------------------------------------------------------------------
  // Shielding
  // -----------------------------------------------------------------------
  ConstructShielding(description, envelope);
  ConstructDumpShielding(description, envelope);
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxECal::ConstructSupportAssembly(dd4hep::Detector& description)
{
  Material opppDetSupportMaterial = description.material(
      description.constant<std::string>("OPPPDetSupportMaterial"));

  double OPPPDetTopPlateX  = description.constant<double>("OPPPDetTopPlateX");
  double OPPPDetTopPlateY  = description.constant<double>("OPPPDetTopPlateY");
  double ECalTopPlateZ     = description.constant<double>("ECalTopPlateZ");
  double ECalY             = description.constant<double>("ECalY");
  double OPPPColdPlateY    = description.constant<double>("OPPPColdPlateY");
  double OPPPStaveLiftY    = description.constant<double>("OPPPStaveLiftY");
  double OPPPDetSupportBallH = description.constant<double>("OPPPDetSupportBallH");

  double trackerh       = OPPPColdPlateY/2. + OPPPStaveLiftY + OPPPDetTopPlateY + OPPPDetSupportBallH;
  double ECalSupportBallH = trackerh - ECalY/2. - OPPPDetTopPlateY;

  Box    solidECalTopPlate(OPPPDetTopPlateX/2., OPPPDetTopPlateY/2., ECalTopPlateZ/2.);
  Volume logicECalTopPlate("logicECalTopPlate", solidECalTopPlate, opppDetSupportMaterial);

  Sphere solidECalSupportBall(0., ECalSupportBallH/2., 0., M_PI, 0., 2.*M_PI);
  Volume logicECalSupportBall("logicECalSupportBall", solidECalSupportBall, opppDetSupportMaterial);

  Assembly ecalSupportAssembly("ECalSupportAssembly");
  ecalSupportAssembly.placeVolume(logicECalTopPlate, Position(0.,0.,0.));

  double dbz = 0.5*ECalTopPlateZ - 1.05*ECalSupportBallH;
  double dbx = 0.5*OPPPDetTopPlateX - 1.05*ECalSupportBallH;
  double dby = -(OPPPDetTopPlateY + ECalSupportBallH)/2.;
  ecalSupportAssembly.placeVolume(logicECalSupportBall, Position( dbx, dby,  dbz));
  ecalSupportAssembly.placeVolume(logicECalSupportBall, Position(-dbx, dby,  dbz));
  ecalSupportAssembly.placeVolume(logicECalSupportBall, Position( 0.,  dby, -dbz));

  return ecalSupportAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxECal::ConstructCasingAssembly(dd4hep::Detector& description,
                                                   double ecalz)
{
  Material ecalCasingMaterial = description.material(
      description.constant<std::string>("ECalCasingMaterial"));

  double ECalX           = description.constant<double>("ECalX");
  double ECalY           = description.constant<double>("ECalY");
  double ECalCasingTopY  = description.constant<double>("ECalCasingTopY");
  double ECalCasingSideX = description.constant<double>("ECalCasingSideX");
  double ECalCasingBackZ = description.constant<double>("ECalCasingBackZ");
  double ECalCasingGapZ  = description.constant<double>("ECalCasingGapZ");
  double OPPPECalPCBY    = description.constant<double>("OPPPECalPCBY");

  double casingz  = ecalz + ECalCasingGapZ;
  double casingy  = ECalY + ECalCasingTopY + OPPPECalPCBY;
  double sideypos = 0.5*(casingy - ECalY);

  Box    solidECalCasingTop ("solidECalCasingTop",  ECalX/2.,           ECalCasingTopY/2.,             casingz/2.);
  Box    solidECalCasingSide("solidECalCasingSide",  ECalCasingSideX/2., casingy/2.,                   (casingz+ECalCasingBackZ)/2.);
  Box    solidECalCasingBack("solidECalCasingBack",  ECalX/2.,           casingy/2.,                    ECalCasingBackZ/2.);

  Volume logicECalCasingTop ("logicECalCasingTop",  solidECalCasingTop,  ecalCasingMaterial);
  Volume logicECalCasingSide("logicECalCasingSide",  solidECalCasingSide, ecalCasingMaterial);
  Volume logicECalCasingBack("logicECalCasingBack",  solidECalCasingBack, ecalCasingMaterial);

  Assembly ecalCasingAssembly("ECalCasingAssembly");
  ecalCasingAssembly.placeVolume(logicECalCasingTop,
      Position(0., sideypos + (casingy - ECalCasingTopY)/2., 0.));
  ecalCasingAssembly.placeVolume(logicECalCasingSide,
      Position( 0.5*(ECalX + ECalCasingSideX), sideypos, ECalCasingBackZ/2.));
  ecalCasingAssembly.placeVolume(logicECalCasingSide,
      Position(-0.5*(ECalX + ECalCasingSideX), sideypos, ECalCasingBackZ/2.));
  ecalCasingAssembly.placeVolume(logicECalCasingBack,
      Position(0., sideypos, 0.5*(casingz + ECalCasingBackZ)));

  return ecalCasingAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Volume LxECal::ConstructPCB(dd4hep::Detector& description,
                                     double ecalx, double ecalz)
{
  Material envMaterial     = description.material(
      description.constant<std::string>("EnvironmentMaterial"));
  Material ecalPCBMaterial = description.material(
      description.constant<std::string>("OPPPECalPCBMaterial"));

  double OPPPECalPCBY  = description.constant<double>("OPPPECalPCBY");
  double OPPPECalPCBZ  = description.constant<double>("OPPPECalPCBZ");
  bool   OverlapTest   = (description.constant<int>("OverlapTest") != 0);

  Box    solidECalPCBContainer(ecalx/2., OPPPECalPCBY/2., ecalz/2.);
  Volume logicECalPCBContainer("logicECalPCBContainer", solidECalPCBContainer, envMaterial);

  Box    solidECalPCB(ecalx/2., OPPPECalPCBY/2., OPPPECalPCBZ/2.);
  Volume logicECalPCB("logicECalPCB", solidECalPCB, ecalPCBMaterial);

  int    n_layers = (int)description.constant<double>("ECalNLayers") - 1;
  double offsetz  = fECalLayerZ - OPPPECalPCBZ;

  for (int il = 0; il < n_layers; ++il) {
    double lzpos = il*fECalLayerZ + offsetz + (OPPPECalPCBZ - ecalz)/2.;
    PlacedVolume pv = logicECalPCBContainer.placeVolume(logicECalPCB,
        Position(0., 0., lzpos));
    pv.addPhysVolID("OPPPECalPCB", il);
    if (OverlapTest) pv.ptr()->CheckOverlaps();
  }

  return logicECalPCBContainer;
}


// ---------------------------------------------------------------------------
void LxECal::ConstructShielding(dd4hep::Detector& description,
                                 dd4hep::Volume&   motherVol)
{
  Material ecalPipeShieldMaterial = description.material(
      description.constant<std::string>("ECalPipeShieldMaterial"));

  double ECalPipeShieldX       = description.constant<double>("ECalPipeShieldX");
  double ECalPipeShieldY       = description.constant<double>("ECalPipeShieldY");
  double ECalPipeShieldZ       = description.constant<double>("ECalPipeShieldZ");
  double ECalPipeShieldShiftZ  = description.constant<double>("ECalPipeShieldShiftZ");
  double ECalXpos              = description.constant<double>("ECalXpos");
  double ECalCasingSideX       = description.constant<double>("ECalCasingSideX");
  double IPMagnetZpos          = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength     = description.constant<double>("FlashMFieldLength");
  double OPPPDetZtoMagnet      = description.constant<double>("OPPPDetZtoMagnet");
  bool   OverlapTest           = (description.constant<int>("OverlapTest") != 0);
  double VacChambertoOPPPDetZGap = description.constant<double>("VacChambertoOPPPDetZGap");

  Box    solidECalPipeShield(ECalPipeShieldX/2., ECalPipeShieldY/2., ECalPipeShieldZ/2.);
  Volume logicECalPipeShield("logicECalPipeShield", solidECalPipeShield, ecalPipeShieldMaterial);

  double bpshielxpos = ECalXpos - ECalCasingSideX - ECalPipeShieldX/2.;
  double correctionECalPipeShieldShiftZ = ECalPipeShieldShiftZ < VacChambertoOPPPDetZGap ? ECalPipeShieldShiftZ : VacChambertoOPPPDetZGap/2.0;
  double bpshielzpos = IPMagnetZpos + FlashMFieldLength/2. + OPPPDetZtoMagnet
                       - correctionECalPipeShieldShiftZ + ECalPipeShieldZ/2.;

  PlacedVolume pv0 = motherVol.placeVolume(logicECalPipeShield,
      Position( bpshielxpos, 0., bpshielzpos));
  pv0.addPhysVolID("ECalPipeShield", 0);
  if (OverlapTest) pv0.ptr()->CheckOverlaps();

  PlacedVolume pv1 = motherVol.placeVolume(logicECalPipeShield,
      Position(-bpshielxpos, 0., bpshielzpos));
  pv1.addPhysVolID("ECalPipeShield", 1);
  if (OverlapTest) pv1.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxECal::ConstructDumpShielding(dd4hep::Detector& description,
                                     dd4hep::Volume&   motherVol)
{
  Material ecalDumpShieldContainerMaterial = description.material(
      description.constant<std::string>("EnvironmentMaterial"));
  Material ecalDumpShieldMaterialInner = description.material(
      description.constant<std::string>("ECalDumpShieldMaterialInner"));
  Material ecalDumpShieldMaterialOuter = description.material(
      description.constant<std::string>("ECalDumpShieldMaterialOuter"));

  double ECalDumpShieldX      = description.constant<double>("ECalDumpShieldX");
  double ECalDumpShieldY      = description.constant<double>("ECalDumpShieldY");
  double ECalDumpShieldInnerZ = description.constant<double>("ECalDumpShieldInnerZ");
  double ECalDumpShieldOuterZ = description.constant<double>("ECalDumpShieldOuterZ");
  double ECalDumpShieldXpos   = description.constant<double>("ECalDumpShieldXpos");
  double ECalDumpShieldYpos   = description.constant<double>("ECalDumpShieldYpos");
  double ECalDumpShieldZpos   = description.constant<double>("ECalDumpShieldZpos");
  double ECalNLayers          = description.constant<double>("ECalNLayers");
  double IPMagnetZpos         = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength    = description.constant<double>("FlashMFieldLength");
  double OPPPDetZtoMagnet     = description.constant<double>("OPPPDetZtoMagnet");
  double OPPPDetTopPlateZ     = description.constant<double>("OPPPDetTopPlateZ");
  double OPPPTrackerECalZ     = description.constant<double>("OPPPTrackerECalZ");
  double FloorSurfaceYpos     = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest          = (description.constant<int>("OverlapTest") != 0);

  double shield_z = ECalDumpShieldInnerZ + 2.*ECalDumpShieldOuterZ;

  Box    solidECalDumpShieldContainer(ECalDumpShieldX/2., ECalDumpShieldY/2., shield_z/2.);
  Volume logicECalDumpShieldContainer("logicECalDumpShieldContainer",
                                       solidECalDumpShieldContainer, ecalDumpShieldContainerMaterial);

  Box    solidECalDumpShieldInner(ECalDumpShieldX/2., ECalDumpShieldY/2., ECalDumpShieldInnerZ/2.);
  Volume logicECalDumpShieldInner("logicECalDumpShieldInner", solidECalDumpShieldInner, ecalDumpShieldMaterialInner);

  Box    solidECalDumpShieldOuter(ECalDumpShieldX/2., ECalDumpShieldY/2., ECalDumpShieldOuterZ/2.);
  Volume logicECalDumpShieldOuter("logicECalDumpShieldOuter", solidECalDumpShieldOuter, ecalDumpShieldMaterialOuter);

  logicECalDumpShieldContainer.placeVolume(logicECalDumpShieldInner,
      Position(0.,0.,0.)).addPhysVolID("ECalDumpShieldInner",0);

  double dz = 0.5*(ECalDumpShieldInnerZ + ECalDumpShieldOuterZ);
  logicECalDumpShieldContainer.placeVolume(logicECalDumpShieldOuter,
      Position(0.,0.,-dz)).addPhysVolID("ECalDumpShieldOuter",0);
  logicECalDumpShieldContainer.placeVolume(logicECalDumpShieldOuter,
      Position(0.,0., dz)).addPhysVolID("ECalDumpShieldOuter",1);

  double ecalz      = ECalNLayers * fECalLayerZ;
  double trackerzpos= IPMagnetZpos + FlashMFieldLength/2. + OPPPDetZtoMagnet + OPPPDetTopPlateZ/2.;
  double ecaldsXpos = ECalDumpShieldXpos + ECalDumpShieldX/2.;
  double ecaldsYpos = ECalDumpShieldYpos;
  double ecaldsZpos = ecalz + trackerzpos + OPPPDetTopPlateZ/2. + OPPPTrackerECalZ
                      + ECalDumpShieldZpos + shield_z/2.;

  PlacedVolume pvDS = motherVol.placeVolume(logicECalDumpShieldContainer,
      Position(ecaldsXpos, ecaldsYpos, ecaldsZpos));
  pvDS.addPhysVolID("ECalDumpShield", 0);
  if (OverlapTest) pvDS.ptr()->CheckOverlaps();

  // Pedestal
  double supx = 0.9*ECalDumpShieldX;
  double supy = -FloorSurfaceYpos - ECalDumpShieldY/2.;
  double supz = 0.9*shield_z;
  Volume pedestal = LxAux::BuildPedestal(description, "ECalDumpShiel", supx, supy, supz);
  PlacedVolume pvPed = motherVol.placeVolume(pedestal,
      Position(ecaldsXpos, -(supy + ECalDumpShieldY)/2., ecaldsZpos));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------
static Ref_t create_LxECal(dd4hep::Detector& description,
                             xml_h             e,
                             dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxECal ecal;
  ecal.Construct(description, sdet, e, sd);

  return sdet;
}

DECLARE_DETELEMENT(LxECal, create_LxECal)
