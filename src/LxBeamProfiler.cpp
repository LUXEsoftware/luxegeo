//
/// \brief Implementation of the LxBeamProfiler class (DD4hep version)
//

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxAux.h"
#include "LxBeamProfiler.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxBeamProfiler::CreateMaterial(dd4hep::Detector& description)
{
//   if (description.material("FR4").ptr()) return;
//
//   TGeoElement* H  = new TGeoElement("Hydrogen_bp", "H",  1,  1.01);
//   TGeoElement* C  = new TGeoElement("Carbon_bp",   "C",  6,  12.01);
//   TGeoElement* O  = new TGeoElement("Oxygen_bp",   "O",  8,  16.00);
//   TGeoElement* Si = new TGeoElement("Silicon_bp",  "Si", 14, 28.09);
//
//   // Epoxy
//   TGeoMixture* epoxy = new TGeoMixture("BPEpoxy", 3, 1.3);
//   epoxy->AddElement(H,  0.1310);
//   epoxy->AddElement(C,  0.5357);
//   epoxy->AddElement(O,  0.3333);
//   description.manager().AddMaterial(epoxy);
//
//   // Fiberglass (SiO2)
//   TGeoMixture* fiberglass = new TGeoMixture("fiberglass_bp", 2, 2.61);
//   fiberglass->AddElement(Si, 1);
//   fiberglass->AddElement(O,  2);
//   description.manager().AddMaterial(fiberglass);
//
//   // FR4 PCB material
//   TGeoMixture* FR4 = new TGeoMixture("FR4", 2, 1.85);
//   FR4->AddMaterial(epoxy,      0.39);
//   FR4->AddMaterial(fiberglass, 0.61);
//   description.manager().AddMaterial(FR4);
}


// ---------------------------------------------------------------------------
void LxBeamProfiler::Construct(dd4hep::Detector&          description,
                                dd4hep::DetElement&        sdet,
                                xml_h&                     e,
                                dd4hep::SensitiveDetector& sd)
{
  Volume fLogicWorld = description.worldVolume();

  xml_comp_t x_det(e);
  Assembly   envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
//   envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  sd.setType("tracker");

  CreateMaterial(description);

  Material bpContainerMaterial = description.material(
      description.constant<std::string>("BeamProfilerContanerMaterial"));
  Material bpBoxMaterial = description.material(
      description.constant<std::string>("BeamProfilerContanerWallMaterial"));
  Material bpWindowMaterial = description.material(
      description.constant<std::string>("BeamProfilerWindowMaterial"));

  double BeamProfilerContanerX    = description.constant<double>("BeamProfilerContanerX");
  double BeamProfilerContanerY    = description.constant<double>("BeamProfilerContanerY");
  double BeamProfilerContanerZ    = description.constant<double>("BeamProfilerContanerZ");
  double BeamProfilerContanerD    = description.constant<double>("BeamProfilerContanerD");
  double BeamProfilerV408Y        = description.constant<double>("BeamProfilerV408Y");
  double BeamProfilerMotorGap     = description.constant<double>("BeamProfilerMotorGap");
  double BeamProfilerQ545Y        = description.constant<double>("BeamProfilerQ545Y");
  double BeamProfilerQ545X        = description.constant<double>("BeamProfilerQ545X");
  double BeamProfilerV408X        = description.constant<double>("BeamProfilerV408X");
  double BeamProfilerPCBX         = description.constant<double>("BeamProfilerPCBX");
  double BeamProfilerPCBHolderX   = description.constant<double>("BeamProfilerPCBHolderX");
  double BeamProfilerPCBHolderZ   = description.constant<double>("BeamProfilerPCBHolderZ");
  double BeamProfilerWindowX      = description.constant<double>("BeamProfilerWindowX");
  double BeamProfilerWindowY      = description.constant<double>("BeamProfilerWindowY");
  double BeamProfilerWindowZ      = description.constant<double>("BeamProfilerWindowZ");
  double FloorSurfaceYpos         = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest              = (description.constant<int>("OverlapTest") != 0);

  // BeamProfilerZpos is a vector — read as indexed constants
  int nzpos = description.constant<int>("BeamProfilerNZpos");
  std::vector<double> BeamProfilerZpos;
  for (int i = 0; i < nzpos; ++i)
    BeamProfilerZpos.push_back(description.constant<double>(
        "BeamProfilerZpos" + std::to_string(i)));

  double zposmean = std::accumulate(BeamProfilerZpos.begin(), BeamProfilerZpos.end(), 0.0)
                    / static_cast<double>(BeamProfilerZpos.size());

  double dw            = BeamProfilerContanerD;
  double assbeambottom = BeamProfilerV408Y + BeamProfilerMotorGap + BeamProfilerQ545Y/2.0;
  double assypos       = assbeambottom + dw - BeamProfilerContanerY/2.0;
  double assxpos       = 0.5*(BeamProfilerV408X - BeamProfilerQ545X);
  double gapx          = BeamProfilerPCBX - BeamProfilerPCBHolderX;
  double pcbdx         = 0.5*(BeamProfilerQ545X + BeamProfilerPCBX) + gapx - assxpos;

  // -----------------------------------------------------------------------
  // Container box
  // -----------------------------------------------------------------------
  Box    solidBPContainer(BeamProfilerContanerX/2.0, BeamProfilerContanerY/2.0, BeamProfilerContanerZ/2.0);
  Volume logicBPContainer("logicBeamProfilerContainer", solidBPContainer, bpContainerMaterial);

  // Container wall (hollow box)
  Box solidBPBox1(BeamProfilerContanerX/2.0, BeamProfilerContanerY/2.0, BeamProfilerContanerZ/2.0);
  Box solidBPBoxCut(BeamProfilerContanerX/2.0 - dw, BeamProfilerContanerY/2.0 - dw, BeamProfilerContanerZ);
  SubtractionSolid solidBPBox("solidBPBox", solidBPBox1, solidBPBoxCut);
  Volume logicBPBox("logicBeamProfilerBox", solidBPBox, bpBoxMaterial);

  PlacedVolume pvBox = logicBPContainer.placeVolume(logicBPBox, Position(0.0, 0.0, 0.0));
  pvBox.addPhysVolID("BeamProfilerBox", 0);
  if (OverlapTest) pvBox.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Motors assembly
  // -----------------------------------------------------------------------
  Assembly motorAssembly = ConstructMotorsAssembly(description);
  logicBPContainer.placeVolume(motorAssembly, Position(assxpos, assypos, 0.0));

  // -----------------------------------------------------------------------
  // PCB+sensor assembly — placed twice (two Z positions)
  // -----------------------------------------------------------------------
  PlacedVolume pvSensor;
  Assembly pcbAssembly = ConstructPCBAssembly(description, sd, pvSensor);

  // Sensor DetElement template — cloned once per pcbAssembly placement
  DetElement sensorDE_template("sensor", 0);
  sensorDE_template.setPlacement(pvSensor);

  // DetElement for the container
  DetElement containerDE(sdet, "BeamProfilerContainer", 0);

  for (int ip = 0; ip < static_cast<int>(BeamProfilerZpos.size()); ++ip) {
    double pcbz = -BeamProfilerPCBHolderZ
                  + (ip > 0 ? BeamProfilerZpos[ip] - BeamProfilerZpos[0] : 0.0);

    Position trsenspcb(-pcbdx, assypos, pcbz);
    PlacedVolume pvpcb = logicBPContainer.placeVolume(pcbAssembly, trsenspcb);

    DetElement moduleDE(containerDE, "module_" + std::to_string(ip), ip);
    pvpcb.addPhysVolID("module", ip);
    moduleDE.setPlacement(pvpcb);

    DetElement sensorDE = sensorDE_template.clone("sensor", 0);
    moduleDE.add(sensorDE);
  }

  // -----------------------------------------------------------------------
  // Front and rear panels with windows
  // -----------------------------------------------------------------------
  double bpfrontx = BeamProfilerContanerX - 2.0*dw;
  double bpfronty = BeamProfilerContanerY - 2.0*dw;

  Box solidBPBoxFront1(bpfrontx/2.0, bpfronty/2.0, dw/2.0);
  Box solidBPBoxFrontCut(BeamProfilerWindowX/2.0, BeamProfilerWindowY/2.0, dw);

  SubtractionSolid solidBPBoxFront("solidBPBoxFront", solidBPBoxFront1, solidBPBoxFrontCut,
      Position(-pcbdx, assypos, 0.0));
  Volume logicBPBoxFront("logicBPBoxFront", solidBPBoxFront, bpBoxMaterial);

  double frontzpos = 0.5*(BeamProfilerContanerZ - dw) - BeamProfilerWindowZ;
  PlacedVolume pvFront0 = logicBPContainer.placeVolume(logicBPBoxFront, Position(0.0, 0.0, -frontzpos));
  pvFront0.addPhysVolID("BeamProfilerBoxFront", 0);
  if (OverlapTest) pvFront0.ptr()->CheckOverlaps();
  PlacedVolume pvFront1 = logicBPContainer.placeVolume(logicBPBoxFront, Position(0.0, 0.0,  frontzpos));
  pvFront1.addPhysVolID("BeamProfilerBoxFront", 1);
  if (OverlapTest) pvFront1.ptr()->CheckOverlaps();

  // Windows
  Box    solidBPBoxWindow(1.1*BeamProfilerWindowX/2.0, 1.1*BeamProfilerWindowY/2.0, BeamProfilerWindowZ/2.0);
  Volume logicBPBoxWindow("logicBPBoxWindow", solidBPBoxWindow, bpWindowMaterial);

  double windowzpos = 0.5*(BeamProfilerContanerZ - BeamProfilerWindowZ);
  PlacedVolume pvWnd0 = logicBPContainer.placeVolume(logicBPBoxWindow,
      Position(-pcbdx, assypos, -windowzpos));
  pvWnd0.addPhysVolID("BeamProfilerBoxWindow", 0);
  if (OverlapTest) pvWnd0.ptr()->CheckOverlaps();
  PlacedVolume pvWnd1 = logicBPContainer.placeVolume(logicBPBoxWindow,
      Position(-pcbdx, assypos,  windowzpos));
  pvWnd1.addPhysVolID("BeamProfilerBoxWindow", 1);
  if (OverlapTest) pvWnd1.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // PCB holder assembly
  // -----------------------------------------------------------------------
  Assembly pcbHolderAssembly = ConstructPCBSupportAssembly(description);
  logicBPContainer.placeVolume(pcbHolderAssembly, Position(-pcbdx, assypos, 0.0));

  // -----------------------------------------------------------------------
  // Place container into world
  // -----------------------------------------------------------------------
  PlacedVolume pvCont = envelope.placeVolume(logicBPContainer,
      Position(pcbdx, -assypos, zposmean));
  pvCont.addPhysVolID("system", x_det.id());
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();
  containerDE.setPlacement(pvCont);

  // -----------------------------------------------------------------------
  // Support table and pedestal
  // -----------------------------------------------------------------------
  double ypestal  = 0.5*dd4hep::m;
  double ylevel   = assbeambottom + dw;
  double tblhight = -ylevel - ypestal - FloorSurfaceYpos;
  double tblx     = 3.0*BeamProfilerContanerX;
  double tblz     = 2.0*BeamProfilerContanerZ;

  Assembly tablesupport = LxAux::BuildTable(description, "ProfilerTable", tblx, tblhight, tblz, 3);
  LxAux::AddAssemblyVolumes(envelope, tablesupport, Position(0.0, -ylevel, zposmean));

  Volume pedestal = LxAux::BuildPedestal(description, "Profiler", tblx, ypestal, tblz);
  PlacedVolume pvPed = envelope.placeVolume(pedestal,
      Position(0.0, 0.5*ypestal + FloorSurfaceYpos, zposmean));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxBeamProfiler::ConstructPCBAssembly(dd4hep::Detector&          description,
                                                       dd4hep::SensitiveDetector& sd,
                                                       dd4hep::PlacedVolume&      pvSensor)
{
  Material bpPCBMaterial       = description.material("FR4");
  Material bpSensorMaterial    = description.material(
      description.constant<std::string>("BeamProfilerSensorMaterial"));
  Material bpEnvMaterial       = description.material(
      description.constant<std::string>("BeamProfilerContanerMaterial"));
  Material bpSensMetalMaterial = description.material(
      description.constant<std::string>("BeamProfilerSensorMetalization"));

  double BeamProfilerPCBX          = description.constant<double>("BeamProfilerPCBX");
  double BeamProfilerPCBY          = description.constant<double>("BeamProfilerPCBY");
  double BeamProfilerPCBZ          = description.constant<double>("BeamProfilerPCBZ");
  double BeamProfilerPCBCutX       = description.constant<double>("BeamProfilerPCBCutX");
  double BeamProfilerPCBCutY       = description.constant<double>("BeamProfilerPCBCutY");
  double BeamProfilerSensorX       = description.constant<double>("BeamProfilerSensorX");
  double BeamProfilerSensorY       = description.constant<double>("BeamProfilerSensorY");
  double BeamProfilerSensorZ       = description.constant<double>("BeamProfilerSensorZ");
  double BeamProfilerMetalizationZ = description.constant<double>("BeamProfilerMetalizationZ");
  bool   OverlapTest               = (description.constant<int>("OverlapTest") != 0);

  // PCB board with notch cut
  Box solidBPPCB1(BeamProfilerPCBX/2.0, BeamProfilerPCBY/2.0, BeamProfilerPCBZ/2.0);
  Box solidBPPCBCut(BeamProfilerPCBCutX/2.0, BeamProfilerPCBCutY/2.0, BeamProfilerPCBZ);
  SubtractionSolid solidBPPCB("solidBPPCB", solidBPPCB1, solidBPPCBCut);
  Volume logicBPPCB("logicBPPCB", solidBPPCB, bpPCBMaterial);

  // Sensor container (holds sensor + metallization)
  double bpsensorcontainerz = BeamProfilerSensorZ + BeamProfilerMetalizationZ;
  Box    solidBPSensorContainer(BeamProfilerSensorX/2.0, BeamProfilerSensorY/2.0, bpsensorcontainerz/2.0);
  Volume logicBPSensorContainer("logicBPSensorContainer", solidBPSensorContainer, bpEnvMaterial);

  // Sensor (sensitive)
  Box    solidBPSensor(BeamProfilerSensorX/2.0, BeamProfilerSensorY/2.0, BeamProfilerSensorZ/2.0);
  Volume logicBPSensor("logicBPSensor", solidBPSensor, bpSensorMaterial);
  logicBPSensor.setSensitiveDetector(sd);

  // Place sensor into container — this is the placement we expose for DetElement cloning
  pvSensor = logicBPSensorContainer.placeVolume(logicBPSensor,
      Position(0.0, 0.0, BeamProfilerMetalizationZ/2.0));
  // Note: no addPhysVolID here — sensor has no additional ID field beyond module
  if (OverlapTest) pvSensor.ptr()->CheckOverlaps();

  // Metallization layer
  Box    solidBPSensorMetal(BeamProfilerSensorX/2.0, BeamProfilerSensorY/2.0, BeamProfilerMetalizationZ/2.0);
  Volume logicBPSensorMetal("logicBPSensorMetal", solidBPSensorMetal, bpSensMetalMaterial);
  PlacedVolume pvMetal = logicBPSensorContainer.placeVolume(logicBPSensorMetal,
      Position(0.0, 0.0, -0.5*(bpsensorcontainerz - BeamProfilerMetalizationZ)));
  pvMetal.addPhysVolID("BeamProfilerSensorMetal", 0);
  if (OverlapTest) pvMetal.ptr()->CheckOverlaps();

  // Build assembly
  Assembly bpPCBAssembly("BPPCBAssembly");

  double pcbdz = BeamProfilerPCBZ/2.0;
  bpPCBAssembly.placeVolume(logicBPPCB, Position(0.0, 0.0, -pcbdz));

  double trbpsensorz = -pcbdz - 0.5*(BeamProfilerPCBZ + bpsensorcontainerz);
  bpPCBAssembly.placeVolume(logicBPSensorContainer, Position(0.0, 0.0, trbpsensorz));

  return bpPCBAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxBeamProfiler::ConstructPCBSupportAssembly(dd4hep::Detector& description)
{
  Material bpPCBHoldMaterial = description.material(
      description.constant<std::string>("BeamProfilerPCBHolderMaterial"));

  double BeamProfilerPCBHolderX    = description.constant<double>("BeamProfilerPCBHolderX");
  double BeamProfilerPCBHolderY    = description.constant<double>("BeamProfilerPCBHolderY");
  double BeamProfilerPCBHolderZ    = description.constant<double>("BeamProfilerPCBHolderZ");
  double BeamProfilerPCBHolderCutR = description.constant<double>("BeamProfilerPCBHolderCutR");
  double BeamProfilerQ545X         = description.constant<double>("BeamProfilerQ545X");
  double BeamProfilerPCBX          = description.constant<double>("BeamProfilerPCBX");

  // Main holder with round cut
  Box  solidBPPCBHold1(BeamProfilerPCBHolderX/2.0, BeamProfilerPCBHolderY/2.0, BeamProfilerPCBHolderZ/2.0);
  Tube solidBPPCBHoldCut1(0.0, BeamProfilerPCBHolderCutR, BeamProfilerPCBHolderZ, 0.0, 2.0*M_PI);
  SubtractionSolid solidBPPCBHold("solidBPPCBHold", solidBPPCBHold1, solidBPPCBHoldCut1);
  Volume logicBPPCBHold("logicBPPCBHold", solidBPPCBHold, bpPCBHoldMaterial);

  // Side bracket
  double hdimx = BeamProfilerQ545X + BeamProfilerPCBX - BeamProfilerPCBHolderX;
  double hdimy = 0.9 * BeamProfilerPCBHolderY;
  Box solidBPPCBHoldSide1(hdimx/2.0, hdimy/2.0, BeamProfilerPCBHolderZ/2.0);
  hdimy -= 20.0*dd4hep::mm;
  Box solidBPPCBHoldSideCut(hdimx/2.0, hdimy/2.0, BeamProfilerPCBHolderZ);

  double hdxcut = 10.0*dd4hep::mm;
  SubtractionSolid solidBPPCBHoldSide("solidBPPCBHoldSide", solidBPPCBHoldSide1, solidBPPCBHoldSideCut,
      Position(-hdxcut, 0.0, 0.0));
  Volume logicBPPCBHoldSide("logicBPPCBHoldSide", solidBPPCBHoldSide, bpPCBHoldMaterial);

  Assembly bpPCBSupportAssembly("BPPCBSupportAssembly");

  bpPCBSupportAssembly.placeVolume(logicBPPCBHold,
      Position(0.0, 0.0, -BeamProfilerPCBHolderZ/2.0));

  double trpcbholdsidex = 0.5*(hdimx + BeamProfilerPCBHolderX);
  bpPCBSupportAssembly.placeVolume(logicBPPCBHoldSide,
      Position(trpcbholdsidex, 0.0, -BeamProfilerPCBHolderZ/2.0));

  return bpPCBSupportAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxBeamProfiler::ConstructMotorsAssembly(dd4hep::Detector& description)
{
  Material bpV408Material = description.material(
      description.constant<std::string>("BeamProfilerV408Material"));
  Material bpQ545Material = description.material(
      description.constant<std::string>("BeamProfilerQ545Material"));
  Material bpaMaterial    = description.material(
      description.constant<std::string>("BeamProfilerAngleSupportMaterial"));

  double BeamProfilerV408X          = description.constant<double>("BeamProfilerV408X");
  double BeamProfilerV408Y          = description.constant<double>("BeamProfilerV408Y");
  double BeamProfilerV408Z          = description.constant<double>("BeamProfilerV408Z");
  double BeamProfilerQ545X          = description.constant<double>("BeamProfilerQ545X");
  double BeamProfilerQ545Y          = description.constant<double>("BeamProfilerQ545Y");
  double BeamProfilerQ545Z          = description.constant<double>("BeamProfilerQ545Z");
  double BeamProfilerMotorGap       = description.constant<double>("BeamProfilerMotorGap");
  double BeamProfilerAngleSupportX  = description.constant<double>("BeamProfilerAngleSupportX");
  double BeamProfilerAngleSupportD  = description.constant<double>("BeamProfilerAngleSupportD");

  Box    solidV408(BeamProfilerV408X/2.0, BeamProfilerV408Y/2.0, BeamProfilerV408Z/2.0);
  Volume logicV408("logicV408", solidV408, bpV408Material);

  Box    solidQ545(BeamProfilerQ545X/2.0, BeamProfilerQ545Y/2.0, BeamProfilerQ545Z/2.0);
  Volume logicQ545("logicQ545", solidQ545, bpQ545Material);

  Assembly motorAssembly("MotorAssembly");

  // Q545 motor
  motorAssembly.placeVolume(logicQ545, Position(0.0, 0.0, 0.5*BeamProfilerQ545Z));

  // V408 motor
  double mdx = 0.5*(BeamProfilerV408X - BeamProfilerQ545X);
  double mdy = BeamProfilerMotorGap + 0.5*(BeamProfilerV408Y + BeamProfilerQ545Y);
  double mdz = 0.5*BeamProfilerQ545Z;
  motorAssembly.placeVolume(logicV408, Position(-mdx, -mdy, mdz));

  // Angular support bracket
  double elangledx = BeamProfilerAngleSupportX;
  double elangledy = BeamProfilerMotorGap + BeamProfilerQ545Y;
  double elanglez  = 0.5*(BeamProfilerV408Z - BeamProfilerQ545Z);
  double eladr     = BeamProfilerAngleSupportD;

  Box solidBPAngSup1(elangledx/2.0, elangledy/2.0, elanglez/2.0);
  Box solidBPAngSupCut1(elangledx, elangledy/2.0, elanglez/2.0);
  SubtractionSolid solidBPAngSup("solidBPAngSup", solidBPAngSup1, solidBPAngSupCut1,
      Position(0.0, eladr, eladr));
  Volume logicBPAngSup("logicBPAngSup", solidBPAngSup, bpaMaterial);

  double asdy = 0.5*(BeamProfilerQ545Y - elangledy);
  double asdz = BeamProfilerQ545Z + 0.5*elanglez;
  motorAssembly.placeVolume(logicBPAngSup, Position(0.0, asdy, asdz));

  return motorAssembly;
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxBeamProfiler(dd4hep::Detector& description,
                                    xml_h             e,
                                    dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxBeamProfiler bp;
  bp.Construct(description, sdet, e, sd);

  return sdet;
}

DECLARE_DETELEMENT(LxBeamProfiler, create_LxBeamProfiler)
