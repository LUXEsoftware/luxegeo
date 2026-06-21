//
/// \brief Implementation of the LxDetectorGammaCalo8 class (DD4hep version)
//

#include <cmath>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Readout.h"
#include "DD4hep/Segmentations.h"
#include "TGeoManager.h"

#include "LxGammaCalo.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxDetectorGammaCalo8::Construct(dd4hep::Detector&   description,
                                     dd4hep::DetElement& sdet,
                                     xml_h&              e,
                                     dd4hep::SensitiveDetector& sd)
{
  Volume fLogicWorld = description.worldVolume();

  xml_comp_t x_det(e);
  Assembly   envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
//   envPV.addPhysVolID("GammaMonitor", x_det.id());
  sdet.setPlacement(envPV);

  Material gammaCaloMaterial           = description.material(description.constant<std::string>("GammaCaloMaterial"));
  Material LGFoilMaterial              = description.material(description.constant<std::string>("GammaCaloFoilMaterial"));
  Material GammaMonitorContainerMaterial = description.material(description.constant<std::string>("GammaMonitorDetMaterial"));
  Material CalSupportMaterial          = description.material(description.constant<std::string>("CalSupportMaterial"));
  Material GMSupportMaterial           = description.material(description.constant<std::string>("GammaMonitorSupportMaterial"));
  Material GMHolderMaterial            = description.material(description.constant<std::string>("GammaMonitorHolderMaterial"));

  double LeadGlassX             = description.constant<double>("LeadGlassX");
  double LeadGlassY             = description.constant<double>("LeadGlassY");
  double LeadGlassZ             = description.constant<double>("LeadGlassZ");
  double GammaCaloFoilThickness = description.constant<double>("GammaCaloFoilThickness");
  double GammaCaloZpos          = description.constant<double>("GammaCaloZpos");
  double GammaCaloZ             = description.constant<double>("GammaCaloZ");
  double CalSupportR            = description.constant<double>("CalSupportR");
  double CalSupporthickness     = description.constant<double>("CalSupporthickness");
  double FloorSurfaceYpos       = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  int    NLGblocks = 8;
  double lgthkns   = GammaCaloFoilThickness;
  double lgdz      = LeadGlassZ + lgthkns;
  double lgdy      = LeadGlassY + 2.0*lgthkns;
  double lgdx      = LeadGlassX + 2.0*lgthkns;

  // -----------------------------------------------------------------------
  // Lead glass block (volume A — sensitive)
  // -----------------------------------------------------------------------
  Box    solidLeadGlass(LeadGlassX/2.0, LeadGlassY/2.0, LeadGlassZ/2.0);
  Volume logicLeadGlass("logicLeadGlass", solidLeadGlass, gammaCaloMaterial);
//   logicLeadGlass.setSensitiveDetector(description.sensitiveDetector("GammaMonitor"));
  logicLeadGlass.setSensitiveDetector(sd);

  // -----------------------------------------------------------------------
  // Foil wrapping (volume B — module)
  // -----------------------------------------------------------------------
  Box    solidLGFoil(lgdx/2.0, lgdy/2.0, lgdz/2.0);
  Volume logicLGFoil("logicLGFoil", solidLGFoil, LGFoilMaterial);

  // Place A into B once — outside the loop
  PlacedVolume pvLG = logicLGFoil.placeVolume(logicLeadGlass,
      Position(0.0, 0.0, -0.5*(lgdz - LeadGlassZ)));
//   pvLG.addPhysVolID("sensor", 0);

  // Sensor DetElement template — will be cloned per module
  DetElement sensorDE_template("sensor", 0);
  sensorDE_template.setPlacement(pvLG);

  // -----------------------------------------------------------------------
  // Container box
  // -----------------------------------------------------------------------
  double gmcontx = CalSupportR + lgdy;
  double gmconty = gmcontx + 20.0*dd4hep::cm;

  Box    solidGammaMonitorContainer(gmcontx, gmconty, lgdz/2.0);
  Volume logicGammaMonitorContainer("logicGammaMonitorContainer",
                                     solidGammaMonitorContainer, GammaMonitorContainerMaterial);

  // Cylindrical support ring
  Tube   solidCalSupport(CalSupportR - CalSupporthickness, CalSupportR, lgdz/2.0, 0.0, 2.0*M_PI);
  Volume logicCalSupport("logicCalSupport", solidCalSupport, CalSupportMaterial);
  PlacedVolume pvSupport = logicGammaMonitorContainer.placeVolume(logicCalSupport,
      Position(0.0, 0.0, 0.0));
  pvSupport.addPhysVolID("CalSupport", 0);
  if (OverlapTest) pvSupport.ptr()->CheckOverlaps();

  // DetElement for the container
  DetElement containerDE(sdet, "GammaMonitor", 0);

  // -----------------------------------------------------------------------
  // Place B (foil) N=8 times in a ring — loop
  // -----------------------------------------------------------------------
  double rlg = CalSupportR + lgdy/2.0;
  for (int i = 0; i < NLGblocks; ++i) {
    double LGphi = i * 2.0*M_PI / static_cast<double>(NLGblocks) - M_PI_2;

    // G4RotationMatrix(G4ThreeVector(0,0,1), LGphi): axis-angle around Z → RotationZYX(LGphi,0,0)
    PlacedVolume pvFoil = logicGammaMonitorContainer.placeVolume(logicLGFoil,
        Transform3D(RotationZYX(LGphi, 0.0, 0.0),
                    Position(rlg*std::cos(LGphi), rlg*std::sin(LGphi), 0.0)));
    pvFoil.addPhysVolID("module", i);
    if (OverlapTest) pvFoil.ptr()->CheckOverlaps();

    // Module DetElement
    DetElement moduleDE(containerDE, "module_" + std::to_string(i), i);
    moduleDE.setPlacement(pvFoil);

    // Sensor DetElement — cloned for each module, pointing to same pvLG
    DetElement sensorDE = sensorDE_template.clone("sensor", 0);
    moduleDE.add(sensorDE);
  }

  // -----------------------------------------------------------------------
  // Metal holders
  // -----------------------------------------------------------------------
  double gmholddx = 1.0*dd4hep::cm;
  double gmholdx  = 1.5 * lgdy;
  double gmholdy  = gmconty - std::sqrt(CalSupportR*CalSupportR - gmholdx*gmholdx);

  Box          solidGammaMonitorHolder0(gmholddx/2.0, gmholdy/2.0, GammaCaloZ/2.0 * 0.5);
  Tube         solidCalSupportCut(0.0, CalSupportR, GammaCaloZ/2.0 * 0.5, 0.0, 2.0*M_PI);
  EllipticalTube solidCalSupportCut1(0.3*gmholdy, 0.3*GammaCaloZ/2.0 * 0.5, gmholddx);

  SubtractionSolid solidGammaMonitorHolder1("solidGammaMonitorHolder1",
      solidGammaMonitorHolder0, solidCalSupportCut,
      Position(-gmholdx, gmconty - gmholdy/2.0, 0.0));

  // G4RotationMatrix(G4ThreeVector(0,1,0), pi/2): axis-angle around Y → RotationZYX(0, pi/2, 0)
  SubtractionSolid solidGammaMonitorHolder("solidGammaMonitorHolder",
      solidGammaMonitorHolder1, solidCalSupportCut1,
      Transform3D(RotationZYX(0.0, M_PI/2.0, 0.0), Position(0.0, 0.0, 0.0)));

  Volume logicGammaMonitorHolder("logicGammaMonitorHolder", solidGammaMonitorHolder, GMHolderMaterial);

  // Commented: alternative placements inside container
  // logicGammaMonitorContainer.placeVolume(logicGammaMonitorHolder,
  //     Transform3D(RotationZYX(0, M_PI, 0),
  //                 Position(-gmholdx, FloorSurfaceYpos + gmsupy + gmholdy/2.0, GammaCaloZpos)));
  // logicGammaMonitorContainer.placeVolume(logicGammaMonitorHolder,
  //     Position(gmholdx, FloorSurfaceYpos + gmsupy + gmholdy/2.0, GammaCaloZpos));

  PlacedVolume pvHolder0 = logicGammaMonitorContainer.placeVolume(logicGammaMonitorHolder,
      // G4RotationMatrix(G4ThreeVector(0,1,0), pi): axis-angle around Y → RotationZYX(0, pi, 0)
      Transform3D(RotationZYX(0.0, M_PI, 0.0),
                  Position(-gmholdx, -(gmconty - gmholdy/2.0), 0.0)));
  pvHolder0.addPhysVolID("GammaMonitorHolder", 0);
  if (OverlapTest) pvHolder0.ptr()->CheckOverlaps();

  PlacedVolume pvHolder1 = logicGammaMonitorContainer.placeVolume(logicGammaMonitorHolder,
      Position(gmholdx, -(gmconty - gmholdy/2.0), 0.0));
  pvHolder1.addPhysVolID("GammaMonitorHolder", 1);
  if (OverlapTest) pvHolder1.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Place container into world
  // -----------------------------------------------------------------------
  PlacedVolume pvCont = envelope.placeVolume(logicGammaMonitorContainer,
      Position(0.0, 0.0, GammaCaloZpos));
  pvCont.addPhysVolID("system", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();
  containerDE.setPlacement(pvCont);

  // -----------------------------------------------------------------------
  // Support pedestal
  // -----------------------------------------------------------------------
  double gmsupx = 3.0 * rlg;
  double gmsupy = -FloorSurfaceYpos - gmconty;
  double gmsupz = 0.5 * GammaCaloZ;

  Box    solidGammaMonitorSupport(gmsupx/2.0, gmsupy/2.0, gmsupz/2.0);
  Volume logicGammaMonitorSupport("logicGammaMonitorSupport", solidGammaMonitorSupport, GMSupportMaterial);
  PlacedVolume pvSup = envelope.placeVolume(logicGammaMonitorSupport,
      Position(0.0, FloorSurfaceYpos + gmsupy/2.0, GammaCaloZpos));
  pvSup.addPhysVolID("GammaMonitorSupport", 0);
  if (OverlapTest) pvSup.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Active sub-detectors
  // -----------------------------------------------------------------------
//  ConstructGammaBeamDump(description, envelope);         // commented — inactive
  ConstructGammaBeamDumpMagnetized(description, envelope);
//  ConstructGammaBeamDumpShielding(description, envelope); // commented — inactive
  ConstructGammaBeamDumpShieldingTight(description, envelope);
}


// ---------------------------------------------------------------------------
void LxDetectorGammaCalo8::ConstructGammaBeamDump(dd4hep::Detector& description,
                                                   dd4hep::Volume&   motherVol)
{
  Material beamDumpMaterial    = description.material(description.constant<std::string>("GammaBeamDumpMaterial"));
  Material dumpPlugMaterial    = description.material(description.constant<std::string>("GammaBeamDumpPlugMaterial"));
  Material environmentMaterial = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material beamDumpWrapMaterial= description.material(description.constant<std::string>("GammaBeamDumpWrapMaterial"));

  double GammaBeamDumpR       = description.constant<double>("GammaBeamDumpR");
  double GammaBeamDumpZ       = description.constant<double>("GammaBeamDumpZ");
  double GammaBeamDumpWrapR   = description.constant<double>("GammaBeamDumpWrapR");
  double GammaBeamDumpPlugRin = description.constant<double>("GammaBeamDumpPlugRin");
  double GammaBeamDumpPlugZ   = description.constant<double>("GammaBeamDumpPlugZ");
  double GammaBeamDumpZpos    = description.constant<double>("GammaBeamDumpZpos");
  bool   OverlapTest          = (description.constant<int>("OverlapTest") != 0);

  double gammabdz = GammaBeamDumpZ + GammaBeamDumpPlugZ;
  double gbdposz  = (gammabdz - GammaBeamDumpZ) / 2.0;

  Tube   solidGammaBeamDumpContainer(0.0, GammaBeamDumpWrapR, gammabdz/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDumpContainer("logicGammaBeamDumpContainer",
                                      solidGammaBeamDumpContainer, environmentMaterial);

  Tube   solidBeamDump(0.0, GammaBeamDumpR, GammaBeamDumpZ/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDump("logicGammaBeamDump", solidBeamDump, beamDumpMaterial);

  Tube   solidBeamDumpWrap(GammaBeamDumpR, GammaBeamDumpWrapR, GammaBeamDumpZ/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDumpWrap("logicGammaBeamDumpWrap", solidBeamDumpWrap, beamDumpWrapMaterial);

  Tube   solidBeamDumpPlug(GammaBeamDumpPlugRin, GammaBeamDumpR, GammaBeamDumpPlugZ/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamDumpPlug("logicBeamDumpPlug", solidBeamDumpPlug, dumpPlugMaterial);

  Tube   solidBeamDumpPlugWrap(GammaBeamDumpR, GammaBeamDumpWrapR, GammaBeamDumpPlugZ/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamDumpPlugWrap("logicBeamDumpPlugWrap", solidBeamDumpPlugWrap, beamDumpWrapMaterial);

  PlacedVolume pvDump = logicGammaBeamDumpContainer.placeVolume(logicGammaBeamDump,
      Position(0.0, 0.0, gbdposz));
  pvDump.addPhysVolID("GammaBeamDump", 0);
  if (OverlapTest) pvDump.ptr()->CheckOverlaps();

  PlacedVolume pvWrap = logicGammaBeamDumpContainer.placeVolume(logicGammaBeamDumpWrap,
      Position(0.0, 0.0, gbdposz));
  pvWrap.addPhysVolID("GammaBeamDumpWrap", 0);
  if (OverlapTest) pvWrap.ptr()->CheckOverlaps();

  double plugzpos = gbdposz - (GammaBeamDumpPlugZ + GammaBeamDumpZ)/2.0;
  PlacedVolume pvPlug = logicGammaBeamDumpContainer.placeVolume(logicBeamDumpPlug,
      Position(0.0, 0.0, plugzpos));
  pvPlug.addPhysVolID("GammaBeamDumpPlug", 0);
  if (OverlapTest) pvPlug.ptr()->CheckOverlaps();

  PlacedVolume pvPlugWrap = logicGammaBeamDumpContainer.placeVolume(logicBeamDumpPlugWrap,
      Position(0.0, 0.0, plugzpos));
  pvPlugWrap.addPhysVolID("GammaBeamDumpPlugWrap", 0);
  if (OverlapTest) pvPlugWrap.ptr()->CheckOverlaps();

  // Commented: dump support (see original for geometry)
  // ...

  PlacedVolume pvCont = motherVol.placeVolume(logicGammaBeamDumpContainer,
      Position(0.0, 0.0, GammaBeamDumpZpos));
  pvCont.addPhysVolID("GammaBeamDumpAssembly", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxDetectorGammaCalo8::ConstructGammaBeamDumpMagnetized(dd4hep::Detector& description,
                                                             dd4hep::Volume&   motherVol)
{
  Material beamDumpMaterial       = description.material(description.constant<std::string>("GammaBeamDumpMaterial"));
  Material dumpPlugMaterial       = description.material(description.constant<std::string>("GammaBeamDumpPlugMaterial"));
  Material environmentMaterial    = description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material beamDumpWrapMaterial   = description.material(description.constant<std::string>("GammaBeamDumpWrapMaterial"));
  Material beamDumpMagnetMaterial = description.material(description.constant<std::string>("GammaBeamDumpMagnetMaterial"));

  double GammaBeamDumpR           = description.constant<double>("GammaBeamDumpR");
  double GammaBeamDumpZ           = description.constant<double>("GammaBeamDumpZ");
  double GammaBeamDumpWrapR       = description.constant<double>("GammaBeamDumpWrapR");
  double GammaBeamDumpPlugRin     = description.constant<double>("GammaBeamDumpPlugRin");
  double GammaBeamDumpPlugZ       = description.constant<double>("GammaBeamDumpPlugZ");
  double GammaBeamDumpZpos        = description.constant<double>("GammaBeamDumpZpos");
  double GammaBeamDumpBSMRAirGapZ = description.constant<double>("GammaBeamDumpBSMRAirGapZ");
  double GammaBeamDumpBSMRAirR    = description.constant<double>("GammaBeamDumpBSMRAirR");
  double GammaBeamDumpMagnetGapZ  = description.constant<double>("GammaBeamDumpMagnetGapZ");
  double GammaBeamDumpMagnetThick = description.constant<double>("GammaBeamDumpMagnetThick");
  bool   OverlapTest              = (description.constant<int>("OverlapTest") != 0);

  double airZ    = GammaBeamDumpZ - GammaBeamDumpBSMRAirGapZ;
  double magnetZ = GammaBeamDumpZ - GammaBeamDumpMagnetGapZ;
  double magnetR = GammaBeamDumpR + GammaBeamDumpMagnetThick;

  double gammabdz = GammaBeamDumpZ + GammaBeamDumpPlugZ;
  double gbdposz  = (gammabdz - GammaBeamDumpZ) / 2.0;

  Tube   solidGammaBeamDumpContainer(0.0, GammaBeamDumpWrapR, gammabdz/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDumpContainer("logicGammaBeamDumpContainer",
                                      solidGammaBeamDumpContainer, environmentMaterial);

  // Dump body with BSM air gap
  Tube   solidGammaBeamDump(0.0, GammaBeamDumpR, GammaBeamDumpZ/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDump("logicGammaBeamDump", solidGammaBeamDump, beamDumpMaterial);

  Tube   solidGammaBeamDumpBSMAir(0.0, GammaBeamDumpBSMRAirR, airZ/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDumpBSMAir("logicGammaBeamDumpBSMAir", solidGammaBeamDumpBSMAir, environmentMaterial);

  double bsmairzpos = 0.5*(GammaBeamDumpZ - airZ);
  PlacedVolume pvAir = logicGammaBeamDump.placeVolume(logicGammaBeamDumpBSMAir,
      Position(0.0, 0.0, bsmairzpos));
  pvAir.addPhysVolID("GammaBeamDumpBSMAir", 0);
  if (OverlapTest) pvAir.ptr()->CheckOverlaps();

  // Wrap + magnet ring
  Tube   solidBeamDumpWrap(GammaBeamDumpR, GammaBeamDumpWrapR, GammaBeamDumpZ/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDumpWrap("logicGammaBeamDumpWrap", solidBeamDumpWrap, beamDumpWrapMaterial);

  Tube   solidGammaBeamDumpMagnet(GammaBeamDumpR, magnetR, magnetZ/2.0, 0.0, 2.0*M_PI);
  Volume logicGammaBeamDumpMagnet("logicGammaBeamDumpMagnet", solidGammaBeamDumpMagnet, beamDumpMagnetMaterial);

  double bsmmagzpos = 0.5*(GammaBeamDumpZ - magnetZ);
  PlacedVolume pvMagnet = logicGammaBeamDumpWrap.placeVolume(logicGammaBeamDumpMagnet,
      Position(0.0, 0.0, bsmmagzpos));
  pvMagnet.addPhysVolID("GammaBeamDumpMagnet", 0);
  if (OverlapTest) pvMagnet.ptr()->CheckOverlaps();

  // Plug
  Tube   solidBeamDumpPlug(GammaBeamDumpPlugRin, GammaBeamDumpR, GammaBeamDumpPlugZ/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamDumpPlug("logicBeamDumpPlug", solidBeamDumpPlug, dumpPlugMaterial);

  Tube   solidBeamDumpPlugWrap(GammaBeamDumpR, GammaBeamDumpWrapR, GammaBeamDumpPlugZ/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamDumpPlugWrap("logicBeamDumpPlugWrap", solidBeamDumpPlugWrap, beamDumpWrapMaterial);

  double plugzpos = gbdposz - (GammaBeamDumpPlugZ + GammaBeamDumpZ)/2.0;

  PlacedVolume pvDump = logicGammaBeamDumpContainer.placeVolume(logicGammaBeamDump,
      Position(0.0, 0.0, gbdposz));
  pvDump.addPhysVolID("GammaBeamDump", 0);
  if (OverlapTest) pvDump.ptr()->CheckOverlaps();

  PlacedVolume pvWrap = logicGammaBeamDumpContainer.placeVolume(logicGammaBeamDumpWrap,
      Position(0.0, 0.0, gbdposz));
  pvWrap.addPhysVolID("GammaBeamDumpWrap", 0);
  if (OverlapTest) pvWrap.ptr()->CheckOverlaps();

  PlacedVolume pvPlug = logicGammaBeamDumpContainer.placeVolume(logicBeamDumpPlug,
      Position(0.0, 0.0, plugzpos));
  pvPlug.addPhysVolID("GammaBeamDumpPlug", 0);
  if (OverlapTest) pvPlug.ptr()->CheckOverlaps();

  PlacedVolume pvPlugWrap = logicGammaBeamDumpContainer.placeVolume(logicBeamDumpPlugWrap,
      Position(0.0, 0.0, plugzpos));
  pvPlugWrap.addPhysVolID("GammaBeamDumpPlugWrap", 0);
  if (OverlapTest) pvPlugWrap.ptr()->CheckOverlaps();

  PlacedVolume pvCont = motherVol.placeVolume(logicGammaBeamDumpContainer,
      Position(0.0, 0.0, GammaBeamDumpZpos));
  pvCont.addPhysVolID("GammaBeamDumpAssembly", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();

  // SetGammaBeamDumpBField — omitted (Geant4-specific field infrastructure)
}


// ---------------------------------------------------------------------------
void LxDetectorGammaCalo8::ConstructGammaBeamDumpShielding(dd4hep::Detector& description,
                                                            dd4hep::Volume&   motherVol)
{
  Material shieldingMaterial     = description.material(description.constant<std::string>("GammaBeamDumpShieldingMaterial"));
  Material shieldingPlugMaterial = description.material(description.constant<std::string>("GammaBeamDumpShieldingPlugMaterial"));

  double GammaBeamDumpWrapR            = description.constant<double>("GammaBeamDumpWrapR");
  double GammaBeamDumpShieldingX       = description.constant<double>("GammaBeamDumpShieldingX");
  double GammaBeamDumpShieldingY       = description.constant<double>("GammaBeamDumpShieldingY");
  double GammaBeamDumpShieldingZ       = description.constant<double>("GammaBeamDumpShieldingZ");
  double GammaBeamDumpShieldingGapX    = description.constant<double>("GammaBeamDumpShieldingGapX");
  double GammaBeamDumpShieldingGapY    = description.constant<double>("GammaBeamDumpShieldingGapY");
  double GammaBeamDumpShieldingPlugZ   = description.constant<double>("GammaBeamDumpShieldingPlugZ");
  double GammaBeamDumpZpos             = description.constant<double>("GammaBeamDumpZpos");
  double GammaBeamDumpZ                = description.constant<double>("GammaBeamDumpZ");
  double GammaBeamDumpPlugZ            = description.constant<double>("GammaBeamDumpPlugZ");
  double FloorSurfaceYpos              = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest                   = (description.constant<int>("OverlapTest") != 0);

  if (GammaBeamDumpWrapR > GammaBeamDumpShieldingGapX/2.0 ||
      GammaBeamDumpWrapR > GammaBeamDumpShieldingGapY/2.0) {
    throw std::runtime_error("LxDetectorGammaCalo8::ConstructGammaBeamDumpShielding: "
                             "gamma beam dump is bigger than the space in the shielding!");
  }

  Box solidGBDShielding0(GammaBeamDumpShieldingX/2.0,
                         GammaBeamDumpShieldingY/2.0,
                         GammaBeamDumpShieldingZ/2.0);
  Box solidGBDShieldingCut(GammaBeamDumpShieldingGapX/2.0,
                            GammaBeamDumpShieldingGapY/2.0,
                            GammaBeamDumpShieldingZ);

  SubtractionSolid solidGBDShielding("solidGBDShielding", solidGBDShielding0, solidGBDShieldingCut,
      Position(0.0, -GammaBeamDumpShieldingY/2.0 - FloorSurfaceYpos, 0.0));
  Volume logicGBDShielding("logicGBDShielding", solidGBDShielding, shieldingMaterial);

  Box    solidGBDShieldingPlug(GammaBeamDumpShieldingGapX/2.0,
                               GammaBeamDumpShieldingGapY/2.0,
                               GammaBeamDumpShieldingPlugZ/2.0);
  Volume logicGBDShieldingPlug("logicGBDShieldingPlug", solidGBDShieldingPlug, shieldingPlugMaterial);

  double gbdszpos  = GammaBeamDumpZpos + 0.5*(GammaBeamDumpShieldingZ
                     - GammaBeamDumpZ - GammaBeamDumpPlugZ);
  double splugzpos = gbdszpos + 0.5*(GammaBeamDumpShieldingPlugZ - GammaBeamDumpShieldingZ)
                     + GammaBeamDumpZ + GammaBeamDumpPlugZ;

  PlacedVolume pvShield = motherVol.placeVolume(logicGBDShielding,
      Position(0.0, FloorSurfaceYpos + GammaBeamDumpShieldingY/2.0, gbdszpos));
  pvShield.addPhysVolID("GammaBeamDumpShielding", 0);
  if (OverlapTest) pvShield.ptr()->CheckOverlaps();

  PlacedVolume pvPlug = motherVol.placeVolume(logicGBDShieldingPlug,
      Position(0.0, 0.0, splugzpos));
  pvPlug.addPhysVolID("GammaBeamDumpShieldingPlug", 0);
  if (OverlapTest) pvPlug.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxDetectorGammaCalo8::ConstructGammaBeamDumpShieldingTight(dd4hep::Detector& description,
                                                                  dd4hep::Volume&   motherVol)
{
  Material shieldingMaterial = description.material(
      description.constant<std::string>("GammaBeamDumpShieldingMaterial"));

  double GammaBeamDumpShieldingX  = description.constant<double>("GammaBeamDumpShieldingX");
  double GammaBeamDumpShieldingY  = description.constant<double>("GammaBeamDumpShieldingY");
  double GammaBeamDumpShieldingZ  = description.constant<double>("GammaBeamDumpShieldingZ");
  double GammaBeamDumpWrapR       = description.constant<double>("GammaBeamDumpWrapR");
  double GammaBeamDumpZpos        = description.constant<double>("GammaBeamDumpZpos");
  double GammaBeamDumpZ           = description.constant<double>("GammaBeamDumpZ");
  double GammaBeamDumpPlugZ       = description.constant<double>("GammaBeamDumpPlugZ");
  double FloorSurfaceYpos         = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest              = (description.constant<int>("OverlapTest") != 0);

  Box  solidGBDShielding0(GammaBeamDumpShieldingX/2.0,
                          GammaBeamDumpShieldingY/2.0,
                          GammaBeamDumpShieldingZ/2.0);
  Tube solidGBDShieldingCut(0.0, GammaBeamDumpWrapR, GammaBeamDumpShieldingZ, 0.0, 2.0*M_PI);

  SubtractionSolid solidGBDShielding("solidGBDShielding", solidGBDShielding0, solidGBDShieldingCut,
      Position(0.0, -GammaBeamDumpShieldingY/2.0 - FloorSurfaceYpos, 0.0));
  Volume logicGBDShielding("logicGBDShielding", solidGBDShielding, shieldingMaterial);

  double gbdszpos = GammaBeamDumpZpos + 0.5*(GammaBeamDumpShieldingZ
                    - GammaBeamDumpZ - GammaBeamDumpPlugZ);

  PlacedVolume pvShield = motherVol.placeVolume(logicGBDShielding,
      Position(0.0, FloorSurfaceYpos + GammaBeamDumpShieldingY/2.0, gbdszpos));
  pvShield.addPhysVolID("GammaBeamDumpShielding", 0);
  if (OverlapTest) pvShield.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxDetectorGammaCalo8(dd4hep::Detector& description,
                                          xml_h             e,
                                          dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  sd.setType("calorimeter");
  LxDetectorGammaCalo8 det;
  det.Construct(description, sdet, e, sd);

  return sdet;
}

DECLARE_DETELEMENT(LxDetectorGammaCalo8, create_LxDetectorGammaCalo8)
