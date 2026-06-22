//
/// \brief Implementation of the LxDetectorComptonFluka class (DD4hep version)
//

#include <cmath>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxAux.h"
#include "LxGammaLanex.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxDetectorComptonFluka::Construct(dd4hep::Detector&          description,
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

  Material beamPipeMaterial    = description.material(description.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial      = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));
  Material CollimatorMaterial  = description.material(description.constant<std::string>("CollimatorMaterial"));
  Material LysoCalMaterial     = description.material(description.constant<std::string>("ComptonLysoMaterial"));
  Material lysoPedestaMaterial = description.material(description.constant<std::string>("LYSOPedestaMaterial"));

  double CollimatorRin          = description.constant<double>("CollimatorRin");
  double CollimatorRout         = description.constant<double>("CollimatorRout");
  double CollimatorZ            = description.constant<double>("CollimatorZ");
  double CollimatorZpos         = description.constant<double>("CollimatorZpos");
  double ComptonLysoX           = description.constant<double>("ComptonLysoX");
  double ComptonLysoY           = description.constant<double>("ComptonLysoY");
  double ComptonLysoZ           = description.constant<double>("ComptonLysoZ");
  double ComptonLysoXpos        = description.constant<double>("ComptonLysoXpos");
  double ComptonLysoZpos        = description.constant<double>("ComptonLysoZpos");
  double BPipeRLG               = description.constant<double>("BPipeRLG");
  double BPipeR                 = description.constant<double>("BPipeR");
  double BPipeThickness         = description.constant<double>("BPipeThickness");
  double HICSDumpFrontXPos      = description.constant<double>("HICSDumpFrontXPos");
  double HICSDumpFrontZPos      = description.constant<double>("HICSDumpFrontZPos");
  double HICSDumpZ              = description.constant<double>("HICSDumpZ");
  double IPMagnetZpos           = description.constant<double>("IPMagnetZpos");
  double GMagnetZpos            = description.constant<double>("GMagnetZpos");
  double FlashMFieldLength      = description.constant<double>("FlashMFieldLength");
  double GTargetZpos            = description.constant<double>("GTargetZpos");
  double GTargetContainerZ      = description.constant<double>("GTargetContainerZ");
  double ComptonElectronBeamtoStageY = description.constant<double>("ComptonElectronBeamtoStageY");
  double OPPPDetTopPlateY       = description.constant<double>("OPPPDetTopPlateY");
  double LYSODetSupportBallH    = description.constant<double>("LYSODetSupportBallH");
  double LYSODetTopPlateX       = description.constant<double>("LYSODetTopPlateX");
  double LYSODetTopPlateZ       = description.constant<double>("LYSODetTopPlateZ");
  double LYSODetBottomSupportZ  = description.constant<double>("LYSODetBottomSupportZ");
  double LYSOBasePlateX         = description.constant<double>("LYSOBasePlateX");
  double LYSOBasePlateZ         = description.constant<double>("LYSOBasePlateZ");
  double FloorSurfaceYpos       = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  double fHICSDumpAngle = std::atan2(HICSDumpFrontXPos, HICSDumpFrontZPos - IPMagnetZpos);

  // -----------------------------------------------------------------------
  // Collimator (2nd only — 1st is commented out in original)
  // -----------------------------------------------------------------------
  Tube   solidCollimator(CollimatorRin, CollimatorRout, CollimatorZ/2.0, 0.0, 2.0*M_PI);
  Volume logicCollimator("logicCollimator", solidCollimator, CollimatorMaterial);

  // Commented: 1st collimator placement
  // envelope.placeVolume(logicCollimator, Position(0.0, 0.0, CollimatorZpos)).addPhysVolID("Collimator", 0);

  PlacedVolume pvColl = envelope.placeVolume(logicCollimator,
      Position(0.0, 0.0, CollimatorZpos + CollimatorZ + 50.0*dd4hep::cm));
  pvColl.addPhysVolID("Collimator", 1);
  if (OverlapTest) pvColl.ptr()->CheckOverlaps();

  Tube   solidCollimatorVac(0.0, CollimatorRin, CollimatorZ/2.0, 0.0, 2.0*M_PI);
  Volume logicCollimatorVac("logicCollimatorVac", solidCollimatorVac, vacuumMaterial);

  // Commented: 1st collimator vacuum
  // envelope.placeVolume(logicCollimatorVac, Position(0.0, 0.0, CollimatorZpos)).addPhysVolID("CollimatorVac", 0);

  PlacedVolume pvCollVac = envelope.placeVolume(logicCollimatorVac,
      Position(0.0, 0.0, CollimatorZpos + CollimatorZ + 50.0*dd4hep::cm));
  pvCollVac.addPhysVolID("CollimatorVac", 1);
  if (OverlapTest) pvCollVac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // LYSO calorimeters (sensitive, volume A — placed twice)
  // -----------------------------------------------------------------------
  Box    solidComptonLyso(ComptonLysoX/2.0, ComptonLysoY/2.0, ComptonLysoZ/2.0);
  Volume logicComptonLysoCal("logicComptonLysoCal", solidComptonLyso, LysoCalMaterial);
  logicComptonLysoCal.setSensitiveDetector(sd);

  // Place logicComptonLysoCal twice directly — no intermediate wrapper volume,
  // each placement gets its own module DetElement
  PlacedVolume pvLyso0 = envelope.placeVolume(logicComptonLysoCal,
      Position(ComptonLysoXpos, 0.0, ComptonLysoZpos));
  pvLyso0.addPhysVolID("module", 0);
  if (OverlapTest) pvLyso0.ptr()->CheckOverlaps();

  DetElement lysoDE0(sdet, "LysoCal_0", 0);
  lysoDE0.setPlacement(pvLyso0);

  PlacedVolume pvLyso1 = envelope.placeVolume(logicComptonLysoCal,
      Position(-ComptonLysoXpos, 0.0, ComptonLysoZpos));
  pvLyso1.addPhysVolID("module", 1);
  if (OverlapTest) pvLyso1.ptr()->CheckOverlaps();

  DetElement lysoDE1(sdet, "LysoCal_1", 1);
  lysoDE1.setPlacement(pvLyso1);

  // -----------------------------------------------------------------------
  // Beam pipe section: target → 1st collimator (CutTubs)
  // -----------------------------------------------------------------------
  double frontZpos    = HICSDumpFrontZPos + HICSDumpFrontXPos * std::tan(fHICSDumpAngle);
  double lpipe_gd     = frontZpos - GTargetZpos - GTargetContainerZ/2.0;
  double lpipe_pos_1st= frontZpos - lpipe_gd/2.0;

  // Commented: straight tube version
  // Tube solidBeamPipeGammaT1stC(BPipeRLG-BPipeThickness, BPipeRLG, lpipe_gd/2.0, 0.0, 2.0*M_PI);

  // Active: CutTubs — low normal (0,0,-1), high normal (-sin,0,cos)
  CutTube solidBeamPipeGammaT1stC(BPipeRLG - BPipeThickness, BPipeRLG, lpipe_gd/2.0, 0.0, 2.0*M_PI,
      0.0, 0.0, -1.0,
      -std::sin(fHICSDumpAngle), 0.0, std::cos(fHICSDumpAngle));
  Volume logicBeamPipeGammaT1stC("logicBeamPipeGammaT1stC", solidBeamPipeGammaT1stC, beamPipeMaterial);
  PlacedVolume pvP1 = envelope.placeVolume(logicBeamPipeGammaT1stC,
      Position(0.0, 0.0, lpipe_pos_1st));
  pvP1.addPhysVolID("BeamPipeGammaT1stC", 0);
  if (OverlapTest) pvP1.ptr()->CheckOverlaps();

  CutTube solidBeamPipeGammaT1stCVac(0.0, BPipeR - BPipeThickness, lpipe_gd/2.0, 0.0, 2.0*M_PI,
      0.0, 0.0, -1.0,
      -std::sin(fHICSDumpAngle), 0.0, std::cos(fHICSDumpAngle));
  Volume logicBeamPipeGammaT1stCVac("logicBeamPipeGammaT1stCVac", solidBeamPipeGammaT1stCVac, vacuumMaterial);
  PlacedVolume pvP1Vac = envelope.placeVolume(logicBeamPipeGammaT1stCVac,
      Position(0.0, 0.0, lpipe_pos_1st));
  pvP1Vac.addPhysVolID("BeamPipeGammaT1stCVac", 0);
  if (OverlapTest) pvP1Vac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Beam pipe section: HICS dump → 2nd collimator (CutTubs)
  // -----------------------------------------------------------------------
  double pos_2ndColl = CollimatorZpos + CollimatorZ + 50.0*dd4hep::cm;
  double rearZpos    = HICSDumpFrontZPos + HICSDumpZ / std::cos(fHICSDumpAngle)
                       + HICSDumpFrontXPos * std::tan(fHICSDumpAngle);
  double lpipe_gd2nd    = pos_2ndColl - rearZpos - CollimatorZ/2.0;
  double lpipe_pos_2nd  = pos_2ndColl - CollimatorZ/2.0 - lpipe_gd2nd/2.0;

  // Commented: straight tube version
  // Tube solidBeamPipeGamma1stC2ndC(BPipeRLG-BPipeThickness, BPipeRLG, lpipe_gd2nd/2.0, 0.0, 2.0*M_PI);

  // Active: CutTubs — low normal (sin,0,-cos), high normal (0,0,1)
  CutTube solidBeamPipeGamma1stC2ndC(BPipeRLG - BPipeThickness, BPipeRLG, lpipe_gd2nd/2.0, 0.0, 2.0*M_PI,
       std::sin(fHICSDumpAngle), 0.0, -std::cos(fHICSDumpAngle),
       0.0, 0.0, 1.0);
  Volume logicBeamPipeGamma1stC2ndC("logicBeamPipeGamma1stC2ndC", solidBeamPipeGamma1stC2ndC, beamPipeMaterial);
  PlacedVolume pvP2 = envelope.placeVolume(logicBeamPipeGamma1stC2ndC,
      Position(0.0, 0.0, lpipe_pos_2nd));
  pvP2.addPhysVolID("BeamPipeGamma1stC2ndC", 0);
  if (OverlapTest) pvP2.ptr()->CheckOverlaps();

  CutTube solidBeamPipeGamma1stC2ndCVac(0.0, BPipeR - BPipeThickness, lpipe_gd2nd/2.0, 0.0, 2.0*M_PI,
       std::sin(fHICSDumpAngle), 0.0, -std::cos(fHICSDumpAngle),
       0.0, 0.0, 1.0);
  Volume logicBeamPipeGamma1stC2ndCVac("logicBeamPipeGamma1stC2ndCVac", solidBeamPipeGamma1stC2ndCVac, vacuumMaterial);
  PlacedVolume pvP2Vac = envelope.placeVolume(logicBeamPipeGamma1stC2ndCVac,
      Position(0.0, 0.0, lpipe_pos_2nd));
  pvP2Vac.addPhysVolID("BeamPipeGamma1stC2ndCVac", 0);
  if (OverlapTest) pvP2Vac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Beam pipe section: 2nd collimator → gamma magnet field (round)
  // -----------------------------------------------------------------------
  double lpipe_gd3d   = GMagnetZpos - FlashMFieldLength/2.0 - lpipe_pos_2nd - CollimatorZ - lpipe_gd2nd/2.0;
  double lpipe_pos3d  = GMagnetZpos - FlashMFieldLength/2.0 - lpipe_gd3d/2.0;

  Tube   solidBeamPipeGamma2ndCMagF(BPipeRLG - BPipeThickness, BPipeRLG, lpipe_gd3d/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeGamma2ndCMagF("logicBeamPipeGamma2ndCMagF", solidBeamPipeGamma2ndCMagF, beamPipeMaterial);
  PlacedVolume pvP3 = envelope.placeVolume(logicBeamPipeGamma2ndCMagF,
      Position(0.0, 0.0, lpipe_pos3d));
  pvP3.addPhysVolID("BeamPipeGamma2ndCMagF", 0);
  if (OverlapTest) pvP3.ptr()->CheckOverlaps();

  Tube   solidBeamPipeGamma2ndCMagFVac(0.0, BPipeR - BPipeThickness, lpipe_gd3d/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeGamma2ndCMagFVac("logicBeamPipeGamma2ndCMagFVac", solidBeamPipeGamma2ndCMagFVac, vacuumMaterial);
  PlacedVolume pvP3Vac = envelope.placeVolume(logicBeamPipeGamma2ndCMagFVac,
      Position(0.0, 0.0, lpipe_pos3d));
  pvP3Vac.addPhysVolID("BeamPipeGamma2ndCMagFVac", 0);
  if (OverlapTest) pvP3Vac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Support assembly for LYSO detectors
  // -----------------------------------------------------------------------
  double supporty;
  Assembly lysoSupportAssembly = ConstructSupportAssembly(description, supporty);

  double chery   = ComptonElectronBeamtoStageY;
  double detxpos = (LYSODetTopPlateX - ComptonLysoX)/2.0 + ComptonLysoXpos;
  double supypos = -chery + OPPPDetTopPlateY/2.0 + LYSODetSupportBallH;
  double detzpos = ComptonLysoZpos - (ComptonLysoZ - LYSODetTopPlateZ)/2.0;

  // +X side
  Position trsupport(detxpos, supypos, detzpos);
  LxAux::AddAssemblyVolumes(envelope, lysoSupportAssembly, trsupport);

  // -X side
  trsupport = Position(-detxpos, supypos, detzpos);
  LxAux::AddAssemblyVolumes(envelope, lysoSupportAssembly, trsupport);

  // -----------------------------------------------------------------------
  // Small pedestal for LYSO detector
  // -----------------------------------------------------------------------
  double pedy = chery - ComptonLysoY/2.0 - LYSODetSupportBallH - OPPPDetTopPlateY;

  Box    solidLYSOPedestal(ComptonLysoX/2.0, pedy/2.0, ComptonLysoZ/2.0);
  Volume logicLYSOPedestal("logicLYSOPedestal", solidLYSOPedestal, LysoCalMaterial);

  PlacedVolume pvPed0 = envelope.placeVolume(logicLYSOPedestal,
      Position(ComptonLysoXpos, -(ComptonLysoY + pedy)/2.0, ComptonLysoZpos));
  pvPed0.addPhysVolID("LysoPedestal", 0);
  if (OverlapTest) pvPed0.ptr()->CheckOverlaps();

  PlacedVolume pvPed1 = envelope.placeVolume(logicLYSOPedestal,
      Position(-ComptonLysoXpos, -(ComptonLysoY + pedy)/2.0, ComptonLysoZpos));
  pvPed1.addPhysVolID("LysoPedestal", 1);
  if (OverlapTest) pvPed1.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Table and pedestal support
  // -----------------------------------------------------------------------
  double ypestal    = 0.5*dd4hep::m;
  double ylevel     = supporty + pedy + ComptonLysoY/2.0;
  double tblhight   = -ylevel - ypestal - FloorSurfaceYpos;
  double tablezpos  = ComptonLysoZpos - (ComptonLysoZ - LYSOBasePlateZ)/2.0;

  Assembly tablesupport = LxAux::BuildTable(description, "GammaSpectrTable",
                                            LYSOBasePlateX, tblhight, LYSOBasePlateZ, 4);
  LxAux::AddAssemblyVolumes(envelope, tablesupport, Position(0.0, -ylevel, tablezpos));

  Volume pedestal = LxAux::BuildPedestal(description, "GammaSpectr",
                                         LYSOBasePlateX, ypestal, LYSOBasePlateZ);
  PlacedVolume pvPedestal = envelope.placeVolume(pedestal,
      Position(0.0, 0.5*ypestal + FloorSurfaceYpos, tablezpos));
  if (OverlapTest) pvPedestal.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Shielding
  // -----------------------------------------------------------------------
  ConstructComptShielding(description, envelope);
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxDetectorComptonFluka::ConstructSupportAssembly(
    dd4hep::Detector& description,
    double&           sphight)
{
  Material lysoDetSupportMaterial = description.material(
      description.constant<std::string>("OPPPDetSupportMaterial"));

  double LYSODetBottomSupportX  = description.constant<double>("LYSODetBottomSupportX");
  double LYSODetBottomSupportZ  = description.constant<double>("LYSODetBottomSupportZ");
  double LYSODetTopPlateX       = description.constant<double>("LYSODetTopPlateX");
  double LYSODetTopPlateZ       = description.constant<double>("LYSODetTopPlateZ");
  double LYSODetSupportBallH    = description.constant<double>("LYSODetSupportBallH");
  double OPPPDetBottomSupportY  = description.constant<double>("OPPPDetBottomSupportY");
  double OPPPDetTopPlateY       = description.constant<double>("OPPPDetTopPlateY");

  Box    solidLYSOBottomSupport(LYSODetBottomSupportX/2.0, OPPPDetBottomSupportY/2.0, LYSODetBottomSupportZ/2.0);
  Volume logicLYSOBottomSupport("logicLYSOBottomSupport", solidLYSOBottomSupport, lysoDetSupportMaterial);

  Box    solidLYSOTopPlate(LYSODetTopPlateX/2.0, OPPPDetTopPlateY/2.0, LYSODetTopPlateZ/2.0);
  Volume logicLYSOTopPlate("logicLYSOTopPlate", solidLYSOTopPlate, lysoDetSupportMaterial);

  // G4Sphere(rmin, rmax, startPhi, dPhi, startTheta, dTheta)
  // → dd4hep::Sphere(rmin, rmax, startTheta, endTheta, startPhi, endPhi)
  Sphere solidLYSODetSupportBall(0.0, LYSODetSupportBallH/2.0, 0.0, M_PI, 0.0, 2.0*M_PI);
  Volume logicLYSODetSupportBall("logicLYSODetSupportBall", solidLYSODetSupportBall, lysoDetSupportMaterial);

  Assembly lysoSupportAssembly("LYSOSupportAssembly");

  // Top plate at origin
  lysoSupportAssembly.placeVolume(logicLYSOTopPlate, Position(0.0, 0.0, 0.0));

  // Bottom support
  double bsy = -(OPPPDetBottomSupportY + OPPPDetTopPlateY)/2.0 - LYSODetSupportBallH;
  double bsz = (LYSODetBottomSupportZ - LYSODetTopPlateZ)/2.0;
  lysoSupportAssembly.placeVolume(logicLYSOBottomSupport, Position(0.0, bsy, bsz));

  // Support balls (3 positions)
  double dbz = 0.5*LYSODetTopPlateZ - 1.05*LYSODetSupportBallH;
  double dbx = 0.5*LYSODetTopPlateX - 1.05*LYSODetSupportBallH;
  double dby = bsy + (OPPPDetBottomSupportY + LYSODetSupportBallH)/2.0;

  lysoSupportAssembly.placeVolume(logicLYSODetSupportBall, Position( dbx, dby,  dbz));
  lysoSupportAssembly.placeVolume(logicLYSODetSupportBall, Position(-dbx, dby,  dbz));
  lysoSupportAssembly.placeVolume(logicLYSODetSupportBall, Position( 0.0, dby, -dbz));

  // Hexapod
  double hexhight = 0.0;
  Assembly hexAssembly = LxAux::BuildHexapod(description, "LYSO", hexhight);
  Position hexpos(0.0, bsy - (OPPPDetBottomSupportY + hexhight)/2.0, bsz);
  LxAux::AddAssemblyVolumes(lysoSupportAssembly, hexAssembly, hexpos);

  sphight = (OPPPDetTopPlateY + hexhight)/2.0 - hexpos.Y();
  return lysoSupportAssembly;
}


// ---------------------------------------------------------------------------
void LxDetectorComptonFluka::ConstructComptShielding(dd4hep::Detector& description,
                                                     dd4hep::Volume&   motherVol)
{
  Material concreteMaterial      = description.material("ShieldingConcrete");
  Material shieldingMaterial     = description.material(description.constant<std::string>("ComptShieldingMaterial"));
  Material shieldingMaterial2    = description.material(description.constant<std::string>("ComptShieldingMaterial2"));
  Material shieldingPlateMaterial= description.material(description.constant<std::string>("ComptShieldingPlateMaterial"));

  double ComptShieldingX    = description.constant<double>("ComptShieldingX");
  double ComptShieldingY    = description.constant<double>("ComptShieldingY");
  double ComptShieldingZ    = description.constant<double>("ComptShieldingZ");
  double ComptShieldingZpos = description.constant<double>("ComptShieldingZpos");
  double BPipeRLG           = description.constant<double>("BPipeRLG");
  double FloorSurfaceYpos   = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest        = (description.constant<int>("OverlapTest") != 0);

  double floordy = -(ComptShieldingY/2.0 + FloorSurfaceYpos);

  // -----------------------------------------------------------------------
  // Main shielding block with beam pipe hole
  // -----------------------------------------------------------------------
  Box  solidShielding0(ComptShieldingX/2.0, ComptShieldingY/2.0, ComptShieldingZ/2.0);
  Tube solidShieldingH(0.0, BPipeRLG, 0.55*ComptShieldingZ, 0.0, 2.0*M_PI);

  SubtractionSolid solidShielding1("solidComptShielding", solidShielding0, solidShieldingH,
      Position(0.0, floordy, 0.0));

  // Commented: 3-layer version
  // Volume logicShielding("logicComptShielding", solidShielding1, shieldingMaterial);
  // Volume logicShielding2("logicComptShielding2", solidShielding1, shieldingMaterial2);
  // double zs = ComptShieldingZ;
  // motherVol.placeVolume(logicShielding,  Position(0, -floordy, ComptShieldingZpos)).addPhysVolID("ComptShieldingFe", 0);
  // motherVol.placeVolume(logicShielding2, Position(0, -floordy, ComptShieldingZpos - zs)).addPhysVolID("ComptShielding1Al", 0);
  // motherVol.placeVolume(logicShielding2, Position(0, -floordy, ComptShieldingZpos + zs)).addPhysVolID("ComptShielding2Al", 0);

  // -----------------------------------------------------------------------
  // Lead plate embedded in front of shielding
  // -----------------------------------------------------------------------
  double ShieldingPlateX    = ComptShieldingX;
  double ShieldingPlateY    = 10.0*dd4hep::cm;
  double ShieldingPlateZ    = 30.0*dd4hep::cm;
  double ShieldingPlatePosZ = ComptShieldingZpos - ComptShieldingZ/2.0 + ShieldingPlateZ/2.0;

  Box  solidShieldingPlate0(ShieldingPlateX/2.0, ShieldingPlateY/2.0, ShieldingPlateZ/2.0);
  Tube solidShieldingT(0.0, BPipeRLG, 0.55*ShieldingPlateZ, 0.0, 2.0*M_PI);

  SubtractionSolid solidShieldingPlate("solidShieldingPlate", solidShieldingPlate0, solidShieldingT,
      Position(0.0, 0.0, 0.0));
  Volume logicShieldingPlate("logicComptShieldingPlate", solidShieldingPlate, shieldingPlateMaterial);

  // -----------------------------------------------------------------------
  // 1-layer concrete shielding with plate cut-out
  // -----------------------------------------------------------------------
  SubtractionSolid solidShieldingP("solidShieldingP", solidShielding1, solidShieldingPlate0,
      Position(0.0, floordy, -ComptShieldingZ/2.0 + ShieldingPlateZ/2.0));
  Volume logicShielding("logicComptShielding", solidShieldingP, concreteMaterial);

  PlacedVolume pvShield = motherVol.placeVolume(logicShielding,
      Position(0.0, -floordy, ComptShieldingZpos));
  pvShield.addPhysVolID("ComptShieldingConcrete", 0);
  if (OverlapTest) pvShield.ptr()->CheckOverlaps();

  PlacedVolume pvPlate = motherVol.placeVolume(logicShieldingPlate,
      Position(0.0, 0.0, ShieldingPlatePosZ));
  pvPlate.addPhysVolID("ComptShieldingPlate", 0);
  if (OverlapTest) pvPlate.ptr()->CheckOverlaps();

  // Commented: pipe through shielding
  // double lshildpipe   = ...;
  // double shildpipezpos= ...;
  // Tube solidShieldingPipe(BPipeR-BPipeThickness, BPipeR, lshildpipe/2.0, 0.0, 2.0*M_PI);
  // ...

  // Commented: concrete support for shielding (not finished)
  // Box solidShieldingSupp(ComptShieldingX/2.0, ComptShieldingSuppY/2.0, 3*ComptShieldingZ/2.0);
  // ...
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxDetectorComptonFluka(dd4hep::Detector& description,
                                            xml_h             e,
                                            dd4hep::SensitiveDetector sd)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxDetectorComptonFluka det;
  det.Construct(description, sdet, e, sd);

  return sdet;
}

DECLARE_DETELEMENT(LxDetectorComptonFluka, create_LxDetectorComptonFluka)
