//
/// \brief Implementation of the LxHICSBeamDump class (DD4hep version)
//

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "Math/Vector3D.h"

#include "LxAux.h"
#include "LxHICSBeamDump.h"

using namespace dd4hep;

// ---------------------------------------------------------------------------
void LxHICSBeamDump::Construct(dd4hep::Detector&   description,
                                dd4hep::DetElement& sdet,
                                xml_h&              e)
{
  Volume fLogicWorld = description.worldVolume();

  xml_comp_t x_det(e);
  Assembly   envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  Material hicsDumpMaterial   = description.material(description.constant<std::string>("HICSDumpMaterial"));
  Material HDInsertAlMaterial = description.material(description.constant<std::string>("HICSDumpAlInsertMaterial"));
  Material environmentMaterial= description.material(description.constant<std::string>("EnvironmentMaterial"));
  Material vacuumMaterial     = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));

  double CollimatorRin       = description.constant<double>("CollimatorRin");
  double HICSDumpR           = description.constant<double>("HICSDumpR");
  double HICSDumpZ           = description.constant<double>("HICSDumpZ");
  double HICSDumpAlInsertR   = description.constant<double>("HICSDumpAlInsertR");
  double HICSDumpAlInsertZ   = description.constant<double>("HICSDumpAlInsertZ");
  double HICSDumpAirInsertR  = description.constant<double>("HICSDumpAirInsertR");
  double HICSDumpAirInsertZ  = description.constant<double>("HICSDumpAirInsertZ");
  double HICSDumpFrontXPos   = description.constant<double>("HICSDumpFrontXPos");
  double HICSDumpFrontZPos   = description.constant<double>("HICSDumpFrontZPos");
  double HICSDump2TPosX      = description.constant<double>("HICSDump2TPosX");
  double IPMagnetZpos        = description.constant<double>("IPMagnetZpos");
  double HICSShieldingXPos   = description.constant<double>("HICSShieldingXPos");
  double HICSShieldingY      = description.constant<double>("HICSShieldingY");
  double HICSShieldingSideZ  = description.constant<double>("HICSShieldingSideZ");
  double HICSShieldingSupportZ = description.constant<double>("HICSShieldingSupportZ");
  double HICSDetBottomSupportX = description.constant<double>("HICSDetBottomSupportX");
  double ComptonElectronBeamtoStageY = description.constant<double>("ComptonElectronBeamtoStageY");
  double ScintFrameThickness = description.constant<double>("ScintFrameThickness");
  double FlashMFieldLength   = description.constant<double>("FlashMFieldLength");
  double ComptonElBackshift  = description.constant<double>("ComptonElBackshift");
  bool   OverlapTest         = (description.constant<int>("OverlapTest") != 0);

  fHICSDumpAngle = std::atan2(HICSDumpFrontXPos, HICSDumpFrontZPos - IPMagnetZpos);

  // -----------------------------------------------------------------------
  // Gamma beam hole cut (shared by dump, Al insert, air insert)
  // -----------------------------------------------------------------------
  Tube solidHICSDumpGammaBeamHoleCut(0.0, CollimatorRin, HICSDumpZ, 0.0, 2.0*M_PI);

  // -----------------------------------------------------------------------
  // Main dump body
  // -----------------------------------------------------------------------
  Solid solidHICSDump;
  if (HICSDumpR > HICSDumpFrontXPos / std::cos(fHICSDumpAngle)) {
    Tube solidHICSDump0(0.0, HICSDumpR, HICSDumpZ/2.0, 0.0, 2.0*M_PI);
    double dx = std::tan(fHICSDumpAngle) * HICSDumpZ/2.0
                + HICSDumpFrontXPos / std::cos(fHICSDumpAngle);
    // G4RotationMatrix(G4ThreeVector(0,1,0), angle): axis-angle around Y → RotationZYX(0, angle, 0)
    solidHICSDump = SubtractionSolid("solidHICSDump", solidHICSDump0, solidHICSDumpGammaBeamHoleCut,
        Transform3D(RotationZYX(0.0, fHICSDumpAngle, 0.0), Position(dx, 0.0, 0.0)));
  } else {
    solidHICSDump = Tube(0.0, HICSDumpR, HICSDumpZ/2.0, 0.0, 2.0*M_PI);
  }

  // -----------------------------------------------------------------------
  // Al insert
  // -----------------------------------------------------------------------
  Solid solidHICSDumpAlInsert;
  if (HICSDumpAlInsertR > HICSDumpFrontXPos / std::cos(fHICSDumpAngle)) {
    Tube solidHICSDumpAlInsert0(0.0, HICSDumpAlInsertR, HICSDumpAlInsertZ/2.0, 0.0, 2.0*M_PI);
    double dx = std::tan(fHICSDumpAngle) * HICSDumpAlInsertZ/2.0
                + HICSDumpFrontXPos / std::cos(fHICSDumpAngle);
    solidHICSDumpAlInsert = SubtractionSolid("solidHICSDumpAlInsert",
        solidHICSDumpAlInsert0, solidHICSDumpGammaBeamHoleCut,
        Transform3D(RotationZYX(0.0, fHICSDumpAngle, 0.0), Position(dx, 0.0, 0.0)));
  } else {
    solidHICSDumpAlInsert = Tube(0.0, HICSDumpAlInsertR, HICSDumpAlInsertZ/2.0, 0.0, 2.0*M_PI);
  }

  // -----------------------------------------------------------------------
  // Air insert
  // -----------------------------------------------------------------------
  Solid solidHICSDumpAirInsert;
  if (HICSDumpAirInsertR > HICSDumpFrontXPos / std::cos(fHICSDumpAngle)) {
    Tube solidHICSDumpAirInsert0(0.0, HICSDumpAirInsertR, HICSDumpAirInsertZ/2.0, 0.0, 2.0*M_PI);
    double dx = std::tan(fHICSDumpAngle) * HICSDumpAirInsertZ/2.0
                + HICSDumpFrontXPos / std::cos(fHICSDumpAngle);
    solidHICSDumpAirInsert = SubtractionSolid("solidHICSDumpAirInsert",
        solidHICSDumpAirInsert0, solidHICSDumpGammaBeamHoleCut,
        Transform3D(RotationZYX(0.0, fHICSDumpAngle, 0.0), Position(dx, 0.0, 0.0)));
  } else {
    solidHICSDumpAirInsert = Tube(0.0, HICSDumpAirInsertR, HICSDumpAirInsertZ/2.0, 0.0, 2.0*M_PI);
  }

  Volume logicHICSDump        ("logicHICSDump",         solidHICSDump,         hicsDumpMaterial);
  Volume logicHICSDumpAlInsert("logicHICSDumpAlInsert", solidHICSDumpAlInsert, HDInsertAlMaterial);
  Volume logicHICSDumpAirInsert("logicHICSDumpAirInsert",solidHICSDumpAirInsert,environmentMaterial);

  // Air insert inside Al insert
  PlacedVolume pvAir = logicHICSDumpAlInsert.placeVolume(logicHICSDumpAirInsert,
      Position(0.0, 0.0, (HICSDumpAirInsertZ - HICSDumpAlInsertZ)/2.0));
  pvAir.addPhysVolID("HICSDumpAirInsert", 0);
  if (OverlapTest) pvAir.ptr()->CheckOverlaps();

  // Al insert inside dump (copy 0: +X side)
  PlacedVolume pvAl0 = logicHICSDump.placeVolume(logicHICSDumpAlInsert,
      Position(0.0, 0.0, (HICSDumpAlInsertZ - HICSDumpZ)/2.0));
  pvAl0.addPhysVolID("HICSDumpAlInsert", 0);
  if (OverlapTest) pvAl0.ptr()->CheckOverlaps();

  // Al insert inside dump (copy 1: -X side, 2T component)
  PlacedVolume pvAl1 = logicHICSDump.placeVolume(logicHICSDumpAlInsert,
      Position(-HICSDump2TPosX, 0.0, (HICSDumpAlInsertZ - HICSDumpZ)/2.0));
  pvAl1.addPhysVolID("HICSDumpAlInsert", 1);
  if (OverlapTest) pvAl1.ptr()->CheckOverlaps();

  // Place dump into world
  double hicsdx = HICSDumpFrontXPos + 0.5*HICSDumpZ * std::sin(fHICSDumpAngle);
  double hicsdz = HICSDumpFrontZPos + 0.5*HICSDumpZ * std::cos(fHICSDumpAngle);
  PlacedVolume pvDump = envelope.placeVolume(logicHICSDump,
      Transform3D(RotationZYX(0.0, -fHICSDumpAngle, 0.0),
                  Position(-hicsdx, 0.0, hicsdz)));
  pvDump.addPhysVolID("HICSDumpAssembly", 0);
  if (OverlapTest) pvDump.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Vacuum in the gamma beam hole through the dump
  // -----------------------------------------------------------------------
  double frontZpos  = HICSDumpFrontZPos + HICSDumpFrontXPos * std::tan(fHICSDumpAngle);
  double lholeVac   = HICSDumpZ / std::cos(fHICSDumpAngle);
  double lpipe_pos_1st = frontZpos + lholeVac/2.0;

  // G4CutTubs low normal:  ( sin,  0, -cos)
  // G4CutTubs high normal: (-sin,  0,  cos)
  CutTube solidBeamPipeGammaHICSdumpVac(0.0, CollimatorRin, lholeVac/2.0, 0.0, 2.0*M_PI,
       std::sin(fHICSDumpAngle), 0.0, -std::cos(fHICSDumpAngle),
      -std::sin(fHICSDumpAngle), 0.0,  std::cos(fHICSDumpAngle));
  Volume logicBeamPipeGammaHICSdumpVac("logicBeamPipeGammaHICSdumpVac",
                                        solidBeamPipeGammaHICSdumpVac, vacuumMaterial);
  PlacedVolume pvVac = envelope.placeVolume(logicBeamPipeGammaHICSdumpVac,
      Position(0.0, 0.0, lpipe_pos_1st));
  pvVac.addPhysVolID("BeamPipeGammaHICSdumpVac", 0);
  if (OverlapTest) pvVac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // 2T dump component
  // -----------------------------------------------------------------------
  double z2t = HICSDumpZ;
  // G4CutTubs low normal:  (-sin, 0, -cos)
  // G4CutTubs high normal: (  0,  0,   1 )
  CutTube solidHICSDump2T0(0.0, HICSDumpR, z2t/2.0, 0.0, 2.0*M_PI,
      -std::sin(fHICSDumpAngle), 0.0, -std::cos(fHICSDumpAngle),
       0.0, 0.0, 1.0);
  Tube solidHICSDump2TCut(0.0, HICSDumpR, HICSDumpZ, 0.0, 2.0*M_PI);

  // G4RotationMatrix(G4ThreeVector(0,1,0), 0): identity → RotationZYX(0,0,0)
  SubtractionSolid solidHICSDump2T("solidHICSDump2T", solidHICSDump2T0, solidHICSDump2TCut,
      Position(HICSDump2TPosX * std::cos(fHICSDumpAngle), 0.0, 0.0));
  Volume logicHICSDump2T("logicHICSDump2T", solidHICSDump2T, hicsDumpMaterial);

  double zpos2t = hicsdz - HICSDump2TPosX * std::sin(fHICSDumpAngle);
  PlacedVolume pvDump2T = envelope.placeVolume(logicHICSDump2T,
      Transform3D(RotationZYX(0.0, -fHICSDumpAngle, 0.0),
                  Position(-hicsdx - HICSDump2TPosX * std::cos(fHICSDumpAngle), 0.0, zpos2t)));
  pvDump2T.addPhysVolID("HICSDump2T", 0);
  if (OverlapTest) pvDump2T.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Electron shielding and neutron absorber — placed directly into envelope
  // -----------------------------------------------------------------------
  ConstructElectronShielding(description, envelope);
  ConstructNeutronAbsorber(description, envelope);

  // -----------------------------------------------------------------------
  // Shield support assembly
  // -----------------------------------------------------------------------
  double supy;
  Assembly supportAssembly = ConstructSupportAssembly(description, supy);
  double supxpos = HICSShieldingXPos;
  double supypos = -(HICSShieldingY + supy)/2.0;
  double supzpos = HICSDumpFrontZPos - HICSDumpR*std::sin(fHICSDumpAngle) - HICSShieldingSideZ;
  envelope.placeVolume(supportAssembly, Position(supxpos, supypos, supzpos));

  // -----------------------------------------------------------------------
  // Detector support assembly
  // -----------------------------------------------------------------------
  double supporthight;
  Assembly hicsDetSupport = ConstructHICSElDetSupportAssembly(description, envelope, supporthight);

  double detxpos  = HICSDetBottomSupportX/2.0;
  double detypos  = -ComptonElectronBeamtoStageY - ScintFrameThickness;
  double detzpos  = IPMagnetZpos + FlashMFieldLength/2.0 + ComptonElBackshift;
  envelope.placeVolume(hicsDetSupport, Position(-detxpos, detypos, detzpos));
}


// ---------------------------------------------------------------------------
void LxHICSBeamDump::ConstructElectronShielding(dd4hep::Detector& description,
                                                 dd4hep::Volume&   motherVol)
{
  Material hicsShielSMaterial = description.material(description.constant<std::string>("HICSShieldingSideMaterial"));
  Material hicsShielMMaterial = description.material(description.constant<std::string>("HICSShieldingMiddleMaterial"));

  double HICSShieldingX         = description.constant<double>("HICSShieldingX");
  double HICSShieldingY         = description.constant<double>("HICSShieldingY");
  double HICSShieldingSideZ     = description.constant<double>("HICSShieldingSideZ");
  double HICSShieldingMiddleZ   = description.constant<double>("HICSShieldingMiddleZ");
  double HICSShieldingXPos      = description.constant<double>("HICSShieldingXPos");
  double HICSDumpFrontXPos      = description.constant<double>("HICSDumpFrontXPos");
  double HICSDumpFrontZPos      = description.constant<double>("HICSDumpFrontZPos");
  double HICSDumpAlInsertR      = description.constant<double>("HICSDumpAlInsertR");
  double HICSDumpAlInsertZ      = description.constant<double>("HICSDumpAlInsertZ");
  double HICSDumpR              = description.constant<double>("HICSDumpR");
  double HICSDump2TPosX         = description.constant<double>("HICSDump2TPosX");
  double BPipeRLG               = description.constant<double>("BPipeRLG");
  double HICSShieldingGapX      = description.constant<double>("HICSShieldingGapX");
  double HICSShieldingGapY      = description.constant<double>("HICSShieldingGapY");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Front shielding plane (with dump, 2T dump, beam pipe, and gap cutouts)
  // -----------------------------------------------------------------------
  Box  solidHICSShieldingSide1(HICSShieldingX/2.0, HICSShieldingY/2.0, HICSShieldingSideZ/2.0);
  Tube solidHICSShieldSideCut(0.0, HICSDumpAlInsertR, HICSDumpAlInsertZ/2.0, 0.0, 2.0*M_PI);

  double dx = -(HICSShieldingXPos + HICSDumpFrontXPos);
  SubtractionSolid solidHICSShieldingSide0("solidHICSShieldingSide0",
      solidHICSShieldingSide1, solidHICSShieldSideCut, Position(dx, 0.0, 0.0));

  // Commented: alternative beam gap cut
  // double xscut = HICSDumpFrontXPos + 4.0*CollimatorRin;
  // Box solidHICSShieldingSideCut1(xscut/2.0, 2.0*CollimatorRin, HICSShieldingSideZ);
  // SubtractionSolid solidHICSShieldingSide2("solidHICSShieldingSide2",
  //     solidHICSShieldingSide0, solidHICSShieldingSideCut1,
  //     Position(dx + xscut/2.0, 0.0, 0.0));

  Tube solidHICSShieldSide4Bpipe(0.0, BPipeRLG, HICSShieldingSideZ, 0.0, 2.0*M_PI);
  SubtractionSolid solidHICSShieldingSide2("solidHICSShieldingSide2",
      solidHICSShieldingSide0, solidHICSShieldSide4Bpipe,
      Position(-HICSShieldingXPos, 0.0, 0.0));

  Box solidHICSShieldingSideGapCut(HICSShieldingGapX/2.0, HICSShieldingGapY/2.0, HICSShieldingSideZ/2.0);
  SubtractionSolid solidHICSShieldingSide3("solidHICSShieldingSide3",
      solidHICSShieldingSide2, solidHICSShieldingSideGapCut,
      Position(dx - HICSShieldingGapX/2.0, 0.0, -HICSShieldingSideZ/2.0));

  double dx2t = -(HICSShieldingXPos + HICSDumpFrontXPos + HICSDump2TPosX);
  SubtractionSolid solidHICSShieldingSide("solidHICSShieldingSide",
      solidHICSShieldingSide3, solidHICSShieldSideCut, Position(dx2t, 0.0, 0.0));

  Volume logicHICSShieldingSide("logicHICSShieldingSide", solidHICSShieldingSide, hicsShielSMaterial);

  // -----------------------------------------------------------------------
  // Back shielding plane (with angled dump cutouts)
  // -----------------------------------------------------------------------
  Box  solidHICSShieldingMiddle0(HICSShieldingX/2.0, HICSShieldingY/2.0, HICSShieldingMiddleZ/2.0);
  Tube solidHICSSieldMiddleCut  (0.0, HICSDumpR, HICSShieldingMiddleZ, 0.0, 2.0*M_PI);

  double dxm = std::tan(fHICSDumpAngle) *
               (0.5*HICSShieldingMiddleZ - HICSDumpR*std::sin(fHICSDumpAngle));

  // G4RotationMatrix(G4ThreeVector(0,-1,0), angle): axis-angle around -Y by angle
  // = rotation around Y by -angle → RotationZYX(0, -angle, 0)
  SubtractionSolid solidHICSShieldingMiddle1("solidHICSShieldingMiddle1",
      solidHICSShieldingMiddle0, solidHICSSieldMiddleCut,
      Transform3D(RotationZYX(0.0, -fHICSDumpAngle, 0.0),
                  Position(-dxm - HICSShieldingXPos - HICSDumpFrontXPos, 0.0, 0.0)));
  SubtractionSolid solidHICSShieldingMiddle("solidHICSShieldingMiddle",
      solidHICSShieldingMiddle1, solidHICSSieldMiddleCut,
      Transform3D(RotationZYX(0.0, -fHICSDumpAngle, 0.0),
                  Position(-dxm - HICSShieldingXPos - HICSDumpFrontXPos
                           - HICSDump2TPosX/std::cos(fHICSDumpAngle), 0.0, 0.0)));

  Volume logicHICSShieldingMiddle("logicHICSShieldingMiddle", solidHICSShieldingMiddle, hicsShielMMaterial);

  // -----------------------------------------------------------------------
  // Placements
  // -----------------------------------------------------------------------
  double zpos = HICSDumpFrontZPos - HICSDumpR*std::sin(fHICSDumpAngle) - HICSShieldingSideZ/2.0;
  PlacedVolume pvSide = motherVol.placeVolume(logicHICSShieldingSide,
      Position(HICSShieldingXPos, 0.0, zpos));
  pvSide.addPhysVolID("HICSShieldingSide", 0);
  if (OverlapTest) pvSide.ptr()->CheckOverlaps();

  zpos += HICSShieldingSideZ/2.0 + HICSShieldingMiddleZ/2.0;
  PlacedVolume pvMiddle = motherVol.placeVolume(logicHICSShieldingMiddle,
      Position(HICSShieldingXPos, 0.0, zpos));
  pvMiddle.addPhysVolID("HICSShieldingMiddle", 0);
  if (OverlapTest) pvMiddle.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxHICSBeamDump::ConstructNeutronAbsorber(dd4hep::Detector& description,
                                               dd4hep::Volume&   motherVol)
{
  Material hicsNeutronAbsorberMaterial = description.material(
      description.constant<std::string>("HICSNeutronAbsorberMaterial"));

  double HICSNeutronAbsorberX   = description.constant<double>("HICSNeutronAbsorberX");
  double HICSNeutronAbsorberY   = description.constant<double>("HICSNeutronAbsorberY");
  double HICSNeutronAbsorberZ   = description.constant<double>("HICSNeutronAbsorberZ");
  double HICSShieldingXPos      = description.constant<double>("HICSShieldingXPos");
  double HICSShieldingGapX      = description.constant<double>("HICSShieldingGapX");
  double HICSShieldingGapY      = description.constant<double>("HICSShieldingGapY");
  double HICSShieldingSideZ     = description.constant<double>("HICSShieldingSideZ");
  double HICSShieldingSupportY  = description.constant<double>("HICSShieldingSupportY");
  double HICSDumpFrontXPos      = description.constant<double>("HICSDumpFrontXPos");
  double HICSDumpFrontZPos      = description.constant<double>("HICSDumpFrontZPos");
  double HICSDumpAlInsertR      = description.constant<double>("HICSDumpAlInsertR");
  double HICSDump2TPosX         = description.constant<double>("HICSDump2TPosX");
  double HICSDumpR              = description.constant<double>("HICSDumpR");
  double BPipeRLG               = description.constant<double>("BPipeRLG");
  double FloorSurfaceYpos       = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  double dx  = -(HICSShieldingXPos + HICSDumpFrontXPos);
  double dx2t= -(HICSShieldingXPos + HICSDumpFrontXPos + HICSDump2TPosX);

  // -----------------------------------------------------------------------
  // Front absorber plane
  // -----------------------------------------------------------------------
  Box  solidHICSNeutronAbsorberSide1(HICSNeutronAbsorberX/2.0, HICSNeutronAbsorberY/2.0, HICSNeutronAbsorberZ/2.0);
  Tube solidNeutronAbsorberSideCut  (0.0, HICSDumpAlInsertR, HICSNeutronAbsorberZ, 0.0, 2.0*M_PI);

  SubtractionSolid solidNeutronAbsorberSide0("solidNeutronAbsorberSide0",
      solidHICSNeutronAbsorberSide1, solidNeutronAbsorberSideCut, Position(dx, 0.0, 0.0));

  Tube solidNeutronAbsorber4Bpipe(0.0, BPipeRLG, HICSShieldingSideZ, 0.0, 2.0*M_PI);
  SubtractionSolid solidNeutronAbsorberSide2("solidNeutronAbsorberSide2",
      solidNeutronAbsorberSide0, solidNeutronAbsorber4Bpipe,
      Position(-HICSShieldingXPos, 0.0, 0.0));

  Box solidNeutronAbsorberGapCut(HICSShieldingGapX/2.0, HICSShieldingGapY/2.0, HICSNeutronAbsorberZ);
  SubtractionSolid solidNeutronAbsorberSide3("solidNeutronAbsorberSide3",
      solidNeutronAbsorberSide2, solidNeutronAbsorberGapCut,
      Position(dx - HICSShieldingGapX/2.0, 0.0, 0.0));

  SubtractionSolid solidNeutronAbsorberSide("solidNeutronAbsorberSide",
      solidNeutronAbsorberSide3, solidNeutronAbsorberSideCut, Position(dx2t, 0.0, 0.0));

  Volume logicNeutronAbsorberSide("logicNeutronAbsorberSide", solidNeutronAbsorberSide,
                                   hicsNeutronAbsorberMaterial);

  // -----------------------------------------------------------------------
  // Top and bottom absorber parts
  // -----------------------------------------------------------------------
  double supy    = -FloorSurfaceYpos - HICSNeutronAbsorberY/2.0;
  double suptopx = 1.1 * HICSNeutronAbsorberX;
  double suptopy = HICSShieldingSupportY - supy;
  double suptopz = HICSNeutronAbsorberZ;

  Box solidHICSNeutronAbsorberTop1(suptopx/2.0, suptopy/2.0, suptopz/2.0);
  Box solidHICSNeutronAbsorberTopCut1(HICSNeutronAbsorberX/2.0, HICSNeutronAbsorberY, HICSShieldingSideZ);
  SubtractionSolid soliHICSNeutronAbsorberTop("soliHICSNeutronAbsorberTop",
      solidHICSNeutronAbsorberTop1, solidHICSNeutronAbsorberTopCut1,
      Position(0.0, -suptopy/2.0, 0.0));
  Volume logicHICSNeutronAbsorberTop("logicHICSNeutronAbsorberTop", soliHICSNeutronAbsorberTop,
                                      hicsNeutronAbsorberMaterial);

  Box    solidHICSNeutronAbsorberBottom(suptopx/2.0, supy/2.0, suptopz/2.0);
  Volume logicHICSNeutronAbsorberBottom("logicHICSNeutronAbsorberBottom",
                                         solidHICSNeutronAbsorberBottom, hicsNeutronAbsorberMaterial);

  // -----------------------------------------------------------------------
  // Placements
  // -----------------------------------------------------------------------
  double zpos = HICSDumpFrontZPos - HICSDumpR*std::sin(fHICSDumpAngle)
                - HICSShieldingSideZ - HICSNeutronAbsorberZ/2.0;

  PlacedVolume pvSide = motherVol.placeVolume(logicNeutronAbsorberSide,
      Position(HICSShieldingXPos, 0.0, zpos));
  pvSide.addPhysVolID("HICSNeutronAbsorberSide", 0);
  if (OverlapTest) pvSide.ptr()->CheckOverlaps();

  PlacedVolume pvTop = motherVol.placeVolume(logicHICSNeutronAbsorberTop,
      Position(HICSShieldingXPos, 0.5*(suptopy - HICSNeutronAbsorberY), zpos));
  pvTop.addPhysVolID("HICSNeutronAbsorberTop", 0);
  if (OverlapTest) pvTop.ptr()->CheckOverlaps();

  PlacedVolume pvBot = motherVol.placeVolume(logicHICSNeutronAbsorberBottom,
      Position(HICSShieldingXPos, FloorSurfaceYpos + 0.5*supy, zpos));
  pvBot.addPhysVolID("HICSNeutronAbsorberBottom", 0);
  if (OverlapTest) pvBot.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxHICSBeamDump::ConstructSupportAssembly(dd4hep::Detector& description,
                                                           double&           supporthight)
{
  Material hicsShieldSupportMaterial = description.material(
      description.constant<std::string>("HICSShieldingSupportMaterial"));
  Material hicsShieldTopMaterial = description.material(
      description.constant<std::string>("HICSShieldingTopMaterial"));

  double HICSShieldingX         = description.constant<double>("HICSShieldingX");
  double HICSShieldingY         = description.constant<double>("HICSShieldingY");
  double HICSShieldingSupportZ  = description.constant<double>("HICSShieldingSupportZ");
  double HICSShieldingSupportY  = description.constant<double>("HICSShieldingSupportY");
  double HICSShieldingSideZ     = description.constant<double>("HICSShieldingSideZ");
  double HICSShieldingMiddleZ   = description.constant<double>("HICSShieldingMiddleZ");
  double HICSShieldingXPos      = description.constant<double>("HICSShieldingXPos");
  double HICSDumpFrontXPos      = description.constant<double>("HICSDumpFrontXPos");
  double HICSDump2TPosX         = description.constant<double>("HICSDump2TPosX");
  double HICSDumpR              = description.constant<double>("HICSDumpR");
  double HICSDumpZ              = description.constant<double>("HICSDumpZ");
  double CollimatorRout         = description.constant<double>("CollimatorRout");
  double FloorSurfaceYpos       = description.constant<double>("FloorSurfaceYpos");

  double supx  = 1.1 * HICSShieldingX;
  double supy  = -FloorSurfaceYpos - HICSShieldingY/2.0;
  double supz  = HICSShieldingSupportZ;

  Volume pedstal = LxAux::BuildPedestal(description, "HICSShield", supx, supy, supz);

  // Top concrete part
  double suptopx    = supx;
  double suptopy    = HICSShieldingSupportY - supy;
  double suptopz    = supz;
  double suptopcutz = HICSShieldingSideZ + HICSShieldingMiddleZ;

  Box solidHICSSSupportTop1(suptopx/2.0, suptopy/2.0, suptopz/2.0);
  Box solidHICSSSupportTopCut1(HICSShieldingX/2.0, HICSShieldingY, suptopcutz);
  SubtractionSolid soliHICSSSupportTop2("soliHICSSSupportTop2",
      solidHICSSSupportTop1, solidHICSSSupportTopCut1,
      Position(0.0, -suptopy/2.0, -suptopz/2.0));

  double topdumpcutx = 1.2 * (2.0*HICSDumpR + HICSDump2TPosX);
  double topdumpcuty = HICSShieldingY;
  double topdumpcutz = 1.1 * HICSDumpZ;
  Box solidHICSSSupportTopCut2(topdumpcutx/2.0, topdumpcuty, topdumpcutz/2.0);
  double trdumpcutx = -HICSShieldingXPos - HICSDumpFrontXPos - 0.5*HICSDump2TPosX;
  double trdumpcuty = -suptopy/2.0;
  double trdumpcutz = -0.5*(suptopz - topdumpcutz) + HICSShieldingSideZ;
  SubtractionSolid soliHICSSSupportTop3("soliHICSSSupportTop3",
      soliHICSSSupportTop2, solidHICSSSupportTopCut2,
      Position(trdumpcutx, trdumpcuty, trdumpcutz));

  Tube solidHICSSSupportTopCut3(0.0, CollimatorRout, supz, 0.0, 2.0*M_PI);
  double trcolimcutx = -HICSShieldingXPos;
  double trcolimcuty = -0.5*(suptopy - HICSShieldingY);
  SubtractionSolid soliHICSSSupportTop("soliHICSSSupportTop4",
      soliHICSSSupportTop3, solidHICSSSupportTopCut3,
      Position(trcolimcutx, trcolimcuty, 0.0));

  Volume logicHICSSSupportTop("logicHICSSSupportTop", soliHICSSSupportTop, hicsShieldTopMaterial);

  Assembly hicsShieldSupportAssembly("HICSShieldSupport_assembly");
  hicsShieldSupportAssembly.placeVolume(pedstal,
      Position(0.0, 0.0, HICSShieldingSupportZ/2.0));
  hicsShieldSupportAssembly.placeVolume(logicHICSSSupportTop,
      Position(0.0, 0.5*(supy + suptopy), HICSShieldingSupportZ/2.0));

  supporthight = supy;
  return hicsShieldSupportAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxHICSBeamDump::ConstructHICSElDetSupportAssembly(
    dd4hep::Detector& description,
    dd4hep::Volume&   motherVol,
    double&           sphight)
{
  Material opppDetSupportMaterial = description.material(
      description.constant<std::string>("HICSDetSupportMaterial"));

  double HICSDetBottomSupportX  = description.constant<double>("HICSDetBottomSupportX");
  double HICSDetBottomSupportY  = description.constant<double>("HICSDetBottomSupportY");
  double HICSDetBottomSupportZ  = description.constant<double>("HICSDetBottomSupportZ");
  double IPMagnetZpos           = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength      = description.constant<double>("FlashMFieldLength");
  double ComptonElBackshift     = description.constant<double>("ComptonElBackshift");
  double OPPPDetBottomSupportY  = description.constant<double>("OPPPDetBottomSupportY");
  double HICSBasePlateX         = description.constant<double>("HICSBasePlateX");
  double HICSBasePlateY         = description.constant<double>("HICSBasePlateY");
  double HICSBasePlateZ         = description.constant<double>("HICSBasePlateZ");
  double FloorSurfaceYpos       = description.constant<double>("FloorSurfaceYpos");
  double ComptonElectronBeamtoStageY = description.constant<double>("ComptonElectronBeamtoStageY");
  double ScintFrameThickness    = description.constant<double>("ScintFrameThickness");
  double CerenkovChannelLayers  = description.constant<double>("CerenkovChannelLayers");
  double CerenkovStrawZFrequency= description.constant<double>("CerenkovStrawZFrequency");
  double CerenkovStrawInnerRadius        = description.constant<double>("CerenkovStrawInnerRadius");
  double CerenkovStrawGraphiteLayerThickness = description.constant<double>("CerenkovStrawGraphiteLayerThickness");
  double CerenkovStrawAlLayerThickness   = description.constant<double>("CerenkovStrawAlLayerThickness");
  double CerenkovStrawKaptonLayerThickness   = description.constant<double>("CerenkovStrawKaptonLayerThickness");
  double CerenkovStrawPolyurethaneLayerThickness = description.constant<double>("CerenkovStrawPolyurethaneLayerThickness");
  double CerenkovStrawPlateBuffer  = description.constant<double>("CerenkovStrawPlateBuffer");
  double CerenkovBoxThickness      = description.constant<double>("CerenkovBoxThickness");
  double CerenkovShieldingPlateThickness = description.constant<double>("CerenkovShieldingPlateThickness");
  bool   OverlapTest            = (description.constant<int>("OverlapTest") != 0);

  Box    soliHICSScintBottomSupport(HICSDetBottomSupportX/2.0,
                                    HICSDetBottomSupportY/2.0,
                                    HICSDetBottomSupportZ/2.0);
  Volume logicHICSScintBottomSupport("logicHICSScintBottomSupport",
                                      soliHICSScintBottomSupport, opppDetSupportMaterial);

  double detzpos  = IPMagnetZpos + FlashMFieldLength/2.0 + ComptonElBackshift;

  int    ccl   = (int)CerenkovChannelLayers;
  double cszf  = CerenkovStrawZFrequency;
  double sir   = CerenkovStrawInnerRadius;
  double sglt  = CerenkovStrawGraphiteLayerThickness;
  double salt  = CerenkovStrawAlLayerThickness;
  double sklt  = CerenkovStrawKaptonLayerThickness;
  double splt  = CerenkovStrawPolyurethaneLayerThickness;
  double cspb  = CerenkovStrawPlateBuffer;
  double cbt   = CerenkovBoxThickness;
  double cspt  = CerenkovShieldingPlateThickness;

  double scmotherdz  = 2*ComptonElectronBeamtoStageY + 3*ScintFrameThickness;
  double cherenkovdz = 2.0 * ((ccl-1)*cszf/2.0 + sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt);
  double basezpos    = detzpos + HICSDetBottomSupportZ/2.0 - scmotherdz/2.0;

  Assembly hicsSupportAssembly("HICSDetSupport_assembly");

  Position asstr(0.0, 0.0, basezpos - detzpos);
  hicsSupportAssembly.placeVolume(logicHICSScintBottomSupport, asstr);

  // Hexapod
  double hexhight = 0.0;
  Assembly hexAssembly = LxAux::BuildHexapod(description, "HICS", hexhight);
  Position hexpos(0.0, asstr.Y() - (OPPPDetBottomSupportY + hexhight)/2.0, basezpos - detzpos);
  LxAux::AddAssemblyVolumes(hicsSupportAssembly, hexAssembly, hexpos);

  double baseypos = -ComptonElectronBeamtoStageY - hexhight
                   - HICSDetBottomSupportY - HICSBasePlateY/2.0;
  sphight = (HICSDetBottomSupportY + hexhight)/2.0 - hexpos.Y();

  // Table support and pedestal — placed into motherVol directly (not into assembly)
  double ypestal  = 0.5 * dd4hep::m;
  double ylevel   = baseypos + HICSBasePlateY/2.0;
  double tblhight = ylevel - ypestal - FloorSurfaceYpos;
  double tblx     = HICSBasePlateX;
  double tblz     = HICSBasePlateZ;

  Assembly tablesupport = LxAux::BuildTable(description, "HICSDetectorTable",
                                            tblx, tblhight, tblz, 4);
  motherVol.placeVolume(tablesupport, Position(-0.5*tblx, ylevel, basezpos));

  Volume pedestal = LxAux::BuildPedestal(description, "HICSDetector", tblx, ypestal, tblz);
  PlacedVolume pvPed = motherVol.placeVolume(pedestal,
      Position(-0.5*tblx, 0.5*ypestal + FloorSurfaceYpos, basezpos));
  if (OverlapTest) pvPed.ptr()->CheckOverlaps();

  return hicsSupportAssembly;
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxHICSBeamDump(dd4hep::Detector& description,
                                    xml_h             e,
                                    dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxHICSBeamDump hics;
  hics.Construct(description, sdet, e);

  return sdet;
}

DECLARE_DETELEMENT(LxHICSBeamDump, create_LxHICSBeamDump)
