//
/// \brief Implementation of the LxInteractionChamber class (DD4hep version)
//

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxAux.h"
#include "LxInteractionChamber.h"

using namespace dd4hep;


// ---------------------------------------------------------------------------
void LxInteractionChamber::Construct(dd4hep::Detector&   description,
                                     dd4hep::DetElement& sdet,
                                     xml_h&              e)
{
  // World (mother) volume + mandatory DetElement placement via assembly envelope.
  // The envelope serves as the single top-level assembly — no intermediate
  // fIPChamberAssembly needed, all volumes are placed directly here.
  Volume fLogicWorld = description.worldVolume();

  xml_comp_t x_det(e);
  Assembly   envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  // -----------------------------------------------------------------------
  // Read constants  (lxs->X  →  description.constant<T>("X"))
  // -----------------------------------------------------------------------
  Material tauIChamberMaterial = description.material(description.constant<std::string>("TAUIChamberMaterial"));
  Material vacuumMaterial      = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));
  Material beamPipeMaterial    = description.material(description.constant<std::string>("BeamPipeMaterial"));

  double TAUIChamberX             = description.constant<double>("TAUIChamberX");
  double TAUIChamberY             = description.constant<double>("TAUIChamberY");
  double TAUIChamberZ             = description.constant<double>("TAUIChamberZ");
  double TAUIChamberSideThickness = description.constant<double>("TAUIChamberSideThickness");
  double TAUIChamberFrontThickness= description.constant<double>("TAUIChamberFrontThickness");
  double TAUIChamberElPipeRIn     = description.constant<double>("TAUIChamberElPipeRIn");
  double TAUIChamberElPipeROut    = description.constant<double>("TAUIChamberElPipeROut");
  double TAUIChamberXpos          = description.constant<double>("TAUIChamberXpos");
  double TAUIChamberYpos          = description.constant<double>("TAUIChamberYpos");
  double TAUIChamberZpos          = description.constant<double>("TAUIChamberZpos");
  double TAUIChamberBottomX       = description.constant<double>("TAUIChamberBottomX");
  double TAUIChamberBottomY       = description.constant<double>("TAUIChamberBottomY");
  double TAUIChamberBottomZ       = description.constant<double>("TAUIChamberBottomZ");
  double TAUIChamberPipeLength    = description.constant<double>("TAUIChamberPipeLength");
  double TAUIChamberPipeFlangedH  = description.constant<double>("TAUIChamberPipeFlangedH");
  double TAUIChamberPipeFlangedL  = description.constant<double>("TAUIChamberPipeFlangedL");
  double TAUIChamberSupportH      = description.constant<double>("TAUIChamberSupportH");
  double TAUIChamberSupportR      = description.constant<double>("TAUIChamberSupportR");
  double FlashMFieldX             = description.constant<double>("FlashMFieldX");
  double FlashMagneteffY          = description.constant<double>("FlashMagneteffY");
  double IPMAgnetBeamPipeXGap     = description.constant<double>("IPMAgnetBeamPipeXGap");
  double BPipeThickness           = description.constant<double>("BPipeThickness");
  double TAUIChamberBeamBottomY   = description.constant<double>("TAUIChamberBeamBottomY");
  double ICBBoardTableX           = description.constant<double>("ICBBoardTableX");
  double ICBBoardTableZ           = description.constant<double>("ICBBoardTableZ");
  double ICBBoardTableY           = description.constant<double>("ICBBoardTableY");
  double ICBBoardTableLegR        = description.constant<double>("ICBBoardTableLegR");
  double ICBBoardTableGapY        = description.constant<double>("ICBBoardTableGapY");
  double FloorSurfaceYpos         = description.constant<double>("FloorSurfaceYpos");
  bool   OverlapTest              = (description.constant<int>("OverlapTest") != 0);

  // World offset — applied to all top-level placements into envelope,
  // mirrors ConstructIPChamber: tr0(TAUIChamberXpos, TAUIChamberYpos, TAUIChamberZpos)
  Position worldPos(TAUIChamberXpos, TAUIChamberYpos, TAUIChamberZpos);

  // -----------------------------------------------------------------------
  // Container box
  // -----------------------------------------------------------------------
  Box    solidTAUICContainer(TAUIChamberX/2.0, TAUIChamberY/2.0, TAUIChamberZ/2.0);
  Volume logicTAUICContainer("logicTAUICContainer", solidTAUICContainer, vacuumMaterial);

  // -----------------------------------------------------------------------
  // Side panels
  // -----------------------------------------------------------------------
  Box    solidTAUICSide(TAUIChamberSideThickness/2.0, TAUIChamberY/2.0, TAUIChamberZ/2.0);
  Volume logicTAUICSide("logicTAUICSide", solidTAUICSide, tauIChamberMaterial);

  double sideposX = (TAUIChamberX - TAUIChamberSideThickness)/2.0;
  PlacedVolume pvSide0 = logicTAUICContainer.placeVolume(logicTAUICSide, Position( sideposX, 0.0, 0.0));
  pvSide0.addPhysVolID("TAUICSide", 0);
  PlacedVolume pvSide1 = logicTAUICContainer.placeVolume(logicTAUICSide, Position(-sideposX, 0.0, 0.0));
  pvSide1.addPhysVolID("TAUICSide", 1);
  if (OverlapTest) { pvSide0.ptr()->CheckOverlaps(); pvSide1.ptr()->CheckOverlaps(); }

  // -----------------------------------------------------------------------
  // Front and rear panels
  // -----------------------------------------------------------------------
  double beamtocy = TAUIChamberYpos;
  double dxf      = TAUIChamberX - 2.0*TAUIChamberSideThickness;

  Box  solidTAUICFront1(dxf/2.0, TAUIChamberY/2.0, TAUIChamberFrontThickness/2.0);
  Tube solidTAUICFPCut(0.0, TAUIChamberElPipeRIn, TAUIChamberFrontThickness, 0.0, 2.0*M_PI);

  Position cutPos1(-TAUIChamberXpos, -beamtocy, 0.0);
  SubtractionSolid solidTAUICFront("solidTAUICFront", solidTAUICFront1, solidTAUICFPCut, cutPos1);
  Volume logicTAUICFront("logicTAUICFront", solidTAUICFront, tauIChamberMaterial);

  double ipmbphx = FlashMFieldX/2.0 - IPMAgnetBeamPipeXGap;
  Box solidTAUICFPRearCut(ipmbphx - BPipeThickness,
                          FlashMagneteffY/2.0 - BPipeThickness,
                          TAUIChamberFrontThickness);
  SubtractionSolid solidTAUICFPRear("solidTAUICFPRear", solidTAUICFront1, solidTAUICFPRearCut, cutPos1);
  Volume logicTAUICFPRear("logicTAUICFPRear", solidTAUICFPRear, tauIChamberMaterial);

  double frontposZ = (TAUIChamberZ - TAUIChamberFrontThickness)/2.0;
  PlacedVolume pvRear  = logicTAUICContainer.placeVolume(logicTAUICFPRear,  Position(0.0, 0.0,  frontposZ));
  PlacedVolume pvFront = logicTAUICContainer.placeVolume(logicTAUICFront,   Position(0.0, 0.0, -frontposZ));
  pvRear.addPhysVolID("TAUICRear",  0);
  pvFront.addPhysVolID("TAUICFront", 0);
  if (OverlapTest) { pvRear.ptr()->CheckOverlaps(); pvFront.ptr()->CheckOverlaps(); }

  // -----------------------------------------------------------------------
  // IP volume
  // -----------------------------------------------------------------------
  double IPVolumeR = description.constant<double>("IPVolumeR");
  double IPVolumeZ = description.constant<double>("IPVolumeZ");
  Tube   solidIPVolume(0.0, IPVolumeR, IPVolumeZ/2.0, 0.0, 2.0*M_PI);
  Volume logicIPVolume("logicIPVolume", solidIPVolume, vacuumMaterial);
  PlacedVolume pvIPVol = logicTAUICContainer.placeVolume(logicIPVolume,
      Position(-TAUIChamberXpos, -TAUIChamberYpos, -TAUIChamberZpos));
  pvIPVol.addPhysVolID("IPVolume", 0);
  if (OverlapTest) pvIPVol.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Mirrors
  // -----------------------------------------------------------------------
  ConstructMirrors(description, logicTAUICContainer);

  // -----------------------------------------------------------------------
  // Top and bottom panels
  // -----------------------------------------------------------------------
  Box    solidTAUICBottom(TAUIChamberBottomX/2.0, TAUIChamberBottomY/2.0, TAUIChamberBottomZ/2.0);
  Volume logicTAUICBottom("logicTAUICBottom", solidTAUICBottom, tauIChamberMaterial);
  Box    solidTAUICTop(TAUIChamberBottomX/2.0, TAUIChamberBottomY/2.0, TAUIChamberBottomZ/2.0);
  Volume logicTAUICTop("logicTAUICTop", solidTAUICTop, tauIChamberMaterial);

  // -----------------------------------------------------------------------
  // In/out pipes
  // -----------------------------------------------------------------------
  Tube solidTAUIChamberBPipe(TAUIChamberElPipeRIn, TAUIChamberElPipeROut,
                             0.5*TAUIChamberPipeLength, 0.0, 2.0*M_PI);
  Volume logicTAUIChamberBPipe("logicTAUIChamberBPipe", solidTAUIChamberBPipe, beamPipeMaterial);

  Tube solidTAUIChamberBPipeVac(0.0, TAUIChamberElPipeRIn,
                                0.5*TAUIChamberPipeLength, 0.0, 2.0*M_PI);
  Volume logicTAUIChamberBPipeVac("logicTAUIChamberBPipeVac", solidTAUIChamberBPipeVac, vacuumMaterial);

  Tube solidTAUIChamberBPipeFlange(TAUIChamberElPipeROut,
                                   TAUIChamberElPipeROut + TAUIChamberPipeFlangedH,
                                   0.5*TAUIChamberPipeFlangedL, 0.0, 2.0*M_PI);
  Volume logicTAUIChamberBPipeFlange("logicTAUIChamberBPipeFlange", solidTAUIChamberBPipeFlange, beamPipeMaterial);

  Box solidTAUIChamberBPipeOut1(ipmbphx, FlashMagneteffY/2.0, 0.5*TAUIChamberPipeLength);
  Box solidTAUIChamberBPipeOutCut(ipmbphx - BPipeThickness,
                                  FlashMagneteffY/2.0 - BPipeThickness,
                                  TAUIChamberPipeLength);
  SubtractionSolid solidTAUIChamberBPipeOut("solidTAUIChamberBPipeOut",
                                            solidTAUIChamberBPipeOut1, solidTAUIChamberBPipeOutCut);
  Volume logicTAUIChamberBPipeOut("logicTAUIChamberBPipeOut", solidTAUIChamberBPipeOut, beamPipeMaterial);

  Box    solidTAUIChamberBPipeOutVac(ipmbphx - BPipeThickness,
                                     FlashMagneteffY/2.0 - BPipeThickness,
                                     0.5*TAUIChamberPipeLength);
  Volume logicTAUIChamberBPipeOutVac("logicTAUIChamberBPipeOutVac", solidTAUIChamberBPipeOutVac, vacuumMaterial);

  // -----------------------------------------------------------------------
  // Support legs (Polycone)
  // -----------------------------------------------------------------------
  const int nsplane = 4;
  double zsf[nsplane]    = {0.0, 0.97, 0.98, 1.0};
  double zsplane[nsplane];
  std::transform(zsf, zsf+nsplane, zsplane,
                 std::bind(std::multiplies<double>(), std::placeholders::_1, TAUIChamberSupportH));
  double rr = TAUIChamberSupportR;
  double rsouter[nsplane] = {rr,  rr,  0.5*rr, 0.5*rr};
  double rsinner[nsplane] = {0.0, 0.0, 0.0,    0.0   };

  std::vector<double> zv(zsplane, zsplane+nsplane);
  std::vector<double> riv(rsinner, rsinner+nsplane);
  std::vector<double> rov(rsouter, rsouter+nsplane);
  Polycone solidTAUICSupport(0.0, 2.0*M_PI, riv, rov, zv);
  Volume   logicTAUICSupport("logicTAUICSupport", solidTAUICSupport, tauIChamberMaterial);

  // -----------------------------------------------------------------------
  // Place everything directly into envelope (+ worldPos offset)
  // -----------------------------------------------------------------------

  // Support legs — 4 corners
  // G4RotationMatrix(0.0, 0.5*M_PI, 0.0) is ZXZ Euler: rotation around X by pi/2 → Z maps to Y
  RotationZYX srot(0.0, 0.0, -M_PI/2.0);
  for (int ii = 0; ii < 4; ++ii) {
    double sdx = std::pow(-1.0,  ii    & 1) * (TAUIChamberX/2.0 - TAUIChamberSupportR);
    double sdz = std::pow(-1.0, (ii>>1)& 1) * (TAUIChamberZ/2.0 - TAUIChamberSupportR);
    Position spos(sdx + TAUIChamberXpos,
                  -TAUIChamberSupportH - TAUIChamberY/2.0 - TAUIChamberBottomY + TAUIChamberYpos,
                  sdz + TAUIChamberZpos);
    PlacedVolume pvS = envelope.placeVolume(logicTAUICSupport, Transform3D(srot, spos));
    pvS.addPhysVolID("TAUICSupport", ii);
    if (OverlapTest) pvS.ptr()->CheckOverlaps();
  }

  // Container (at worldPos)
  PlacedVolume pvCont = envelope.placeVolume(logicTAUICContainer, worldPos);
  pvCont.addPhysVolID("TAUICContainer", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();

  // Top
  PlacedVolume pvTop = envelope.placeVolume(logicTAUICTop,
      Position(TAUIChamberXpos,
               (TAUIChamberY + TAUIChamberBottomY)/2.0 + TAUIChamberYpos,
               TAUIChamberZpos));
  pvTop.addPhysVolID("TAUICTop", 0);
  if (OverlapTest) pvTop.ptr()->CheckOverlaps();

  // Bottom
  PlacedVolume pvBot = envelope.placeVolume(logicTAUICBottom,
      Position(TAUIChamberXpos,
               -(TAUIChamberY + TAUIChamberBottomY)/2.0 + TAUIChamberYpos,
               TAUIChamberZpos));
  pvBot.addPhysVolID("TAUICBottom", 0);
  if (OverlapTest) pvBot.ptr()->CheckOverlaps();

  // Out pipe end (+Z side)
  Position ipcbppos(TAUIChamberXpos - TAUIChamberXpos,
                    TAUIChamberYpos - TAUIChamberYpos,
                    (TAUIChamberZ + TAUIChamberPipeLength)/2.0 + TAUIChamberZpos);
  // simplifies to:
  ipcbppos = Position(0.0, 0.0, (TAUIChamberZ + TAUIChamberPipeLength)/2.0 + TAUIChamberZpos);
  PlacedVolume pvBPOut    = envelope.placeVolume(logicTAUIChamberBPipeOut,    ipcbppos);
  PlacedVolume pvBPOutVac = envelope.placeVolume(logicTAUIChamberBPipeOutVac, ipcbppos);
  pvBPOut.addPhysVolID("TAUIChamberBPipeOut",    0);
  pvBPOutVac.addPhysVolID("TAUIChamberBPipeOutVac", 0);
  if (OverlapTest) { pvBPOut.ptr()->CheckOverlaps(); pvBPOutVac.ptr()->CheckOverlaps(); }

  // In pipe end (-Z side)
  ipcbppos = Position(0.0, 0.0, -(TAUIChamberZ + TAUIChamberPipeLength)/2.0 + TAUIChamberZpos);
  PlacedVolume pvBP    = envelope.placeVolume(logicTAUIChamberBPipe,    ipcbppos);
  PlacedVolume pvBPVac = envelope.placeVolume(logicTAUIChamberBPipeVac, ipcbppos);
  pvBP.addPhysVolID("TAUIChamberBPipe",    0);
  pvBPVac.addPhysVolID("TAUIChamberBPipeVac", 0);
  if (OverlapTest) { pvBP.ptr()->CheckOverlaps(); pvBPVac.ptr()->CheckOverlaps(); }

  // Flange (-Z side, further out)
  ipcbppos = Position(0.0, 0.0,
                      -((TAUIChamberZ - TAUIChamberPipeFlangedL)/2.0 + TAUIChamberPipeLength)
                      + TAUIChamberZpos);
  PlacedVolume pvFlange = envelope.placeVolume(logicTAUIChamberBPipeFlange, ipcbppos);
  pvFlange.addPhysVolID("TAUIChamberBPipeFlange", 0);
  if (OverlapTest) pvFlange.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Support for IC internal breadboard  (LxAux::BuildTable equivalent)
  // -----------------------------------------------------------------------
  double tblhight = -TAUIChamberBeamBottomY - ICBBoardTableGapY - FloorSurfaceYpos;
  Assembly icbboardtable = LxAux::BuildTable(description, "ICBBoardTable",
                                             ICBBoardTableX, tblhight, ICBBoardTableZ,
                                             4, ICBBoardTableY, ICBBoardTableLegR);
  double ylevel = TAUIChamberY/2.0 + TAUIChamberBottomY + ICBBoardTableGapY;
  Position trbbtable(TAUIChamberXpos, -ylevel + TAUIChamberYpos, TAUIChamberZpos);
//   LxAux::AddAssemblyVolumes(envelope, icbboardtable, trbbtable);
  PlacedVolume pvICBBoardTable = envelope.placeVolume(icbboardtable, trbbtable);

  Assembly icbboardSupport = LxAux::BuildTable(description, "ICBBoardSupport",
                                               ICBBoardTableX, ICBBoardTableGapY, ICBBoardTableZ,
                                               4, 0.0);
  double yslevel = TAUIChamberY/2.0 + TAUIChamberBottomY;
  Position trbbsupport(TAUIChamberXpos, -yslevel + TAUIChamberYpos, TAUIChamberZpos);
//   LxAux::AddAssemblyVolumes(envelope, icbboardSupport, trbbsupport);
  PlacedVolume pvICBBoardSupport = envelope.placeVolume(icbboardSupport, trbbsupport);
}


// ---------------------------------------------------------------------------
void LxInteractionChamber::ConstructMirrors(dd4hep::Detector& description,
                                            dd4hep::Volume&   logicTAUICContainer)
{
  Material tauICBBoardMaterial = description.material(
      description.constant<std::string>("TAUIChamberBBoardMaterial"));

  double TAUIChamberBBoardX    = description.constant<double>("TAUIChamberBBoardX");
  double TAUIChamberBBoardY    = description.constant<double>("TAUIChamberBBoardY");
  double TAUIChamberBBoardZ    = description.constant<double>("TAUIChamberBBoardZ");
  double TAUIChamberBBoardGapY = description.constant<double>("TAUIChamberBBoardGapY");
  double TAUIChamberY          = description.constant<double>("TAUIChamberY");
  double TAUIChamberYpos       = description.constant<double>("TAUIChamberYpos");
  double TAUIChamberXpos       = description.constant<double>("TAUIChamberXpos");
  double TAUIChamberZpos       = description.constant<double>("TAUIChamberZpos");
  bool   OverlapTest           = (description.constant<int>("OverlapTest") != 0);

  // Breadboard
  Box    solidTAUICBBoard(TAUIChamberBBoardX/2.0, TAUIChamberBBoardY/2.0, TAUIChamberBBoardZ/2.0);
  Volume logicTAUICBBoard("logicTAUICBBoard", solidTAUICBBoard, tauICBBoardMaterial);
  double bbposy = 0.5*(TAUIChamberY - TAUIChamberBBoardY) - TAUIChamberBBoardGapY;
  PlacedVolume pvBB = logicTAUICContainer.placeVolume(logicTAUICBBoard,
                          Position(0.0, -bbposy, 0.0));
  pvBB.addPhysVolID("ICBBoard", 0);
  if (OverlapTest) pvBB.ptr()->CheckOverlaps();

  // Mirrors
  int nmirrors = description.constant<int>("TAUIChamberNMirrors");
  double mh = -TAUIChamberYpos + TAUIChamberY/2.0
              - TAUIChamberBBoardY - TAUIChamberBBoardGapY;

  for (int mi = 0; mi < nmirrors; ++mi) {
    std::string idx   = std::to_string(mi);
    double r          = description.constant<double>("TAUIChamberMirrorR_"    + idx);
    double d          = description.constant<double>("TAUIChamberMirrorD_"    + idx);
    double mxpos      = description.constant<double>("TAUIChamberMirrorXpos_" + idx);
    double mzpos      = description.constant<double>("TAUIChamberMirrorZpos_" + idx);
    double mtheta     = description.constant<double>("TAUIChamberMirrorTheta_"+ idx);
    std::string mname = description.constant<std::string>("TAUIChamberMirrorName_" + idx);

    Assembly mirrorAsm = ConstructMirrorAssembly(description, r, d, mh, mname);

    Position mpos(mxpos - TAUIChamberXpos, -TAUIChamberYpos, mzpos - TAUIChamberZpos);
    RotationZYX mrot(0.0, mtheta, 0.0);
    PlacedVolume pvM = logicTAUICContainer.placeVolume(mirrorAsm, Transform3D(mrot, mpos));
    pvM.addPhysVolID("ICMirror", mi);
    if (OverlapTest) pvM.ptr()->CheckOverlaps();
  }
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxInteractionChamber::ConstructMirrorAssembly(
    dd4hep::Detector&  description,
    double r, double d, double h,
    const std::string& mname)
{
  Material tauICMirrorMaterial     = description.material(description.constant<std::string>("TAUIChamberMirrorMaterial"));
  Material tauICMirrorRingMaterial = description.material(description.constant<std::string>("TAUIChamberMirrorRingMaterial"));
  Material tauICMirrorHoldMaterial = description.material(description.constant<std::string>("TAUIChamberMirrorHolderMaterial"));

  double rh = description.constant<double>("TAUIChamberMirrorRingH");

  std::string lgc = "logic";

  Tube   solidICMirror(0.0, r-rh, d/2.0, 0.0, 2.0*M_PI);
  Volume logicICMirror(lgc+"ICMirror"+mname, solidICMirror, tauICMirrorMaterial);

  Tube   solidICMirrorRing(r-rh, r, d/2.0, 0.0, 4.0*M_PI/3.0);
  Volume logicICMirrorRing(lgc+"ICMirrorRing"+mname, solidICMirrorRing, tauICMirrorRingMaterial);

  Tube   solidICMirrorHold(0.0, d/2.0, (h-r)/2.0, 0.0, 2.0*M_PI);
  Volume logicICMirrorHold(lgc+"ICMirrorHold"+mname, solidICMirrorHold, tauICMirrorHoldMaterial);

  Assembly mirrorAssembly(std::string("mirrorAssembly_") + mname);

  mirrorAssembly.placeVolume(logicICMirror, Position(0.0, 0.0, 0.0));

  mirrorAssembly.placeVolume(logicICMirrorRing,
      Transform3D(RotationZYX(4.0*M_PI/3.0, 0.0, 0.0), Position(0.0, 0.0, 0.0)));

  // G4RotationMatrix(G4ThreeVector(1,0,0), pi/2): axis-angle, rotation around X by pi/2
  mirrorAssembly.placeVolume(logicICMirrorHold,
      Transform3D(RotationZYX(0.0, 0.0, M_PI/2.0), Position(0.0, -(h+r)/2.0, 0.0)));

  return mirrorAssembly;
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxInteractionChamber(dd4hep::Detector& description,
                                         xml_h             e,
                                         dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxInteractionChamber lxic;
  lxic.Construct(description, sdet, e);

  return sdet;
}

DECLARE_DETELEMENT(LxInteractionChamber, create_LxInteractionChamber)
