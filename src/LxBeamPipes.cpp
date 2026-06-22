//
/// \brief Implementation of the LxBeamPipes class (DD4hep version)
//

#include <cmath>
#include <stdexcept>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoManager.h"

#include "LxBeamPipes.h"

using namespace dd4hep;


LxBeamPipes::LxBeamPipes(dd4hep::Detector& description, dd4hep::Volume& motherVol)
  : fDescription(description), fMotherVol(motherVol)
{
  fFunctionMap["BeamPipeToDump"]      = &LxBeamPipes::ConstructBeamPipeToDump;
  fFunctionMap["BeamPipeToIP"]        = &LxBeamPipes::ConstructBeamPipeToIP;
  fFunctionMap["GammaVacuumChamber"]  = &LxBeamPipes::ConstructGammaVacuumChamber;
  fFunctionMap["BeamPipeTM"]          = &LxBeamPipes::ConstructBeamPipeTM;
  fFunctionMap["BeamPipeInc"]         = &LxBeamPipes::ConstructBeamPipeInc;
  fFunctionMap["BeamPipeOPPPDGT"]     = &LxBeamPipes::ConstructBeamPipeOPPPDGT;
}


void LxBeamPipes::Construct(const std::string& pipeName)
{
  // Check TGeoManager — skip if this section was already placed
  const std::string asmName = pipeName + "_beampipe";
  TGeoVolume* existing = fDescription.manager().GetVolume(asmName.c_str());
  if (existing) {
    printout(INFO, "LxBeamPipes::Construct",
             "Beam pipe section '%s' already placed, skipping.", pipeName.c_str());
    return;
  }

  auto it = fFunctionMap.find(pipeName);
  if (it == fFunctionMap.end())
    throw std::runtime_error("LxBeamPipes: unknown pipe section '" + pipeName + "'");

  (this->*(it->second))();
}


// ---------------------------------------------------------------------------
void LxBeamPipes::ConstructBeamPipeToDump()
{
  Material beamPipeMaterial     = fDescription.material(fDescription.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial       = fDescription.material(fDescription.constant<std::string>("BeamPipeVacuumMaterial"));
  Material beamPipeDumpMaterial = fDescription.material(fDescription.constant<std::string>("BeamDumpPipeMaterial"));

  double BPipeR             = fDescription.constant<double>("BPipeR");
  double BPipeThickness     = fDescription.constant<double>("BPipeThickness");
  double BeamPipeSplitterZ  = fDescription.constant<double>("BeamPipeSplitterZ");
  double DumpMagnetZpos     = fDescription.constant<double>("DumpMagnetZpos");
  double BeamDumpZpos       = fDescription.constant<double>("BeamDumpZpos");
  double BeamDumpXpos       = fDescription.constant<double>("BeamDumpXpos");
  double BeamDumpZ          = fDescription.constant<double>("BeamDumpZ");
  double BeamDumpR          = fDescription.constant<double>("BeamDumpR");
  double BeamDumpAngle      = fDescription.constant<double>("BeamDumpAngle");
  double magnetz            = fDescription.constant<double>("FlashMFieldLength");
  bool   OverlapTest        = (fDescription.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Beam splitter box
  // -----------------------------------------------------------------------
  Box solidBeamSplitContainer(2.0*(BPipeR + BPipeThickness), BPipeR, BeamPipeSplitterZ/2.0);
  Volume logicBeamSplitContainer("logicBeamSplitContainer", solidBeamSplitContainer, vacuumMaterial);

  Box solidBeamSplitOuter(2.0*(BPipeR + BPipeThickness), BPipeR, BeamPipeSplitterZ/2.0);
  Box solidBeamSplitInner(2.0*BPipeR + BPipeThickness,
                          BPipeR - BPipeThickness,
                          BeamPipeSplitterZ/2.0 - BPipeThickness);
  SubtractionSolid solidBeamSplitXY("solidBeamSplitXY", solidBeamSplitOuter, solidBeamSplitInner);

  Tube solidSHoleIn(0.0, BPipeR - BPipeThickness, BPipeThickness, 0.0, 2.0*M_PI);

  // Commented: round hole for income pipe
  // G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0,0,1), 0),
  //     G4ThreeVector(BPipeR+2*BPipeThickness, 0, -(BeamPipeSplitterZ-BPipeThickness)/2));
  // SubtractionSolid solidBeamSplit1("solidBeamSplit1", solidBeamSplitXY, solidSHoleIn,
  //     Position(BPipeR+2.0*BPipeThickness, 0.0, -(BeamPipeSplitterZ-BPipeThickness)/2.0));

  // Rectangular hole for income pipe
  Box solidInCut(2.0*BPipeR + BPipeThickness, BPipeR - BPipeThickness, BPipeThickness);
  SubtractionSolid solidBeamSplit1("solidBeamSplit1", solidBeamSplitXY, solidInCut,
      Position(0.0, 0.0, -(BeamPipeSplitterZ - BPipeThickness)/2.0));

  // Round hole for +X exit pipe
  SubtractionSolid solidBeamSplit2("solidBeamSplit2", solidBeamSplit1, solidSHoleIn,
      Position(BPipeR + 2.0*BPipeThickness, 0.0, (BeamPipeSplitterZ - BPipeThickness)/2.0));

  // Round hole for -X exit pipe
  SubtractionSolid solidBeamSplit("solidBeamSplit", solidBeamSplit2, solidSHoleIn,
      Position(-(BPipeR + 2.0*BPipeThickness), 0.0, (BeamPipeSplitterZ - BPipeThickness)/2.0));

  Volume logicBeamSplit("logicBeamSplit", solidBeamSplit, beamPipeMaterial);

  PlacedVolume pvSplit = logicBeamSplitContainer.placeVolume(logicBeamSplit, Position(0.0, 0.0, 0.0));
  pvSplit.addPhysVolID("BeamSplit", 0);
  if (OverlapTest) pvSplit.ptr()->CheckOverlaps();

  // Commented placement (original alternative):
  // fMotherVol.placeVolume(logicBeamSplitContainer,
  //     Position(-(BPipeR+2*BPipeThickness), 0,
  //              DumpMagnetZpos + 0.5*magnetz + BeamPipeSplitterZ/2.0));

  // G4RotationMatrix(G4ThreeVector(0,0,1), pi/2): axis-angle around Z → RotationZYX(pi/2, 0, 0)
  PlacedVolume pvSplitCont = fMotherVol.placeVolume(logicBeamSplitContainer,
      Transform3D(RotationZYX(M_PI/2.0, 0.0, 0.0),
                  Position(0.0, -(BPipeR + 2.0*BPipeThickness),
                           DumpMagnetZpos + 0.5*magnetz + BeamPipeSplitterZ/2.0)));
  pvSplitCont.addPhysVolID("BeamSplitContainer", 0);
  if (OverlapTest) pvSplitCont.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Straight pipe section (beam pipe after magnet, before dump)
  // -----------------------------------------------------------------------
  double lpipe = BeamDumpZpos - DumpMagnetZpos - 0.5*magnetz
                 - 0.5*BeamDumpZ*std::cos(BeamDumpAngle)
                 - BeamDumpR*std::sin(BeamDumpAngle) - BeamPipeSplitterZ;

  Tube   solidStrightPipe(BPipeR - BPipeThickness, BPipeR, lpipe/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeMB("logicBeamPipeMB", solidStrightPipe, beamPipeMaterial);
  PlacedVolume pvMB = fMotherVol.placeVolume(logicBeamPipeMB,
      Position(0.0, 0.0, DumpMagnetZpos + 0.5*magnetz + BeamPipeSplitterZ + lpipe/2.0));
  pvMB.addPhysVolID("BeamPipeMB", 0);
  if (OverlapTest) pvMB.ptr()->CheckOverlaps();

  Tube   solidStrightPipeVac(0.0, BPipeR - BPipeThickness, lpipe/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeMBVac("logicBeamPipeMBVac", solidStrightPipeVac, vacuumMaterial);
  PlacedVolume pvMBVac = fMotherVol.placeVolume(logicBeamPipeMBVac,
      Position(0.0, 0.0, DumpMagnetZpos + 0.5*magnetz + BeamPipeSplitterZ + lpipe/2.0));
  pvMBVac.addPhysVolID("BeamPipeMBVac", 0);
  if (OverlapTest) pvMBVac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Angled pipe section to beam dump
  // -----------------------------------------------------------------------
  double dxdump    = BeamDumpXpos - 0.5*BeamDumpZ * std::tan(BeamDumpAngle);
  double dxsplit   = 2.0*(BPipeR + 2.0*BPipeThickness);
  double zpipedump = BeamDumpZpos - DumpMagnetZpos - 0.5*magnetz
                     - 0.5*BeamDumpZ*std::cos(BeamDumpAngle) - BeamPipeSplitterZ;
  double lpipedump = std::sqrt(std::pow(zpipedump, 2.0) + std::pow(dxdump - dxsplit, 2.0));
  double dumpipeangle = std::atan2(dxdump - dxsplit, zpipedump);

  // Commented alternative (XZ plane tilt):
  // CutTube solidPipeToDump(BPipeR-BPipeThickness, BPipeR, lpipedump/2.0, 0.0, 2.0*M_PI,
  //     ROOT::Math::XYZVector(-(dxdump-dxsplit)/lpipedump, 0.0, -zpipedump/lpipedump),
  //     ROOT::Math::XYZVector(std::sin(dumpipeangle-BeamDumpAngle), 0.0, std::cos(dumpipeangle-BeamDumpAngle)));

  // YZ plane tilt (active version)
  // dd4hep::CutTube(rmin, rmax, hz, startPhi, endPhi, lx, ly, lz, tx, ty, tz)
  // low normal:  (0, -(dxdump-dxsplit)/lpipedump, -zpipedump/lpipedump)
  // high normal: (0,  sin(dumpipeangle-BeamDumpAngle), cos(dumpipeangle-BeamDumpAngle))
  CutTube solidPipeToDump(BPipeR - BPipeThickness, BPipeR, lpipedump/2.0, 0.0, 2.0*M_PI,
      0.0, -(dxdump - dxsplit)/lpipedump, -zpipedump/lpipedump,
      0.0,  std::sin(dumpipeangle - BeamDumpAngle), std::cos(dumpipeangle - BeamDumpAngle));
  Volume logicBeamPipeMD("logicBeamPipeMD", solidPipeToDump, beamPipeDumpMaterial);

  // Commented alternative rotation (around Y):
  // fMotherVol.placeVolume(logicBeamPipeMD,
  //     Transform3D(RotationZYX(0.0, std::atan((dxdump-dxsplit)/zpipedump), 0.0),
  //                 Position(-(dxdump+dxsplit)/2.0, 0.0,
  //                          DumpMagnetZpos+0.5*magnetz+BeamPipeSplitterZ+zpipedump/2.0)));

  // G4RotationMatrix(G4ThreeVector(-1,0,0), angle): axis-angle around -X → RotationZYX(0, 0, -angle)
  PlacedVolume pvMD = fMotherVol.placeVolume(logicBeamPipeMD,
      Transform3D(RotationZYX(0.0, 0.0, std::atan((dxdump - dxsplit)/zpipedump)),
                  Position(0.0, -(dxdump + dxsplit)/2.0,
                           DumpMagnetZpos + 0.5*magnetz + BeamPipeSplitterZ + zpipedump/2.0)));
  pvMD.addPhysVolID("BeamPipeMD", 0);
  if (OverlapTest) pvMD.ptr()->CheckOverlaps();

  // Commented alternative rotation (around Y):
  // fMotherVol.placeVolume(logicBeamPipeMDVac,
  //     Transform3D(RotationZYX(0.0, std::atan((dxdump-dxsplit)/zpipedump), 0.0),
  //                 Position(-(dxdump+dxsplit)/2.0, 0.0,
  //                          DumpMagnetZpos+0.5*magnetz+BeamPipeSplitterZ+zpipedump/2.0)));

  CutTube solidPipeToDumpVac(0.0, BPipeR - BPipeThickness, lpipedump/2.0, 0.0, 2.0*M_PI,
      0.0, -(dxdump - dxsplit)/lpipedump, -zpipedump/lpipedump,
      0.0,  std::sin(dumpipeangle - BeamDumpAngle), std::cos(dumpipeangle - BeamDumpAngle));
  Volume logicBeamPipeMDVac("logicBeamPipeMDVac", solidPipeToDumpVac, vacuumMaterial);

  PlacedVolume pvMDVac = fMotherVol.placeVolume(logicBeamPipeMDVac,
      Transform3D(RotationZYX(0.0, 0.0, std::atan((dxdump - dxsplit)/zpipedump)),
                  Position(0.0, -(dxdump + dxsplit)/2.0,
                           DumpMagnetZpos + 0.5*magnetz + BeamPipeSplitterZ + zpipedump/2.0)));
  pvMDVac.addPhysVolID("BeamPipeMDVac", 0);
  if (OverlapTest) pvMDVac.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxBeamPipes::ConstructBeamPipeToIP()
{
  Material beamPipeMaterial = fDescription.material(fDescription.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial   = fDescription.material(fDescription.constant<std::string>("BeamPipeVacuumMaterial"));

  double BPipeR               = fDescription.constant<double>("BPipeR");
  double BPipeThickness       = fDescription.constant<double>("BPipeThickness");
  double magnetz              = fDescription.constant<double>("FlashMFieldLength");
  double IPContainerZ         = fDescription.constant<double>("IPContainerZ");
  double ShieldingZ           = fDescription.constant<double>("ShieldingZ");
  double TAUIChamberZpos      = fDescription.constant<double>("TAUIChamberZpos");
  double TAUIChamberElPipeRIn = fDescription.constant<double>("TAUIChamberElPipeRIn");
  double TAUIChamberElPipeROut= fDescription.constant<double>("TAUIChamberElPipeROut");
  double IPMagnetZpos         = fDescription.constant<double>("IPMagnetZpos");
  double FlashMFieldX         = fDescription.constant<double>("FlashMFieldX");
  double FlashMagneteffY      = fDescription.constant<double>("FlashMagneteffY");
  double IPMAgnetBeamPipeXGap = fDescription.constant<double>("IPMAgnetBeamPipeXGap");
  bool   OverlapTest          = (fDescription.constant<int>("OverlapTest") != 0);

  // -----------------------------------------------------------------------
  // Pipe from shielding to IP chamber (conical, expanding toward IP)
  // -----------------------------------------------------------------------

  double BeamDumpR = fDescription.constant<double>("BeamDumpR");
  double BeamDumpAngle = fDescription.constant<double>("BeamDumpAngle");
  double BeamDumpZ = fDescription.constant<double>("BeamDumpZ");
  double ShieldingDeepMargine = fDescription.constant<double>("ShieldingDeepMargine");
  double ShieldingDeepZ = fDescription.constant<double>("ShieldingDeepZ");
  double BeamDumpZpos = fDescription.constant<double>("BeamDumpZpos");

  double ddd = 0.5*(2.0*BeamDumpR * sin(BeamDumpAngle) + BeamDumpZ * cos(BeamDumpAngle))
                 + ShieldingDeepMargine;
  double zdeeppos =  0.5 * (ShieldingDeepZ - ShieldingZ);
  double ShieldingZpos = BeamDumpZpos - ShieldingDeepZ + 0.5*ShieldingZ + ddd;
  double lpipe_sip = -ShieldingZpos - IPContainerZ/2.0 - ShieldingZ/2.0 + TAUIChamberZpos;

  // Commented: straight tube version
  // Tube solidBeamPipeSIP(BPipeR-BPipeThickness, BPipeR, lpipe_sip/2.0, 0.0, 2.0*M_PI);

  // Conical pipe (Cons): rmin1/rmax1 at -z face, rmin2/rmax2 at +z face
  Cone solidBeamPipeSIP(lpipe_sip/2.0,
                        BPipeR - BPipeThickness, BPipeR,
                        TAUIChamberElPipeRIn, TAUIChamberElPipeROut);
  Volume logicBeamPipeSIP("logicBeamPipeSIP", solidBeamPipeSIP, beamPipeMaterial);
  PlacedVolume pvSIP = fMotherVol.placeVolume(logicBeamPipeSIP,
      Position(0.0, 0.0, -(lpipe_sip + IPContainerZ)/2.0 + TAUIChamberZpos));
  pvSIP.addPhysVolID("BeamPipeSIP", 0);
  if (OverlapTest) pvSIP.ptr()->CheckOverlaps();

  // Commented: straight vacuum version
  // Tube solidBeamPipeSIPVac(0.0, BPipeR-BPipeThickness, lpipe_sip/2.0, 0.0, 2.0*M_PI);

  Cone solidBeamPipeSIPVac(lpipe_sip/2.0,
                           0.0, BPipeR - BPipeThickness,
                           0.0, TAUIChamberElPipeRIn);
  Volume logicBeamPipeSIPVac("logicBeamPipeSIPVac", solidBeamPipeSIPVac, vacuumMaterial);
  PlacedVolume pvSIPVac = fMotherVol.placeVolume(logicBeamPipeSIPVac,
      Position(0.0, 0.0, -(lpipe_sip + IPContainerZ)/2.0 + TAUIChamberZpos));
  pvSIPVac.addPhysVolID("BeamPipeSIPVac", 0);
  if (OverlapTest) pvSIPVac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Rectangular pipe from IP chamber to IP magnet
  // -----------------------------------------------------------------------
  double lpipe_ipm = IPMagnetZpos - magnetz/2.0 - IPContainerZ/2.0 - TAUIChamberZpos;
  double ipmbphx   = FlashMFieldX/2.0 - IPMAgnetBeamPipeXGap;

  // Commented: round tube version
  // Tube solidBeamPipeIPM(BPipeR-BPipeThickness, BPipeR, lpipe_ipm/2.0, 0.0, 2.0*M_PI);

  Box solidBeamPipeIPM1(ipmbphx, FlashMagneteffY/2.0, lpipe_ipm/2.0);
  Box solidBeamPipeIPMCut(ipmbphx - BPipeThickness,
                          FlashMagneteffY/2.0 - BPipeThickness,
                          lpipe_ipm);
  SubtractionSolid solidBeamPipeIPM("solidBeamPipeIPM", solidBeamPipeIPM1, solidBeamPipeIPMCut);
  Volume logicBeamPipeIPM("logicBeamPipeIPM", solidBeamPipeIPM, beamPipeMaterial);
  PlacedVolume pvIPM = fMotherVol.placeVolume(logicBeamPipeIPM,
      Position(0.0, 0.0, (lpipe_ipm + IPContainerZ)/2.0 + TAUIChamberZpos));
  pvIPM.addPhysVolID("BeamPipeIPM", 0);
  if (OverlapTest) pvIPM.ptr()->CheckOverlaps();

  // Commented: round vacuum version
  // Tube solidBeamPipeIPMVac(0.0, BPipeR-BPipeThickness, lpipe_ipm/2.0, 0.0, 2.0*M_PI);

  Box solidBeamPipeIPMVac(ipmbphx - BPipeThickness,
                          FlashMagneteffY/2.0 - BPipeThickness,
                          lpipe_ipm/2.0);
  Volume logicBeamPipeIPMVac("logicBeamPipeIPMVac", solidBeamPipeIPMVac, vacuumMaterial);
  PlacedVolume pvIPMVac = fMotherVol.placeVolume(logicBeamPipeIPMVac,
      Position(0.0, 0.0, (lpipe_ipm + IPContainerZ)/2.0 + TAUIChamberZpos));
  pvIPMVac.addPhysVolID("BeamPipeIPMVac", 0);
  if (OverlapTest) pvIPMVac.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
void LxBeamPipes::ConstructGammaVacuumChamber()
{
  Material beamPipeMaterial   = fDescription.material(fDescription.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial     = fDescription.material(fDescription.constant<std::string>("BeamPipeVacuumMaterial"));
  Material bpipeWindowMaterial= fDescription.material(fDescription.constant<std::string>("BeamPipeWindowMaterial"));
  Material bpipeLidMaterial   = fDescription.material(fDescription.constant<std::string>("BeamPipeLidMaterial"));

  double BPipeThickness          = fDescription.constant<double>("BPipeThickness");
  double ComptonLysoZpos         = fDescription.constant<double>("ComptonLysoZpos");
  double ComptonLysoZ            = fDescription.constant<double>("ComptonLysoZ");
  double ComptonLysoX            = fDescription.constant<double>("ComptonLysoX");
  double ComptonLysoXpos         = fDescription.constant<double>("ComptonLysoXpos");
  double GMagnetZpos             = fDescription.constant<double>("GMagnetZpos");
  double GamVacChamberGap        = fDescription.constant<double>("GamVacChamberGap");
  double BeamPipe2ProfilerZ      = fDescription.constant<double>("BeamPipe2ProfilerZ");
  double BeamPipeLidThickness    = fDescription.constant<double>("BeamPipeLidThickness");
  bool   OverlapTest             = (fDescription.constant<int>("OverlapTest") != 0);

  // Variables shared with GammaMagnet geometry
  double gmpipex    = fDescription.constant<double>("QBeamPipeContainerX");
  double gmpiptby   = fDescription.constant<double>("GChamberTBWallThickness");
  double gmpipsidex = fDescription.constant<double>("GammaBPipeWindowThickness");
  double gfieldz    = fDescription.constant<double>("FlashMFieldLength") / 2.0;
  double gfieldy    = fDescription.constant<double>("FlashMagneteffY")   / 2.0;
  double gfieldx    = fDescription.constant<double>("FlashMFieldX")      / 2.0;
  double GammaMagnetBeamPipeXGap = fDescription.constant<double>("GammaMagnetBeamPipeXGap");

  double wpipecontainerx = gfieldx - gmpipex/2.0 - GammaMagnetBeamPipeXGap;
  double swy             = gfieldy - gmpiptby;

  double vclength = ComptonLysoZpos - GMagnetZpos - ComptonLysoZ/2.0 - gfieldz - GamVacChamberGap;
  double magnetx  = wpipecontainerx + gmpipex/2.0;
  double detx     = ComptonLysoX + gmpipsidex;

  // -----------------------------------------------------------------------
  // Trapezoidal vacuum chamber container
  // -----------------------------------------------------------------------
  // G4Trd(dx1, dx2, dy1, dy2, dz): half-lengths; dx1/dy1 at -z, dx2/dy2 at +z
  // dd4hep::Trapezoid(dx1, dx2, dy1, dy2, dz) — same convention
  Trapezoid solidGVCContainer(magnetx, detx, gfieldy, gfieldy, vclength/2.0);
  Volume    logicGVCContainer("logicGVCContainer", solidGVCContainer, vacuumMaterial);

  // Top and bottom trapezoidal walls
  Trapezoid solidGamChamTopBot(magnetx, detx, gmpiptby/2.0, gmpiptby/2.0, vclength/2.0);
  Volume    logicGamChamTopBot("logicGamChamTopBot", solidGamChamTopBot, beamPipeMaterial);

  PlacedVolume pvTop = logicGVCContainer.placeVolume(logicGamChamTopBot,
      Position(0.0, gfieldy - gmpiptby/2.0, 0.0));
  pvTop.addPhysVolID("GamChamTop", 0);
  if (OverlapTest) pvTop.ptr()->CheckOverlaps();

  PlacedVolume pvBot = logicGVCContainer.placeVolume(logicGamChamTopBot,
      Position(0.0, -(gfieldy - gmpiptby/2.0), 0.0));
  pvBot.addPhysVolID("GamChamBottom", 1);
  if (OverlapTest) pvBot.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Exit window (flat, at +Z face, with beam pipe cutout)
  // -----------------------------------------------------------------------
  double wtheta   = std::atan2(detx - magnetx, vclength);
  double vcwndz   = gmpipsidex;
  double vcwallx  = detx - vcwndz * std::tan(wtheta);
  double BPipeRProf = ComptonLysoXpos - ComptonLysoX/2.0;

  Box  solidVCWindow1(vcwallx, swy, vcwndz/2.0);
  Tube solidCVBeamPipeHole(0.0, BPipeRProf - BPipeThickness, vcwndz, 0.0, 2.0*M_PI);
  SubtractionSolid solidVCWindow("solidVCWindow", solidVCWindow1, solidCVBeamPipeHole);
  Volume logicGVCWindow("logicGVCWindow", solidVCWindow, bpipeWindowMaterial);

  PlacedVolume pvWnd = logicGVCContainer.placeVolume(logicGVCWindow,
      Position(0.0, 0.0, (vclength - vcwndz)/2.0));
  pvWnd.addPhysVolID("GVCWindow", 0);
  if (OverlapTest) pvWnd.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Side walls (G4Para → dd4hep::ParallelWorldPlacement not available;
  // use EightPointSolid or Para — DD4hep wraps ROOT TGeoPara)
  // G4Para(dx, dy, dz, alpha, theta, phi):
  //   alpha=0, theta=wtheta, phi=0 → TGeoPara(dx, dy, dz, alpha_deg, theta_deg, phi_deg)
  // -----------------------------------------------------------------------
  double swdx = 0.5*gmpipsidex / std::cos(wtheta);

  TGeoPara* tgeoPara = new TGeoPara("solidGamChamSide",
    swdx/2.0 / dd4hep::cm,      // ROOT works in cm
    swy      / dd4hep::cm,
    (vclength - vcwndz)/2.0 / dd4hep::cm,
    0.0,                          // alpha in degrees
    wtheta * 180.0 / M_PI,        // theta in degrees
    0.0);                          // phi in degrees
  Solid solidGamChamSide(tgeoPara);
  Volume logicGamChamSideWall("logicGamChamSideWall", solidGamChamSide, beamPipeMaterial);
// //   // dd4hep::Para(dx, dy, dz, alpha, theta, phi) — angles in radians
// //   Para   solidGamChamSide(swdx/2.0, swy, (vclength - vcwndz)/2.0, 0.0, wtheta, 0.0);
// //   Volume logicGamChamSideWall("logicGamChamSideWall", solidGamChamSide, beamPipeMaterial);

  PlacedVolume pvSide0 = logicGVCContainer.placeVolume(logicGamChamSideWall,
      Position(0.5*(vclength - vcwndz)*std::tan(wtheta) - swdx/2.0 + magnetx,
               0.0, -vcwndz/2.0));
  pvSide0.addPhysVolID("GamChamSide", 0);
  if (OverlapTest) pvSide0.ptr()->CheckOverlaps();

  // G4RotationMatrix(G4ThreeVector(0,0,1), pi): axis-angle around Z → RotationZYX(pi, 0, 0)
  PlacedVolume pvSide1 = logicGVCContainer.placeVolume(logicGamChamSideWall,
      Transform3D(RotationZYX(M_PI, 0.0, 0.0),
                  Position(-0.5*(vclength - vcwndz)*std::tan(wtheta) + swdx/2.0 - magnetx,
                           0.0, -vcwndz/2.0)));
  pvSide1.addPhysVolID("GamChamSide", 1);
  if (OverlapTest) pvSide1.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Place vacuum chamber container into world
  // -----------------------------------------------------------------------
  double chpos = GMagnetZpos + vclength/2.0 + gfieldz;
  PlacedVolume pvCont = fMotherVol.placeVolume(logicGVCContainer,
      Position(0.0, 0.0, chpos));
  pvCont.addPhysVolID("GMagnetWideChamber", 0);
  if (OverlapTest) pvCont.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Short beam pipe joining chamber to photon detectors
  // -----------------------------------------------------------------------
  double lbpipe = BeamPipe2ProfilerZ;
  double bppos  = chpos + (vclength + lbpipe)/2.0;

  Tube   solidBeamPipeVCG(BPipeRProf - BPipeThickness, BPipeRProf, 0.5*lbpipe, 0.0, 2.0*M_PI);
  Volume logicBeamPipeVCG("logicBeamPipeVCG", solidBeamPipeVCG, beamPipeMaterial);
  PlacedVolume pvVCG = fMotherVol.placeVolume(logicBeamPipeVCG, Position(0.0, 0.0, bppos));
  pvVCG.addPhysVolID("BeamPipeVCGProf", 0);
  if (OverlapTest) pvVCG.ptr()->CheckOverlaps();

  Tube   solidBeamPipeVCGVac(0.0, BPipeRProf - BPipeThickness, 0.5*lbpipe, 0.0, 2.0*M_PI);
  Volume logicBeamPipeVCGVac("logicBeamPipeVCGVac", solidBeamPipeVCGVac, vacuumMaterial);
  PlacedVolume pvVCGVac = fMotherVol.placeVolume(logicBeamPipeVCGVac, Position(0.0, 0.0, bppos));
  pvVCGVac.addPhysVolID("BeamPipeVCGProfVac", 0);
  if (OverlapTest) pvVCGVac.ptr()->CheckOverlaps();

  // Lid at end of beam pipe
  Tube   solidBeamPipeVCGLid(0.0, BPipeRProf, 0.5*BeamPipeLidThickness, 0.0, 2.0*M_PI);
  Volume logicBeamPipeVCGLid("logicBeamPipeVCGLid", solidBeamPipeVCGLid, bpipeLidMaterial);
  PlacedVolume pvLid = fMotherVol.placeVolume(logicBeamPipeVCGLid,
      Position(0.0, 0.0, bppos + (lbpipe + BeamPipeLidThickness)/2.0));
  pvLid.addPhysVolID("BeamPipeVCGProfLid", 0);
  if (OverlapTest) pvLid.ptr()->CheckOverlaps();
}


void LxBeamPipes::ConstructBeamPipeTM()
{
  Material beamPipeMaterial = fDescription.material(
      fDescription.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial   = fDescription.material(
      fDescription.constant<std::string>("BeamPipeVacuumMaterial"));

  double BPipeR          = fDescription.constant<double>("BPipeR");
  double BPipeThickness  = fDescription.constant<double>("BPipeThickness");
  double DumpMagnetZpos  = fDescription.constant<double>("DumpMagnetZpos");
  double BTargetZpos     = fDescription.constant<double>("BTargetZpos");
  double TargetChamberZ  = fDescription.constant<double>("TargetChamberZ");
  double FlashMFieldLength = fDescription.constant<double>("FlashMFieldLength");
  bool   OverlapTest     = (fDescription.constant<int>("OverlapTest") != 0);

//   double BTargetChamberMagnetGapZ = fDescription.constant<double>("BTargetChamberMagnetGapZ");
//   double dumpMagnetZ = fDescription.constant<double>("FlashMagnetCoilZ");
//   G4double dumpMagnetZ = lxs->FlashMagnetCoilZ;
//
//   double fZposAbs = DumpMagnetZpos - 0.5*(dumpMagnetZ + TargetChamberZ) - BTargetChamberMagnetGapZ;
//   BTargetZpos = fZposAbs;
//  <constant name="BTargetZpos" value="DumpMagnetZpos - 0.5*(DumpMagnetZ + TargetChamberZ) - BTargetChamberMagnetGapZ"/>


  double lpipe_tm = DumpMagnetZpos - BTargetZpos
                    - 0.5*(TargetChamberZ + FlashMFieldLength);

  Tube   solidBeamPipeTM(BPipeR - BPipeThickness, BPipeR, lpipe_tm/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeTM("logicBeamPipeTM", solidBeamPipeTM, beamPipeMaterial);
  PlacedVolume pvTM = fMotherVol.placeVolume(logicBeamPipeTM,
      Position(0.0, 0.0, BTargetZpos + (TargetChamberZ + lpipe_tm)/2.0));
  pvTM.addPhysVolID("BeamPipeTM", 0);
  if (OverlapTest) pvTM.ptr()->CheckOverlaps();

  Tube   solidBeamPipeTMVac(0.0, BPipeR - BPipeThickness, lpipe_tm/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeTMVac("logicBeamPipeTMVac", solidBeamPipeTMVac, vacuumMaterial);
  PlacedVolume pvTMVac = fMotherVol.placeVolume(logicBeamPipeTMVac,
      Position(0.0, 0.0, BTargetZpos + (TargetChamberZ + lpipe_tm)/2.0));
  pvTMVac.addPhysVolID("BeamPipeTMVac", 0);
  if (OverlapTest) pvTMVac.ptr()->CheckOverlaps();
}


void LxBeamPipes::ConstructBeamPipeInc()
{
  Material beamPipeMaterial = fDescription.material(
      fDescription.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial   = fDescription.material(
      fDescription.constant<std::string>("BeamPipeVacuumMaterial"));

  double BPipeR         = fDescription.constant<double>("BPipeR");
  double BPipeThickness = fDescription.constant<double>("BPipeThickness");
  double BTargetZpos    = fDescription.constant<double>("BTargetZpos");
  double TargetChamberZ = fDescription.constant<double>("TargetChamberZ");
  double WorldSizeZ     = fDescription.constant<double>("world_z");
  bool   OverlapTest    = (fDescription.constant<int>("OverlapTest") != 0);

  double lpipe_inc = BTargetZpos + WorldSizeZ/2.0 - TargetChamberZ/2.0;

  Tube   solidBeamPipeInc(BPipeR - BPipeThickness, BPipeR, lpipe_inc/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeInc("logicBeamPipeInc", solidBeamPipeInc, beamPipeMaterial);
  PlacedVolume pvInc = fMotherVol.placeVolume(logicBeamPipeInc,
      Position(0.0, 0.0, BTargetZpos - (TargetChamberZ + lpipe_inc)/2.0));
  pvInc.addPhysVolID("BeamPipeInc", 0);
  if (OverlapTest) pvInc.ptr()->CheckOverlaps();

  Tube   solidBeamPipeIncVac(0.0, BPipeR - BPipeThickness, lpipe_inc/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeIncVac("logicBeamPipeIncVac", solidBeamPipeIncVac, vacuumMaterial);
  PlacedVolume pvIncVac = fMotherVol.placeVolume(logicBeamPipeIncVac,
      Position(0.0, 0.0, BTargetZpos - (TargetChamberZ + lpipe_inc)/2.0));
  pvIncVac.addPhysVolID("BeamPipeIncVac", 0);
  if (OverlapTest) pvIncVac.ptr()->CheckOverlaps();
}


void LxBeamPipes::ConstructBeamPipeOPPPDGT()
{
  Material beamPipeMaterial = fDescription.material(
      fDescription.constant<std::string>("BeamPipeMaterial"));
  Material vacuumMaterial   = fDescription.material(
      fDescription.constant<std::string>("BeamPipeVacuumMaterial"));

  double BPipeR          = fDescription.constant<double>("BPipeR");
  double BPipeThickness  = fDescription.constant<double>("BPipeThickness");
  double FlashMFieldLength = fDescription.constant<double>("FlashMFieldLength");
  double IPMagnetZpos    = fDescription.constant<double>("IPMagnetZpos");
  double OPPPDetZtoMagnet= fDescription.constant<double>("OPPPDetZtoMagnet");
  double GTargetZpos     = fDescription.constant<double>("GTargetZpos");
  double TargetChamberZ  = fDescription.constant<double>("TargetChamberZ");
  bool   OverlapTest     = (fDescription.constant<int>("OverlapTest") != 0);

  // Commented: alternative magnet field length
  // double dumpMagnetZ = fDescription.constant<double>("TypMBFieldLength");
  double dumpMagnetZ  = FlashMFieldLength;
  double opppdetzpos  = IPMagnetZpos + dumpMagnetZ/2.0 + OPPPDetZtoMagnet;
  double lpipe_dt     = GTargetZpos - opppdetzpos - 0.5*TargetChamberZ;

  Tube   solidBeamPipeOPPPDGT(BPipeR - BPipeThickness, BPipeR, lpipe_dt/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeOPPPDGT("logicBeamPipeOPPPDGT", solidBeamPipeOPPPDGT, beamPipeMaterial);
  PlacedVolume pvDGT = fMotherVol.placeVolume(logicBeamPipeOPPPDGT,
      Position(0.0, 0.0, GTargetZpos - 0.5*(TargetChamberZ + lpipe_dt)));
  pvDGT.addPhysVolID("BeamPipeOPPPDGT", 0);
  if (OverlapTest) pvDGT.ptr()->CheckOverlaps();

  Tube   solidBeamPipeOPPPDGTVac(0.0, BPipeR - BPipeThickness, lpipe_dt/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeOPPPDGTVac("logicBeamPipeOPPPDGTVac", solidBeamPipeOPPPDGTVac, vacuumMaterial);
  PlacedVolume pvDGTVac = fMotherVol.placeVolume(logicBeamPipeOPPPDGTVac,
      Position(0.0, 0.0, GTargetZpos - 0.5*(TargetChamberZ + lpipe_dt)));
  pvDGTVac.addPhysVolID("BeamPipeOPPPDGTVac", 0);
  if (OverlapTest) pvDGTVac.ptr()->CheckOverlaps();
}

