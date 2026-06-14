//
/// \brief Implementation of the LxDetChamberAssembly class (DD4hep version)
//

#include <cmath>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TGeoArb8.h"   // TGeoTrap
#include "TGeoPara.h"

#include "LxAux.h"
#include "LxDetChamberAssembly.h"

using namespace dd4hep;


LxDetChamberAssembly::LxDetChamberAssembly(const std::string& typeName)
  : fChamberType(typeName)
{}


// ---------------------------------------------------------------------------
void LxDetChamberAssembly::Construct(dd4hep::Detector&   description,
                                     dd4hep::DetElement& sdet,
                                     xml_h&              e)
{
  Volume       fLogicWorld = description.worldVolume();
  xml_comp_t   x_det(e);
  Assembly     envelope(x_det.nameStr() + "_assembly");
  PlacedVolume envPV = fLogicWorld.placeVolume(envelope, Transform3D());
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);

  double FlashMagnetCoilZ = description.constant<double>("FlashMagnetCoilZ");
  double FlashMagnetZ     = description.constant<double>("FlashMagnetZ");
  double OPPPDetZtoMagnet = description.constant<double>("OPPPDetZtoMagnet");
  double VacChambertoOPPPDetZGap = description.constant<double>("VacChambertoOPPPDetZGap");
  double IPMagnetZpos     = description.constant<double>("IPMagnetZpos");
  double FlashMFieldLength= description.constant<double>("FlashMFieldLength");
  bool   OverlapTest      = (description.constant<int>("OverlapTest") != 0);

  double gaptomag = 0.5 * (FlashMagnetCoilZ - FlashMagnetZ);
  double bpipel   = OPPPDetZtoMagnet - VacChambertoOPPPDetZGap - gaptomag;
  double zpos     = IPMagnetZpos + FlashMFieldLength/2.0 + bpipel/2.0 + gaptomag;

  Assembly chamberAssembly = ConstructVacuumChamber(description);

  PlacedVolume pvAsm = envelope.placeVolume(chamberAssembly, Position(0.0, 0.0, zpos));
  pvAsm.addPhysVolID("DetChamber", 0);
  if (OverlapTest) pvAsm.ptr()->CheckOverlaps();
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxDetChamberAssembly::ConstructVacuumChamber(dd4hep::Detector& description)
{
  Material vcWindowMaterial       = description.material(description.constant<std::string>("VacChamberWindowMaterial"));
  Material vacuumMaterial         = description.material(description.constant<std::string>("BeamPipeVacuumMaterial"));
  Material beamPipeMaterial       = description.material(description.constant<std::string>("BeamPipeMaterial"));

  // For <constant name="..." value="..." type="string"/> automatic substitution does not work, so this is work around:
  auto mat = [&](const std::string& name) {return description.constant<std::string>(name);};
  Material vacChamberMaterial     = description.material(description.constant<std::string>(mat("VacChamberMaterial")));
  Material vacChamberSideMaterial = description.material(description.constant<std::string>(mat("VacChamberSideMaterial")));

  double FlashMagnetCoilZ              = description.constant<double>("FlashMagnetCoilZ");
  double FlashMagnetZ                  = description.constant<double>("FlashMagnetZ");
  double FlashMagneteffY               = description.constant<double>("FlashMagneteffY");
  double OPPPDetZtoMagnet              = description.constant<double>("OPPPDetZtoMagnet");
  double VacChambertoOPPPDetZGap       = description.constant<double>("VacChambertoOPPPDetZGap");
  double VacChamberMagXZWidth          = description.constant<double>("VacChamberMagXZWidth");
  double VacChamberDetXZWidth          = description.constant<double>("VacChamberDetXZWidth");
  double VacChamberSideWallThickness   = description.constant<double>("VacChamberSideWallThickness");
  double VacChamberXZWallThickness     = description.constant<double>("VacChamberXZWallThickness");
  double VacChamberHight               = description.constant<double>("VacChamberHight");
  double VacChamberFrontThickness      = description.constant<double>("VacChamberFrontThickness");
  double VacChamberWindowPanelZ        = description.constant<double>("VacChamberWindowPanelZ");
  double VacChamberWindowPanelX        = description.constant<double>("VacChamberWindowPanelX");
  double VacChamberWindowPanelY        = description.constant<double>("VacChamberWindowPanelY");
  double VacChamberWindowPanelCutX     = description.constant<double>("VacChamberWindowPanelCutX");
  double VacChamberWindowPanelCutY     = description.constant<double>("VacChamberWindowPanelCutY");
  double VacChamberWindowPanelCutXpos  = description.constant<double>("VacChamberWindowPanelCutXpos");
  double VacChamberWindowThickness     = description.constant<double>("VacChamberWindowThickness");
  double VacChamberWindowPanelXpos     = description.constant<double>("VacChamberWindowPanelXpos");
  double VacChamberFrontCutX           = description.constant<double>("VacChamberFrontCutX");
  double VacChamberFrontCutY           = description.constant<double>("VacChamberFrontCutY");
  double VacChamberFrontCutXpos        = description.constant<double>("VacChamberFrontCutXpos");
  double VacChamberReinforceBarX       = description.constant<double>("VacChamberReinforceBarX");
  double VacChamberReinforceBarZ       = description.constant<double>("VacChamberReinforceBarZ");
  double VacChamberReinforceFilletXZ   = description.constant<double>("VacChamberReinforceFilletXZ");
  double BPipeR                        = description.constant<double>("BPipeR");
  double BPipeThickness                = description.constant<double>("BPipeThickness");
  double VacChamberFrameLX             = description.constant<double>("VacChamberFrameLX");
  double VacChamberFrameLY             = description.constant<double>("VacChamberFrameLY");
  double VacChamberFrameHX             = description.constant<double>("VacChamberFrameHX");
  double VacChamberFrameHY             = description.constant<double>("VacChamberFrameHY");
  double VacChamberFrameHZ             = description.constant<double>("VacChamberFrameHZ");
  bool   OverlapTest                   = (description.constant<int>("OverlapTest") != 0);

  double gaptomag  = 0.5 * (FlashMagnetCoilZ - FlashMagnetZ);
  double vclength  = OPPPDetZtoMagnet - VacChambertoOPPPDetZGap - gaptomag;
  double magnetx   = 0.25 * VacChamberMagXZWidth;
  double detx      = 0.5  * VacChamberDetXZWidth + VacChamberSideWallThickness;
  double dwall     = 0.5  * VacChamberXZWallThickness;
  double vctheta   = std::atan2(detx - magnetx, vclength + gaptomag);
  double vcminx    = magnetx + gaptomag * std::tan(vctheta);
  double thetaside = std::atan2(2.0*(detx - magnetx), vclength + gaptomag);
  double fronttotz = VacChamberFrontThickness + VacChamberWindowPanelZ;
  double vcwalll   = vclength - fronttotz;
  double vcwallx   = detx - 0.5*fronttotz * std::tan(thetaside);
  double lbpipe    = VacChambertoOPPPDetZGap;

  // -----------------------------------------------------------------------
  // Container: G4Trap
  // G4Trap(dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2)
  // "half-trapezia" form: here dy1==dy2, dx1==dx2 (no tilt in y), alpha==0
  // Equivalent TGeoTrap(dz, theta_deg, phi_deg, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2)
  //   dz    = vclength/2
  //   theta = vctheta (polar angle of axis from Z — in degrees for ROOT)
  //   phi   = 0
  //   h1,h2 = VacChamberHight/2 (same top and bottom)
  //   bl1   = vcminx (half-width at -z)
  //   tl1   = vcminx
  //   bl2   = detx   (half-width at +z)
  //   tl2   = detx
  //   alpha = 0
  TGeoTrap* tgeoVCContainer = new TGeoTrap("solidVCContainer",
      vclength/2.0     / dd4hep::cm,
      vctheta * 180.0 / M_PI,   // theta in degrees
      0.0,                       // phi
      VacChamberHight/2.0 / dd4hep::cm,
      vcminx / dd4hep::cm,
      vcminx / dd4hep::cm,
      0.0,                       // alpha1
      VacChamberHight/2.0 / dd4hep::cm,
      detx   / dd4hep::cm,
      detx   / dd4hep::cm,
      0.0);                      // alpha2
  Solid  solidVCContainer(tgeoVCContainer);
  Volume logicVCContainer("logicVCContainer", solidVCContainer, vacuumMaterial);

  // -----------------------------------------------------------------------
  // Top/bottom (XZ) walls: also G4Trap, same shape as container but thinner
  TGeoTrap* tgeoVCXZWall = new TGeoTrap("solidVCXZWall",
      vcwalll/2.0      / dd4hep::cm,
      vctheta * 180.0 / M_PI,
      0.0,
      dwall / dd4hep::cm,
      vcminx  / dd4hep::cm,
      vcminx  / dd4hep::cm,
      0.0,
      dwall / dd4hep::cm,
      vcwallx / dd4hep::cm,
      vcwallx / dd4hep::cm,
      0.0);
  Solid  solidVCXZWall(tgeoVCXZWall);
  Volume logicVCXZWall("logicVCXZWall", solidVCXZWall, vacChamberMaterial);

  Position wtrans((vcwallx - detx)/2.0, -(VacChamberHight/2.0 - dwall), -fronttotz/2.0);
  PlacedVolume pvWall0 = logicVCContainer.placeVolume(logicVCXZWall, wtrans);
  pvWall0.addPhysVolID("VCXZWall", 0);
  if (OverlapTest) pvWall0.ptr()->CheckOverlaps();
  PlacedVolume pvWall1 = logicVCContainer.placeVolume(logicVCXZWall,
      Position(wtrans.X(), VacChamberHight/2.0 - dwall, wtrans.Z()));
  pvWall1.addPhysVolID("VCXZWall", 1);
  if (OverlapTest) pvWall1.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Side wall: G4Para
  // G4Para(dx, dy, dz, alpha, theta, phi) — theta is tilt of Z axis in XZ plane
  // TGeoPara(dx, dy, dz, alpha_deg, theta_deg, phi_deg)
  double fpnldx = VacChamberSideWallThickness / std::cos(thetaside);
  double fpnldy = VacChamberHight - 2.0 * VacChamberXZWallThickness;
  TGeoPara* tgeoSideWall = new TGeoPara("solidSideWall",
      fpnldx/2.0  / dd4hep::cm,
      fpnldy/2.0  / dd4hep::cm,
      vcwalll/2.0 / dd4hep::cm,
      0.0,
      thetaside * 180.0 / M_PI,
      0.0);
  Solid  solidSideWall(tgeoSideWall);
  Volume logicSideWall("logicSideWall", solidSideWall, vacChamberSideMaterial);

  Position swtrans((vcwallx - detx - fpnldx)/2.0 + 0.5*(vcwallx + vcminx), 0.0, -fronttotz/2.0);
  PlacedVolume pvSide = logicVCContainer.placeVolume(logicSideWall, swtrans);
  pvSide.addPhysVolID("VCSideWall", 0);
  if (OverlapTest) pvSide.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Front panel (with beam pipe hole + rectangular cut + round fillets)
  double frontpnlz = VacChamberFrontThickness;
  Box  solidVCFront1(vcwallx, VacChamberHight/2.0, frontpnlz/2.0);
  Tube solidCVBeamPipeHole(0.0, BPipeR - BPipeThickness, frontpnlz, 0.0, 2.0*M_PI);

  SubtractionSolid solidVCFront2("solidVCFront2", solidVCFront1, solidCVBeamPipeHole,
      Position(-vcwallx, 0.0, 0.0));

  Tube solidCVFrontCutR(0.0, VacChamberFrontCutY/2.0, frontpnlz, 0.0, 2.0*M_PI);
  Box  solidVCFrontCut1(VacChamberFrontCutX/2.0, VacChamberFrontCutY/2.0, frontpnlz);
  double reccutposx = VacChamberFrontCutXpos + VacChamberFrontCutX/2.0 - vcwallx;

  SubtractionSolid solidVCFront3("solidVCFront3", solidVCFront2, solidVCFrontCut1,
      Position(reccutposx, 0.0, 0.0));
  SubtractionSolid solidVCFront4("solidVCFront4", solidVCFront3, solidCVFrontCutR,
      Position(reccutposx - VacChamberFrontCutX/2.0, 0.0, 0.0));
  SubtractionSolid solidVCFront("solidVCFront", solidVCFront4, solidCVFrontCutR,
      Position(reccutposx + VacChamberFrontCutX/2.0, 0.0, 0.0));

  Volume logicVCFront("logicVCFront", solidVCFront, vacChamberMaterial);
  Position wndtrans(vcwallx - (vcminx + detx)/2.0, 0.0, (vcwalll - fronttotz + frontpnlz)/2.0);
  PlacedVolume pvFront = logicVCContainer.placeVolume(logicVCFront, wndtrans);
  pvFront.addPhysVolID("VCFront", 0);
  if (OverlapTest) pvFront.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Reinforcement bar (with beam pipe hole)
  double reinfy = VacChamberHight - 2.0 * VacChamberXZWallThickness;
  Box solidVCReinforceBar1(VacChamberReinforceBarX/2.0, reinfy/2.0, VacChamberReinforceBarZ/2.0);
  SubtractionSolid solidVCReinforceBar("solidVCReinforceBar", solidVCReinforceBar1, solidCVBeamPipeHole,
      Position(-VacChamberReinforceBarX/2.0, 0.0, 0.0));
  Volume logicVCReinforceBar("logicVCReinforceBar", solidVCReinforceBar, vacChamberMaterial);

  Position trreinbar(VacChamberReinforceBarX/2.0 - (vcminx + detx)/2.0, 0.0,
                     (vcwalll - fronttotz - VacChamberReinforceBarZ)/2.0);
  PlacedVolume pvReinBar = logicVCContainer.placeVolume(logicVCReinforceBar, trreinbar);
  pvReinBar.addPhysVolID("VCReinforceBar", 0);
  if (OverlapTest) pvReinBar.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Reinforcement fillet
  Box  solidVCReinforceFillet1(VacChamberReinforceFilletXZ/2.0, reinfy/2.0, VacChamberReinforceFilletXZ/2.0);
  Tube solidCVReinforceFilletCut(0.0, VacChamberReinforceFilletXZ, reinfy, 0.0, 2.0*M_PI);

  // G4RotationMatrix(G4ThreeVector(1,0,0), pi/2): axis-angle around X by pi/2 → RotationZYX(0,0,pi/2)
  SubtractionSolid solidVCReinforceFillet2("solidVCReinforceFillet2",
      solidVCReinforceFillet1, solidCVReinforceFilletCut,
      Transform3D(RotationZYX(0.0, 0.0, M_PI/2.0),
                  Position(VacChamberReinforceFilletXZ/2.0, 0.0, -VacChamberReinforceFilletXZ/2.0)));

  Box solidVCReiFiletCut(VacChamberReinforceFilletXZ/2.0, VacChamberFrontCutY/2.0, VacChamberReinforceFilletXZ);
  double fltcutposx = VacChamberFrontCutXpos - VacChamberReinforceBarX;
  SubtractionSolid solidVCReinforceFillet3("solidVCReinforceFillet3",
      solidVCReinforceFillet2, solidVCReiFiletCut,
      Position(fltcutposx, 0.0, 0.0));
  SubtractionSolid solidVCReinforceFillet("solidVCReinforceFillet",
      solidVCReinforceFillet3, solidCVFrontCutR,
      Position(fltcutposx - VacChamberReinforceFilletXZ/2.0, 0.0, 0.0));

  Volume logicVCReinforceFillet("logicVCReinforceFillet", solidVCReinforceFillet, vacChamberMaterial);
  Position trfilet(VacChamberReinforceBarX + VacChamberReinforceFilletXZ/2.0 - (vcminx + detx)/2.0,
                   0.0,
                   (vcwalll - fronttotz - VacChamberReinforceFilletXZ)/2.0);
  PlacedVolume pvFillet = logicVCContainer.placeVolume(logicVCReinforceFillet, trfilet);
  pvFillet.addPhysVolID("VCReinforceFillet", 0);
  if (OverlapTest) pvFillet.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Window panel (with detector cut + beam pipe hole)
  Box  solidVCWindowPanel1(VacChamberWindowPanelX/2.0, VacChamberWindowPanelY/2.0, VacChamberWindowPanelZ/2.0);
  Box  solidVCWindowCut   (VacChamberWindowPanelCutX/2.0, VacChamberWindowPanelCutY/2.0, VacChamberWindowPanelZ/2.0);
  double wincutposx = VacChamberWindowPanelCutXpos
                      + (VacChamberWindowPanelCutX - VacChamberWindowPanelX)/2.0;
  SubtractionSolid solidVCWindowPanel2("solidVCWindowPanel2", solidVCWindowPanel1, solidVCWindowCut,
      Position(wincutposx, 0.0, -VacChamberWindowThickness));

  Tube solidCVWBeamPipeCut(0.0, BPipeR, VacChamberWindowPanelZ, 0.0, 2.0*M_PI);
  double winpnlpipecutx = -VacChamberWindowPanelX/2.0 - VacChamberWindowPanelXpos;
  SubtractionSolid solidVCWindowPanel("solidVCWindowPanel", solidVCWindowPanel2, solidCVWBeamPipeCut,
      Position(winpnlpipecutx, 0.0, 0.0));
  Volume logicVCWindowPanel("logicVCWindowPanel", solidVCWindowPanel, vcWindowMaterial);

  double winpnlposx = VacChamberWindowPanelXpos + VacChamberWindowPanelX/2.0 - (vcminx + detx)/2.0;
  Position wndpnltrans(winpnlposx, 0.0,
                       frontpnlz + (vcwalll - fronttotz + VacChamberWindowPanelZ)/2.0);
  // Note: wndpnltrans.z() == frontpnlz + (vcwalll-fronttotz+VacChamberWindowPanelZ)/2.0
  PlacedVolume pvWndPanel = logicVCContainer.placeVolume(logicVCWindowPanel, wndpnltrans);
  pvWndPanel.addPhysVolID("VCWindowPanel", 0);
  if (OverlapTest) pvWndPanel.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Short beam pipe between window panels
  // G4Tubs half-circle: startPhi=0, dPhi=pi — maps to DD4hep Tube(rmin,rmax,hz,startPhi,endPhi)
  Tube   solidVCWindBeamPipe(BPipeR - BPipeThickness, BPipeR,
                             VacChamberWindowPanelZ/2.0, 0.0, M_PI);
  Volume logicVCWindowBeamPipe("logicVCWindowBeamPipe", solidVCWindBeamPipe, beamPipeMaterial);

  // G4RotationMatrix(G4ThreeVector(0,0,1), pi/2): axis-angle around Z by pi/2 → RotationZYX(pi/2,0,0)
  Position wndpipetrans(-(vcminx + detx)/2.0, 0.0, wndpnltrans.Z());
  PlacedVolume pvWndPipe = logicVCContainer.placeVolume(logicVCWindowBeamPipe,
      Transform3D(RotationZYX(-M_PI/2.0, 0.0, 0.0), wndpipetrans));
  pvWndPipe.addPhysVolID("VCWindowBeamPipe", 0);
  if (OverlapTest) pvWndPipe.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Short pipe joining vacuum chamber to photon detector pipe
  Tube   solidBeamPipeVCG(BPipeR - BPipeThickness, BPipeR, lbpipe/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeVCG("logicBeamPipeVCG", solidBeamPipeVCG, beamPipeMaterial);

  Tube   solidBeamPipeVCGVac(0.0, BPipeR - BPipeThickness, lbpipe/2.0, 0.0, 2.0*M_PI);
  Volume logicBeamPipeVCGVac("logicBeamPipeVCGVac", solidBeamPipeVCGVac, vacuumMaterial);

  // -----------------------------------------------------------------------
  // Trapezoidal joining section from field volume to vacuum chamber (-Z end)
  double jclength  = 0.8 * gaptomag;
  double maptheta  = std::atan2(2.0*(vcminx - magnetx), jclength);
  double hx1       = 2.0 * magnetx;
  double hx2       = 2.0 * magnetx + jclength * std::tan(maptheta);
  double jvch1     = FlashMagneteffY;
  double jvch2     = VacChamberHight;

  Trapezoid solidVCMagFieldJoin(hx1, hx2, jvch1/2.0, jvch2/2.0, jclength/2.0);
  Volume    logicVCMagFieldJoin("logicVCMagFieldJoin", solidVCMagFieldJoin, beamPipeMaterial);

  hx1  -= BPipeThickness;
  hx2  -= VacChamberSideWallThickness / std::cos(thetaside);
  jvch1-= 2.0 * BPipeThickness;
  jvch2-= 2.0 * VacChamberXZWallThickness;
  Trapezoid solidVCMagFieldJoinVac(hx1, hx2, jvch1/2.0, jvch2/2.0, jclength/2.0);
  Volume    logicVCMagFieldJoinVac("logicVCMagFieldJoinVac", solidVCMagFieldJoinVac, vacuumMaterial);

  PlacedVolume pvJoinVac = logicVCMagFieldJoin.placeVolume(logicVCMagFieldJoinVac,
                                                            Position(0.0, 0.0, 0.0));
  pvJoinVac.addPhysVolID("VCMagFieldJoinVac", 0);
  if (OverlapTest) pvJoinVac.ptr()->CheckOverlaps();

  // Rectangular extension to fill remaining gap
  double hx3      = 2.0 * magnetx;
  double jclength1= gaptomag - jclength;
  Box    solidVCMagFieldJoin1(hx3, FlashMagneteffY/2.0, jclength1/2.0);
  Volume logicVCMagFieldJoin1("logicVCMagFieldJoin1", solidVCMagFieldJoin1, beamPipeMaterial);

  hx3 -= BPipeThickness;
  Box    solidVCMagFieldJoin1Vac(hx3, FlashMagneteffY/2.0 - BPipeThickness, jclength1/2.0);
  Volume logicVCMagFieldJoin1Vac("logicVCMagFieldJoin1Vac", solidVCMagFieldJoin1Vac, vacuumMaterial);

  PlacedVolume pvJoin1Vac = logicVCMagFieldJoin1.placeVolume(logicVCMagFieldJoin1Vac,
                                                              Position(0.0, 0.0, 0.0));
  pvJoin1Vac.addPhysVolID("VCMagFieldJoin1Vac", 0);
  if (OverlapTest) pvJoin1Vac.ptr()->CheckOverlaps();

  // -----------------------------------------------------------------------
  // Build the full chamber assembly
  // -----------------------------------------------------------------------
  Assembly fChamberAssembly(fChamberType + "_assembly");

  // Container placed at +X (right side)
  double pipetransx = 0.5*vclength*std::tan(vctheta) + vcminx;
  fChamberAssembly.placeVolume(logicVCContainer, Position( pipetransx, 0.0, 0.0));

  // Container placed at -X (left side), rotated pi around Z
  // G4RotationMatrix(G4ThreeVector(0,0,1), pi): axis-angle around Z → RotationZYX(pi,0,0)
  fChamberAssembly.placeVolume(logicVCContainer,
      Transform3D(RotationZYX(M_PI, 0.0, 0.0), Position(-pipetransx, 0.0, 0.0)));

  // Short pipe to photon detectors (+Z end)
  Position vcbpg(0.0, 0.0, (vclength + lbpipe)/2.0);
  fChamberAssembly.placeVolume(logicBeamPipeVCG,    vcbpg);
  fChamberAssembly.placeVolume(logicBeamPipeVCGVac, vcbpg);

  // Trapezoidal join to field volume (-Z end)
  fChamberAssembly.placeVolume(logicVCMagFieldJoin,
      Position(0.0, 0.0, -(vclength + jclength)/2.0));

  // Rectangular extension further toward -Z
  fChamberAssembly.placeVolume(logicVCMagFieldJoin1,
      Position(0.0, 0.0, -jclength - (vclength + jclength1)/2.0));

  // Frame assembly — top (+Y) and bottom (-Y)
  Assembly frameAssembly = ConstructFrameAssembly(description, vclength);

  Position trframe(0.0, VacChamberHight/2.0, -VacChamberWindowPanelZ/2.0);
  LxAux::AddAssemblyVolumes(fChamberAssembly, frameAssembly, trframe);

  // G4RotationMatrix(G4ThreeVector(0,0,1), pi) → RotationZYX(pi,0,0)
  LxAux::AddAssemblyVolumes(fChamberAssembly, frameAssembly,
                            Position(0.0, -VacChamberHight/2.0, -VacChamberWindowPanelZ/2.0),
                            RotationZYX(M_PI, 0.0, 0.0));

  return fChamberAssembly;
}


// ---------------------------------------------------------------------------
dd4hep::Assembly LxDetChamberAssembly::ConstructFrameAssembly(dd4hep::Detector& description,
                                                               double vcLength)
{
  auto mat = [&](const std::string& name) {return description.constant<std::string>(name);};
  Material vacChamberMaterial     = description.material(description.constant<std::string>(mat("VacChamberMaterial")));

  double VacChamberWindowPanelZ = description.constant<double>("VacChamberWindowPanelZ");
  double VacChamberFrameLX      = description.constant<double>("VacChamberFrameLX");
  double VacChamberFrameLY      = description.constant<double>("VacChamberFrameLY");
  double VacChamberFrameHX      = description.constant<double>("VacChamberFrameHX");
  double VacChamberFrameHY      = description.constant<double>("VacChamberFrameHY");
  double VacChamberFrameHZ      = description.constant<double>("VacChamberFrameHZ");

  double frameLZ = vcLength - VacChamberWindowPanelZ;

  Box    solidVCFrameLBar(VacChamberFrameLX/2.0, VacChamberFrameLY/2.0, frameLZ/2.0);
  Volume logicVCFrameLBar("logicVCFrameLBar", solidVCFrameLBar, vacChamberMaterial);

  Box    solidVCFrameHBar(VacChamberFrameHX/2.0, VacChamberFrameHY/2.0, VacChamberFrameHZ/2.0);
  Volume logicVCFrameHBar("logicVCFrameHBar", solidVCFrameHBar, vacChamberMaterial);

  Assembly frameAssembly("frameAssembly");

  double trlbarx = (VacChamberFrameHX + VacChamberFrameLX)/2.0;
  frameAssembly.placeVolume(logicVCFrameLBar,
      Position( trlbarx, VacChamberFrameLY/2.0, 0.0));
  frameAssembly.placeVolume(logicVCFrameLBar,
      Position(-trlbarx, VacChamberFrameLY/2.0, 0.0));
  frameAssembly.placeVolume(logicVCFrameHBar,
      Position(0.0, VacChamberFrameHY/2.0, 0.0));

  return frameAssembly;
}


// ---------------------------------------------------------------------------
// DD4hep plugin entry point
// ---------------------------------------------------------------------------

static Ref_t create_LxDetChamberAssembly(dd4hep::Detector& description,
                                          xml_h             e,
                                          dd4hep::SensitiveDetector /* sd */)
{
  xml_comp_t  x_det(e);
  std::string detName = x_det.nameStr();
  int         detID   = x_det.id();

  dd4hep::DetElement sdet(detName, detID);

  LxDetChamberAssembly vc("DetectorVacuumChamber");
  vc.Construct(description, sdet, e);

  return sdet;
}

DECLARE_DETELEMENT(LxDetChamberAssembly, create_LxDetChamberAssembly)
