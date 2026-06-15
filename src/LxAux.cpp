//
/// \brief Implementation of the LxAux class (DD4hep version)
//

#include <cmath>
#include <vector>
#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#include "LxAux.h"

using namespace dd4hep;


dd4hep::Assembly LxAux::BuildTable(dd4hep::Detector&  description,
                                   const std::string& tname,
                                   double xs, double ys, double zs,
                                   int nleg,
                                   double ytop,
                                   double rleg)
{
  Material tableTopPlateMaterial = description.material(
      description.constant<std::string>("OPPPBasePlateMaterial"));

  double syplate = ytop;
  if (syplate < 0.0)
    syplate = description.constant<double>("OPPPBasePlateY");

  Assembly tableAssembly(tname);

  std::string lgc = "logic";

  // -----------------------------------------------------------------------
  // Top plate
  // -----------------------------------------------------------------------
  if (syplate > 1.0e-15) {
    std::string topname = tname + "Top";
    Box    solidTopPlate(xs/2.0, syplate/2.0, zs/2.0);
    Volume logicTopPlate(lgc+topname, solidTopPlate, tableTopPlateMaterial);
    tableAssembly.placeVolume(logicTopPlate, Position(0.0, -syplate/2.0, 0.0));
  }

  // -----------------------------------------------------------------------
  // Legs
  // -----------------------------------------------------------------------
  double rout = rleg;
  if (rout < 0.0)
    rout = description.constant<double>("OPPPBasePlateR");

  double blegz = ys - syplate;
  Tube   solidPlateLeg(0.0, rout, blegz/2.0, 0.0, 2.0*M_PI);
  Volume logicPlateLeg(lgc + tname + "Leg", solidPlateLeg, tableTopPlateMaterial);

  int corrnleg = nleg;
  if (nleg > 4)                      corrnleg = 4;
  if (xs < 6.0*rout && corrnleg > 2) corrnleg = 2;
  if (zs < 6.0*rout)                 corrnleg = 1;

  double dx = xs/2.0 - 1.5*rout;
  double dy = syplate + blegz/2.0;
  double dz = zs/2.0 - 1.5*rout;

  std::vector<Position> lmove;
  switch (corrnleg) {
    case 1: lmove.push_back(Position(0.0,  -dy,  0.0)); break;
    case 2: lmove.push_back(Position(0.0,  -dy, -dz));
            lmove.push_back(Position(0.0,  -dy,  dz)); break;
    case 3: lmove.push_back(Position(-dx,  -dy, -dz));
            lmove.push_back(Position( dx,  -dy, -dz));
            lmove.push_back(Position(0.0,  -dy,  dz)); break;
    case 4: lmove.push_back(Position(-dx,  -dy, -dz));
            lmove.push_back(Position(-dx,  -dy,  dz));
            lmove.push_back(Position( dx,  -dy, -dz));
            lmove.push_back(Position( dx,  -dy,  dz)); break;
  }

  // G4RotationMatrix(G4ThreeVector(-1,0,0), pi/2): rotation around -X by pi/2,
  // maps Z -> -Y, standing the tube (aligned along Z) vertically along Y.
  // DD4hep RotationZYX(rz,ry,rx): rx = -pi/2 achieves the same.
  RotationZYX legRot(0.0, 0.0, -M_PI/2.0);

  for (const auto& trv : lmove)
    tableAssembly.placeVolume(logicPlateLeg, Transform3D(legRot, trv));

  return tableAssembly;
}


void LxAux::AddAssemblyVolumes(dd4hep::Assembly&          dstAssembly,
                               dd4hep::Assembly&          srcAssembly,
                               const dd4hep::Position&    translation,
                               const dd4hep::RotationZYX& rotation)
{
  // Re-place all daughters of srcAssembly into dstAssembly with an additional
  // outer transformation (rotation, translation) applied.
  //
  // Composition rule:
  //   finalRot = outerRot * innerRot
  //   finalPos = outerRot * innerPos + outerTranslation
  //
  // The inner position is in the srcAssembly local frame and must be rotated
  // into the destination frame before adding the outer translation.
  // This is what standard Transform3D multiplication does, so outerTrf * innerTrf
  // is actually correct — the earlier bug was the unit conversion: ROOT TGeoMatrix
  // stores translations in cm, DD4hep works in mm, so the factor must be dd4hep::cm
  // (not dd4hep::mm which would be a no-op).

  ROOT::Math::Rotation3D outerRot(rotation);

  TGeoVolume* srcVol = srcAssembly.ptr();
  if (!srcVol) {
    printout(ERROR, "LxAux::AddAssemblyVolumes", "srcAssembly has null TGeoVolume pointer");
    return;
  }

  int ndaughters = srcVol->GetNdaughters();
  for (int i = 0; i < ndaughters; ++i) {
    TGeoNode*   node = srcVol->GetNode(i);
    TGeoVolume* vol  = node->GetVolume();
    TGeoMatrix* mat  = node->GetMatrix();

    const Double_t* trans = mat->GetTranslation();
    const Double_t* rot   = mat->GetRotationMatrix();

    ROOT::Math::Rotation3D innerRot(
        rot[0], rot[1], rot[2],
        rot[3], rot[4], rot[5],
        rot[6], rot[7], rot[8]);

    // ROOT stores translations in cm; convert to DD4hep internal units (mm)
    Position innerPos(trans[0]*dd4hep::cm,
                      trans[1]*dd4hep::cm,
                      trans[2]*dd4hep::cm);

    // Compose outer and inner transforms
    ROOT::Math::Rotation3D finalRot = outerRot * innerRot;
    Position               finalPos = outerRot * innerPos + translation;

    Volume dd4vol(vol);
    dstAssembly.placeVolume(dd4vol, Transform3D(finalRot, finalPos));
  }
}


dd4hep::Volume LxAux::BuildPedestal(dd4hep::Detector&  description,
                                    const std::string& pname,
                                    double xs, double ys, double zs)
{
  Material pedestalMaterial = description.material("ShieldingConcrete");
  std::string volname = pname + "Pedestal";
  Box    solidPedestal(xs/2.0, ys/2.0, zs/2.0);
  Volume logicPedestal("logic" + volname, solidPedestal, pedestalMaterial);
  return logicPedestal;
}


dd4hep::Assembly LxAux::BuildHexapod(dd4hep::Detector&  description,
                                     const std::string& pname,
                                     double&            hexhight)
{
  // Hexapod support assembly, centred in Y.
  // hexhight is first used to shorten the middle section from OPPPHexapodY,
  // then updated to the actual height of the assembly returned.
  Material opppHexapodMaterial = description.material(
      description.constant<std::string>("OPPPHexapodMaterial"));

  double OPPPHexapodUpR        = description.constant<double>("OPPPHexapodUpR");
  double OPPPHexapodUpH        = description.constant<double>("OPPPHexapodUpH");
  double OPPPHexapodDownR      = description.constant<double>("OPPPHexapodDownR");
  double OPPPHexapodDownH      = description.constant<double>("OPPPHexapodDownH");
  double OPPPHexapodY          = description.constant<double>("OPPPHexapodY");

  std::string lgc = "logic";
  std::string sld = "solid";

  std::string topvolname = pname + "HexTop";
  std::string botvolname = pname + "HexBottom";
  std::string midvolname = pname + "HexMiddle";

  Tube   solidOPPPHexTop(0.0, OPPPHexapodUpR,   OPPPHexapodUpH/2.0,   0.0, 2.0*M_PI);
  Volume logicOPPPHexTop(lgc + topvolname, solidOPPPHexTop, opppHexapodMaterial);

  Tube   solidOPPPHexBottom(0.0, OPPPHexapodDownR, OPPPHexapodDownH/2.0, 0.0, 2.0*M_PI);
  Volume logicOPPPHexBottom(lgc + botvolname, solidOPPPHexBottom, opppHexapodMaterial);

  double hexmh    = OPPPHexapodY - OPPPHexapodUpH - OPPPHexapodDownH - hexhight;
  double hexmrup  = 0.9 * OPPPHexapodUpR;
  double hexmrdwn = 0.9 * OPPPHexapodDownR;
  hexhight = OPPPHexapodUpH + OPPPHexapodDownH + hexmh;

  // Commented: cone version
  // Cone solidOPPPHexMiddle(hexmh/2.0,
  //     hexmrdwn - OPPPHexapodThickness, hexmrdwn,
  //     hexmrup  - OPPPHexapodThickness, hexmrup);

  // Polycone leg profile
  double hsphi = M_PI / 6.0;
  Position lupposl(hexmrup  * std::cos(-hsphi/2.0),  hexhight/2.0 - OPPPHexapodUpH,   hexmrup  * std::sin(-hsphi/2.0));
  Position ldwposl(hexmrdwn * std::cos( hsphi/2.0), -hexhight/2.0 + OPPPHexapodDownH, hexmrdwn * std::sin( hsphi/2.0));
  ROOT::Math::XYZVector lposdifl(lupposl.X() - ldwposl.X(),
                                  lupposl.Y() - ldwposl.Y(),
                                  lupposl.Z() - ldwposl.Z());
  double hexsidel = lposdifl.R();

  std::vector<double> rohexleg{0.0, 7.5,  7.5,  25.0, 25.0, 7.5,  7.5,  0.0};
  std::vector<double> rihexleg{0.0, 5.0,  5.0,  22.0, 22.0, 5.0,  5.0,  0.0};
  std::vector<double>  lhexleg{0.0, 0.04, 0.14, 0.16, 0.84, 0.86, 0.96, 1.0};
  std::for_each(rohexleg.begin(), rohexleg.end(), [=](double& x){ x *= dd4hep::mm; });
  std::for_each(rihexleg.begin(), rihexleg.end(), [=](double& x){ x *= dd4hep::mm; });
  // Scale lhexleg by hexsidel and convert units to mm (values above are in mm already)
  std::for_each(lhexleg.begin(), lhexleg.end(), [=](double& x){ x *= hexsidel; });

  Polycone solidOPPPHexMiddle(0.0, 2.0*M_PI, rihexleg, rohexleg, lhexleg);
  Volume   logicOPPPHexMiddle(lgc + midvolname, solidOPPPHexMiddle, opppHexapodMaterial);

  Assembly oppHexapodAssembly(pname + "HexapodAssembly");

  // Top disc — G4RotationMatrix(G4ThreeVector(-1,0,0), pi/2): axis-angle around -X by pi/2
  // → RotationZYX(0, 0, -pi/2)
  RotationZYX discRot(0.0, 0.0, -M_PI/2.0);
  oppHexapodAssembly.placeVolume(logicOPPPHexTop,
      Transform3D(discRot, Position(0.0,  (hexhight - OPPPHexapodUpH)/2.0,   0.0)));
  oppHexapodAssembly.placeVolume(logicOPPPHexBottom,
      Transform3D(discRot, Position(0.0, -(hexhight - OPPPHexapodDownH)/2.0, 0.0)));

  // Six legs: two groups of three, each rotated 120 degrees apart
  // Geant4: rotateY(theta) then rotateZ(phi) on default matrix
  // = intrinsic Y then Z = fixed-frame Z first then Y = RotationZYX(phi, theta, 0)
  // This did not work and was replaced explicit rotations as in original G4 code.
  for (int nl = 0; nl < 3; ++nl) {
    double phinl = nl * 2.0*M_PI/3.0;
    Position luppos(hexmrup  * std::cos(phinl - hsphi/2.0),  hexhight/2.0 - OPPPHexapodUpH,   hexmrup  * std::sin(phinl - hsphi/2.0));
    Position ldwpos(hexmrdwn * std::cos(phinl + hsphi/2.0), -hexhight/2.0 + OPPPHexapodDownH, hexmrdwn * std::sin(phinl + hsphi/2.0));
    ROOT::Math::XYZVector lposdif(luppos.X()-ldwpos.X(), luppos.Y()-ldwpos.Y(), luppos.Z()-ldwpos.Z());
    double theta = lposdif.Theta();
    double phi   = lposdif.Phi();
    RotationZYX rleg = RotationZYX() * RotationZ(phi) * RotationY(theta);
    oppHexapodAssembly.placeVolume(logicOPPPHexMiddle, Transform3D(rleg, ldwpos));
  }

  for (int nl = 0; nl < 3; ++nl) {
    double phinl = nl * 2.0*M_PI/3.0 - 2.0*hsphi;
    Position luppos(hexmrup  * std::cos(phinl + hsphi/2.0),  hexhight/2.0 - OPPPHexapodUpH,   hexmrup  * std::sin(phinl + hsphi/2.0));
    Position ldwpos(hexmrdwn * std::cos(phinl - hsphi/2.0), -hexhight/2.0 + OPPPHexapodDownH, hexmrdwn * std::sin(phinl - hsphi/2.0));
    ROOT::Math::XYZVector lposdif(luppos.X()-ldwpos.X(), luppos.Y()-ldwpos.Y(), luppos.Z()-ldwpos.Z());
    double theta = lposdif.Theta();
    double phi   = lposdif.Phi();
    RotationZYX rleg = RotationZYX() * RotationZ(phi) * RotationY(theta);
    oppHexapodAssembly.placeVolume(logicOPPPHexMiddle, Transform3D(rleg, ldwpos));
  }

  return oppHexapodAssembly;
}
