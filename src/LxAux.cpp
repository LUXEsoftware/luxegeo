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
