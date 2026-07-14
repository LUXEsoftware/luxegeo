//
/// \brief Implementation of the Magnetic field models classes
//

#include <algorithm>
#include <functional>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <typeinfo>

#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DD4hep/Factories.h"
#include "DD4hep/DD4hepUnits.h"
#include "XML/XML.h"

#include "LxBFields.h"


FieldData2D::FieldData2D(const std::string fname) : FieldDistribution(), fFName(fname), fNsort(10)
{
  LoadFieldData(fFName);
  frMin = fFieldData[0];
  frMax = fFieldData[0];
  for (const auto &rr : fFieldData) {
    if (frMin.X() > rr.X()) frMin.SetX(rr.X());
    if (frMin.Y() > rr.Y()) frMin.SetY(rr.Y());
    if (frMax.X() < rr.X()) frMax.SetX(rr.X());
    if (frMax.Y() < rr.Y()) frMax.SetY(rr.Y());
  }
}



void FieldData2D::LoadFieldData(const std::string fname)
{
  std::fstream  fdataf;

  fdataf.open(fname, std::ios::in);
  if (!fdataf.is_open()) {
    throw std::runtime_error("FieldData2D::LoadFieldData(): failed to open magnetic field data file " + fname);
  }

  double x, z, bf;
  while (!fdataf.eof()) {
      fdataf >> x >> z >> bf;
      if (fdataf.eof()) break;
      fFieldData.push_back(dd4hep::Position(x*dd4hep::mm, z*dd4hep::mm, bf*dd4hep::tesla));
//       std::cout << x << " " << z << "  " << bf << std::endl;
  }
  fdataf.close();
  std::cout << fFieldData.size() << " Field data points were loaded from the file " << fname << std::endl;
}



double FieldData2D::DataField(const double xx, const double yy) const
{
  if (!(xx > frMin.X() && xx < frMax.X() && yy > frMin.Y() && yy < frMax.Y()) ) return 0.0;
  
  std::vector<double> drv;
  std::for_each(fFieldData.begin(), fFieldData.end(), 
                [&](const dd4hep::Position &vv){drv.push_back( pow(vv.X()-xx, 2.0) + pow(vv.Y()-yy, 2.0 ));} );
  std::vector<size_t> vind(drv.size());
  size_t ind = 0;
  std::generate(vind.begin(), vind.end(), [&](){return ind++;});
  point_cmp<double> pcmp(&drv);
  std::partial_sort (vind.begin(), vind.begin()+fNsort, vind.end(), pcmp);
   
  size_t ind3 = 2;
  while (CheckCollinear(fFieldData[vind[0]], fFieldData[vind[1]], fFieldData[vind[ind3]]) 
//          && (!CheckInside(fFieldData[vind[0]], fFieldData[vind[1]], fFieldData[vind[ind3]], xx, yy)) 
        ) {
    ind3++;
    if (ind3 > fNsort) {
      throw std::runtime_error("FieldData2D::DataField(): number of elements for partial sorting is not sufficient");
    }
  }

//   std::cout << "DataField: " << xx << "  " << yy << "  " << drv[vind[0]] 
//             << "  " << drv[vind[1]] << "  " << drv[vind[ind3]] << std::endl;
  return FindZ(fFieldData.at(vind[0]), fFieldData.at(vind[1]), fFieldData.at(vind[ind3]), xx, yy);
}



double FieldData2D::FindZ(const dd4hep::Position &r1, const dd4hep::Position &r2, const dd4hep::Position &r3, 
                               const double xx, const double yy) const
{
// Makes linear interpolation finding Z coordinate for a given (x,y) 
// as a point on the surface defined by three 3D points r1, r2, r3.
  dd4hep::Position dr1 = r2 -r1;
  dd4hep::Position dr2 = r3 -r1;

  dd4hep::Position normv = dr1.Cross(dr2).Unit();

  return r1.Z() - (normv.X()*(xx-r1.X()) + normv.Y()*(yy-r1.Y())) / normv.Z();
}



bool FieldData2D::CheckCollinear(const dd4hep::Position &r1, const dd4hep::Position &r2, const dd4hep::Position &r3) const
{
// Check if 2D projection to xy-plane of three 3D points are on one line.
  const double eps = 1.0e-15;
  dd4hep::Position dr1 = r2 -r1;
  dd4hep::Position dr2 = r3 -r1;
  dr1.SetZ(0.0);  
  dr2.SetZ(0.0);
  return ((dr1.Cross(dr2).Mag2()) < eps);
}


/////////////////////////////////////////////////////////////////////
LxBField::LxBField(const VFieldComponent *bx, const VFieldComponent *by, 
                   const VFieldComponent *bz, const dd4hep::Position pos):
  dd4hep::CartesianField::Object(), fPos(pos)
{
  fFieldXYZ.push_back(bx);
  fFieldXYZ.push_back(by);
  fFieldXYZ.push_back(bz);
}



FieldDistribution* LxBFieldAux::CreateFieldDistribution(const std::string &fmodel, const std::string &params)
{
  FieldDistribution *bfd = 0;
  std::stringstream paramstr(params.data());
  std::string vunits;
  if (fmodel == "const") {
    double cmin, cmax, vu;
    paramstr >> cmin >> cmax >> vunits;
    vu = dd4hep::_toDouble(vunits);
    bfd = new FieldConst(cmin*vu, cmax*vu);
  }
  if (fmodel == "f_fd") {
    double cmin, cmax, tmin, tmax, vu;
    paramstr >> cmin >> cmax >> tmin >> tmax >> vunits;
    vu = dd4hep::_toDouble(vunits);
    bfd = new FieldFD(cmin*vu, cmax*vu, tmin*vu, tmax*vu);
  }
  if (fmodel == "f_err") {
    double cmin, cmax, tmin, tmax, vu;
    paramstr >> cmin >> cmax >> tmin >> tmax >> vunits;
    vu = dd4hep::_toDouble(vunits);
    bfd = new FieldErrF(cmin*vu, cmax*vu, tmin*vu, tmax*vu);
  }
  if (fmodel == "cylinder") {
    double rmax, vu;
    paramstr >> rmax >> vunits;
    vu = dd4hep::_toDouble(vunits);
    bfd = new FieldCylinder2D(rmax*vu);
  }
  return bfd;
}


////////////////////////////////////////////////////////////////////
LxDipoleFields::LxDipoleFields(dd4hep::Detector& description, const std::string& name)
{
//   LoadFieldsConfigurations(description, name);
}



void LxDipoleFields::PrintBFieldModel()
{
  std::cout << "======= Following setting are used to create magnetic fields =======\n";
  for (const auto &mitr : fBFieldModelsInfo) {
    const auto &tinfo = mitr.second;
    std::cout << "MagnetID: " << mitr.first << std::endl;
    for (const auto &vitr : tinfo) {
      std::cout << "     " << std::get<0>(vitr) << "  " << std::get<1>(vitr)
             << "  " << std::get<2>(vitr) << "  " << std::get<3>(vitr) << std::endl;
    }
  }
  std::cout << "====================================================================\n";
}



LxBField* LxDipoleFields::GetLxDipoleField(dd4hep::Detector& description, const std::string& name)
{
  LoadFieldsConfigurations(description, name);
  using FnMap = std::function<LxBField* (dd4hep::Detector&)>;
  std::map<std::string, FnMap> magMap;
  magMap["Brems"]  = [this](dd4hep::Detector& dct)->LxBField* { return AddBremsMagField(dct); };
  magMap["IP"]     = [this](dd4hep::Detector& dct)->LxBField* { return AddIPMagField(dct);    };
  magMap["Gamma"]  = [this](dd4hep::Detector& dct)->LxBField* { return AddGammaMagField(dct); };

  return magMap[name](description);
}



void LxDipoleFields::LoadFieldsConfigurations(dd4hep::Detector& description, const std::string& name)
{
  // This map can be reduced to vector, but the names in c++ code and in mac settings for the magnets are different. This maps one to another.
  std::map<std::string, std::string> fieldNames {{"Brems", "DumpMag"}, {"IP", "IPMag"}, {"Gamma", "GammaMag"}};
  std::vector<std::string> fieldComponetNames {"Bx", "By", "Bz"};
  std::vector<std::string> fieldDistribConstNames {"x", "y", "z"};

  auto mag = fieldNames.find(name);
  if (mag == fieldNames.end()) {
    dd4hep::except("LxDipoleFields::LoadFieldsConfigurations","No configuration for the magnet %s", name.c_str());
  }

  for (const auto &fcval : fieldComponetNames) {
    std::string readCName = mag->second + "Field" + fcval;
    std::string cparams = description.constant<std::string>(readCName);
    fBFieldModelsInfo[name].push_back(std::make_tuple(fcval, "bvalue", "bvalue", cparams));

    for (const auto &fdistr : fieldDistribConstNames) {
      std::string readDName = mag->second + "Distrib" + fcval + fdistr;
      std::string newValue = description.constant<std::string>(readDName);
      std::istringstream istr(newValue);
      std::string fmodel;
      istr >> fmodel;
      std::string dparams(newValue.c_str() + istr.tellg());
      fBFieldModelsInfo[name].push_back(std::make_tuple(fcval, fdistr, fmodel, dparams));
    }
  }

  PrintBFieldModel();
}


LxBField* LxDipoleFields::AddBremsMagField(dd4hep::Detector& description)
{
  double DumpMagnetXpos = description.constant<double>("DumpMagnetXpos");
  double DumpMagnetYpos = description.constant<double>("DumpMagnetYpos");
  double DumpMagnetZpos = description.constant<double>("DumpMagnetZpos");
  double DumpMagFieldY = description.constant<double>("DumpMagFieldY");

  std::string magid = "Brems";

  dd4hep::Position trm(DumpMagnetXpos, DumpMagnetYpos, DumpMagnetZpos);
  LxBField *ipfield = ComposeFieldObject(magid, trm, DumpMagFieldY);

  return ipfield;
}



LxBField* LxDipoleFields::AddGammaMagField(dd4hep::Detector& description)
{
  double GammaMagnetXpos = description.constant<double>("GMagnetXpos");
  double GammaMagnetYpos = description.constant<double>("GMagnetYpos");
  double GammaMagnetZpos = description.constant<double>("GMagnetZpos");
  double GammaMagFieldY = description.constant<double>("GMagFieldY");

  std::string magid = "Gamma";

  dd4hep::Position trm(GammaMagnetXpos, GammaMagnetYpos, GammaMagnetZpos);
  LxBField *ipfield = ComposeFieldObject(magid, trm, GammaMagFieldY);

  return ipfield;
}



LxBField* LxDipoleFields::AddIPMagField(dd4hep::Detector& description)
{
  double IPMagnetXpos = description.constant<double>("IPMagnetXpos");
  double IPMagnetYpos = description.constant<double>("IPMagnetYpos");
  double IPMagnetZpos = description.constant<double>("IPMagnetZpos");
  double IPMagFieldY = description.constant<double>("IPMagFieldY");

  std::string magid = "IP";

  dd4hep::Position trm(IPMagnetXpos, IPMagnetYpos, IPMagnetZpos);
  LxBField *ipfield = ComposeFieldObject(magid, trm, IPMagFieldY);

  return ipfield;
}



LxBField* LxDipoleFields::ComposeFieldObject(const std::string magid, const dd4hep::Position magpos, const double bval)
{
  double bxval, byval, bzval;
  bxval = byval = bzval = 0.0;

  if (magid == "IP") byval = bval;
  if (magid == "Brems") bxval = bval;
  if (magid == "Gamma") byval = bval;

  typedef VFieldComponentDistrib<FieldDistribution, FieldDistribution, FieldDistribution>  LxTFieldComponent;

  std::vector<FieldDistribution*> bcj(9,0);
  const auto msit = fBFieldModelsInfo.find(magid);
  if (msit != fBFieldModelsInfo.end()) {
    const auto &msetv = msit->second;
    for (auto svitr = msetv.cbegin(); svitr != msetv.cend(); ++svitr) {
      const std::string &fcomponent = std::get<0>(*svitr);
      const std::string &coordinate = std::get<1>(*svitr);
      const std::string &dmodel = std::get<2>(*svitr);
      const std::string &params = std::get<3>(*svitr);

      if (fcomponent == "Bx") {
        if (coordinate == "x") {
          bcj[0] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "y") {
          bcj[1] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "z") {
          bcj[2] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "bvalue") {
          bxval = dd4hep::_toDouble(params.data());
        }
      }

      if (fcomponent == "By") {
        if (coordinate == "x") {
          bcj[3] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "y") {
          bcj[4] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "z") {
          bcj[5] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "bvalue") {
          byval = dd4hep::_toDouble(params.data());
        }
      }

      if (fcomponent == "Bz") {
        if (coordinate == "x") {
          bcj[6] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "y") {
          bcj[7] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "z") {
          bcj[8] = LxBFieldAux::CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "bvalue") {
          bzval = dd4hep::_toDouble(params.data());
        }
      }

    }
  }

  std::for_each(bcj.begin(), bcj.end(), [](FieldDistribution* &bd) {if (!bd) bd = new FieldDistribution(); });

  LxTFieldComponent *bx = new LxTFieldComponent(bcj[0], bcj[1], bcj[2], bxval);
  LxTFieldComponent *by = new LxTFieldComponent(bcj[3], bcj[4], bcj[5], byval);
  LxTFieldComponent *bz = new LxTFieldComponent(bcj[6], bcj[7], bcj[8], bzval);

  LxBField *bfield = new LxBField(bx, by, bz, magpos);

  std::cout << "=========  " << magid << " magnet settings" << "  =========" << std::endl;
  std::cout << "Bx: x: " << typeid(*(bcj[0])).name()
          << "  y: " << typeid(*(bcj[1])).name()
          << "  z: " << typeid(*(bcj[2])).name() << "  Value: " << bxval/dd4hep::tesla << "T" << std::endl;
  std::cout << "By: x: " << typeid(*(bcj[3])).name()
          << "  y: " << typeid(*(bcj[4])).name()
          << "  z: " << typeid(*(bcj[5])).name() << "  Value: " << byval/dd4hep::tesla << "T"  << std::endl;
  std::cout << "Bz: x: " << typeid(*(bcj[6])).name()
          << "  y: " << typeid(*(bcj[7])).name()
          << "  z: " << typeid(*(bcj[8])).name() << "  Value: " << bzval/dd4hep::tesla << "T"  << std::endl;
  std::cout << "====================================" << std::endl;

 return bfield;
}



/////////////////////////////////////////////////////////////
// Sample the field along z at x = y = 0, through the CartesianField handle.
static void dumpFieldAlongZ(const dd4hep::CartesianField& field,
                            int    n_points = 500,
                            double z_min    = -8.0 * dd4hep::m,
                            double z_max    = 20.0 * dd4hep::m) {
  if (!field.isValid() || n_points < 1) return;

  const double dz = (n_points > 1) ? (z_max - z_min) / (n_points - 1) : 0.0;

  std::printf("# %10s %14s %14s %14s %14s\n",
              "z[m]", "Bx[T]", "By[T]", "Bz[T]", "|B|[T]");

  for (int i = 0; i < n_points; ++i) {
    const double z      = z_min + i * dz;      // DD4hep units (mm)
    const double pos[3] = {0.0, 0.0, z};
    double       b[3]   = {0.0, 0.0, 0.0};     // zero first: value() does NOT clear

    field.value(pos, b);                       // handle → assigned object's fieldComponents

    const double bx = b[0] / dd4hep::tesla;
    const double by = b[1] / dd4hep::tesla;
    const double bz = b[2] / dd4hep::tesla;
    const double bmag = std::sqrt(bx*bx + by*by + bz*bz);

//     std::printf("  %10.4f %14.6e %14.6e %14.6e %14.6e\n",
//                 z / dd4hep::m, bx, by, bz, bmag);
    std::printf("  %10.4f %14.6e\n",
                z / dd4hep::m, bmag);

  }
}


// Debug helper: sample the field along z at x = y = 0.
// Calls the object directly, so it tests only your class's fieldComponents().
static void dumpFieldAlongZ(LxBField* obj,
                            int    n_points = 100,
                            double z_min    = -20.0 * dd4hep::m,
                            double z_max    = 20.0 * dd4hep::m) {
  if (!obj || n_points < 1) return;

  const double dz = (n_points > 1) ? (z_max - z_min) / (n_points - 1) : 0.0;

  std::printf("# %10s %14s %14s %14s %14s\n",
              "z[m]", "Bx[T]", "By[T]", "Bz[T]", "|B|[T]");

  for (int i = 0; i < n_points; ++i) {
    const double z      = z_min + i * dz;       // DD4hep units (mm)
    const double pos[3] = {0.0, 0.0, z};
    double       b[3]   = {0.0, 0.0, 0.0};      // zero before sampling

    obj->fieldComponents(pos, b);               // your class fills b (DD4hep units)

    const double bx = b[0] / dd4hep::tesla;
    const double by = b[1] / dd4hep::tesla;
    const double bz = b[2] / dd4hep::tesla;
    const double bmag = std::sqrt(bx*bx + by*by + bz*bz);

    std::printf("  %10.4f %14.6e %14.6e %14.6e %14.6e\n",
                z / dd4hep::m, bx, by, bz, bmag);
  }
}



static dd4hep::CartesianField create_LxCustomField(dd4hep::Detector& description, xml_h e)
{
  xml_comp_t x_field(e);
  std::string mname =  x_field.nameStr();
  std::cout << "Creating filed for magnet " << mname << std::endl;
  LxDipoleFields df(description, mname);
  LxBField* obj = df.GetLxDipoleField(description, mname);
  obj->field_type = dd4hep::CartesianField::MAGNETIC;
//   dumpFieldAlongZ(obj);

  dd4hep::CartesianField field;
  field.assign(obj, x_field.nameStr(), x_field.typeStr());
//   dumpFieldAlongZ(field);

  return field;
}

//============================================================================
// IMPORTANT — DD4hep field accumulation contract (undocumented in the base API)
//
// This method MUST ADD its contribution to MagField[] (use +=), it must NOT
// overwrite it (=).  The double* out-parameter looks like "fill this in", but
// DD4hep's semantics are different:
//
//   dd4hep::OverlayedField (DDCore/src/Fields.cpp) holds a list of registered
//   magnetic fields.  On each query it zeroes the output array ONCE, then calls
//   every field's fieldComponents() on that SAME array, in sequence, without
//   re-zeroing between calls and without summing into a temporary.  Each field
//   is therefore expected to accumulate into whatever is already there.
//
//   -> With '=', every field overwrites the previous one's contribution.  A
//      region-limited field returns 0 outside its region, so the fields called
//      after this one clobber our value with 0.  Net effect: only the LAST
//      registered field is ever visible.  (This is a silent failure: it
//      compiles, runs, and produces a plausible-looking single field.)
//   -> With '+=', region-limited fields add 0 outside their region and leave
//      earlier contributions intact, so spatially-separated fields superpose
//      correctly.
//
// This contract is not stated on CartesianField::Object::fieldComponents.  It is
// only observable (a) in the OverlayedField loop above, and (b) in that every
// built-in field in DDCore/src/FieldTypes.cpp (Constant, Solenoid, Dipole,
// Multipole) uses '+='.
//
// Corollary of the same rule: because nothing re-zeroes between fields, any code
// that calls fieldComponents() directly on a single field (e.g. a debug dump)
// must zero the output array itself before the call.  fieldComponents() adds; it
// never clears.

//   dd4hep::OverlayedField::combinedMagnetic() (DDCore/src/Fields.cpp) zeroes the
//   output array ONCE, then calls calculate_combined_field(), which loops
//   `for (const auto& i : v) i.value(pos, field);` over all registered fields on
//   that SAME array -- no re-zeroing between calls, no temporary, no external sum.
//   CartesianField::value() forwards straight to fieldComponents() without
//   clearing.  Each field is therefore expected to accumulate into field[].
//============================================================================

DECLARE_XMLELEMENT(LxDipoleField, create_LxCustomField)


