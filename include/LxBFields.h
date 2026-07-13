
#ifndef LxBFields_h
#define LxBFields_h 1

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DD4hep/Objects.h"   // dd4hep::Position ( == ROOT::Math::XYZVector )


////////////////////////////////////////////////////////////////////
class FieldDistribution
{
public:
  FieldDistribution() {};
  virtual ~FieldDistribution() {};
  virtual double GetField(const double *xx = 0) { return 1.0; }
  virtual bool IsLimited() const {return false;}
};


class FieldConst: public FieldDistribution
{
public:
  FieldConst(): FieldDistribution(), fxmin(0.0), fxmax(0.0) {};
  FieldConst(const double xmin, const double xmax) : FieldDistribution(), fxmin(xmin), fxmax(xmax) {};
  virtual ~FieldConst() {};
  double GetField(const double *xx) { if ((xx[0] >= fxmin) && (xx[0] <= fxmax)) return 1.0; return 0.0; }
  bool IsLimited() const {return true;}

protected:
  double fxmin, fxmax;
};


class FieldFD: public FieldDistribution
{
public:
  FieldFD(): FieldDistribution(), fxmin(0.0), fxmax(0.0), ftmin(1.0), ftmax(1.0) {};
  FieldFD(const double xmin, const double xmax, const double tmin, const double tmax) :
     FieldDistribution(), fxmin(xmin), fxmax(xmax), ftmin(tmin), ftmax(tmax) {};
  virtual ~FieldFD() {};
  double GetField(const double *xx) {
    const double x = xx[0];
    return 1.0 / ( (1.0 + exp((fxmin-x)/ftmin)) * (1.0 + exp((x-fxmax)/ftmax)) );
  }
  bool IsLimited() const {return true;}

protected:
  double fxmin, fxmax, ftmin, ftmax;
};


class FieldErrF: public FieldDistribution
{
public:
  FieldErrF(): FieldDistribution(), fxmin(0.0), fxmax(0.0), ftmin(1.0), ftmax(1.0) {};
  FieldErrF(const double xmin, const double xmax, const double tmin, const double tmax) :
     FieldDistribution(), fxmin(xmin), fxmax(xmax), ftmin(tmin), ftmax(tmax) {};
  virtual ~FieldErrF() {};
  double GetField(const double *xx) {
    const double x = xx[0];
    return (1.0 + erf((x-fxmin)/ftmin)) * (1.0 + erf((fxmax-x)/ftmax));
  }
  bool IsLimited() const {return true;}

protected:
  double fxmin, fxmax, ftmin, ftmax;
};



class FieldData2D: public FieldDistribution
{
public:
  FieldData2D(): FieldDistribution(), fNsort(10) {};
  FieldData2D(const std::string fname);

  virtual ~FieldData2D() {};
  double GetField(const double *const xx) { return  DataField(xx[0], xx[1]); }
  bool IsLimited() const {return true;}

protected:
  void LoadFieldData(const std::string fname);
  double FindZ(const dd4hep::Position &r1, const dd4hep::Position &r2, const dd4hep::Position &r3,
                               const double xx, const double yy) const;
  bool CheckCollinear(const dd4hep::Position &r1, const dd4hep::Position &r2, const dd4hep::Position &r3) const;
  double DataField(const double xx, const double yy) const;

  template <class T> class point_cmp
  {
    public:
    point_cmp(const std::vector<T> *ptr): p(ptr) {};
    bool operator() (size_t i, size_t j) const { if ((i < p->size()) && (j < p->size()) ) return (p->at(i) < p->at(j));
                                                 return false; }
    const std::vector<T> *p;
  };

protected:
  std::string fFName;
  dd4hep::Position frMin, frMax;

  std::vector<dd4hep::Position> fFieldData;
  size_t fNsort;
};



class FieldCylinder2D: public FieldDistribution
{
public:
  FieldCylinder2D(): FieldDistribution(), fR(0.0) {};
  FieldCylinder2D(const double rmax): FieldDistribution(), fR(rmax) {};

  virtual ~FieldCylinder2D() {};
  double GetField(const double *const xx) { return  sqrt(xx[0]*xx[0] + xx[1]*xx[1]) > fR ? 0.0 : 1.0; }
  bool IsLimited() const {return true;}

protected:
  double fR;
};


/////////////////////////////////////////////////////////////
class VFieldComponent
{
public:
  VFieldComponent() : fBField(0.0) {};
  VFieldComponent(const double bf) : fBField(bf) {};
  virtual ~VFieldComponent() {};
  virtual double GetField(const double *const x = 0)  const { return fBField; }
  virtual bool IsLimited() const {return false;}
  virtual bool IsZLimited() const {return false;}

protected:
  double  fBField;
};


template <class FM2, class FM1, int C0, int C1, int C2>
class VFieldComponentDistrib2D : public VFieldComponent
{
public:
  VFieldComponentDistrib2D() : fFDistrib2(0), fFDistrib1(0) {};
  VFieldComponentDistrib2D(FM2 *fcc, FM1 *fc, const double bf) : VFieldComponent(bf), fFDistrib2(fcc), fFDistrib1(fc) {};
  virtual ~VFieldComponentDistrib2D() {
     if (fFDistrib2) delete fFDistrib2;
     if (fFDistrib1) delete fFDistrib1;
     fFDistrib2 = 0; fFDistrib1 = 0;
   };
  double GetField(const double *const x)  const
  {
    double xy[2] = {x[C0], x[C1]};
    double xx[1] = {x[C2]};
    return fBField * fFDistrib2->GetField(xy) * fFDistrib1->GetField(xx);
  }
  bool IsLimited() const {return fFDistrib2->IsLimited() && fFDistrib1->IsLimited();}
  bool IsZLimited() const {return (C2==2) ? fFDistrib1->IsLimited() : fFDistrib2->IsLimited();}

  FM2 *fFDistrib2;
  FM1 *fFDistrib1;
};


template <class FMX, class FMY, class FMZ>
class VFieldComponentDistrib : public VFieldComponent
{
public:
  VFieldComponentDistrib() : fFDistribX(0), fFDistribY(0), fFDistribZ(0) {};
  VFieldComponentDistrib(FMX *fcx, FMY *fcy, FMZ *fcz, const double bf) : VFieldComponent(bf),
                         fFDistribX(fcx), fFDistribY(fcy), fFDistribZ(fcz) {};
  virtual ~VFieldComponentDistrib()
   { if (fFDistribX) delete fFDistribX;
     if (fFDistribY) delete fFDistribY;
     if (fFDistribZ) delete fFDistribZ;
     fFDistribX = 0; fFDistribY = 0; fFDistribZ = 0;
   };
  double GetField(const double *const x) const {
    return fBField * fFDistribX->GetField(x) * fFDistribY->GetField(&x[1]) * fFDistribZ->GetField(&x[2]);
  }
  bool IsLimited() const {return fFDistribX->IsLimited() && fFDistribY->IsLimited() && fFDistribZ->IsLimited();}
  bool IsZLimited() const {return fFDistribZ->IsLimited();}

  FMX *fFDistribX;
  FMY *fFDistribY;
  FMZ *fFDistribZ;
};



///////////////////////////////////////////////////////
class LxBField: public dd4hep::CartesianField::Object
{
public:
  LxBField() : dd4hep::CartesianField::Object(), fFieldXYZ(3), fPos(dd4hep::Position()) {};
  LxBField(const VFieldComponent *bx, const VFieldComponent *by, const VFieldComponent *bz, const dd4hep::Position pos);
  ~LxBField() {for (auto &vv : fFieldXYZ) { if (vv) delete vv; } fFieldXYZ.clear(); };

  void fieldComponents(const double* ppos, double* MagField) override {
    dd4hep::Position rr(ppos[0], ppos[1], ppos[2]);
    dd4hep::Position rloc = rr - fPos;
    double rc[3];
    rloc.GetCoordinates(rc);                 // rc = { x, y, z }
    for (int ii = 0; ii < 3; ++ii) {
      MagField[ii] += fFieldXYZ[ii]->GetField(rc);
    }
  }

  bool IsZLimited() const {
    bool zl = true;
    for (const auto &bf : fFieldXYZ) { zl = zl && bf->IsZLimited(); }
     return zl;
 }

protected:

  std::vector<const VFieldComponent*>  fFieldXYZ;
  dd4hep::Position fPos;
};



///////////////////////////////////////////////////////
class LxBFieldAux
{
public:
  LxBFieldAux() {};
//   static void AddFieldToLogVolumes(G4FieldManager* fieldMgr, const std::vector<std::string> vname);
  static FieldDistribution* CreateFieldDistribution(const std::string &fmodel, const std::string &params);
//   static void ConfigureFieldManager(G4FieldManager *fieldMgr, G4MagneticField *mfield);
};



///////////////////////////////////////////////////////
class LxDipoleFields
{
public:
  LxDipoleFields(dd4hep::Detector& description, const std::string& name);

  LxBField* GetLxDipoleField(dd4hep::Detector& description, const std::string& name);
  void LoadFieldsConfigurations(dd4hep::Detector& description, const std::string& name);
  LxBField* ComposeFieldObject(const std::string magid, const dd4hep::Position magpos, const double bval);
  LxBField* AddBremsMagField(dd4hep::Detector& description);
  LxBField* AddIPMagField(dd4hep::Detector& description);
  LxBField* AddGammaMagField(dd4hep::Detector& description);
  void PrintBFieldModel();

protected:
  std::map<std::string, std::vector<std::tuple<std::string, std::string, std::string, std::string> > > fBFieldModelsInfo;

};

#endif
