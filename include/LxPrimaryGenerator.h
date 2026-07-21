//
/// \brief Definition of the LxPrimaryGenerator class
//
/// DDG4 port of the Geant4 PrimaryGeneratorAction. One Geant4GeneratorAction covers all
/// beam modes (gaussian, mono, MC from text/HDF5/ROOT-ntuple). The old G4UImessenger is
/// replaced by declareProperty() knobs with ctor defaults, so a bare run needs no
/// steering; overrides are optional. Per-event MC-primary info is published as an
/// LxPrimaryInfo event extension in every mode; a separate output event action consumes
/// it (that replaces the old SavePrimaryTrack()/EventInfo path).

#ifndef LxPrimaryGenerator_h
#define LxPrimaryGenerator_h 1

#include <string>
#include <vector>
#include <cstddef>

#include "DDG4/Geant4GeneratorAction.h"

#include "globals.hh"
#include "G4ThreeVector.hh"

#include "LxPrimaryInfo.h"

class G4Event;
class G4ParticleGun;
class PrimarySpectra;
class LuxeMCGenerator;

namespace dd4hep {
namespace sim {

class LxPrimaryGenerator : public dd4hep::sim::Geant4GeneratorAction
{
public:
  LxPrimaryGenerator(dd4hep::sim::Geant4Context* ctxt, const std::string& nam);
  virtual ~LxPrimaryGenerator();

  /// DDG4 generator callback (was GeneratePrimaries in the G4 version)
  virtual void operator()(G4Event* anEvent) override;

public:
  void ApplyProperties();
  void SetDefaultKinematic();
  void SetBeamType(G4String val);
  void SetSpectraFile(G4String val);
  void SetMCParticleFile(G4String val, const G4bool list = false);

  void EmitInteraction(G4Event* /*anEvent*/);
  void GeneratePrimaries(G4Event* anEvent);
  void GenerateGaussian(G4Event* anEvent);
  void GenerateMono(G4Event* anEvent);
  void GeneratefromMC(G4Event* anEvent);

  const G4ParticleGun* GetParticleGun() const {return fParticleGun;}

  G4double GetX0() const {return fx0;}
  G4double GetY0() const {return fy0;}
  G4double GetZ0() const {return fz0;}
  G4double GetSigmaX() const {return fsigmax;}
  G4double GetBetaX()  const {return fbetax;}
  G4double GetEmittanceX() const {return femittancex;}
  G4double GetBeamFocus() const {return fbeamfocus;}

  G4double GetSigmaY() const {return fsigmay;}
  G4double GetBetaY()  const {return fbetay;}
  G4double GetEmittanceY() const {return femittancey;}
  tBeamType GetBeamType() const {return fBeamType;}
  G4String GetMCfile() const {return fMCfile;}
  const PrimarySpectra *GetSpectraSettings() const {return fSpectra;}

  void SetSigmaX(const G4double sigma);
  void SetSigmaY(const G4double sigma);
  void SetSigmaZ(const G4double sigma);
  void SetBeamPosZ(const G4double posz);
  void SetBeamFocus(const G4double zfocus) { fbeamfocus = zfocus; }
  void SetMCWeightScale(const G4int scale) { fMCWeightScale = scale; }
  void SetBeamPosition(const G4ThreeVector pos){ fx0=pos.x(); fy0=pos.y(); fz0=pos.z(); fHasBeamPos=true; }
  void SetSkipEvents(const size_t nev) {fSkipEvents = nev; }

  G4double GetSigmaZ() const {return fsigmaz;}
  size_t GetSkipEvents() const {return fSkipEvents;}

  void AddMCSelectParticle(const G4int pdgid);

protected:
  /// One-time setup on the first event: apply property values through the setters
  /// (this replaces what the messenger used to do at PreInit/Idle).
  void Prepare();
  /// Build and attach the per-event LxPrimaryInfo extension (every mode).
  void AttachPrimaryInfo();

  G4double TestHitTarget(const std::vector <double> &pp, const double *vtx);
  void InitReader();

  /// dd4hep(compact) -> Geant4(CLHEP) unit rebasing. Identity on Geant4-unit builds.
  static double toG4Length(double v);
  static double toG4Energy(double v);

private:
  G4ParticleGun*         fParticleGun;
  tBeamType              fBeamType;
  G4double               fx0, fy0, fz0;
  G4double               fsigmax, fsigmay, fsigmaz;
  G4double               femittancex, femittancey;
  G4double               fbetax, fbetay;
  G4double               fbeamfocus;

  static LuxeMCGenerator   *lxgen;
  PrimarySpectra        *fSpectra;
  G4String               fMCfile;
  G4bool                 fMCList;
  G4int                  fnfixparticles;

  G4int                  fMCWeightScale;
  LxPrimaryInfo          fMCInfo;         ///< persistent scratch (was EventInfo* fMCEventInfo)

  std::vector<G4int>     fselectpdgvec;
  G4bool                 fselectpdg;

  size_t                 fSkipEvents;

  G4int                  fMCPrimeParentId;
  G4int                  fMCPrimeNDescendant;
  G4double               fMCPrimeXi;
  G4int                  fMCPhysProc;

  // ---- geometry constants, read from the compact <define> in the ctor and
  //      rebased to Geant4/CLHEP units (were fDetector->GetzstartAbs() and
  //      LXSetUp::Instance()->GTargetX / ->GTargetZpos) --------------------------
  G4double               fzstartAbs;
  G4double               fGTargetX;
  G4double               fGTargetZpos;

  // ---- steering properties (mirror the old /lxphoton/gun/ messenger commands).
  //      Length/position values arrive in dd4hep units and are rebased in Prepare().
  bool                   fPrepared;
  std::string            fBeamTypeProp;
  std::string            fSpectraFileProp;
  std::string            fMCFileProp;
  std::string            fMCFileListProp;
  double                 fSigmaXProp;
  double                 fSigmaYProp;
  double                 fSigmaZProp;
  double                 fPosZProp;
  std::vector<double>    fPositionProp;
  std::vector<int>       fSelectPdgProp;
  int                    fSkipEventsProp;
  double                 fEnergyProp;        // NaN -> keep 16.5 GeV
  double                 fEnergy;            // because of units in dd4hep and G4
  std::string            fParticleName;

  G4ThreeVector          fDir;
  bool                   fHasBeamPos;
  std::vector<double>    fDirectionProp;
  // "Mask" is required by ddsim, not by us: DD4hepSimulation.run() unconditionally
  // does gen.Mask = <index> on every userInputPlugin generator. The base
  // Geant4GeneratorAction does not declare it, so we must -- otherwise ddsim throws
  // "Cannot set Primary.Mask". The value (ddsim overwrites it) is the interaction
  // mask keying the primary-interaction/merger chain; inert for our direct-gun path.
  int                    fMask;
};
}}
#endif
