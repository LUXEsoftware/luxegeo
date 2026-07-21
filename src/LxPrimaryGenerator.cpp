//
/// \brief Implementation of the LxPrimaryGenerator class (DDG4 port)
//

#include <algorithm>
#include <string>
#include <cmath>
#include <limits>

#include "LxPrimaryGenerator.h"

#include "DDG4/Factories.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/InstanceCount.h"
#include "DDG4/Geant4Primary.h"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4AutoLock.hh"
#include "globals.hh"

#include "PrimarySpectra.h"
#include "LuxeTestGenerator.h"

// using dd4hep::sim::Geant4Context;
// using dd4hep::sim::Geant4GeneratorAction;

namespace dd4hep {
namespace sim {
namespace { G4Mutex MCReadMutex = G4MUTEX_INITIALIZER; }

LuxeMCGenerator *LxPrimaryGenerator::lxgen = 0;


double LxPrimaryGenerator::toG4Length(double v) { return v / dd4hep::mm  * CLHEP::mm;  }
double LxPrimaryGenerator::toG4Energy(double v) { return v / dd4hep::GeV * CLHEP::GeV; }


LxPrimaryGenerator::LxPrimaryGenerator(Geant4Context* ctxt, const std::string& nam)
: Geant4GeneratorAction(ctxt, nam),
  fParticleGun(0), fBeamType(beamGauss),
  fx0(0.0), fy0(0.0), fz0(0.0), fsigmax(5.0*CLHEP::um), fsigmay(5.0*CLHEP::um), fsigmaz(24.0*CLHEP::um),
  femittancex(0.0), femittancey(0.0), fbetax(0.0), fbetay(0.0), fbeamfocus(0.0), fSpectra(0),
  fMCfile(""), fMCList(false), fnfixparticles(0),
  fselectpdgvec(0), fselectpdg(false), fSkipEvents(0), fMCPrimeParentId(-1), fMCPrimeNDescendant(-1),
  fMCPrimeXi(-1.0), fMCPhysProc(0),
  fzstartAbs(0.0), fGTargetX(0.0), fGTargetZpos(0.0),
  fPrepared(false),
  fPositionProp(), fSelectPdgProp(),
  fEnergy(16.5*CLHEP::GeV), fParticleName("e-"), fDir(0.,0.,1.), fHasBeamPos(false), fDirectionProp()
{
  declareProperty("BeamType",      fBeamTypeProp = "gaussian");
  declareProperty("SpectraFile",   fSpectraFileProp = "");
  declareProperty("MCFile",        fMCFileProp = "");
  declareProperty("MCFileList",    fMCFileListProp = "");
  declareProperty("SigmaX",        fSigmaXProp = 5.0 *dd4hep::um);
  declareProperty("SigmaY",        fSigmaYProp = 5.0 *dd4hep::um);
  declareProperty("SigmaZ",        fSigmaZProp = 24.0 *dd4hep::um);
  declareProperty("PosZ",          fPosZProp = std::numeric_limits<double>::quiet_NaN());
  declareProperty("Position",      fPositionProp);
  declareProperty("MCWeightScale", fMCWeightScale = 1);
  declareProperty("SelectPdg",     fSelectPdgProp);
  declareProperty("SkipEvents",    fSkipEventsProp = 0);
  declareProperty("Particle",      fParticleName = "e-");
  declareProperty("Energy",        fEnergyProp = std::numeric_limits<double>::quiet_NaN());
  declareProperty("Direction",     fDirectionProp);

  // "Mask" is required by ddsim, not by us: DD4hepSimulation.run() unconditionally
  // does gen.Mask = <index> on every userInputPlugin generator. The base
  // Geant4GeneratorAction does not declare it, so we must -- otherwise ddsim throws
  // "Cannot set Primary.Mask". The value (ddsim overwrites it) is the interaction
  // mask keying the primary-interaction/merger chain; inert for our direct-gun path.
  declareProperty("Mask",          fMask = 4);

//   G4int n_particle = 1;
//   fParticleGun = new G4ParticleGun(n_particle);

  // Geometry numbers formerly from DetectorConstruction / LXSetUp, now compact <define>
  // constants. detectorDescription() is valid at construction time; read once + rebase.
  dd4hep::Detector& description = context()->detectorDescription();
  double BTargetZpos   = toG4Length(description.constant<double>("BTargetZpos"));
  double BTargetZ      = toG4Length(description.constant<double>("BTargetZ"));
  fzstartAbs = BTargetZpos-0.5*BTargetZ;
  fGTargetX    = toG4Length(description.constant<double>("GTargetX"));
  fGTargetZpos = toG4Length(description.constant<double>("GTargetZpos"));

  // fMCInfo default-constructs with weight 1.0, trackId -1 (was new EventInfo(1.0,-1)).

  dd4hep::InstanceCount::increment(this);
  std::cout << "=================  LxPrimaryGenerator::LxPrimaryGenerator done =================\n";
}



LxPrimaryGenerator::~LxPrimaryGenerator()
{
  if (fParticleGun) delete fParticleGun;
  if (fSpectra) delete fSpectra;
  G4AutoLock mcfilemutex(&MCReadMutex);
  if (lxgen) {
    delete lxgen;
    lxgen = 0;
  }
  mcfilemutex.unlock();
  dd4hep::InstanceCount::decrement(this);
}



void LxPrimaryGenerator::ApplyProperties()
{
  SetBeamType(fBeamTypeProp);

  if (!fSpectraFileProp.empty())  SetSpectraFile(fSpectraFileProp);

  if (!std::isnan(fSigmaXProp))   fsigmax = toG4Length(fSigmaXProp);
  if (!std::isnan(fSigmaYProp))   fsigmay = toG4Length(fSigmaYProp);

  if (!std::isnan(fPosZProp))     SetBeamPosZ(toG4Length(fPosZProp));
  if (fPositionProp.size() == 3)
    SetBeamPosition(G4ThreeVector(toG4Length(fPositionProp[0]),
                                  toG4Length(fPositionProp[1]),
                                  toG4Length(fPositionProp[2])));

  if (fDirectionProp.size() == 3) {
    G4ThreeVector d(fDirectionProp[0], fDirectionProp[1], fDirectionProp[2]);
    if (d.mag2() > 0.0) fDir = d.unit();     // dimensionless -> just normalize
  }

  if (!fMCFileListProp.empty())   SetMCParticleFile(fMCFileListProp, true);
  else if (!fMCFileProp.empty())  SetMCParticleFile(fMCFileProp, false);

  for (G4int pdgid : fSelectPdgProp) AddMCSelectParticle(pdgid);

  if (fSkipEventsProp > 0)        fSkipEvents = static_cast<size_t>(fSkipEventsProp);

  if (!std::isnan(fEnergyProp)) fEnergy       = toG4Energy(fEnergyProp);   // dd4hep -> CLHEP
}



void LxPrimaryGenerator::SetDefaultKinematic()
{
  // Particle & energy: no property in this generator -> defaults ("the rest as is").
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(fParticleName);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(fEnergy);

  // z: property (PosZ/Position) wins; otherwise the drift default.
  if (!fHasBeamPos) fz0 = fzstartAbs - 20.0*cm;

  // Beam optics from kinetic_energy/mass and the (property-aware) sigmas.
  G4double lf = 1.0 + fParticleGun->GetParticleEnergy() / particle->GetPDGMass();
  femittancex = 1.4e-3 * mm / lf;
  femittancey = 1.4e-3 * mm / lf;
  fbetax = fsigmax*fsigmax/femittancex;
  fbetay = fsigmay*fsigmay/femittancey;

  // Push members into the gun.
  fParticleGun->SetParticleMomentumDirection(fDir);
  fParticleGun->SetParticlePosition(G4ThreeVector(fx0, fy0, fz0));
}



void LxPrimaryGenerator::Prepare()
{
  if (fPrepared) return;

  // Gun created here, not in the ctor (needs a ready G4ParticleTable).
  if (!fParticleGun) fParticleGun = new G4ParticleGun(1);

  ApplyProperties();        // 1) properties -> member variables (no gun access)
  SetDefaultKinematic();    // 2) members -> gun, plus derived beam optics

  fPrepared = true;
}



void LxPrimaryGenerator::operator()(G4Event* anEvent)
{
  Prepare();                   // one-time, on the first event
  GeneratePrimaries(anEvent);  // configures fParticleGun (mode dispatch); NO GeneratePrimaryVertex now
  AttachPrimaryInfo();
}


void LxPrimaryGenerator::EmitInteraction(G4Event* /*anEvent*/)
{
  // Mirror Geant4InputAction::operator(): publish a Geant4PrimaryInteraction into the
  // per-event Geant4PrimaryEvent (created upstream by Geant4GenerationInit), so that
  // Geant4PrimaryHandler builds the G4 primaries WITH a backing Geant4Particle and
  // Geant4ParticleHandler::begin() finds them. We do NOT fire the gun into the G4Event.
  Geant4Event&        evt  = context()->event();
  Geant4PrimaryEvent* prim = evt.extension<Geant4PrimaryEvent>();

  Geant4PrimaryInteraction* inter = prim->get(fMask);
  if ( !inter ) {
    inter = new Geant4PrimaryInteraction();
    inter->mask = fMask;
    prim->add(fMask, inter);
  }

  const int nparts = fParticleGun->GetNumberOfParticles();
  if (nparts <= 0) return;                 // filtered-out event -> empty interaction (like primaries.empty())

  G4ParticleDefinition* pdef = fParticleGun->GetParticleDefinition();
  const G4ThreeVector   pos  = fParticleGun->GetParticlePosition();
  const G4ThreeVector   dir  = fParticleGun->GetParticleMomentumDirection().unit();
  const G4double        ekin = fParticleGun->GetParticleEnergy();          // CLHEP (MeV)
  const G4double        mass = pdef->GetPDGMass();                         // CLHEP (MeV)
  const G4double        pmag = std::sqrt(ekin*ekin + 2.0*ekin*mass);
  const G4ThreeVector   mom  = pmag * dir;
  const G4double        tim  = fParticleGun->GetParticleTime();

  // NOTE(units): Geant4Particle/Geant4Vertex are stored in Geant4/CLHEP units (MeV, mm, ns)
  // -- the DDG4 vertex dump is labelled [mm]/[ns]. Gun values are already CLHEP, so no rebase.
  // If momenta/positions come out off by a constant factor, rebase HERE (see end note).

  Geant4Vertex* v = new Geant4Vertex();
  v->mask = fMask;
  v->x = pos.x();  v->y = pos.y();  v->z = pos.z();
  v->time = tim;

  for (int i = 0; i < nparts; ++i) {
    const int pid = i;
    Geant4Particle* p = new Geant4Particle(pid);
    p->mask      = fMask;
    p->pdgID     = pdef->GetPDGEncoding();
    p->charge    = (char) lround(3.0 * pdef->GetPDGCharge() / eplus);   // stored as 3x charge
    p->mass      = mass;
    p->genStatus = 1;                                                   // stable final state
    p->status   |= G4PARTICLE_GEN_STABLE;                              // == PropertyMask::set
    p->vsx = pos.x();  p->vsy = pos.y();  p->vsz = pos.z();
    p->psx = mom.x();  p->psy = mom.y();  p->psz = mom.z();
    p->time = tim;

    v->out.insert(pid);                          // primary emanates from this vertex (the crucial link)
    inter->particles.emplace(pid, p);
  }
  inter->vertices[fMask].emplace_back(v);
}


void LxPrimaryGenerator::AttachPrimaryInfo()
{
  LxPrimaryInfo* info = new LxPrimaryInfo();
  info->fMode        = fBeamType;
  info->fWeight      = fMCInfo.GetWeight();
  info->fMCTrackId   = fMCInfo.GetMCTrackId();
  info->fParentId    = fMCPrimeParentId;
  info->fNDescendant = fMCPrimeNDescendant;
  info->fXi          = fMCPrimeXi;
  info->fPhysProc    = fMCPhysProc;
  info->fNGenerated  = fParticleGun->GetNumberOfParticles();

  G4ParticleDefinition* pd = fParticleGun->GetParticleDefinition();
  info->fPdg  = pd ? pd->GetPDGEncoding() : 0;
  info->fEkin = fParticleGun->GetParticleEnergy();
  info->fTime = fParticleGun->GetParticleTime();
  const G4ThreeVector vtx = fParticleGun->GetParticlePosition();
  const G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  info->fVertex[0]    = vtx.x();  info->fVertex[1]    = vtx.y();  info->fVertex[2]    = vtx.z();
  info->fDirection[0] = dir.x();  info->fDirection[1] = dir.y();  info->fDirection[2] = dir.z();

  // DDG4 owns the extension and deletes it at end of event.
  context()->event().addExtension<LxPrimaryInfo>(info);
}



void LxPrimaryGenerator::SetBeamType(G4String val)
{
  if (val == "gaussian") {
    fBeamType = beamGauss;
  } else if (val == "mono") {
    fBeamType = beamMono;
  } else if (val == "mc") {
    fBeamType = beamMC;
  } else if (val == "mchdf5") {
    fBeamType = beamMCh5;
  } else if (val == "mctupleg4") {
    fBeamType = beamMCTupleG4;
  } else {
    G4cout << "LxPrimaryGenerator::SetBeamType: <" << val << ">"
           << " is not defined. Using default: gaussian."
           << G4endl;
  }
}



void LxPrimaryGenerator::SetSigmaX(const G4double sigma)
{
  fsigmax = sigma;
  fbetax = fsigmax*fsigmax/femittancex;
  std::cout << "LxPrimaryGenerator: Set beam sigmaX at IP to : " << G4BestUnit(fsigmax, "Length") << std::endl;
}



void LxPrimaryGenerator::SetSigmaY(const G4double sigma)
{
  fsigmay = sigma;
  fbetay = fsigmay*fsigmay/femittancey;
  std::cout << "LxPrimaryGenerator: Set beam sigmaY at IP to : " << G4BestUnit(fsigmay, "Length") << std::endl;
}


void LxPrimaryGenerator::SetSigmaZ(const G4double sigma)
{
  fsigmaz = sigma;
  std::cout << "LxPrimaryGenerator: Set beam sigmaZ at IP to : " << G4BestUnit(fsigmaz, "Length") << std::endl;
}


void LxPrimaryGenerator::SetBeamPosZ(const G4double posz)
{
  fz0 = posz;
  std::cout << "LxPrimaryGenerator: Set beam Z position : " << G4BestUnit(fz0, "Length") << std::endl;
  fHasBeamPos = true;
}



void LxPrimaryGenerator::SetSpectraFile(G4String val)
{
  if (fSpectra) delete fSpectra;
  fSpectra = new PrimarySpectra();
  int ndata = fSpectra->LoadData(val);
  fSpectra->SetScale(CLHEP::GeV);
  G4cout << "LxPrimaryGenerator: Primary spectra with " << ndata
         << " data points was loaded from the file " << val << G4endl;
}



void LxPrimaryGenerator::SetMCParticleFile(G4String val, const G4bool list)
{
  G4AutoLock mcfilemutex(&MCReadMutex);
  if (lxgen) { delete lxgen; lxgen = 0; }
  mcfilemutex.unlock();

  fMCfile = val;
  fMCList = list;
  if (fMCList) G4cout << "File with the list of files with primary MC particles " << fMCfile << G4endl;
  else         G4cout << "File with primary MC particles " << fMCfile << G4endl;
}



void LxPrimaryGenerator::AddMCSelectParticle(const G4int pdgid)
{
  if (std::find(fselectpdgvec.begin(), fselectpdgvec.end(), pdgid) == fselectpdgvec.end() ) {
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle(pdgid);
    if (particle) {
      fselectpdgvec.push_back(pdgid);
      fselectpdg = true;
      G4cout << "Particle with PDG_ID = " << pdgid << " was added to the list of selected MC patricles."  << G4endl;
    } else {
      G4cout << "Attempt to add particle with PDG_ID = " << pdgid
             << " to the list of selected MC patricles failed. PDG_ID is not valid! Ignore."  << G4endl;
    }
  }
}



void LxPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  if (fSpectra) fParticleGun->SetParticleEnergy(fSpectra->GetRandom());

  if (fBeamType == beamGauss) {
    GenerateGaussian(anEvent);
  } else if (fBeamType == beamMono) {
    GenerateMono(anEvent);
  }
    else if (fBeamType == beamMC || fBeamType == beamMCh5 || fBeamType == beamMCTupleG4) {
    GeneratefromMC(anEvent);
  }
}



void LxPrimaryGenerator::GenerateGaussian(G4Event* anEvent)
{
  G4double z0 = G4RandGauss::shoot(fz0, fsigmaz);
  G4double zdrift = z0 - fbeamfocus;  // This is needed to have correct drift distance for x, y distribution.

  G4double sigmax = fsigmax * sqrt(1.0 + pow(zdrift/fbetax, 2.0));
  G4double x0 = G4RandGauss::shoot(fx0, sigmax);
  G4double meandx = x0*zdrift / (zdrift*zdrift + fbetax*fbetax);
  G4double sigmadx = sqrt( femittancex*fbetax / (zdrift*zdrift + fbetax*fbetax) );
  G4double dx0 = G4RandGauss::shoot(meandx, sigmadx);

  G4double sigmay = fsigmay * sqrt(1.0 + pow(zdrift/fbetay, 2.0));
  G4double y0 = G4RandGauss::shoot(fy0, sigmay);
  G4double meandy = y0*zdrift / (zdrift*zdrift + fbetay*fbetay);
  G4double sigmady = sqrt( femittancey*fbetay / (zdrift*zdrift + fbetay*fbetay) );
  G4double dy0 = G4RandGauss::shoot(meandy, sigmady);

  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  G4double E = fParticleGun->GetParticleEnergy();

  G4double pz = sqrt( (E*E - mass*mass)/ (dx0*dx0 + dy0*dy0 + 1.0) );

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx0*pz, dy0*pz, pz));
//   fParticleGun->GeneratePrimaryVertex(anEvent);
  EmitInteraction(anEvent);
}



void LxPrimaryGenerator::GenerateMono(G4Event* anEvent)
{
//   G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
//   G4double E = fParticleGun->GetParticleEnergy();
//  G4double pz = sqrt(E*E - mass*mass);
//   G4double pz = sqrt(E*E + 2.0*E*mass);  // E is kinetic energy

//   fParticleGun->SetParticlePosition(G4ThreeVector(fx0, fy0, fz0));
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, pz));
//   std::cout << " ===== Generating mono, E/pdg/position: " <<
//   fParticleGun->GeneratePrimaryVertex(anEvent);
  EmitInteraction(anEvent);
}



void LxPrimaryGenerator::InitReader()
{
  if (fBeamType == beamMC) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    LuxeTestGenerator *gipr = new LuxeTestGenerator();
    if (fMCList) gipr->SetFileList(fMCfile);
    else gipr->AddEventFile(fMCfile);
    std::cout << "Processing file: " << fMCfile << std::endl;
//    gipr->SetFileType("out", 9, 9);
    gipr->SetFileType("out", 10, 10);
    lxgen = gipr;
    std::cout << "MC file type is set "  << std::endl;
    mcfilemutex.unlock();

  } else if (fBeamType == beamMCh5) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    lxgen = new LxHDF5Reader(fMCfile);
    mcfilemutex.unlock();

  } else if (fBeamType == beamMCTupleG4) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    lxgen = new LxNTupleReader(fMCfile);
    mcfilemutex.unlock();

  } else {
    G4String msgstr("Value of fBeamType: ");
    msgstr += std::to_string(fBeamType) + G4String(" is not supported\n");
    G4Exception("LxPrimaryGenerator::", "InitReader()", FatalException, msgstr.c_str());
  }

  if (fSkipEvents > 0) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    lxgen->SkipEvents(fSkipEvents);
    std::cout << "Skipping " << fSkipEvents << " events form the file " << fMCfile << std::endl;
    fSkipEvents = 0;   // This is to ensure that only one thread does it
    mcfilemutex.unlock();
  }
}



void LxPrimaryGenerator::GeneratefromMC(G4Event* anEvent)
{
  if (!lxgen) {
    InitReader();
  }

  if (fnfixparticles > 0) {
    --fnfixparticles;
//     fParticleGun->GeneratePrimaryVertex(anEvent);
    EmitInteraction(anEvent);
    // MC primary info is published by AttachPrimaryInfo() after this call; the
    // persistent fMCInfo / fMCPrime* still hold this split particle's values.
    return;
  }

  std::vector < std::vector <double> > ptcls;

  G4AutoLock mcfilemutex(&MCReadMutex);

  int nscat = lxgen->GetEventFromFile(ptcls);
  mcfilemutex.unlock();

  if (nscat > 1) {
      G4String msgstr("Error reading particle from a file! More than one were read, it is not supported!\n");
      G4Exception("LxPrimaryGenerator::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
  }
  if (nscat <= 0) {
    if (nscat == -1) {
      G4cout << "All particle from the file " << fMCfile << " are processed." << G4endl;
      fParticleGun->SetNumberOfParticles(0);
      ptcls.push_back(std::vector<double>(10, 90.0));
      G4RunManager::GetRunManager()->AbortRun(true);
    } else {
      G4String msgstr("Error reading MC particle from the file!\n");
      G4Exception("LxPrimaryGenerator::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
    }
  }

  std::vector <double> pdata = ptcls[0];
  int pid = static_cast<int>(pdata[7]);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(pid);
  if (!particle) {
    if (pid == 90)  {
      particle = particleTable->FindParticle("geantino");
      G4String msgstr("Initial particle with pdg_id=");
      msgstr += std::to_string(pid) + G4String(" is in the file! Ignore it\n");
      G4Exception("LxPrimaryGenerator::", "GeneratefromMC(Event)", JustWarning, msgstr.c_str());
    }
    else {
      G4String msgstr("Error setting initial particle from a file!\n");
      G4Exception("LxPrimaryGenerator::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
    }
  }

  if (pdata.size()>12) {
    fMCPrimeParentId = pdata[10];
    fMCPrimeNDescendant = pdata[11];
    fMCPrimeXi = pdata[12];
  }
  if (pdata.size() > 14)  {
    fMCPhysProc = pdata[14];
  }

  std::vector <G4double> pp(4);
  G4double  vtx[3];
  G4double  wghtf, wght = 1.0;

  wghtf = pdata[8];
  fnfixparticles = 1;
  if (wghtf >= fMCWeightScale)   fnfixparticles = fMCWeightScale;
  else if (wghtf > 1.0)  fnfixparticles = static_cast<G4int>(floor(wghtf));
  wght = wghtf/static_cast<G4double>(fnfixparticles);
  fMCInfo.SetWeight(wght);
  fMCInfo.SetMCTrackId(pdata[9]);

  std::copy(pdata.begin()+4, pdata.begin()+7, pp.begin());
  std::copy(pdata.begin()+1, pdata.begin()+4, vtx);
  std::for_each(vtx, vtx+3, [=](G4double &x){x *= CLHEP::um;});
  pp[3] = pdata[0];

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(pp[3]*CLHEP::GeV - particle->GetPDGMass());
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pp[0], pp[1], pp[2]));
  fParticleGun->SetParticlePosition(G4ThreeVector(vtx[0], vtx[1], vtx[2]));
  if (pdata.size() > 13)  fParticleGun->SetParticleTime(pdata[13]);

  if (    (!fselectpdg && (pid != 90))
       || (fselectpdg && std::find(fselectpdgvec.begin(), fselectpdgvec.end(), pid) != fselectpdgvec.end()) ) {
    if (true) {
//    if (TestHitTarget(pp, vtx) < 1.2) {
      fParticleGun->SetNumberOfParticles(1);
      fnfixparticles -= 1;
    } else {
      fParticleGun->SetNumberOfParticles(0);
      fnfixparticles = 0;
    }
  } else {
    fParticleGun->SetNumberOfParticles(0);
    fnfixparticles = 0;
  }

//   fParticleGun->GeneratePrimaryVertex(anEvent);
  EmitInteraction(anEvent);
}



G4double LxPrimaryGenerator::TestHitTarget(const std::vector <double> &pp, const double *vtx)
{
   G4ThreeVector pv = G4ThreeVector(pp[0], pp[1], pp[2]);
   G4ThreeVector rv = G4ThreeVector(vtx[0], vtx[1], vtx[2]);
   G4ThreeVector rt = rv + pv.unit() * (fGTargetZpos - vtx[2]);
   return std::fabs(2.0*rt.x()/fGTargetX);
}

}}

// DECLARE_GEANT4ACTION(LxPrimaryGenerator)
DECLARE_GEANT4ACTION_NS(dd4hep::sim, LxPrimaryGenerator)
