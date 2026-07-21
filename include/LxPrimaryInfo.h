//
/// \brief Definition of the LxPrimaryInfo class
//
/// Per-event Monte-Carlo primary information, carried as a DDG4 event extension.
///
/// Attached in EVERY beam mode via context()->event().addExtension<LxPrimaryInfo>(...).
/// DDG4 owns the object and deletes it at the end of the event -- do NOT delete it.
/// The output event action retrieves it with event().extension<LxPrimaryInfo>() and,
/// depending on fMode, writes the relevant fields to the ROOT tree.
///
/// This class also absorbs the role of the old EventInfo: fWeight and fMCTrackId with
/// the same GetWeight()/SetWeight()/GetMCTrackId()/SetMCTrackId() accessors. In this
/// generator the primary particle fully defines the initial state of the event, so no
/// separate event-level object is needed.
///
/// Kinematics are stored in Geant4/CLHEP units, exactly as held by the G4ParticleGun.

#ifndef LxPrimaryInfo_h
#define LxPrimaryInfo_h 1

/// Beam / generation mode. Kept as a global unscoped enum with the original
/// enumerator names so both the generator and the output action use the bare names.
enum tBeamType : int { beamGauss, beamMono, beamMC, beamMCh5, beamMCTupleG4, beamUnknown = -1 };


class LxPrimaryInfo
{
public:
  LxPrimaryInfo() = default;

  // ---- former EventInfo interface -------------------------------------------
  double GetWeight() const     { return fWeight; }
  void   SetWeight(double w)   { fWeight = w; }
  long   GetMCTrackId() const  { return fMCTrackId; }
  void   SetMCTrackId(long id) { fMCTrackId = id; }

  // ---- data (plain per-event record) ----------------------------------------
  tBeamType fMode         { beamUnknown };
  double    fWeight       { 1.0 };
  long      fMCTrackId    { -1 };

  // MC-file provenance (meaningful for MC* modes; defaults otherwise)
  long      fParentId     { -1 };
  int       fNDescendant  { -1 };
  double    fXi           { -1.0 };
  int       fPhysProc     { 0 };

  // primary kinematics as launched (Geant4/CLHEP units)
  int       fPdg          { 0 };
  double    fEkin         { 0.0 };
  double    fTime         { 0.0 };
  double    fVertex[3]    { 0.0, 0.0, 0.0 };
  double    fDirection[3] { 0.0, 0.0, 0.0 };

  // gun multiplicity for this event: 0 means the primary was filtered out
  // (the old code skipped SavePrimaryTrack() in that case).
  int       fNGenerated   { 0 };
};

#endif
