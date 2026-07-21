
#ifndef __LUXETESTGENERATOR_H
#define __LUXETESTGENERATOR_H

#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "H5Cpp.h"

class FileData;
class LxHDF5Reader;
class G4String;
class G4Tubs;


class LuxeMCGenerator
{
public:
  LuxeMCGenerator() {};
  virtual ~LuxeMCGenerator() {};
  virtual int GetEventFromFile (std::vector < std::vector <double> > &ptcls) = 0;
  virtual int SkipEvents(const size_t nevskip) = 0;
};



class LuxeTestGenerator : public LuxeMCGenerator
{
public:
  LuxeTestGenerator ();
  virtual ~LuxeTestGenerator ();

  void AddEventFile (const std::string fname);
  void SetFileType (const std::string ftype, const int n_col = 17, const int n_skip = 0);
  void SetFileList(const std::vector<std::string> fnlist);
  void SetFileList(const std::string fnamelist);
  virtual int GetEventFromFile (std::vector < std::vector <double> > &ptcls);
  virtual int SkipEvents(const size_t nevskip);

protected:
  int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);

protected:

  FileData     *datf;
};



class FileData
{
public:
  FileData();
  FileData(const std::string fname);
  FileData(const std::vector<std::string> &fnamelist);
  virtual ~FileData();
  void AddFile(const std::string fname);
  void SetNColumns (const int n) { n_columns = n; }
  void SetSkipLines (const int n) { skip_lines = n; }
  void SetDebug (const int n) { debugl = n; }

  int GetNColumns () const { return n_columns; }
  int GetData(std::vector < std::vector <double> > &particles);
  void ListFiles(void) const;
  std::string GetCurrentFileName(void) const { return fcname; }
  int SetFileType(std::string ftype);

protected:
  int ReadRecord(std::vector < std::vector <double> > &particles);
  int ReadRecordOut(std::vector < std::vector <double> > &particles);

protected:
  std::vector<std::string> fname_list;
  std::string              fcname;
  std::fstream             fdata;
  int                      fid;
  int                      n_columns, skip_lines;
  int                     (FileData::*freadfn)(std::vector < std::vector <double> > &particles);
  int                      debugl;
};



class LxHDF5Reader : public LuxeMCGenerator
{
public:
  LxHDF5Reader(const std::string fnamemc);
  virtual ~LxHDF5Reader();
  int DumpData(const hsize_t start, const hsize_t count);
  virtual int GetEventFromFile (std::vector < std::vector <double> > &ptcls);
  virtual int SkipEvents(const size_t nevskip);

protected:
  void OpenFile();
  void OpenGroupDatasets (const int pdg);
  bool NextGroup();

protected:

  std::string   fFileName;
  H5::H5File   *fFile;
  
  std::map<int, std::string> fParticleGroup;
  std::map<int, hsize_t> fParticleDataSize;
  std::map<std::string, H5::DataSet> fDsMap;
  std::map<std::string, H5::DataSpace> fDpMap;
  std::map<int, std::vector<std::string> > FDsNames;
  int fCurrentGroup;
  hsize_t fCurrentOffset;
};



class LxNTupleReader : public LuxeMCGenerator
{
public:
  LxNTupleReader(const std::string fnamemc);
  virtual ~LxNTupleReader();
  virtual int GetEventFromFile (std::vector < std::vector <double> > &ptcls);
  virtual int SkipEvents(const size_t nevskip);

protected:
  void OpenFile();
  G4int FindNextParticle();
  void FillOutput(std::vector < std::vector <double> > &ptcls, const G4int vppos);
//   G4int Transform(G4ThreeVector &rr, G4ThreeVector &pp);

protected:

  G4String          fFileName;
  G4int             fTreeID;
  G4String          fTreeNme;
  G4int             fReadEvent;
  G4int             fEPCounter;

  std::map <G4int, std::vector<G4double> > fvolTrackDMap;
  std::map <G4int, std::vector<G4int> >    fvolTrackIMap;
  G4double                                 fWeight;
  G4int                                    fEventID;

  std::vector<G4int>  fDetIDv;
  G4bool              fDoTransformation;
  G4Tubs             *fDumpSolid;

};

#endif

