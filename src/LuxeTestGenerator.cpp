

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <map>

#include "H5Cpp.h"

#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RootAnalysisReader.hh"
using G4AnalysisReader = G4RootAnalysisReader;

// #include "LXSetUp.hh"
#include "LuxeTestGenerator.h"


LuxeTestGenerator::LuxeTestGenerator()
  : datf(0)  {

}

LuxeTestGenerator::~LuxeTestGenerator() {

  if (datf) delete datf;
}



void LuxeTestGenerator::AddEventFile (const std::string fname)
{
  if (!datf) {
    datf = new FileData();  
  }
  datf->AddFile(fname);
}


void LuxeTestGenerator::SetFileType (const std::string ftype, const int n_col, const int n_skip)
{
  if (!datf) {
    datf = new FileData();
  }
  datf->SetFileType(ftype);
  datf->SetNColumns(n_col);
  datf->SetSkipLines(n_skip);
}


int LuxeTestGenerator::GetEventFromFile (std::vector < std::vector <double> > &ptcls)
{
  if (!datf) return -5;
  int nrdat = datf->GetData(ptcls);
  return nrdat;
}



int LuxeTestGenerator::SkipEvents(const size_t nevskip)
{
// poor implementation, never needed
  if (!datf) return -5;
  std::vector < std::vector <double> > ptcls;
  size_t nev(nevskip);
  int nrdat(1);
  while (nev-- && nrdat > 0) nrdat = datf->GetData(ptcls);
  return nrdat;
}



void LuxeTestGenerator::SetFileList(const std::vector<std::string> fnlist)
{
  for (std::vector<std::string>::const_iterator itr = fnlist.begin(); itr != fnlist.end(); ++itr) {
    AddEventFile(*itr);  
  }
}



void LuxeTestGenerator::SetFileList(const std::string fnamelist)
{
  std::vector<std::string>  flist;
  ProcessList(fnamelist, flist);
  SetFileList(flist);
}



int LuxeTestGenerator::ProcessList(const std::string &fnamelist, std::vector<std::string> &flist)
{
  std::fstream  fdata;
  fdata.open(fnamelist, std::ios::in);
  if (!fdata.is_open()) {
    throw std::runtime_error(std::string("Error reding data from the file ") + fnamelist);
  }
  
  unsigned long lid = 0;
  while (!fdata.eof()) {
    std::string  ffname;
    fdata >> ffname;
    if (!fdata.fail()) { 
//       std::cout << "File name " << ffname << " is read from the list file" << std::endl;
      flist.push_back(ffname);
    }
    else if (fdata.eof()) { break; }
    else {
      std::cout << "ProcessList(..)  :  Error reading data from the file " << fnamelist 
                << ",  line: " << lid << ". Exit." << std::endl;
      fdata.close();          
      return -2;
    }
    ++lid;
  }
  
  fdata.close();

  return 0;
}  
  


//////////////////////////////////////////////////////////////////////////

FileData::FileData() : fcname(""), fid(-1), n_columns(17), skip_lines(0), freadfn(&FileData::ReadRecord), debugl(0)
{
}


FileData::FileData(const std::string fname) : fcname(""), fid(-1), n_columns(17), skip_lines(0), 
   freadfn(&FileData::ReadRecord), debugl(0)
{
  fname_list.push_back(fname);
}


FileData::FileData(const std::vector<std::string> &fnamelist) :
  fname_list(fnamelist), fcname(""), fid(-1), n_columns(17), skip_lines(0), 
  freadfn(&FileData::ReadRecord), debugl(0)
{
}


FileData::~FileData()
{
  if (fdata.is_open()) fdata.close();
}


int FileData::SetFileType(std::string ftype)
{
  if(!ftype.compare("hepevt")) {
    freadfn = &FileData::ReadRecord;
    n_columns = 18;
    skip_lines = 0;
  } else if (!ftype.compare("out")) { 
      freadfn = &FileData::ReadRecordOut;
      n_columns = 9;
      skip_lines = 9;
  } else {
    std::cout << "FileData::SetFileType: Error: File type " << ftype << " is not supported!\n";    
    return -1;
  }
  return 0;
}


void FileData::AddFile(const std::string fname)
{
  if (!fdata.is_open()) {
    fname_list.push_back(fname);
  } else {
    std::cout << "FileData::AddFile :  Cannot add file while reading data\n";
  }
}


int FileData::GetData(std::vector < std::vector <double> > &particles)
{
// It would be better probably to add exceptions for error handling not just return values
  int  ndata(0);  
  if (fdata.is_open() ) {
    ndata = (this->*freadfn)(particles);
    if (ndata >= 0 ) {
      return ndata;
    } else {  
      if ( fdata.eof() ) { fdata.close();  }
      else  { fdata.close(); return ndata; }
    }
  }

  if (!fdata.is_open()) {
    if ( ++fid < 0 || fname_list.size() > static_cast<unsigned int>(fid) ) {
      fdata.open(fname_list[fid], std::ios::in);
      fcname = fname_list[fid];
      if (debugl > 0) { std::cout << "FileData::GetData: opening file: " << fcname << std::endl; }
    } else {
      std::cout << "FileData::GetData(..) :  All data from " << fid << " file(s) have been read\n";
      fcname = "";
      --fid;
      return -1;
    }
  }
    
  if (!fdata.is_open()) {
    std::cerr << "FileData::GetData(..) :  Error open file " << fname_list[fid] << ". Exit.\n";
    fcname = "";
    return -2;
  }

  fdata.clear();
  std::string tmpstr;
  for (int ii = 0; ii < skip_lines; ++ii) {std::getline(fdata, tmpstr); }
  
  ndata = (this->*freadfn)(particles);
  if (ndata < 0) fdata.close();
      
  return ndata;
}


int FileData::ReadRecord(std::vector < std::vector <double> > &particles)
{
  std::vector <double> vdat(n_columns);
  int evid(0), npart(0);
   
  particles.clear(); 
  fdata >> evid >> npart;
  if (fdata.fail() || fdata.eof()) { 
    return -1; 
  }
  for (int ii = 0; ii < npart; ++ii) {
    for (int jj = 0; jj < n_columns; ++jj) {
      double xx;  
      fdata >> xx;
      if (!fdata.fail()) { vdat[jj] = xx; }
      else {
        std::cout << "FileData::ReadRecord(..)  :  Error reading data from the file " << fname_list[fid] 
                  << ",  event: " << evid << "  particle: " << ii << "  column: " << jj << std::endl;
        particles.clear();
        return -2;
      }
    }
    particles.push_back(vdat);
  }
  if (debugl > 1) {
    std::cout << "FileData::ReadRecord : Info: file: " << fcname << "  event: " << evid 
              << "  number of particles: " << npart << std::endl;
    for (auto itr = particles.begin(); itr != particles.end(); ++itr) {
      std::for_each(itr->begin(), itr->end(), [](const double x) {std::cout << x << "  ";} );
      std::cout << std::endl;
    }
  }
  return npart;
}


int FileData::ReadRecordOut(std::vector < std::vector <double> > &particles)
{
  std::vector <double> vdat(n_columns);
  int evid(0), npart(0);
   
  particles.clear(); 
  for (int jj = 0; jj < n_columns; ++jj) {
    double xx;  
    fdata >> xx;
//     std::cout << xx << "  ";
    if (!fdata.fail()) { vdat[jj] = xx; }
    else if (fdata.eof()) { return -1; }
    else {
      std::cout << "FileData::ReadRecordOut(..)  :  Error reading data from the file " << fname_list[fid] 
                << "  column: " << jj << std::endl;
      particles.clear();
      return -2;
    }
  }
  particles.push_back(vdat);
// std::cout << std::endl;

  if (debugl > 1) {
    std::cout << "FileData::ReadRecordOut : Info: file: " << fcname << "  event: " << evid << std::endl;
    for (auto itr = particles.begin(); itr != particles.end(); ++itr) {
      std::for_each(itr->begin(), itr->end(), [](const double x) {std::cout << x << "  ";} );
      std::cout << std::endl;
    }
  }
  return ++npart;
}


void FileData::ListFiles(void) const
{
  std::for_each(fname_list.begin(), fname_list.end(), [](const std::string ss) {std::cout << ss << std::endl;} );
}



//////////////////////////////////////////////////////////////////////////

LxHDF5Reader::LxHDF5Reader(const std::string fnamemc): 
fFileName(fnamemc), fFile(0), fCurrentGroup(0), fCurrentOffset(-1)
{
  OpenFile();
  fParticleGroup[11] = "/final-state/electron";
  fParticleGroup[-11] = "/final-state/positron";
  fParticleGroup[22] = "/final-state/photon";
  
  fParticleDataSize[11] = -1;
  fParticleDataSize[-11] = -1;
  fParticleDataSize[22] = -1;

  std::vector<std::string> dsnamesel{"momentum", "position", "weight", "id", "parent_id", "n_gamma"};
  std::vector<std::string> dsnamesph{"momentum", "position", "weight", "id", "parent_id", "n_pos", "a0_at_creation"};
  std::vector<std::string> dsnamespo{"momentum", "position", "weight", "id", "parent_id", "n_gamma"};
  FDsNames[11] = dsnamesel;
  FDsNames[-11] = dsnamespo;
  FDsNames[22] = dsnamesph;

  fCurrentGroup = 11;
  OpenGroupDatasets(fCurrentGroup);
}



LxHDF5Reader::~LxHDF5Reader()
{
  if (fFile) {
    fFile->close();
    delete fFile;
  }
}



void LxHDF5Reader::OpenFile()
{
  try {
    fFile = new H5::H5File(fFileName.c_str(), H5F_ACC_RDONLY);

  } catch( H5::FileIException &error ) {
    fFile = 0;
    error.printErrorStack();
  }
}



int LxHDF5Reader::SkipEvents(const size_t nevskip)
{
  if (nevskip < fParticleDataSize[fCurrentGroup]-fCurrentOffset) {
   fCurrentOffset += nevskip;
   return nevskip + 1;
  } else {
    size_t ngskip = nevskip - (fParticleDataSize[fCurrentGroup] - fCurrentOffset);
    if (NextGroup()) {
      OpenGroupDatasets(fCurrentGroup);
      return SkipEvents(ngskip);
    }
  }
  fCurrentOffset = fParticleDataSize[fCurrentGroup];
  return -1;
}



bool LxHDF5Reader::NextGroup()
{
  if (fCurrentGroup == 11) {
    fCurrentGroup = 22;
  } else if (fCurrentGroup == 22) {
    fCurrentGroup = -11;
  } else if (fCurrentGroup == -11) {
    return false;
  }
  return true;
}


int LxHDF5Reader::GetEventFromFile (std::vector < std::vector <double> > &ptcls)
{
  ptcls.clear();

  while (fCurrentOffset >= fParticleDataSize[fCurrentGroup]) {
    if (NextGroup())
      OpenGroupDatasets(fCurrentGroup);
    else
      return -1;
  }

  const std::vector<std::string>  &grdsnames = FDsNames[fCurrentGroup];

  std::map<std::string, double*> pdata;
  pdata[grdsnames[0]] = new double[4];
  pdata[grdsnames[1]] = new double[4];
  std::for_each(grdsnames.begin()+2, grdsnames.end(), [&pdata](const std::string ss) {pdata[ss] = new double;} );

  try {
    for (const auto &dsname: grdsnames) {
      H5::DataSpace &dspace = fDpMap[dsname];
      int rank = dspace.getSimpleExtentNdims();

//       std::cout << dsname << "  rank: " << rank << std::endl;
//       if (rank != 1) {
//       std::cout << "rank != 1 is not supported\n";
//       }

      hsize_t count = 1;
      hsize_t memdim[1];
      memdim[0] = count;
      H5::DataSpace memspace(rank, memdim);

      hsize_t  pcount[1];    // size of the hyperslab in the file
      hsize_t  poffset[1];          // hyperslab offset in the file
      poffset[0] = fCurrentOffset;
      pcount[0]  = count;
      dspace.selectHyperslab(H5S_SELECT_SET, pcount, poffset);

      H5T_class_t tpclass = fDsMap[dsname].getTypeClass();
      if (tpclass == H5T_INTEGER) {
        int intbuf;
        fDsMap[dsname].read(&intbuf, H5::PredType::NATIVE_INT, memspace, dspace);
        *pdata[dsname] = intbuf;
      } else if (tpclass == H5T_FLOAT) {
        fDsMap[dsname].read(pdata[dsname], H5::PredType::NATIVE_DOUBLE, memspace, dspace);
      } else if (tpclass == H5T_ARRAY) {
        hsize_t dims[1] = {4};
        fDsMap[dsname].read(pdata[dsname], H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, dims), memspace, dspace);
      }

      memspace.close();
    }

    ++fCurrentOffset;
  }
  catch (H5::DataSetIException &error) {
    error.printErrorStack();
  }
  catch( H5::DataSpaceIException &error ) {
    error.printErrorStack();
  }

  std::vector<double> vdat(13);
  vdat[0] = pdata[grdsnames[0]][0];
  std::copy(pdata[grdsnames[0]]+1, pdata[grdsnames[0]]+4, vdat.begin()+4);

//   std::for_each(pdata[FDsNames[1]], pdata[FDsNames[1]]+4, [=](double &x){x *= 1.0e6;}); // m to um
  std::for_each(pdata[grdsnames[1]], pdata[grdsnames[1]]+4, [=](double &x){x *= 1.0e3;}); // m to um
  std::copy(pdata[grdsnames[1]]+1, pdata[grdsnames[1]]+4, vdat.begin()+1);

  vdat[7] = fCurrentGroup;            // pdg
  vdat[8] = pdata[grdsnames[2]][0];   // weight
  vdat[9] = pdata[grdsnames[3]][0];   // MC particle id
  vdat[10] = pdata[grdsnames[4]][0];  // MC parent_particle id
  vdat[11] = pdata[grdsnames[5]][0];  // MC n_secondaries
  vdat[12] = grdsnames.size()>6 ? pdata[grdsnames[6]][0] : 0.0;  // a0_at_creation (xi)

  ptcls.push_back(vdat);

  std::for_each (pdata.begin(), pdata.end(), [](std::pair<std::string, double*> pp) {delete [] pp.second; });

  return 1;
}


void LxHDF5Reader::OpenGroupDatasets (const int pdg)
{
  try {
    for (auto dsitr = fDsMap.begin(); dsitr != fDsMap.end(); ++dsitr) {
      if (dsitr->second.getId() > 0) {
        dsitr->second.close();
        fDpMap[dsitr->first].close();
      }
    }
    fDsMap.clear();
    fDpMap.clear();

    const std::vector<std::string>  &grdsnames = FDsNames[pdg];
    for (const auto &dsname: grdsnames) {

      std::string dspath = fParticleGroup[pdg] + "/" + dsname;
      fDsMap[dsname] = fFile->openDataSet(dspath);
      H5D_space_status_t status;
      fDsMap[dsname].getSpaceStatus(status);
      if (status == H5D_SPACE_STATUS_NOT_ALLOCATED) {
        fParticleDataSize[pdg] = 0;
        break;
      }	
      fDpMap[dsname] = fDsMap[dsname].getSpace();

      if (dsname == grdsnames[0]) {
        hsize_t npcount[1];
        fDpMap[dsname].getSimpleExtentDims(npcount, NULL);
        fParticleDataSize[pdg] = npcount[0];
      }
    }
    fCurrentOffset = 0;
  }
  catch (H5::DataSetIException &error) {
    error.printErrorStack();
  }
  catch( H5::DataSpaceIException &error ) {
    error.printErrorStack();
  }
}



int LxHDF5Reader::DumpData(const hsize_t start, const hsize_t count)
{
  if (!fFile) {
    OpenFile();
  }

  std::string epdsetname = "/final-state/electron/momentum";
  try {
//     H5::Exception::dontPrint();

    H5::DataSet epdset = fFile->openDataSet(epdsetname);
    H5::DataSpace epdspace = epdset.getSpace();
    int rank = epdspace.getSimpleExtentNdims();
    std::cout << epdsetname << "  rank: " << rank << std::endl;
    if (rank != 1) {
      std::cout << "rank != 1 is not supported\n";
    }

    hsize_t nposcount[1];
//     int ndims = epdspace.getSimpleExtentDims(nposcount, NULL);
    epdspace.getSimpleExtentDims(nposcount, NULL);
    size_t npos = nposcount[0];
    std::cout << "Number of positrons: " << npos << std::endl;

    H5::DataType epdt =	epdset.getDataType();
    std::cout << "DataType: " << epdt.getClass() << std::endl;

    H5::ArrayType epdart = epdset.getArrayType();
    std::cout << "AttayType: " << epdart.getClass() << std::endl;

    if (start > npos) {
      std::cout << "Start entry " << start << " is greater than number of positrons." << std::endl;
      return 0;
    }

    hsize_t memdim[1];
    memdim[0] = count;
    H5::DataSpace memspace(rank, memdim);

    hsize_t  pcount[1];    // size of the hyperslab in the file
    hsize_t  poffset[1];          // hyperslab offset in the file
    double *x = new double[4*count];
    poffset[0] = start;
    pcount[0]  = count;
    epdspace.selectHyperslab(H5S_SELECT_SET, pcount, poffset);

    epdset.read(x, epdt, memspace, epdspace);

    std::cout << "Showing read data:\n";
    for (size_t ii = 0; ii < count; ++ii) {
      std::cout << std::setprecision(16) << x[ii*4] << "   "
                << std::setprecision(16) << x[ii*4 + 1] << "   "
                << std::setprecision(16) << x[ii*4 + 2] << "   "
                << std::setprecision(16) << x[ii*4 + 3] << std::endl;
    }

    memspace.close();
    epdspace.close();
    epdset.close();

    delete [] x;

  }
  catch( H5::FileIException &error ) {
    error.printErrorStack();
    return -1;
  }
  catch (H5::DataSetIException &error) {
    error.printErrorStack();
    return -1;
  }
  catch( H5::DataSpaceIException &error ) {
    error.printErrorStack();
    return -1;
  }
  catch( H5::DataTypeIException &error ) {
    error.printErrorStack();
    return -1;
  }

  return 0;
}


//////////////////////////////////////////////////////////////////////////

LxNTupleReader::LxNTupleReader(const std::string fname) :
  fFileName(fname), fTreeID(-1), fTreeNme(), fReadEvent(-1), fEPCounter(-1), fWeight(1.0), fEventID(-1.0),
  fDetIDv(0), fDoTransformation(false), fDumpSolid(0)
{
  OpenFile();
}



LxNTupleReader::~LxNTupleReader()
{}



void LxNTupleReader::OpenFile()
{
  fTreeNme = "Tracks";
  auto analysisReader = G4AnalysisReader::Instance();
  fTreeID = analysisReader->GetNtuple(fTreeNme, fFileName);

  if (fTreeID < 0) {
    G4String msgstr("Error reading NTuple ");
    msgstr += fTreeNme + " from the file " + fFileName + "\n";
    G4Exception("LxNTupleReader::", "OpenFile()", FatalException, msgstr.c_str());
  }

  analysisReader->SetNtupleIColumn(fTreeID, "eventid",  fEventID);                    //0
  analysisReader->SetNtupleIColumn(fTreeID, "trackid",  fvolTrackIMap[16]);           //1
  analysisReader->SetNtupleIColumn(fTreeID, "detid",    fvolTrackIMap[17]);           //2
  analysisReader->SetNtupleIColumn(fTreeID, "pdg",      fvolTrackIMap[18]);           //3
  analysisReader->SetNtupleIColumn(fTreeID, "physproc", fvolTrackIMap[19]);           //4
  analysisReader->SetNtupleDColumn(fTreeID, "E", fvolTrackDMap[0]);                   //5
  analysisReader->SetNtupleDColumn(fTreeID, "x", fvolTrackDMap[1]);                   //6
  analysisReader->SetNtupleDColumn(fTreeID, "y", fvolTrackDMap[2]);                   //7
  analysisReader->SetNtupleDColumn(fTreeID, "z", fvolTrackDMap[3]);                   //8
  analysisReader->SetNtupleDColumn(fTreeID, "t", fvolTrackDMap[4]);                   //9
  analysisReader->SetNtupleDColumn(fTreeID, "vtxx", fvolTrackDMap[5]);                //10
  analysisReader->SetNtupleDColumn(fTreeID, "vtxy", fvolTrackDMap[6]);                //11
  analysisReader->SetNtupleDColumn(fTreeID, "vtxz", fvolTrackDMap[7]);                //12
  analysisReader->SetNtupleDColumn(fTreeID, "px",   fvolTrackDMap[8]);                //13
  analysisReader->SetNtupleDColumn(fTreeID, "py",   fvolTrackDMap[9]);                //14
  analysisReader->SetNtupleDColumn(fTreeID, "pz",   fvolTrackDMap[10]);               //15
  analysisReader->SetNtupleDColumn(fTreeID, "theta",  fvolTrackDMap[11]);             //16
  analysisReader->SetNtupleDColumn(fTreeID, "phi",    fvolTrackDMap[12]);             //17
  analysisReader->SetNtupleDColumn(fTreeID, "xlocal", fvolTrackDMap[13]);             //18
  analysisReader->SetNtupleDColumn(fTreeID, "ylocal", fvolTrackDMap[14]);             //19
  analysisReader->SetNtupleDColumn(fTreeID, "zlocal", fvolTrackDMap[15]);             //20
  analysisReader->SetNtupleDColumn(fTreeID, "weight", fWeight);                       //21
  analysisReader->SetNtupleIColumn(fTreeID, "ptrackid", fvolTrackIMap[20]);           //22
  analysisReader->SetNtupleIColumn(fTreeID, "nsecondary", fvolTrackIMap[21]);         //23
  analysisReader->SetNtupleDColumn(fTreeID, "esecondary", fvolTrackDMap[22]);         //24
//   analysisReader->SetNtupleIColumn(fTreeID, "runid");                                 //25
}



int LxNTupleReader::GetEventFromFile (std::vector < std::vector <double> > &ptcls)
{
  G4int vpos = FindNextParticle();
  if (vpos >= 0) {
    FillOutput(ptcls, vpos);
    return 1;
  } else {
    auto analysisReader = G4AnalysisReader::Instance();
    while ( analysisReader->GetNtupleRow(fTreeID) ) {
      ++fReadEvent;
      fEPCounter = -1;
      G4int vposi = FindNextParticle();
      if (vposi >= 0) {
        FillOutput(ptcls, vposi);
        return 1;
      }
    }
  }
  return -1;
}



int LxNTupleReader::SkipEvents(const size_t nevskip)
{
  size_t nev(nevskip);
  G4int vpos(0);

  while (nev-- && vpos >= 0) {
    vpos = FindNextParticle();
    if (vpos < 0) {
      auto analysisReader = G4AnalysisReader::Instance();
      while ( analysisReader->GetNtupleRow(fTreeID) ) {
        ++fReadEvent;
        fEPCounter = -1;
        vpos = FindNextParticle();
        if (vpos >= 0) break;
      }
    }
  }
  return nevskip - nev;
}



G4int LxNTupleReader::FindNextParticle()
{
  std::vector<G4int> &detidv = fvolTrackIMap[17];
//   std::vector<G4int> &trkidv = fvolTrackIMap[16];
  if (fDetIDv.empty()) {
    if ( ++fEPCounter < detidv.size() ) {
      return fEPCounter;
    } else {
      return -1;
    }
  } else {
    G4String msgstr("Error generating primary particle from NTuple, selection is not supported!\n");
    G4Exception("LxNTupleReader::", "FindNextParticle()", FatalException, msgstr.c_str());
  }
  return -1;
}



void LxNTupleReader::FillOutput(std::vector < std::vector <double> > &ptcls, const G4int vppos)
{
  ptcls.clear();

  G4int pdg = fvolTrackIMap[18][vppos];
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(pdg);
  if (!particle) {
    G4String msgstr("Error generating primary particle from NTuple, particle PDG not found in the table!\n");
    G4Exception("LxNTupleReader::", "FillOutput()", FatalException, msgstr.c_str());
  }
  G4double pmass = particle->GetPDGMass();

  G4ThreeVector rr(fvolTrackDMap[1][vppos], fvolTrackDMap[2][vppos], fvolTrackDMap[3][vppos]);
  G4ThreeVector pp(fvolTrackDMap[8][vppos], fvolTrackDMap[9][vppos], fvolTrackDMap[10][vppos]);
  rr *= mm;
//   if (fDoTransformation) Transform(rr, pp);
  rr *= 1.0/um;

  std::vector<double> vdat(15);
  vdat[0] = (fvolTrackDMap[0][vppos]*GeV + pmass) /GeV;   // E
  vdat[1] = rr.x();               // x
  vdat[2] = rr.y();               // y
  vdat[3] = rr.z();               // z

  vdat[4] = pp.x();     // px
  vdat[5] = pp.y();     // py
  vdat[6] = pp.z();     // pz

  vdat[7] = pdg;                           // pdg
  vdat[8] = fWeight;                       // weight
  vdat[9] = fvolTrackIMap[16][vppos];      // MC particle id
  vdat[10] = fEventID;                     // MC parent_particle id
  vdat[11] = fvolTrackIMap[21][vppos];     // MC n_secondaries
  vdat[12] = 0.0;                          // a0_at_creation (xi)
  vdat[13] = fvolTrackDMap[4][vppos];      // t, time
  vdat[14] = fvolTrackIMap[19][vppos];     // Phys process

  ptcls.push_back(vdat);
}


/*
G4int LxNTupleReader::Transform(G4ThreeVector &rr, G4ThreeVector &pp)
{
  G4ThreeVector dir(-1.0*pp.unit());
  G4ThreeVector ppos(rr);
  LXSetUp *lxs = LXSetUp::Instance();
  G4double rmax = lxs->HICSDumpR;
  G4double dz = lxs->HICSDumpZ;
  G4double distToDump = 100.0 *mm;

  if (!fDumpSolid) fDumpSolid = new G4Tubs ("AuxTransforSolid", 0.0, rmax, dz/2.0, 0.0, 2.0*M_PI);

  // move position 10cm in front of the beam dump (solid cylinder in its ref. frame), this is where the particles are generated
  ppos.setZ(-dz/2.0 - distToDump);

//   G4cout << "ppos: " << ppos.x() << "  " << ppos.y() << "  " << ppos.z()
//          << "  dir: " << dir.x() << "  " << dir.y() << "  " << dir.z() << G4endl;
  G4double ld = fDumpSolid->DistanceToIn(ppos, dir);
//   G4cout << "Distance: " << ld << G4endl;
//   {G4ThreeVector pspos = ppos + dir * ld;
//   G4cout << "pspos: " << pspos.x() << "  " << pspos.y() << "  " << pspos.z() << G4endl;}

  G4ThreeVector pspos(ppos);
  if (ld != kInfinity) {
    pspos += dir * ld;
  } else {
    G4String msgstr("Transformation of primary particle failed! Using initial position.\n");
    G4Exception("LxNTupleReader::", "Transform()", JustWarning, msgstr.c_str());
  }

//   G4cout << "pspos: " << pspos.x() << "  " << pspos.y() << "  " << pspos.z() << G4endl;
  G4double alf = atan2(lxs->HICSDumpFrontXPos, lxs->HICSDumpFrontZPos - lxs->IPMagnetZpos);
  pspos.setZ(pspos.z() + dz/2.0 + lxs->HICSDumpFrontZPos/cos(alf));
  rr = pspos.rotateY(-alf);
  pp.rotateY(-alf);

  return 0;
}*/
