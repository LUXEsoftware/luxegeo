
def lxGenerator(dd4hepSimulation):
    from DDG4 import GeneratorAction, Kernel
    from dd4hep import units as u
    gen = GeneratorAction(Kernel(), "LxPrimaryGenerator/Primary", True)

    #gen.BeamType = "mc"
    #gen.MCFile   = "test_data4.out"

    gen.SpectraFile = "spectra_test_compt_electron.txt"
    #gen.SpectraFile = "spectra_bw_positron.txt"

    gen.BeamType = "mono"
    gen.Position = [0, 0, 0]
    #gen.Position = [0, 0, -7.4*u.m]
    gen.Energy = 16.5*u.GeV

    gen.Particle = "e-";

    #gen.BeamType = "mchdf5"
    #gen.MCFile = "e0gpc_2.0_0936_particles.h5"

    gen.enableUI()
    return gen

SIM.inputConfig.userInputPlugin = [lxGenerator]

  #declareProperty("BeamType",      fBeamTypeProp = "gaussian");
  #declareProperty("SpectraFile",   fSpectraFileProp = "");
  #declareProperty("MCFile",        fMCFileProp = "");
  #declareProperty("MCFileList",    fMCFileListProp = "");
  #declareProperty("SigmaX",        fSigmaXProp = 5.0 *dd4hep::um);
  #declareProperty("SigmaY",        fSigmaYProp = 5.0 *dd4hep::um);
  #declareProperty("SigmaZ",        fSigmaZProp = 24.0 *dd4hep::um);
  #declareProperty("PosZ",          fPosZProp = std::numeric_limits<double>::quiet_NaN());
  #declareProperty("Position",      fPositionProp);
  #declareProperty("MCWeightScale", fMCWeightScale = 1);
  #declareProperty("SelectPdg",     fSelectPdgProp);
  #declareProperty("SkipEvents",    fSkipEventsProp = 0);
  #declareProperty("Particle",      fParticleName = "e-");
  #declareProperty("Energy",        fEnergyProp = std::numeric_limits<double>::quiet_NaN());
  #declareProperty("Direction",     fDirectionProp);

