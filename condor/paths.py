#Location of the main Hgamma directory
HGAMMA_DIR          = "/users/mrogul/Work/Hgamma/Hgamma"
#Where to store selection trees
SELECTION_DIR       = HGAMMA_DIR+"/results/selection/"
#Where to store skimmed NanoAODs
SKIM_DIR            = "/STORE/matej/Hgamma_skims/v9/"
#Where to store condor logs for selection jobs
SELECTION_JOB_DIR   = HGAMMA_DIR+"/condor/selection_jobs"
#Where to store condor logs for skimming jobs
SKIM_JOB_DIR        = HGAMMA_DIR+"/condor/skim_jobs"
#Location of the directory where timber-env in located (not the location of timber-env itself!)
TIMBERENV_DIR       = "/users/mrogul/Work/Hgamma/"
#Location of the CMSSW
CMSSW_DIR           = "/users/mrogul/Work/Hgamma/CMSSW_12_3_0/"
#Where to store condor logs for template jobs
TEMPLATE_JOB_DIR    = HGAMMA_DIR+"/condor/template_jobs/"

#Where CR selection trees are stored
SELECTION_CR_DIR    = HGAMMA_DIR+"/results/selection_CR/"
#Where to store condor logs for selection CR jobs
SELECTION_CR_JOB_DIR   = HGAMMA_DIR+"/condor/selection_CR_jobs"
#Where to store condor logs for template CR jobs
TEMPLATE_CR_JOB_DIR = HGAMMA_DIR+"/condor/template_CR_jobs/"
#Where to store skimmed CR NanoAODs
SKIM_CR_DIR         = "/STORE/matej/Zbb_skims/v9_py3_with_photon/"