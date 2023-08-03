#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/m/mrogulji/Hgamma/CMSSW_12_3_0/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/m/mrogulji/Hgamma/
source timber-env/bin/activate

cd /afs/cern.ch/work/m/mrogulji/Hgamma/Hgamma/triggerStudies
echo "$(pwd)"

echo "triggerStudy_CR_alt.py $*"
python triggerStudy_CR_alt.py $*