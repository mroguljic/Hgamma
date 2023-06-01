# Hgamma
Initial setup:
```
cmsrel CMSSW_12_3_0
cd CMSSW_12_3_0
cmsenv
cd ..

python3 -m virtualenv timber-env
echo -e "\ncd CMSSW_12_3_0\neval \`scram runtime -sh\`\ncd ..\n" >> timber-env/bin/activate
source timber-env/bin/activate
git clone -b Zbb_branch_py3 git@github.com:mroguljic/TIMBER.git
cd TIMBER
mkdir bin
cd bin
git clone git@github.com:fmtlib/fmt.git
cd ../..

#Boost library path (the boost version as well!) may change depending on scram arch and the CMSSW version so this may need to be modified by hand 
#Copy the whole multi-line string to activate script
cat <<EOT >> timber-env/bin/activate
nexport BOOSTPATH=/cvmfs/cms.cern.ch/slc7_amd64_gcc10/external/boost/1.75.0/lib
if grep -q '\${BOOSTPATH}' <<< '\${LD_LIBRARY_PATH}'
then
  echo 'BOOSTPATH already on LD_LIBRARY_PATH'
else
  export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${BOOSTPATH}
  echo 'BOOSTPATH added to PATH'
fi
EOT

source timber-env/bin/activate
cd TIMBER
source setup.sh
```

From then on:
```
source timber-env/bin/activate
```
