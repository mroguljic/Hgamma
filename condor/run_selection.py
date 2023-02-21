#!/usr/bin/env python

import os, sys, re
from templates import *
from run_skim import createDirIfNotExist
import ROOT as r
from paths import SELECTION_DIR, SELECTION_JOB_DIR, SKIM_DIR

def split_jobs(files, njobs):
    for i in range(0, len(files), njobs):
        yield files[i:i + njobs]

def removeEmptyFiles(inputPath):
    with open(inputPath, "r") as reader:
        nonEmptyInput   = []
        for line in reader:
            filePath    = line.strip()
            rFile       = r.TFile.Open(filePath)
            ttree       = rFile.Get("Events")
            try:
                treeSize    = ttree.GetEntriesFast()
            except:
                print("Couldn't get ttree in file {0}".format(filePath))
                treeSize = 0
            if(treeSize!=0):
                nonEmptyInput.append(line)
    with open(inputPath, "w+") as writer:
        writer.writelines(nonEmptyInput)

def checkForCongregateResult(outDir,sample):
    nomFile = os.path.join(outDir+"_nom.root")
    if not os.path.exists(nomFile):
        return False
    return True

def checkIfAlreadyProcessed(fileBase,outDir,sample):
    if("ZGamma" in sample or "WGamma" in sample or "Hgamma" in sample):
        varFlag = True
    else:
        varFlag = False

    if varFlag:
        variations  = ["nom","jerUp","jerDown","jesUp","jesDown","jmsUp","jmsDown","jmrUp","jmrDown"]
    else:
        variations  = ["nom"]

    for variation in variations:
        fileName    = os.path.join(outDir,fileBase.replace(".root","_{0}.root".format(variation)))
        if not os.path.exists(fileName):
            return False

    return True

def create_jobs(config,year="2016",jobs_dir="",out_dir="",nFiles=1,checkInput=False,submitFlag=False):
    submissionCmds     = []
    for sample, sample_cfg in config.items():

        sampleJobs_dir  = os.path.join(jobs_dir,sample)
        sampleOut_dir   = os.path.join(out_dir, sample)
        #Create dir to store jobs and dir to store output
        createDirIfNotExist(os.path.join(sampleJobs_dir, 'input'))
        createDirIfNotExist(os.path.join(sampleJobs_dir, 'output'))
        createDirIfNotExist(sampleOut_dir)
        if("ZGamma" in sample or "WGamma" in sample or "Hgamma" in sample):
            exeScript = selection_template.replace("JOB_DIR",sampleJobs_dir)
        else:
            exeScript = selection_template_data.replace("JOB_DIR",sampleJobs_dir)
        nPerJob   = nFiles
        open(os.path.join(sampleJobs_dir, 'input', 'run_{}.sh'.format(sample)), 'w').write(exeScript)
        os.system("chmod +x {0}".format(os.path.join(sampleJobs_dir, 'input', 'run_{}.sh'.format(sample))))

        condor_script = re.sub('EXEC',os.path.join(sampleJobs_dir, 'input', 'run_{}.sh'.format(sample)), selection_condor)
        condor_script = re.sub('ARGFILE',os.path.join(sampleJobs_dir, 'input', 'args_{}.txt'.format(sample)), condor_script)
        condor_script = re.sub('OUTPUT',os.path.join(sampleJobs_dir, 'output'), condor_script)
        open(os.path.join(sampleJobs_dir, 'input', 'condor_{}.condor'.format(sample)), 'w').write(condor_script)

        #Split input files
        skimDirectory   = os.path.join(SKIM_DIR, year,sample)
        skimFiles       = [os.path.join(skimDirectory, f) for f in os.listdir(skimDirectory) if os.path.isfile(os.path.join(skimDirectory, f))]
        job_list        = split_jobs(skimFiles, nPerJob)

        #Create file with arguments to the python script
        argsFile        = open(os.path.join(sampleJobs_dir, 'input', 'args_{}.txt'.format(sample)), 'w')

        nTot            = 0
        nToRun          = 0

        congregateFlag  = checkForCongregateResult(sampleOut_dir,sample)

        for n, l  in enumerate(list(job_list)):
            inputPath       = os.path.join(sampleJobs_dir, 'input', 'input_{}.txt'.format(n))
            outputPath      = os.path.join(sampleOut_dir,'{0}_{1}.root'.format(sample,n))
            open(inputPath, 'w').writelines("{}\n".format(root_file) for root_file in l)
            processedFlag   = checkIfAlreadyProcessed('{0}_{1}.root'.format(sample,n),sampleOut_dir,sample)
            
            if processedFlag==False and congregateFlag==False:
                argsFile.write("-i {0} -o {1} -p {2} -y {3}\n".format(inputPath,outputPath,sample,year))
                nToRun+=1

                if checkInput:
                    removeEmptyFiles(inputPath)

            nTot+=1


        print("{0}:\t{1}/{2} files processed".format(sample,nTot-nToRun,nTot))
        if(congregateFlag):
            print("Hadded result exists")

        if(nToRun!=0 and congregateFlag==False):
            submissionCmds.append("condor_submit {0}".format(os.path.join(sampleJobs_dir, 'input', 'condor_{}.condor'.format(sample))))
    #For some reason the last submission does not correctly take arguments
    #Needs to be submitted to condor manually from shell
    for cmd in submissionCmds:
        print(cmd)
        if(submitFlag):
            os.system(cmd)

    return len(submissionCmds)

def haddResults(outDir):
    os.chdir(outDir)
    directories= [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
    for d in directories:
        if("ZGamma" in d or "WGamma" in d or "Hgamma" in d):
           variations  = ["nom","jerUp","jerDown","jesUp","jesDown","jmsUp","jmsDown","jmrUp","jmrDown"]
        else:
           variations  = ["nom"]
        for variation in variations:
            cmd = "hadd {0}_{1}.root {0}/*{1}*root".format(d,variation)
            if not os.path.exists("{0}_{1}.root".format(d,variation)):
                print(cmd)
                os.system(cmd)

        os.system("rm -r {0}".format(d))
            

def runYear(year,submit=False):
    import json
    datasets    = "selection_configs/{0}.json".format(year)
    outDir      = "{0}/{1}".format(SELECTION_DIR,year)
    jobsDir     = "{0}/{1}".format(SELECTION_JOB_DIR,year)
    filesPerJob = 50

    createDirIfNotExist(outDir)

    if os.listdir(outDir):
        print("Output directory {0} is not empty".format(outDir))
        delFlag = input("Delete files in {0}? [y/N]".format(outDir))
        if delFlag=="y":
            os.system("rm -r {0}/*".format(outDir))


    checkInputFlag = input("Check input for empty skim files? [y/N]")

    with open(datasets) as config_file:
        config  = json.load(config_file)
        nJobs   = create_jobs(config, year=year,out_dir=outDir,jobs_dir=jobsDir,nFiles=filesPerJob,checkInput=checkInputFlag,submitFlag=submit)

    if(nJobs==0):
        print("All files processed, hadding results")
        haddResults(outDir)


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-y', '--year', help='Dataset year',default="2016")
    args   = parser.parse_args()

    if(args.year=="RunII"):
        years=["2016APV","2016","2017","2018"]
    else:
        years = [args.year]
    
    for year in years:
        runYear(year,submit=False)
                    

            

if __name__ == "__main__":
    main()

