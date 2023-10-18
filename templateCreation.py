import os
import sys

def createDirIfNotExist(path):
    if not os.path.exists(path):
        print("CREATING DIR: ", path)
        os.makedirs(path)

iDir    = sys.argv[1]
sample  = sys.argv[2]
outDir  = sys.argv[3]
wpUp   = sys.argv[4]
wpLo   = sys.argv[5]

if("2016/" in iDir):
    year="2016"
if("2016APV/" in iDir):
    year="2016APV"
if("2017/" in iDir):
    year="2017"
if("2018/" in iDir):
    year="2018"

if "_CR" in outDir:
    CRflag = True
else:
    CRflag = False


#variations = ["nom","jesUp","jesDown","jerUp","jerDown","jmsUp","jmsDown","jmrUp","jmrDown","pnetUp","pnetDown","photonEsUp","photonEsDown","photonErUp","photonErDown"]

#no pnet
variations       = ["nom","jesUp","jesDown","jerUp","jerDown","jmsUp","jmsDown","jmrUp","jmrDown","photonEsUp","photonEsDown","photonErUp","photonErDown"]
varProcessesSR   = ["Hgamma", "ZGamma", "WGamma"] 
for variation in variations:
    if(variation!="nom"):
        #only run variations on some processes
        if(CRflag and ("JetHT" in sample or "TTbar" in sample or "QCD" in sample)):
            continue
        if(not CRflag and not ("Hgamma" in sample or "ZGamma" in sample or "WGamma" in sample)):
            continue

    if("photon" in variation and CRflag):
        #Don't run photon variations in no-photon CR
        continue

    inputTag = variation

    inputFile = "{0}/{1}_{2}.root".format(iDir,sample,inputTag)
    outputFile = os.path.join(outDir,"nonScaled/",sample)
    outputFile = outputFile+".root"
    createDirIfNotExist(outDir+"/nonScaled")
    createDirIfNotExist(outDir+"/scaled")
    if(variation=="nom"):
        mode="RECREATE"
    else:
        mode="UPDATE"
    if CRflag:
        cmd = "python templateMaker_CR.py -i {0} -o {1} -y {2} -p {3} -v {4} -m {5} -w {6} -w {7}".format(inputFile,outputFile,year,sample,variation,mode,wpUp,wpLo)
    else:
        cmd = "python templateMaker.py -i {0} -o {1} -y {2} -p {3} -v {4} -m {5} -w {6} -w {7}".format(inputFile,outputFile,year,sample,variation,mode,wpUp,wpLo)
    print(cmd)
    os.system(cmd)
