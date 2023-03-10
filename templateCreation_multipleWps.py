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


variations = ["nom","jesUp","jesDown","jerUp","jerDown","jmsUp","jmsDown","jmrUp","jmrDown","pnetUp","pnetDown"]
for variation in variations:
    if(variation!="nom"):
        if not("Hgamma" in sample or "ZGamma" in sample or "WGamma" in sample):
            #Running variations only on some processes
            continue
        if(("pnet" in variation) and not ("ZGamma" in sample or "Hgamma" in sample)):
            #Xbb SF variation only on processes with X->bb
            continue           

    inputTag = variation
    if "pnet" in variation:
        inputTag = "nom"

    inputFile = "{0}/{1}_{2}.root".format(iDir,sample,inputTag)
    outputFile = os.path.join(outDir,"nonScaled/",sample)
    outputFile = outputFile+".root"
    createDirIfNotExist(outDir+"/nonScaled")
    createDirIfNotExist(outDir+"/scaled")
    if(variation=="nom"):
        mode="RECREATE"
    else:
        mode="UPDATE"
    cmd = "python templateMaker_multipleWps.py -i {0} -o {1} -y {2} -p {3} -v {4} -m {5} -w {6} -w {7}".format(inputFile,outputFile,year,sample,variation,mode,wpUp,wpLo)
    print(cmd)
    os.system(cmd)
