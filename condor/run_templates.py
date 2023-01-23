import os
import glob
import sys
from paths import SELECTION_DIR, TEMPLATE_JOB_DIR
from templates import templates_template, selection_condor
from pathlib import Path
import re
import stat

#python run_templates.py tight

TEMPLATE_DIR = SELECTION_DIR.replace("selection","templates")

wp   = sys.argv[1]
args = ""
template_jobs_wp_dir = os.path.join(TEMPLATE_JOB_DIR,wp)
Path(template_jobs_wp_dir).mkdir(exist_ok=True, parents=True)

template_jobs_log_dir = os.path.join(template_jobs_wp_dir,"output")
Path(template_jobs_log_dir).mkdir(exist_ok=True, parents=True)

#ParticleNet pnet0
wp_tight  	= {"2016APV":0.9883,"2016":0.9883,"2017":0.9870,"2018":0.9880}
wp_medium 	= {"2016APV":0.9737,"2016":0.9735,"2017":0.9714,"2018":0.9734}
wp_loose	= {"2016APV":0.9088,"2016":0.9137,"2017":0.9105,"2018":0.9172}


wp_vals   	= {"tight":wp_tight, "medium":wp_medium, "loose":wp_loose}
#for year in ["2016","2016APV","2017","2018"]:
for year in ["2017"]:
	evtSelDir = "{0}/{1}/".format(SELECTION_DIR,year)
	tplDir    = "{0}/{1}/{2}".format(TEMPLATE_DIR,wp,year)
	nomFiles  = glob.glob('{0}/*nom.root'.format(evtSelDir))
	for nomFile in nomFiles:
		sample  = nomFile.split("/")[-1].replace("_nom.root","")
		argLine = "{0} {1} {2} {3}\n".format(evtSelDir,sample,tplDir,wp_vals[wp][year])
		args+=argLine
f = open("{0}/args.txt".format(template_jobs_wp_dir),"w")
f.write(args)
f.close()


open(os.path.join(template_jobs_wp_dir, 'run_script.sh'), 'w').write(templates_template)
os.system("chmod +x {0}".format(os.path.join(template_jobs_wp_dir, 'run_script.sh')))

condor_script = re.sub('EXEC',os.path.join(template_jobs_wp_dir, 'run_script.sh'), selection_condor)
condor_script = re.sub('ARGFILE',os.path.join(template_jobs_wp_dir, 'args.txt'), condor_script)
condor_script = re.sub('OUTPUT',os.path.join(template_jobs_wp_dir, 'output'), condor_script)
open(os.path.join(template_jobs_wp_dir, 'condor_submit.condor'), 'w').write(condor_script)

cmd_to_run    = "condor_submit {0}/condor_submit.condor".format(template_jobs_wp_dir)
print(cmd_to_run)