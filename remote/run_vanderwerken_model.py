import subprocess
import time
import sys
import os
import shutil
import argparse
import yaml


#
MODEL_DIR = "/n/airoldifs2/lab/gbasse/parallelMCMC/models/"
REMOTE_DIR = "/n/airoldifs2/lab/gbasse/parallelMCMC/remote/"
RESULT_ADAPTIVE_NAME = "adaptive.rdata"
RESULT_BRIDGE_NAME = "bridge.rdata"
CONFIG_NAME = "params.yml"

#
CONFIG_FILE = REMOTE_DIR + CONFIG_NAME
fin = open(CONFIG_FILE)
params = yaml.load(fin)
fin.close()
MODEL_NAME = params["model"]


#
RESULTS_DIR = REMOTE_DIR +  "results/"
SCRIPT_DIR = REMOTE_DIR + "scripts/"
TEMPLATES_DIR = SCRIPT_DIR + "templates/"
ADAPTIVE_NAME = "auto_adaptive-easy-vhd.sbatch"
BRIDGE_NAME = "auto_bridge-easy-vhd.sbatch"
PLOT_NAME = "test_partitioning.R"
SCRIPT_ADAPTIVE_PATH = SCRIPT_DIR + ADAPTIVE_NAME
SCRIPT_BRIDGE_PATH = SCRIPT_DIR + BRIDGE_NAME
SCRIPT_PLOT_PATH = REMOTE_DIR + PLOT_NAME
MODEL = MODEL_DIR + MODEL_NAME


# step 0: setting up directories

def get_current_max_experiment():
    res = os.listdir(RESULTS_DIR)
    res_exp = filter(lambda x: "experiment" in x, res)
    if len(res_exp) == 0:
        return 0
    max_ind = max(map(lambda x: int(x.split("_")[1]), res_exp))
    return max_ind

def create_new_result_directory():
    current_max_ind = get_current_max_experiment()
    new_ind = str(current_max_ind + 1)
    if int(new_ind) <= 9:
        new_ind = "0" + new_ind
    os.mkdir(RESULTS_DIR + "experiment_" + new_ind)
    return RESULTS_DIR + "experiment_" + new_ind + "/"


# First, we start by filling the templates of jobs to submit

## start with the adaptive part
def setup_adaptive_script():
    print("Writing first submission script")
    fin = open(TEMPLATES_DIR + "adaptive-easy-vhd.sbatch", "r")
    fout = open(SCRIPT_ADAPTIVE_PATH, "w")
    script_content = fin.read() % (MODEL, RESULT_ADAPTIVE)
    fin.close()
    fout.write(script_content)
    fout.close()
    shutil.copyfile(SCRIPT_ADAPTIVE_PATH, os.path.join(RESULTS_DIR, ADAPTIVE_NAME))


## now the bridge part
def setup_bridge_script():
    print("Writing second submission script")
    fin = open(TEMPLATES_DIR + "bridge-easy-vhd.sbatch", "r")
    fout = open(SCRIPT_BRIDGE_PATH, "w")
    script_content = fin.read() % (MODEL, RESULT_ADAPTIVE, RESULT_BRIDGE) 
    fin.close()
    fout.write(script_content)
    fout.close()
    shutil.copyfile(SCRIPT_BRIDGE_PATH, os.path.join(RESULTS_DIR, BRIDGE_NAME))


## finally, copy the parameter file
def setup_parameter_file():
    shutil.copyfile(CONFIG_FILE, os.path.join(RESULTS_DIR,  CONFIG_NAME))

#sys.exit()

# Job Submission

# First we submit the first job
def submit_adaptive():
    p = subprocess.Popen(["sbatch", SCRIPT_ADAPTIVE_PATH], stdout=subprocess.PIPE)
    message = p.communicate()[0]
    jobid=message.rstrip("\n").split(" ")[-1]
    print message.rstrip("\n")
    return(jobid)

#sys.exit()

# Wait until first job is finished
def get_status(jobid):
    p = subprocess.Popen(["sacct", "-j", "%s"%jobid, "--format", 'State'], stdout=subprocess.PIPE)
    message = p.communicate()[0]
    status = message.split('\n')[-2].strip(' ')
    return status

def wait_before_bridge(jobid):
    status = get_status(jobid)
    running = False
    while status != 'COMPLETED':
        if status == 'FAILED':
            print "Job failed, terminating..."
            sys.exit() # early exit
        if status == 'PENDING':
            print "Job is still pending..."
        if status == 'RUNNING':
            if not running:
                print "Job is currently running, please wait..."
                running = True
        time.sleep(10)
        status = get_status(jobid)

# We can now submit the second job since the first one was successfully completed
def submit_bridge():
    p = subprocess.Popen(["sbatch", SCRIPT_BRIDGE_PATH], stdout=subprocess.PIPE)
    message = p.communicate()[0]
    jobid=message.rstrip("\n").split(" ")[-1]
    print message.rstrip("\n")

def run_plot():
    p = subprocess.Popen(["Rscript", SCRIPT_PLOT_PATH, MODEL, RESULT_ADAPTIVE, RESULTS_DIR], stdout=subprocess.PIPE)
    message = p.communicate()[0]

def get_options_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--experiment", dest="experiment",
                        help="Which experiment directory to use. By default create a new directory")
    parser.add_argument("-p", "--plot", action="store_true", default=False, dest="plot_p", help="Only run the plotting script")
    parser.add_argument("-b", "--bridge", action="store_true", default=False, dest="bridge_p", help="Only run the bridge script")
    parser.add_argument("-a", "--adaptive", action="store_true", default=False, dest="adaptive_p", 
                        help="Only run the adaptive script")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    # parses arguments
    args = get_options_args()

    # setting up the new results directory
    if(args.experiment is None):
        RESULTS_DIR = create_new_result_directory()
    else:
        RESULTS_DIR = os.path.join(RESULTS_DIR, args.experiment)
    
    RESULT_ADAPTIVE = os.path.join(RESULTS_DIR, RESULT_ADAPTIVE_NAME)
    RESULT_BRIDGE = os.path.join(RESULTS_DIR, RESULT_BRIDGE_NAME)
    
    if args.plot_p or args.bridge_p or args.adaptive_p:
        if args.adaptive_p:
            setup_adaptive_script()
            setup_parameter_file()
            jobid = submit_adaptive()
            if args.bridge_p or args.plot_p:
                wait_before_bridge(jobid)
        if args.bridge_p:
            setup_bridge_script()
            submit_bridge()
        if args.plot_p:
            run_plot()
    else:
        # setting the scripts
        setup_adaptive_script()
        setup_bridge_script()
        setup_parameter_file()

        # submitting jobs
        jobid = submit_adaptive()
        wait_before_bridge(jobid)
        submit_bridge()
        run_plot()
    
