#!/usr/bin/env python

import os
import sys
import shutil
from subprocess import *
import fileinput
from numpy import *


################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################


globalDir = os.getcwd()

# Add line to the top of a file;
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


def create_folder(runId):
    # prepare folders
    subdir = "/equil/umbreref_600.0"
    simFolder = "%s/%s/"%(globalDir, subdir)
    rundir = simFolder + "/run%02d/"%(runId)
    if ( os.path.exists(rundir) is not True ):
        os.makedirs(rundir)
    else:
        shutil.rmtree(rundir)
    	os.makedirs(rundir)

    # create qbias file
    # The first run plumed file and the following runs are treated differently, because we don't want to have moving bias again for the following runs
    if (runId == 0):

        fh = open('%s/input.colvars'%(rundir), 'w')
        cv_tmp = fileinput.input('input.colvars_ini.dat')

        for line in cv_tmp:
            fh.write(line)
        fh.close()
    else:
        fh = open('%s/input.colvars'%(rundir), 'w')
        cv_tmp = fileinput.input('input.colvars_again.dat')

        for line in cv_tmp:
            fh.write(line)
        fh.close()

    for line in cv_tmp:
        fh.write(line)
    fh.close()
#
#   -------------------------------------------     #
#
def submit(staId, endId):
    subdir = "./equil/umbreref_600.0"
    simFolder = "%s/%s/"%(globalDir, subdir)
    workDir = os.getcwd()

    # create pbs file
    pbsFile = simFolder + "myjob_eofe.slurm"
    pbs_tmp = fileinput.input('./myjob_eofe.slurm')
    pf = open(pbsFile, 'w')
    # header
    head = '''#!/bin/bash

#SBATCH --job-name=single_nucl_expsol
#SBATCH --array=%d-%d%%1
''' %(staId, endId)
    pf.write(head)
    # tmp
    for line in pbs_tmp:
        pf.write(line)

    pf.close()

    cmd = "cd %s; sbatch myjob_eofe.slurm; cd %s"%(simFolder, workDir)
    q = Popen(cmd, shell=True, stdout=PIPE)
    q.communicate()


if __name__ == '__main__':
    nsta = int(sys.argv[1])
    nend = int(sys.argv[2])
    for runId in my_le_range(nsta, nend, 1):
        create_folder(runId)

    submit(nsta, nend)
