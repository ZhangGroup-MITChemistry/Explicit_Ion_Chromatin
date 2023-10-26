#!/usr/bin/env python

import sys
from numpy import *
from subprocess import *
import fileinput
import shutil

#
#   Lammps input with parallel temp restarting extracted from previous log file
#
def genLammpsInput(restartFolder, outFolder, runId):

    lammps_template = '../../../template.lammps'
    tempList = [300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410]
    
    # We need to decide whether it is the first time run or the restarting run;
    if (runId == 0):

        fh = open(outFolder + "/proteinDna_sim.in", "w")
    
        for line in open(lammps_template).readlines()[0:]:
            items = line.split()
            if len(items)>=2 and items[0] == "read_data":
                fh.write(line)
                fh.write("read_dump start_job_temp${T}.dcd 0 x y z box no format molfile dcd /home/xclin/lib/vmd.old/plugins/LINUXAMD64/molfile\n")
            else:
		# Write everything, including the minimization line in;
		fh.write(line)
    
        fh.close()
	
    else:
        # If it is a restarting run;
        lammps_log_file = restartFolder + "/log.lammps"
        replicaIndex = open(lammps_log_file).readlines()


        # we have to go through this trouble because data file is saved every 1000 steps, but log file is saved every 100 steps
        il = -1
        dumpStep = 5000
        while True:
            step = int(replicaIndex[il].split()[0])
            if mod(step, dumpStep) == 0:
                break
            else:
                il -= 1
        lastLine = replicaIndex[il]

        lastTimeStep = int(lastLine[0])


        tempIndex = lastLine.split()[1::]

        fh = open(outFolder + "/proteinDna_sim.in", "w")
        fh.write("variable              tempIndex   world ")
        for item in tempIndex:
            fh.write(item + " ")
        fh.write("\n")

        # This format will ignore the last line (python range() format)
        for line in open(lammps_template).readlines()[0:-3]:
            items = line.split()
            if len(items)>=2 and items[0] == "read_data" and runId>=1:
                # We need to keep the read_data line;
                fh.write(line)
                for temperature in tempList:
                    cmd = '/pool001/xclin/bin/miniconda3/bin/mdconvert --force -i %d -o DUMP_Extract_temp%d.dcd ../run%02d/DUMP_FILE_temp%d.dcd'%(-1, temperature, runId-1, temperature)
                    q = Popen(cmd, shell=True, stdout=PIPE)
                    q.communicate()
                    src='../run%02d/DUMP_Extract_temp%d.atom'%(runId-1, temperature)
                    dst='./DUMP_FILE_temp%d.atom'%(temperature)
                    shutil.copyfile(src, dst)
#                fh.write("read_dump DUMP_Extract_temp${T}.dcd 0 x y z box yes format molfile dcd /home/xclin/lib/vmd.old/plugins/LINUXAMD64/molfile\n")
#                fh.write("read_restart state_temp${T}.cpt\n")
                    fh.write("read_dump DUMP_FILE_temp${T}.atom 100000 x y z ix iy iz box no format native\n")
            # For the restarting input files, we don't need minimization and velocity regeneration;
            # Actually, we needed velocity regeneration because we are using read_dump here;
            elif runId>=1 and len(items)>=2 and items[0] == "minimize":
                pass
#            elif runId>=1 and len(items)>=2 and items[0] == "velocity":
#                pass
            else:
                fh.write(line)

        fh.write("temper 100000 100 ${T} myfix 3847 58382 ${tempIndex} rigidbodyfix rigidBody\n")
        fh.write("write_restart state_temp${T}.cpt\n")
        fh.write("write_dump all atom DUMP_Extract_temp${T}.atom modify image yes scale no sort id\n")        
	fh.close()


if __name__ == "__main__":
    genLammpsInput(sys.argv[1], sys.argv[2], int(sys.argv[3]))

