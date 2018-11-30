#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 14:00:31 2018

@author: engels
"""
class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

def err( msg ):
    print("")
    print( bcolors.FAIL + "ERROR! " + bcolors.ENDC + msg)
    print("")

def warn( msg ):
    print("")
    print( bcolors.WARNING + "WARNING! " + bcolors.ENDC + msg)
    print("")

print("----------------------------------------")
print("%sIRENE%s submission preflight  " %(bcolors.OKGREEN, bcolors.ENDC))
print("----------------------------------------")

import sys
import os
import wabbit_tools
import numpy as np

# fetch jobfile from call
jobfile = sys.argv[1]


if os.path.isfile( jobfile ):

    # read jobfile
    with open(jobfile) as f:
        # loop over all lines
        for line in f:
            if "ccc_mprun" in line:
                runline = line
            elif "#MSUB -n " in line:
                cpuline = line
            elif "#MSUB -T" in line:
                wtimeline = line

    # now runline contains the actual run:
    #       ccc_mprun ./wabbit suzuki.ini --memory=550.0GB
    runlist = runline.split()

    progfile   = runlist[1]
    paramsfile = runlist[2]
    memory     = float( runlist[3].replace('--memory=','').replace('GB','').replace('gb','') )

    cpulist = cpuline.split()
    ncpu = float(cpulist[2])

    wtimelist = wtimeline.split()
    wtime = float(wtimelist[2])

    core_per_node = 48
    maxmem = ncpu*3.75 #GB

    if not os.path.isfile(paramsfile):
        err("Paramsfile %s not found!" % (paramsfile))
    else:
        print("Paramsfile %s is there." % (paramsfile))


    print("program           = %s%s%s" % (bcolors.OKBLUE, progfile, bcolors.ENDC) )
    print("paramsfile        = %s%s%s" % (bcolors.OKBLUE, paramsfile, bcolors.ENDC) )
    print("memory in call    = %s%f%s GB" % (bcolors.OKBLUE, memory, bcolors.ENDC) )
    print("max memory        = %s%i%s GB" % (bcolors.OKBLUE, maxmem, bcolors.ENDC) )
    print("max memory (safe) = %s%i%s GB" % (bcolors.OKBLUE, maxmem-5.0, bcolors.ENDC) )
    print("ncpu              = %s%i%s" % (bcolors.OKBLUE, ncpu, bcolors.ENDC) )
    print("wtime (jobfile)   = %s%i%s sec (%f hours)" % (bcolors.OKBLUE, wtime, bcolors.ENDC, wtime/3600.0) )
    wtime_ini = wabbit_tools.get_ini_parameter(paramsfile, "Time", "walltime_max", float)
    # hours to seconds
    wtime_ini *= 3600.0
    print("wtime (inifile)   = %s%i%s sec (%f hours)" % (bcolors.OKBLUE, wtime_ini, bcolors.ENDC, wtime_ini/3600.0) )



    if memory >= 0.98*maxmem:
        err("Memory limit exceeded")
    else:
        print("Memory okay.")

    if abs(ncpu/core_per_node - float(round(ncpu/core_per_node))) > 0.0:
        warn('You did not specify N*48 CPUS')
    else:
        print('nCPU is multiple of cpu_per_node')


    if wtime_ini > wtime:
        err("Walltime in ini file greater than walltime in job file!")
    else:
        print('Walltime is okay')

    Jmax = wabbit_tools.get_ini_parameter( paramsfile, 'Blocks', 'max_treelevel', int)
    L = wabbit_tools.get_ini_parameter( paramsfile, 'Domain', 'domain_size', float, vector=True)
    Bs = wabbit_tools.get_ini_parameter( paramsfile, 'Blocks', 'number_block_nodes', int)
    CFL = wabbit_tools.get_ini_parameter( paramsfile, 'Time', 'CFL', float)
    c0 =  wabbit_tools.get_ini_parameter( paramsfile, 'ACM-new', 'c_0', float)
    ceta =  wabbit_tools.get_ini_parameter( paramsfile, 'VPM', 'C_eta', float)

    dx = L[0]*(2**-Jmax)/(Bs-1)

    # this is an assumption:
    umax = 1.0
    u_eigen = abs(umax) + np.sqrt( umax**2 + c0**2 )
    dt_c0 = CFL*dx / u_eigen
    dt_ceta = 0.99*ceta

    print( "dt_c0   = %e" % (dt_c0))
    print( "dt_ceta = %e" % (dt_ceta))

    print('Launching wabbit_tools ini file check now:')
    wabbit_tools.check_parameters_for_stupid_errors( paramsfile )



else:
    err("Jobfile %s not found" % (jobfile))
    raise ValueError( )
