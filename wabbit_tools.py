#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:41:48 2017
@author: engels
"""
import numpy as np

class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

def warn( msg ):
    print( bcolors.WARNING + "WARNING! " + bcolors.ENDC + msg)

def err( msg ):
    print( bcolors.FAIL + "CRITICAL! " + bcolors.ENDC + msg)

def info( msg ):
    print( bcolors.OKBLUE + "Information:  " + bcolors.ENDC + msg)

#%%
def check_parameters_for_stupid_errors( file ):
    """ For a given WABBIT parameter file, check for the most common stupid errors
    the user can commit: Jmax<Jmain, negative time steps, etc.
    """
    import os

    print("We scan %s for stupid errors." % (file) )

    # check if the file exists, at least
    if not os.path.isfile(file):
        raise ValueError("Stupidest error of all: we did not find the INI file.")


    bs = get_ini_parameter( file, 'Blocks', 'number_block_nodes', int)
    if bs % 2 == 0:
        warn('The block size is bs=%i which is an EVEN number.' % (bs) )
    if bs < 3:
        warn('The block size is bs=%i is very small or even negative.' % (bs) )

    # if somebody modifies the standard parameter file, users have to update their
    # ini files they use. this is often forgoten and obnoxious. Hence, if we find
    # value sthat no longer exist, warn the user.
    if exists_ini_parameter( file, "Blocks", "number_data_fields" ) :
        warn('Found deprecated parameter: [Blocks]::number_data_fields')

    if exists_ini_parameter( file, "Physics", "initial_cond" ) :
        warn('Found deprecated parameter: [Physics]::initial_cond')

    if exists_ini_parameter( file, "Dimensionality", "dim" ) :
        warn('Found deprecated parameter: [Dimensionality]::dim')

    if exists_ini_parameter( file, "DomainSize", "Lx" ) :
        warn('Found deprecated parameter: [DomainSize]::Lx')

    if exists_ini_parameter( file, "Time", "time_step_calc" ) :
        warn('Found deprecated parameter: [Time]::time_step_calc')

    ceta = get_ini_parameter( file, 'VPM', 'C_eta', float)
    csponge = get_ini_parameter( file, 'Sponge', 'C_sponge', float)
    penalized = get_ini_parameter( file, 'VPM', 'penalization', bool)
    sponged = get_ini_parameter( file, 'Sponge', 'use_sponge', bool)

    if penalized and sponged:
        if abs(ceta-csponge) > 1.0e-7:
            warn( 'Sponge and penalization parameter are different: C_eta=%e C_sponge=%e' % (ceta,csponge))

    # loop over ini file and check that each non-commented line with a "=" contains the trailing semicolon ";"
    with open(file) as f:
        # loop over all lines
        linenumber = 0
        for line in f:
            # remove trailing & leading spaces
            line = line.strip()
            linenumber += 1
            if line is not "" :
                if line[0] is not "!" and line[0] is not "#" and line[0] is not ";" :
                    if "=" in line and ";" not in line:
                        warn('It appears the line #%i does not contain the semicolon' % (linenumber) )

    restart = get_ini_parameter( file, 'Physics', 'read_from_files', int)
    print("read_from_files=%i" %(restart))

    if restart is 1:
        info("This simulation is being resumed from file")

        infiles = get_ini_parameter( file, 'Physics', 'input_files', str)
        infiles = infiles.split()
        for file in infiles:
            print(file)
            if not os.path.isfile(file):
                raise ValueError("CRUTIAL: read_from_files=1 but infiles NOT found!.")
    else:
        info("This simulation is being started from initial condition (and not from file)")

#%%
def get_ini_parameter( inifile, section, keyword, dtype=float, vector=False ):
    """ From a given ini file, read [Section]::keyword and return the value
        If the value is not found, an error is raised
    """
    import configparser
    import os
    import numpy as np

    # check if the file exists, at least
    if not os.path.isfile(inifile):
        raise ValueError("Stupidest error of all: we did not find the INI file.")

    # initialize parser object
    config = configparser.ConfigParser(allow_no_value=True)
    # read (parse) inifile.
    config.read(inifile)

    # use configparser to find the value
    value_string = config.get( section, keyword, fallback='UNKNOWN')

    # check if that worked
    if value_string is 'UNKNOWN':
        raise ValueError("NOT FOUND! file=%s section=%s keyword=%s" % (inifile, section, keyword) )

    if not vector:
        # configparser returns "0.0;" so remove trailing ";"
        i = value_string.find(';')
        value_string = value_string[:i]
        return dtype(value_string)
    else:
        # you can use the strip() to remove trailing and leading spaces.
        i = value_string.find(';')
        value_string = value_string[:i]
        value_string.strip()
        l = value_string.split()
        return np.asarray( [float(i) for i in l] )



#%%
def exists_ini_parameter( inifile, section, keyword ):
    """ check if a given parameter in the ini file exists or not. can be used to detect
        deprecated entries somebody removed
    """
    found_section = False
    found_parameter = False

    # read jobfile
    with open(inifile) as f:
        # loop over all lines
        for line in f:

            # once found, do not run to next section
            if found_section and line[0] == "[":
                found_section = False

            # until we find the section
            if "["+section+"]" in line:
                found_section = True

            # only if were in the right section the keyword counts
            if found_section and keyword+"=" in line:
                found_parameter = True

    return found_parameter
#%%
def get_inifile_dir( dir ):
    """ For a given simulation dir, return the *.INI file. Warning is issued if the
        choice is not unique. In such cases, we should return the longest one. REASON:
        For insects, we save the wingbeat kinematics in smaller INI files.
    """
    import glob

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir + '/'

    inifile = glob.glob(dir+'*.ini')

    if len(inifile) > 1:
        # if more ethan one ini file is found, try looking for the most appropriate

        print('Warning: we ')
    else:
        return inifile[0]

#%%
def prepare_resuming_backup( inifile, state_vector_prefixes=['ux','uy','uz','p'] ):
    """ we look for the latest *.h5 files
        to resume the simulation, and prepare the INI file accordingly.
        Some errors are caught.
    """
    import numpy as np
    import os
    import glob
    import insect_tools

    # does the ini file exist?
    if not os.path.isfile(inifile):
        raise ValueError("Inifile not found!")

    Tmax = get_ini_parameter(inifile, "Time", "time_max", float)

    # find list of H5 files for first prefix.
    files = glob.glob( state_vector_prefixes[0] + "*.h5" )
    files.sort()

    if not files:
        raise ValueError( "Something is wrong: no h5 files found for resuming" )

    print('Latest file is:           ' + files[-1])
    timestamp = insect_tools.get_timestamp_name( files[-1] )
    t0 = float(timestamp) / 1e6
    print('Latest file is at time:   %f' % (t0))

    d = np.loadtxt('timesteps_info.t')
    t1 = d[-1,0]
    print('Last time stamp in logs is: %f' % (t1))

    # time check when resuming a backup
    if t0 > t1:
        print( "Something is wrong: the latest H5 file is at LATER time than the log files. Is this the right data?" )

    if t0 < 1.0e-6:
        print("Something is wrong: the latest H5 file is almost at t=0. That means no backup has been saved?" )

    if t1 > t0:
        print('Warning: the latest H5 file is younger than the last entry in the log: we will have to compute some times twice.')

    if abs(t1-t0) < 1.0e-4:
        print('Good news: timestamp in H5 file and time in log file match!')

    if t1 >= 0.99*Tmax or t0 >= 0.99*Tmax:
        raise ValueError( "Something is wrong: the run seems to be already finnished!" )

    # check if all required input files exist
    for prefix in state_vector_prefixes:
        if not os.path.isfile( prefix + '_' + timestamp + '.h5'):
            raise ValueError( "file not found!!!! " + prefix + '_' + timestamp + '.h5' )


    # create the string we will put in the ini file
    infiles_string = ""
    for prefix in state_vector_prefixes:
        infiles_string += prefix + '_' + timestamp + '.h5' + ' '

    # remove trailing space:
    infiles_string = infiles_string.strip()
    # add colon
    infiles_string += ';'
    # information (debug)
    print(infiles_string)

    f1 = open( inifile, 'r')
    f2 = open( inifile+'.tmptmp', 'w')
    found, okay1, okay2 = False, False, False

    for line in f1:
        # remove trailing space:
        line_cpy = line.strip()

        if '[Physics]' in line_cpy:
            found = True

        if 'read_from_files=' in line_cpy and found and line_cpy[0] is not ";":
            line = "read_from_files=1;\n"
            okay1 = True

        if 'input_files=' in line_cpy and found and line_cpy[0] is not ";":
            line = "input_files=" + infiles_string + "\n"
            okay2 = True

        f2.write( line )

    f1.close()
    f2.close()

    if okay1 and okay2:
        os.rename( inifile+'.tmptmp', inifile )


#%%
def block_level_distribution_file( file ):
    """ Read a 2D/3D wabbit file and return a list of how many blocks are at the different levels
    """
    import h5py
    import numpy as np

    # open the h5 wabbit file
    fid = h5py.File(file,'r')

    # read treecode table
    b = fid['block_treecode'][:]
    treecode = np.array(b, dtype=float)

    # close file
    fid.close()

    # number of blocks
    Nb = treecode.shape[0]

    # min/max level. required to allocate list!
    jmin, jmax = get_max_min_level( treecode )
    counter = np.zeros(jmax+1)

    # fetch level for each block and count
    for i in range(Nb):
        J = treecode_level(treecode[i,:])
        counter[J] += 1

    return counter

#%%
def read_wabbit_hdf5(file):
    """ Read a wabbit-type HDF5 of block-structured data.
    Return time, x0, dx, box, data, treecode.
    Get number of blocks and blocksize as
    N, Bs = data.shape[0], data.shape[1]
    """
    import h5py
    import numpy as np

    fid = h5py.File(file,'r')
    b = fid['coords_origin'][:]
    x0 = np.array(b, dtype=float)

    b = fid['coords_spacing'][:]
    dx = np.array(b, dtype=float)

    b = fid['blocks'][:]
    data = np.array(b, dtype=float)

    b = fid['block_treecode'][:]
    treecode = np.array(b, dtype=float)

    # get the dataset handle
    dset_id = fid.get('blocks')

    # from the dset handle, read the attributes
    time = dset_id.attrs.get('time')
    iteration = dset_id.attrs.get('iteration')
    box = dset_id.attrs.get('domain-size')

    fid.close()

    jmin, jmax = get_max_min_level( treecode )
    N = data.shape[0]
    Bs = data.shape[1]

    #print("~~~~~~~~~~~~~~~~~~~~~~~~~")
    #print("Reading file %s" % (file) )
    #print("Time=%e it=%i N=%i Bs=%i Jmin=%i Jmax=%i" % (time, iteration, N, Bs, jmin, jmax) )
    #print("~~~~~~~~~~~~~~~~~~~~~~~~~")

    return time, x0, dx, box, data, treecode

#%%
def read_treecode_hdf5(file):
    """ Read a wabbit-type HDF5 of block-structured data.
    same as read_wabbit_hdf5, but reads ONLY the treecode array.
    """
    import h5py
    import numpy as np

    fid = h5py.File(file,'r')

    b = fid['block_treecode'][:]
    treecode = np.array(b, dtype=float)

    return treecode

#%%
def write_wabbit_hdf5( file, time, x0, dx, box, data, treecode, iteration = 0,  ):
    """ Write data from wabbit to an HDF5 file """
    import h5py
    import numpy as np


    Level = np.size(treecode,1)
    if len(data.shape)==4:
        Bs=[0]*3
        # 3d data
        N, Bs[0], Bs[1], Bs[2] = data.shape
        print( "Writing to file=%s max=%e min=%e size=%i %i %i " % (file, np.max(data), np.min(data), Bs[0], Bs[1], Bs[2]) )

    else:
        # 2d data
        Bs=[0]*2
        N, Bs[0], Bs[1] = data.shape
        print("~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("Writing file %s" % (file) )
        print("Time=%e it=%i N=%i Bs[0]=%i Bs[1]=%i Level=%i" % (time, iteration, N, Bs[0], Bs[1],Level) )
        print("~~~~~~~~~~~~~~~~~~~~~~~~~")


    fid = h5py.File( file, 'w')

    fid.create_dataset( 'coords_origin', data=x0, dtype=np.float32 )
    fid.create_dataset( 'coords_spacing', data=dx, dtype=np.float32 )
    fid.create_dataset( 'blocks', data=data, dtype=np.float32 )
    fid.create_dataset( 'block_treecode', data=treecode, dtype=np.float32 )

    fid.close()

    fid = h5py.File(file,'a')
    dset_id = fid.get( 'blocks' )
    dset_id.attrs.create('time', time)
    dset_id.attrs.create('iteration', iteration)
    dset_id.attrs.create('domain-size', box )
    dset_id.attrs.create('total_number_blocks', N )
    fid.close()

#%%


def read_wabbit_hdf5_dir(dir):
    """ Read all h5 files in directory dir.
    Return time, x0, dx, box, data, treecode.
    Use data["phi"][it] to reference quantity phi at iteration it
    """
    import numpy as np
    import re
    import ntpath
    import os

    it=0
    data={'time': [],'x0':[],'dx':[],'treecode':[]}
    # we loop over all files in the given directory
    for file in os.listdir(dir):
        # filter out the good ones (ending with .h5)
        if file.endswith(".h5"):
            # from the file we can get the fieldname
            fieldname=re.split('_',file)[0]
            print(fieldname)
            time, x0, dx, box, field, treecode = read_wabbit_hdf5(os.path.join(dir, file))
            #increase the counter
            data['time'].append(time[0])
            data['x0'].append(x0)
            data['dx'].append(dx)
            data['treecode'].append(treecode)
            if fieldname not in data:
                # add the new field to the dictionary
                data[fieldname]=[]
                data[fieldname].append(field)
            else: # append the field to the existing data field
                data[fieldname].append(field)
            it=it+1
    # the size of the domain
    data['box']=box
    #return time, x0, dx, box, data, treecode
    return data



#%%
def wabbit_error(dir, show=False, norm=2, file=None):
    """ This function is in fact an artefact from earlier times, where we only simulated
    convection/diffusion with WABBIT. Here, the ERROR is computed with respect to an
    exact field (a Gaussian blob) which is set on the grid.
    """
    import numpy as np
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import glob

    if file is None:
        files = glob.glob(dir+'/phi_*.h5')
        files.sort()
        file = files[-1]

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5(file)
    # get number of blocks and blocksize
    N, Bs = data.shape[0], data.shape[1]


    if show:
        plt.close('all')
        plt.figure()
        for i in range(N):
            [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1])
            block = data[i,:,:].copy().transpose()
            y = plt.pcolormesh( Y, X, block, cmap='rainbow' )
            y.set_clim(0,1.0)
            plt.gca().add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], fill=False ))
        plt.colorbar()
        plt.axis('equal')

    # compute error:
    err = []
    exc = []
    for i in range(N):
        [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0] - 0.75, np.arange(Bs)*dx[i,1]+x0[i,1] -0.5 )

        # periodization
        X[X<-0.5] = X[X<-0.5] + 1.0
        Y[Y<-0.5] = Y[Y<-0.5] + 1.0

        exact = np.exp( -((X)**2 + (Y)**2) / 0.01  )
        block = data[i,:,:].copy().transpose()

        err = np.hstack( (err,np.ndarray.flatten(exact - block)) )
        exc = np.hstack( (exc,np.ndarray.flatten(exact)) )

    if show:
        plt.figure()
        c1 = []
        c2 = []
        h = []

        for i in range(N):
            [X,Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0] - 0.75, np.arange(Bs)*dx[i,1]+x0[i,1] -0.5 )
            [X1,Y1] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1] )

            # periodization
            X[X<-0.5] = X[X<-0.5] + 1.0
            Y[Y<-0.5] = Y[Y<-0.5] + 1.0

            exact = np.exp( -((X)**2 + (Y)**2) / 0.01  )
            block = data[i,:,:].copy().transpose()

            err = np.hstack( (err,np.ndarray.flatten(exact - block)) )
            exc = np.hstack( (exc,np.ndarray.flatten(exact)) )

            y = plt.pcolormesh( Y1, X1, exact-block, cmap='rainbow' )
            a = y.get_clim()

            h.append(y)
            c1.append(a[0])
            c2.append(a[1])

            plt.gca().add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], fill=False ))

        for i in range(N):
            h[i].set_clim( (min(c1),max(c2))  )

        plt.colorbar()
        plt.axis('equal')


    return( np.linalg.norm(err, ord=norm) / np.linalg.norm(exc, ord=norm) )

# %%

def key_parameters(dir):
    import configparser
    import glob

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    print(inifile[0])
    config = configparser.ConfigParser()
    config.read(inifile[0])


    adapt_mesh=config.get('Blocks','adapt_mesh',fallback='0')
    adapt_inicond=config.get('Blocks','adapt_inicond',fallback='0')
    eps= config.get('Blocks','eps',fallback='0')

    namestring = adapt_mesh+adapt_inicond+eps

    print(namestring)

# %%

def fetch_dt_dir(dir):

    inifile = get_inifile_dir(dir)
    dt = get_ini_parameter( inifile, 'Time', 'dt_fixed', dtype=float )

    return(dt)



def fetch_Nblocks_dir(dir, return_Bs=False):
    import numpy as np
    import os.path

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    if os.path.isfile(dir+'timesteps_info.t'):
        d = np.loadtxt(dir+'timesteps_info.t')
        if d.shape[1]==5:
            # old data has 5 cols
            N = np.max(d[:,2])
        else:
            # new data we added the col for cpu time (...seufz)
            N = np.max(d[:,3])
    else:
        raise ValueError('timesteps_info.t not found in dir.'+dir)

    if (return_Bs):
        # get blocksize
        Bs = fetch_Bs_dir(dir)
        return(N,Bs)
    else:
        return(N)



def fetch_Nblocks_RHS_dir(dir, return_Bs=False):
    import numpy as np
    import os.path

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    # does the file exist?
    if os.path.isfile(dir+'blocks_per_mpirank_rhs.t'):

        # sometimes, glorious fortran writes ******* eg for iteration counter
        # in which case the loadtxt will fail.
        try:
            # read file
            d = np.loadtxt(dir+'blocks_per_mpirank_rhs.t')

            # unfortunately, we changed the format of this file at some point:
            if d.shape[1] == 10:
                # old format requires computing the number...
                d[:,2] = np.sum( d[:,2:], axis=1 )
                N = np.max(d[:,2])
            else:
                # new format saves Number of blocks
                N = np.max(d[:,2])
        except:
            print("Sometimes, glorious fortran writes ******* eg for iteration counter")
            print("in which case the loadtxt will fail. I think that happended.")
            N = 0

    else:
        raise ValueError('blocks_per_mpirank_rhs.t not found in dir: '+dir)

    if (return_Bs):
        # get blocksize
        Bs = fetch_Bs_dir(dir)
        return(N,Bs)
    else:
        return(N)


def fetch_eps_dir(dir):
    import glob
    import configparser

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    return( get_ini_parameter(inifile[0], 'Blocks', 'eps') )


def fetch_Ceta_dir(dir):
    import glob
    import configparser

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    return( get_ini_parameter(inifile[0], 'VPM', 'C_eta') )


def fetch_Bs_dir(dir):
    import glob
    import configparser

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    return( get_ini_parameter(inifile[0], 'Blocks', 'number_block_nodes') )


def fetch_jmax_dir(dir):
    import glob
    import configparser

    if dir[-1] == '/':
        dir = dir
    else:
        dir = dir+'/'

    inifile = glob.glob(dir+'*.ini')

    if (len(inifile) > 1):
        print('ERROR MORE THAN ONE INI FILE')

    return( get_ini_parameter(inifile[0], 'Blocks', 'max_treelevel') )


def fetch_compression_rate_dir(dir):
    import numpy as np
    if dir[-1] != '/':
        dir = dir+'/'

    d = np.loadtxt(dir+'timesteps_info.t')
    d2 = np.loadtxt(dir+'blocks_per_mpirank_rhs.t')

    # how many blocks was the RHS evaluated on, on all mpiranks?
    nblocks_rhs = d2[:,2]

    # what was the maximum number in time?
    it = np.argmax( nblocks_rhs )

    # what was the maximum level at that time?
#    Jmax = d[it,-1]
    Jmax = fetch_jmax_dir(dir)


    # compute compression rate w.r.t. this level
    compression = nblocks_rhs[it] / ((2**Jmax)**2)

    return( compression )


def convergence_order(N, err):
    """ This is a small function that returns the convergence order, i.e. the least
    squares fit to the log of the two passed lists.
    """
    import numpy as np

    if len(N) != len(err):
        raise ValueError('Convergence order args do not have same length')

    A = np.ones([len(err), 2])
    B = np.ones([len(err), 1])
    # ERR = A*N + B
    for i in range( len(N) ) :
        A[i,0] = np.log(N[i])
        B[i] = np.log(err[i])

    x, residuals, rank, singval  = np.linalg.lstsq(A,B)

    return x[0]

def logfit(N, err):
    """ This is a small function that returns the logfit, i.e. the least
    squares fit to the log of the two passed lists.
    """
    import numpy as np

    if len(N) != len(err):
        raise ValueError('Convergence order args do not have same length')

    A = np.ones([len(err), 2])
    B = np.ones([len(err), 1])
    # ERR = A*N + B
    for i in range( len(N) ) :
        A[i,0] = np.log10(N[i])
        B[i] = np.log10(err[i])

    x, residuals, rank, singval  = np.linalg.lstsq(A,B)

    return x


def plot_wabbit_dir(d, savepng=False):
    import glob

    files = glob.glob(d+'/*.h5')
    files.sort()

    for file in files:
        plot_wabbit_file(file, savepng)


# given a treecode tc, return its level
def treecode_level( tc ):
    for j in range(len(tc)):
            if (tc[j]==-1):
                break

    level = j - 1 + 1 # note one-based level indexing.
    return(level)



# for a treecode list, return max and min level found
def get_max_min_level( treecode ):
    #import numpy as np

    min_level = 99
    max_level = -99
    N = treecode.shape[0]
    for i in range(N):
        tc = treecode[i,:].copy()
        level = treecode_level(tc)
        min_level = min([min_level,level])
        max_level = max([max_level,level])

    return min_level, max_level


# %%
def plot_wabbit_file( file, savepng=False, savepdf=False, cmap='rainbow', caxis=None, caxis_symmetric=False, title=True, mark_blocks=True,
                     gridonly=False, contour=False, ax=None, fig=None, ticks=True, colorbar=True, dpi=300, block_edge_color='k',
                     block_edge_alpha=0.5 , shading='flat',
                     gridonly_coloring='mpirank', flipud=False):

    """ Read a (2D) wabbit file and plot it as a pseudocolor plot.
    Keyword arguments:
        * savepng directly store a image file
        * cmap colormap for data
        * caxis manually specify glbal min / max for color plot
        * block_edge_alpha and block_edge_color defines the transparency
        and color of the blocks
    """
    import numpy as np
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import h5py

    # read procs table, if we want to draw the grid only
    if gridonly:
        fid = h5py.File(file,'r')

        # read procs array from file
        b = fid['procs'][:]
        procs = np.array(b, dtype=float)

        if gridonly_coloring in ['refinement-status', 'refinement_status']:
            b = fid['refinement_status'][:]
            ref_status = np.array(b, dtype=float)

        if gridonly_coloring is 'lgt_id':
            b = fid['lgt_ids'][:]
            lgt_ids = np.array(b, dtype=float)

        fid.close()

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5( file )
    # get number of blocks and blocksize
    N, Bs = data.shape[0], data.shape[1]

    # we need these lists to modify the colorscale, as each block usually gets its own
    # and we would rather like to have a global one.
    h = []
    c1 = []
    c2 = []

    if fig is None:
        fig = plt.gcf()
        fig.clf()
    if ax is None:
        ax = fig.gca()

    # clear axes
    ax.cla()

    # if only the grid is plotted, we use grayscale for the blocks, and for
    # proper scaling we need to know the max/min level in the grid
    jmin, jmax = get_max_min_level( treecode )

    for i in range(N):
        if gridonly_coloring not in ['level', 'white']:
            if not flipud :
                [X, Y] = np.meshgrid( np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1])
            else:
                [X, Y] = np.meshgrid( box[0]-np.arange(Bs)*dx[i,0]+x0[i,0], np.arange(Bs)*dx[i,1]+x0[i,1])

            block = data[i,:,:].copy().transpose()

            if not contour:
                if gridonly:
                    # draw some other qtys (mpirank, lgt_id or refinement-status)
                    if gridonly_coloring in ['mpirank', 'cpu']:
                        block[:,:] = np.round( procs[i] )

                    elif gridonly_coloring in ['refinement-status', 'refinement_status']:
                        block[:,:] = ref_status[i]

                    elif gridonly_coloring is 'lgt_id':
                        block[:,:] = lgt_ids[i]
                        tag = "%i" % (lgt_ids[i])
                        x = Bs/2*dx[i,1]+x0[i,1]
                        if not flipud:
                            y = Bs/2*dx[i,0]+x0[i,0]
                        else:
                            y = box[0] - Bs/2*dx[i,0]+x0[i,0]
                        plt.text( x, y, tag, fontsize=6, horizontalalignment='center',
                                 verticalalignment='center')

                    else:
                        raise ValueError("ERROR! The value for gridonly_coloring is unkown")

                # actual plotting of patch
                hplot = ax.pcolormesh( Y, X, block, cmap=cmap, shading=shading )

            else:
                # contour plot only with actual data
                hplot = ax.contour( Y, X, block, [0.1, 0.2, 0.5, 0.75] )

            # use rasterization for the patch we just draw
            hplot.set_rasterized(True)

            # unfortunately, each patch of pcolor has its own colorbar, so we have to take care
            # that they all use the same.
            h.append(hplot)
            a = hplot.get_clim()
            c1.append(a[0])
            c2.append(a[1])

            if mark_blocks and not gridonly:
                # empty rectangle
                ax.add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0],
                                                fill=False, edgecolor=block_edge_color, alpha=block_edge_alpha ))

            # unfortunately, each patch of pcolor has its own colorbar, so we have to take care
            # that they all use the same.
            if caxis is None:
                if not caxis_symmetric:
                    # automatic colorbar, using min and max throughout all patches
                    for hplots in h:
                        hplots.set_clim( (min(c1),max(c2))  )
                else:
                    # automatic colorbar, but symmetric, using the SMALLER of both absolute values
                    c= min( [abs(min(c1)), max(c2)] )
                    for hplots in h:
                        hplots.set_clim( (-c,c)  )
            else:
                # set fixed (user defined) colorbar for all patches
                for hplots in h:
                    hplots.set_clim( (min(caxis),max(caxis))  )
        else:
            # if we color the blocks simply with grayscale depending on their level
            # well then just draw rectangles. note: you CAN do that with mpirank etc, but
            # then you do not have a colorbar.

            if gridonly_coloring is 'level':
                level = treecode_level( treecode[i,:] )
                color = 0.9 - 0.75*(level-jmin)/(jmax-jmin)
            else:
                color = 1.0
            ax.add_patch( patches.Rectangle( (x0[i,1],x0[i,0]), (Bs-1)*dx[i,1], (Bs-1)*dx[i,0], facecolor=[color,color,color], edgecolor=block_edge_color ))


    if colorbar:
        plt.colorbar(h[0], ax=ax)

    if title:
        plt.title( "t=%f Nb=%i Bs=%i" % (time,N,Bs) )


    if not ticks:
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='left',      # ticks along the bottom edge are off
        top='right',         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off

    plt.axis('tight')
    plt.axes().set_aspect('equal')
    plt.gcf().canvas.draw()

    if not gridonly:
        if savepng:
            plt.savefig( file.replace('h5','png'), dpi=dpi, transparent=True, bbox_inches='tight' )

        if savepdf:
            plt.savefig( file.replace('h5','pdf'), bbox_inches='tight' )
    else:
        if savepng:
            plt.savefig( file.replace('.h5','-grid.png'), dpi=dpi, transparent=True, bbox_inches='tight' )

        if savepdf:
            plt.savefig( file.replace('.h5','-grid.pdf'), bbox_inches='tight' )

#%%
def wabbit_error_vs_flusi(fname_wabbit, fname_flusi, norm=2, dim=2):
    """ Compute the error (in some norm) wrt a flusi field.
    Useful for example for the half-swirl test where no exact solution is available
    at mid-time (the time of maximum distortion)
    NOTE: We require the wabbit-field to be already full (but still in block-data) so run
    ./wabbit-post 2D --sparse-to-dense input_00.h5 output_00.h5
    first
    """
    import numpy as np
    import insect_tools
    import matplotlib.pyplot as plt

    # read in flusi's reference solution
    time_ref, box_ref, origin_ref, data_ref = insect_tools.read_flusi_HDF5( fname_flusi )
    ny = data_ref.shape[1]

    # wabbit field to be analyzed: note has to be full already
    time, x0, dx, box, data, treecode = read_wabbit_hdf5( fname_wabbit )
    Bs = data.shape[1]
    Jflusi = (np.log2(ny/(Bs-1)))
    print("Flusi resolution: %i %i %i so desired level is Jmax=%f" % (data_ref.shape[0], data_ref.shape[2], data_ref.shape[2], Jflusi) )

    if dim==2:
        # squeeze 3D flusi field (where dim0 == 1) to true 2d data
        data_ref = data_ref[0,:,:].copy()
        box_ref = box_ref[1:2].copy()

    # convert wabbit to dense field
    data_dense, box_dense = dense_matrix( x0, dx, data, treecode, dim )

    if data_dense.shape[0] < data_ref.shape[0]:
        # both datasets have different size
        s = int( data_ref.shape[0] / data_dense.shape[0] )
        data_ref = data_ref[::s, ::s].copy()
        raise ValueError("ERROR! Both fields are not a the same resolutionn")

    if data_dense.shape[0] > data_ref.shape[0]:
        raise ValueError("ERROR! The reference solution is not fine enough for the comparison")

    # we need to transpose the flusi data...
    data_ref = data_ref.transpose()

    err = np.ndarray.flatten(data_dense-data_ref)
    exc = np.ndarray.flatten(data_ref)

    err = np.linalg.norm(err, ord=norm) / np.linalg.norm(exc, ord=norm)
    print( "error was e=%e" % (err) )

    return err


#%%
def flusi_error_vs_flusi(fname_flusi1, fname_flusi2, norm=2, dim=2):
    """ compute error given two flusi fields
    """
    import numpy as np
    import insect_tools

    # read in flusi's reference solution
    time_ref, box_ref, origin_ref, data_ref = insect_tools.read_flusi_HDF5( fname_flusi1 )

    time, box, origin, data_dense = insect_tools.read_flusi_HDF5( fname_flusi2 )

    if len(data_ref) is not len(data_dense):
        raise ValueError("ERROR! Both fields are not a the same resolutionn")

    err = np.ndarray.flatten(data_dense-data_ref)
    exc = np.ndarray.flatten(data_ref)

    err = np.linalg.norm(err, ord=norm) / np.linalg.norm(exc, ord=norm)

    print( "error was e=%e" % (err) )

    return err


#%%
def to_dense_grid( fname_in, fname_out):
    """ Convert a WABBIT grid to a full dense grid in a single matrix.
    We asssume here that interpolation has already been performed, i.e. all
    blocks are on the same (finest) level.
    """
    import numpy as np
    import insect_tools
    import matplotlib.pyplot as plt

    # read data
    time, x0, dx, box, data, treecode = read_wabbit_hdf5( fname_in )

    # convert blocks to complete matrix
    field, box = dense_matrix(  x0, dx, data, treecode )

    plt.figure()
    plt.pcolormesh(field)
    plt.axis('equal')
    plt.colorbar()

    # write data to FLUSI-type hdf file
    insect_tools.write_flusi_HDF5( fname_out, time, box, field)

#%%
def compare_two_grids( treecode1, treecode2 ):
    """ Compare two grids. The number returned is the % of blocks from treecode1
    which have also been found in treecode2 """
    import numpy as np

    common_blocks = 0

    for i in range(treecode1.shape[0]):
        # we look for this tree code in the second array
        code1 = treecode1[i,:]

        for j in range(treecode2.shape[0]):
            code2 = treecode2[j,:]
            if np.linalg.norm( code2-code1 ) < 1.0e-13:
                # found code1 in the second array
                common_blocks += 1
                break

    print( "Nblocks1=%i NBlocks2=%i common blocks=%i" % (treecode1.shape[0], treecode2.shape[0], common_blocks) )

    return common_blocks / treecode1.shape[0]


#%%
def overwrite_block_data_with_level(treecode, data):
    """On all blocks of the data array, replace any function values by the level of the block"""

    N = treecode.shape[0]
    for i in range(N):
        level = treecode_level(treecode[i,:])
        data[i,:,:,:] = float( level )

    return data


#%%
def dense_matrix(  x0, dx, data, treecode, dim=2 ):

    import math
    """ Convert a WABBIT grid to a full dense grid in a single matrix.
    We asssume here that interpolation has already been performed, i.e. all
    blocks are on the same (finest) level.
    returns the full matrix and the domain size. Note matrix is periodic and can
    directly be compared to FLUSI-style results (x=L * 0:nx-1/nx)
    """
    # number of blocks
    N = data.shape[0]
    # size of each block
    Bs = data.shape[1]

    # check if all blocks are on the same level or not
    jmin, jmax = get_max_min_level( treecode )
    if jmin != jmax:
        print("ERROR! not an equidistant grid yet...")
        raise ValueError("ERROR! not an equidistant grid yet...")

    # note skipping of redundant points, hence the -1
    if dim==2:
        nx = int( np.sqrt(N)*(Bs-1) )
    else:
        nx = int( math.pow(N,1.0/dim)*(Bs-1)) +1

    # all spacings should be the same - it does not matter which one we use.
    ddx = dx[0,0]

    #print("Number of blocks %i" % (N))
    #print("Spacing %e domain %e" % (ddx, ddx*nx))

    if dim==2:
        # allocate target field
        field = np.zeros([nx,nx])
        print("Dense field resolution %i x %i" % (nx, nx) )
        # domain size
        box = [dx[0,0]*nx, dx[0,1]*nx]
    else:
        # allocate target field
        field = np.zeros([nx,nx,nx])
        print("Dense field resolution %i x %i x %i" % (nx, nx, nx) )
        # domain size
        box = [dx[0,0]*nx, dx[0,1]*nx, dx[0,2]*nx]

    for i in range(N):
        # get starting index of block
        ix0 = int(round(x0[i,0]/dx[i,0]))
        iy0 = int(round(x0[i,1]/dx[i,1]))
        if dim==3:
            iz0 = int(round(x0[i,2]/dx[i,2]))
            # copy block content to data field. Note we skip the last points, which
            # are the redundant nodes.
            field[ ix0:ix0+Bs-1, iy0:iy0+Bs-1, iz0:iz0+Bs-1 ] = data[i,0:-1,0:-1, 0:-1]
        else:
            # copy block content to data field. Note we skip the last points, which
            # are the redundant nodes.
            field[ ix0:ix0+Bs-1, iy0:iy0+Bs-1 ] = data[i,0:-1,0:-1]

    return(field, box)

#%%
def prediction1D( signal1, order=4 ):
    import numpy as np

    N = len(signal1)

    signal_interp = np.zeros( [2*N-1] ) - 7

#    function f_fine = prediction1D(f_coarse)
#    % this is the multiresolution predition operator.
#    % it pushes a signal from a coarser level to the next higher by
#    % interpolation

    signal_interp[0::2] = signal1

    a =   9.0/16.0
    b =  -1.0/16.0

    signal_interp[1]  = (5./16.)*signal1[0] + (15./16.)*signal1[1] -(5./16.)*signal1[2] +(1./16.)*signal1[3]
    signal_interp[-2] = (1./16.)*signal1[-4] -(5./16.)*signal1[-3] +(15./16.)*signal1[-2] +(5./16.)*signal1[-1]

    for k in range(1,N-2):
        signal_interp[2*k+1] = a*signal1[k]+a*signal1[k+1]+b*signal1[k-1]+b*signal1[k+2]

    return signal_interp


#%%
# calculates treecode from the index of the block
# Note: other then in fortran we start counting from 0
def blockindex2treecode(ix, dim, treeN):

    treecode = np.zeros(treeN)
    for d in range(dim):
        # convert block index to binary
        binary = list(format(ix[d],"b"))
        # flip array and convert to numpy
        binary = np.asarray(binary[::-1],dtype=int)
        # sum up treecodes
        lt = np.size(binary)
        treecode[:lt] = treecode[:lt] + binary *2**d

    # flip again befor returning array
    return treecode[::-1]

#%%
def dense_to_wabbit_hdf5(ddata, name , Bs, box_size = None, time = 0, iteration = 0):

    """
    This function creates a <name>_<time>.h5 file with the wabbit
    block structure from a given dense data matrix.
    Therefore the dense data is divided into equal blocks,
    similar as sparse_to_dense option in wabbit-post.
     Input:
         - ddata... 2d/3D array of the data you want to write to a file
         - name ... prefix name of the datafile
         - Bs   ... number of grid points per block
                    is a 2D/3D dimensional array with Bs[0] being the number of
                    grid points in x direction etc.
                    The data size in each dimension has to be dividable by Bs.
                    Optional Input:
                        - box_size... 2D/3D array of the size of your box
                        - time    ... time of the data
                        - iteration ... iteration of the time snappshot
    """
    # concatenate filename in the same style as wabbit does
    fname = name + "_%12.12d" % int(time*1e6) + ".h5"
    Ndim = ddata.ndim
    Nsize = np.asarray(ddata.shape)
    level = 0

    #########################################################
    # do some initial checks on the input data
    # 1) check if the size of the domain is given
    if box_size is None:
        box = np.ones(Ndim)
    else:
        box = box_size

    if (type(Bs) is int):
        Bs = [Bs]*Ndim
    # 2) check if number of lattice points is block decomposable
    # loop over all dimensions
    for d in range(Ndim):
        # check if Block is devidable by Bs
        if (np.remainder(Nsize[d], Bs[d]-1) == 0):
            if(is_power2(Nsize[d]//(Bs[d]-1))):
                level = max(level, int(np.log2(Nsize[d]/(Bs[d]-1))))
            else:
                err("Number of Intervals must be a power of 2!")
        else:
            err("datasize must be multiple of Bs!")
    # 3) check dimension of array:
    if Ndim < 2 or Ndim > 3:
        err("dimensions are wrong")
    #########################################################

    # assume periodicity:
    data = np.zeros(Nsize+1)
    if Ndim == 2:
        data[:-1, :-1] = ddata
        # copy first row and column for periodicity
        data[-1, :] = data[0, :]
        data[:, -1] = data[:, 0]
    else:
        data[:-1, :-1, :-1] = ddata
        # copy for periodicity
        data[-1, :, :] = data[0, :, :]
        data[:, -1, :] = data[:, 0, :]
        data[:, :, -1] = data[:, :, 0]

    # number of intervals in each dimension
    Nintervals = [int(2**level)]*Ndim  # note [val]*3 means [val, val , val]
    Lintervals = box/Nintervals
    x0 = []
    treecode = []
    dx = []
    bdata = []
    if Ndim == 3:
        for ibx in range(Nintervals[0]):
            for iby in range(Nintervals[1]):
                for ibz in range(Nintervals[2]):
                    x0.append([ibx, iby, ibz]*Lintervals)
                    dx.append(Lintervals/(Bs-1))
                    lower = x0[-1]/dx[-1]
                    upper = lower + Bs
                    treecode.append(blockindex2treecode([ibx, iby, ibz], 3, level))
                    bdata.append(data[lower[0]:upper[0], lower[1]:upper[1], lower[2]:upper[2]])
    else:
        for ibx in range(Nintervals[0]):
            for iby in range(Nintervals[1]):
                x0.append([ibx, iby]*Lintervals)
                dx.append(Lintervals/(Bs-1))
                lower = x0[-1]/dx[-1]
                upper = lower + Bs
                treecode.append(blockindex2treecode([ibx, iby], 2, level))
                bdata.append(data[lower[0]:upper[0], lower[1]:upper[1]])

    x0 = np.asarray(x0)
    dx = np.asarray(dx)
    treecode = np.asarray(treecode)
    block_data = np.asarray(bdata)

    write_wabbit_hdf5(fname, time, x0, dx, box, block_data, treecode, iteration )


# %%
def is_power2(num):
    'states if a number is a power of two'
    return num != 0 and ((num & (num - 1)) == 0)




#%%
def dir_dense_matrix_2D( direc, var ):
    #Extract all the wabbit fields (all variables or just one) to a field type matrix of dimensions(n_timesteps x n_x x n_y)
    import glob
    import numpy as np
    #import wabbit_tools as w2p
    if var is None:
        path = direc+'/*.h5'
        files=glob.glob(path)
        files.sort()
    
    else:
        path = direc+'/'+var+'*.h5'
        files=glob.glob(path)
        files.sort()


    noF = len(files) # number of total files
    fin = m=read_wabbit_hdf5(files[0])
    hel = len(dense_matrix(fin[1], fin[2], fin[4], fin[5], dim=2)[0])
    field = np.zeros((noF,hel,hel))
    for file in range(noF):
        m=read_wabbit_hdf5(files[file])
        mm=dense_matrix(m[1],m[2],m[4],m[5], dim=2)
        field[file,:,:]=mm[0]       


    return field

#%%
def space_time_slice(direc,var,row,column):
    import numpy as np
    X_T_MAT = dir_dense_matrix_2D(direc,var)
    if row is None and column is None:
        print('!!!ERROR!!! A row OR a column has to be specified')
        X_T = 0
    elif row is not None and column is not None:
        print('!!!ERROR!!! A row OR a column has to be specified')
        X_T = 0
    elif column is None:
        X_T = np.zeros((np.shape(X_T_MAT)[0],np.shape(X_T_MAT)[2]))
        for line in range(np.shape(X_T)[0]):
            X_T[line,:]= X_T_MAT[line,row,:]
    elif row is None:
        X_T = np.zeros((np.shape(X_T_MAT)[0],np.shape(X_T_MAT)[1]))
        for line in range(np.shape(X_T)[0]):
            X_T[line,:]= X_T_MAT[line,:,column]
    return X_T

#%%










