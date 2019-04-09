import time, os, hashlib
import configobj as cfg
import os.path as op
import scipy
import numpy as np
import subprocess as sp
import logging
from glob import glob
import os.path as opa

THIS_PATH = opa.dirname(opa.abspath(__file__))
EFT_PATH = opa.abspath(THIS_PATH)
metafilpath = opa.join(EFT_PATH,'metafil/')
if opa.isdir(metafilpath):
    import sys
    if metafilpath not in sys.path: 
        sys.path.append(metafilpath)
    import metafil
else:
    raise Exception('Module not found at ' + metafilpath)

# zbEFT file naming convetion
FBSNM = 'PowerSpectra'
FEFT = '1loop'
FLIN = 'Linear'
Trianglepath=opa.abspath(opa.join(EFT_PATH,'/TrianglesConfiguration')) 

# Parameters for config files
DEFAULTCONFIG_CLASS = cfg.ConfigObj({
                'ln10^{10}A_s':3.094,
                'n_s':0.95,
                'h':0.7,
                'omega_b':0.022393,
                'omega_cdm':0.111867,
                'output':'mPk', 
                'P_k_max_h/Mpc':20,
                'root':'class_',
                'headers':'no',
                'format':'camb',
                 'YHe':0.2454,
                 'z_pk':0.55,
                # 'tau_reio':0.05430842,
                # 'N_ur':2.03066666667,
                # 'omega_cdm':0.1201075,
                # 'N_ncdm':1,
                # 'omega_ncdm':0.0006451439
                }) 
DEFAULTCONFIG_zbEFT = cfg.ConfigObj({
                'knl':1,
                'km':1,
                'nbar':3e-4,
                'PathToLinearPowerSpectrum':'class_pk.dat',
                'PathToOutput':'output',
                'ComputePowerSpectrum':'yes',
                'ComputeBispectrum':'no',
                'EpsRel_IntegrBispectrumAP':1e-3,
                'PathToTriangles' : Trianglepath,
                # The below should never be set separately from CLASS!
                'z_pk':DEFAULTCONFIG_CLASS['z_pk'],
                'omega_b':DEFAULTCONFIG_CLASS['omega_b'],
                'omega_cdm':DEFAULTCONFIG_CLASS['omega_cdm'],
                'ln10^{10}A_s':DEFAULTCONFIG_CLASS['ln10^{10}A_s'],
                'n_s':DEFAULTCONFIG_CLASS['n_s'],
                'h':DEFAULTCONFIG_CLASS['h'],
                'aperp':1,
                'apar':1
                })
F2 = 0.67**2
DEFAULTCONFIG_zbEFTw = cfg.ConfigObj({
                'outpath':'./output/',
                'basename':'zbEFT',
                'pid':'', # Enable a process ID stamp to avoid overwriting files
                'CLASS_configf':'class.ini',
                'zbEFT_configf':'zbEFT.ini',
                'zbEFTw_configf':'zbEFTw.ini',
                'logfile':'zbEFT.log',
                'CLASS_path':'../../class/',
                'CLASS_exe':'class',
                'CLASS_pre':'pk_ref.pre',
                'zbEFT_path':'./',
                'zbEFT_exe':'cbird',
                'DM':False, # EFT parameter values are different for DM-only to galaxies
                'kren':0.001, # This needs a k-dependent factor that we set in the code
                'b1':1.,'b2':1.,'b3':1.,'b4':0., # bias parameters
                'b5':-61./315.,'b6':-F2*166./105., 'b7':-F2*2.*(46.+35*F2)/105., # counter-term parameters
                'b8':0.,'b9':1.,'b10':0. # stochastic-term parameters
                })
def cH(Om,a):
    return np.sqrt(Om/a+a**2*(1-Om))
def DgN(Om,a):

    return 5./2*Om*cH(Om,a)/a*scipy.integrate.quad(lambda x:cH(Om,x)**-3,0,a)[0]
def fN(Om,a):
    return (Om*(5*a - 3*DgN(Om,a)))/(2.*(a**3*(1 - Om) + Om)*DgN(Om,a))

def Classtomultipoles(setk,Pkclass,Om,zpk):
    kmulti=np.concatenate([setk,setk,setk])
    f=fN(Om,1./(+zpk))
    P0k=(1 + (2*f)/3. + f**2/5.)*Pkclass
    P2k=((4*f*(7 + 3*f))/21.)*Pkclass
    P4k=((8*f**2)/35.)*Pkclass
    return np.array([kmulti,np.concatenate([P0k,P2k,P4k])])
    
def multipolestoClass(kmulti,Pmulti,Om,zpk):
    kclass=kmulti[:len(kmulti)/3]
    P0k=Pmulti[:len(kmulti)/3]
    P2k=Pmulti[len(kmulti)/3:2*len(kmulti)/3]
    P4k=Pmulti[2*len(kmulti)/3:3*len(kmulti)/3]
    f=fN(Om,1./(+zpk))
    Pkclass=P0k/(1 + (2*f)/3. + f**2/5.)
    return np.array([kclass,Pkclass])
def get_timestamp(prefix='',suffix=''):
    t = (2016,1,1,0,0,0,0,0,0)
    if type(prefix) is list or type(prefix) == np.array: prefix = '_'.join(prefix)
    if type(suffix) is list or type(suffix) == np.array: suffix = '_'.join(suffix)
    return str(prefix)+'_'+str(int(time.time() - time.mktime(t)))+'_'+str(suffix)

def make_hash(config_dict):
    hasher = hashlib.sha1()
    hasher.update(str(config_dict))
    return hasher.hexdigest()[:7]

def combine_hashes(configs):
    return '-'.join([str(make_hash(c)) for c in configs])

def setup_outputs(config_class,config_zbEFT,config_zbEFTw):

    # Create the name for the unique output directory using this structure:
    # path/basename_timestamp_githash_classconfighash_zbeftconfighash_zbeftwconfighash/
    git = metafil.GitEnv()
    basestr = get_timestamp(config_zbEFTw['basename'],
                            [git.get_hash(7), combine_hashes([config_class,config_zbEFT,config_zbEFTw])])
    if config_zbEFTw['pid']: basestr+='_'+config_zbEFTw['pid']
    if basestr in config_zbEFTw['outpath']:
        raise IOError('It seems like you are reusing an ini file from a previous run. '+\
            'Please create a copy and set the outpath parameter to a parent folder!')
    config_zbEFTw['outpath'] = op.abspath(op.join(config_zbEFTw['outpath'],basestr))
    
    # Replace the paths in the config so that they can be easily accessed as absolute    
    # Check what has been provided - prepend folder if none is present
    [
        config_zbEFTw['CLASS_exe'], # Where to find the class executable
        config_zbEFTw['CLASS_pre'],
        config_class['root'], # Where to save the class linear pk
        config_zbEFTw['zbEFT_exe'], # Where to find the zbeft executable
        config_zbEFT['PathToOutput'], # General output directory for zbeft
        config_zbEFTw['logfile'], # What to save the wzbeft log file to
    ] = safe_prepend_folder([
            config_zbEFTw['CLASS_path'],
            config_zbEFTw['CLASS_path'],
            config_zbEFTw['outpath'],
            config_zbEFTw['zbEFT_path'],
            config_zbEFTw['outpath'],
            config_zbEFTw['outpath'],
            ],[
            config_zbEFTw['CLASS_exe'],
            config_zbEFTw['CLASS_pre'],
            config_class['root'],
            config_zbEFTw['zbEFT_exe'],
            '', # PathToOutput should always be set automatically!!
            config_zbEFTw['logfile'],
    ])

    # Finish the CLASS power spectrum file name
    config_zbEFT['PathToLinearPowerSpectrum'] = config_class['root']+'pk.dat' # hardcoded, because this is what class does!
    #config_zbEFT['PathToLinearPowerSpectrum'] = config_class['root']+'pk_cb.dat'
    
    # Make the output folder if it doesn't yet exist
    try:
        if not op.isdir(config_zbEFTw['outpath']):
            os.makedirs(config_zbEFTw['outpath'])
        else:
            os.makedirs(config_zbEFT['PathToOutput'])
    except:
        raise IOError('Cannot create directory and subdirectories: '+config_zbEFTw['outpath'])
    else:
        logging.info('Created new output directory and subdirectories: '+config_zbEFTw['outpath']) 
    
    return basestr

def make_configfiles(configs, filekeys=['CLASS_configf','zbEFT_configf','zbEFTw_configf'], indfmeta=2):

    for i,c in enumerate(configs):
        found=True
        for fk in filekeys:
            if fk not in c.keys(): found=False
        if found: break
    if not found: raise Exception('Configuration for the zbEFT wrapper not found in input!')
    if filekeys[i] not in configs[i].keys(): raise Exception('The order of the inputs must match!')

    for c,f in zip(configs,filekeys):
        c.filename = op.join(configs[indfmeta]['outpath'],configs[indfmeta][f])
        c.write()
        [configs[indfmeta][f]] = safe_prepend_folder(configs[indfmeta]['outpath'],configs[indfmeta][f])

    return 0

def get_config(bigconfig_file=None, bigconfig=None, cat=False):
    """ Splits the configuration into 3: one for CLASS, one for the EFT code, 
        one for the wrapper. Whatever is unset should be set to the default value.
        Inputs
        ------
        EITHER
        bigconfig_file : str
            path to the full config file
        OR
        bigconfig : dict or cfg.ConfigObj
            the full configuration in a dict or ConfigObj
    """

    if bigconfig_file is not None:
        if bigconfig is not None: raise Exception('get_config only takes one input!')
        bigconfig_file = op.abspath(bigconfig_file)
        if not op.isfile(bigconfig_file):
            raise IOError(bigconfig_file + ' file not found!')
        bigconfig = cfg.ConfigObj(bigconfig_file)
        bigconfig.filename = bigconfig_file

    check_config(bigconfig)

    allkeys = DEFAULTCONFIG_CLASS.keys()
    saved = []

    cclass = cfg.ConfigObj(); czbEFT = cfg.ConfigObj(); czbEFTw = cfg.ConfigObj()
    for key in DEFAULTCONFIG_CLASS.keys():
        try:
            cclass[key] = bigconfig[key]
        except KeyError:
            cclass[key] = DEFAULTCONFIG_CLASS[key]
        saved.append(key)

    for key in DEFAULTCONFIG_zbEFT.keys():
        try:
            czbEFT[key] = bigconfig[key]
        except KeyError:
            czbEFT[key] = DEFAULTCONFIG_zbEFT[key]
        saved.append(key)

    for key in DEFAULTCONFIG_zbEFTw.keys():
        try:
            czbEFTw[key] = bigconfig[key]
        except KeyError:
            czbEFTw[key] = DEFAULTCONFIG_zbEFTw[key]
        saved.append(key)

    for key in bigconfig.keys():
        if key not in saved:
            czbEFTw[key] = bigconfig[key]
            saved.append(key)

    if not cat:
        return cclass, czbEFT, czbEFTw
    else:
        return dict(cclass.items() + czbEFT.items() + czbEFTw.items())

def check_config(config):

    # Things I know might go wrong
    if 'z' in config.keys() and 'z_pk' in config.keys():
        if config['z']!=config['z_pk']:
            raise IOError('Your z and z_pk should be the same!!')

def prepend_folder(paths,fnames):
    if not type(fnames) == list: fnames = [fnames]
    if not type(paths) == list: 
        paths = [paths]*len(fnames)
    elif len(paths) != len(fnames):
        raise Exception("Can't combine path and file lists of lengths "+str(len(paths))+" and "+str(len(fnames))+"!")
    for i,[p,fn] in enumerate(zip(paths,fnames)):
        fnames[i] = op.abspath(op.join(p,fn))
    return fnames

def safe_prepend_folder(paths,fnames):
    for i, fn in enumerate(fnames):
        if op.isdir(op.dirname(fn.strip('/'))):
            paths[i] = ''
    return prepend_folder(paths,fnames)

def runcommand(command, logfile, outfile=None, cwd='.'):
    with open(logfile,"wb") as out:
        process = sp.Popen(command, stdout=out, stderr=out, cwd=cwd)
        try:
            process.wait()
        except KeyboardInterrupt as e:
            process.kill()
            raise e
    
    if outfile is not None:
        if len(glob(outfile))==0:
            errmsg = 'Command: '
            errmsg+= ' '.join(command)
            errmsg+= ' failed for unknown reasons.'
            errmsg+= ' The following expected file not found:'
            errmsg+= outfile
            errmsg+= '.'
            raise Exception(errmsg)
    return 0

def read_file(filenamebase, multipole, path='./', verb=0):
    
    filename = op.abspath(op.join(path,filenamebase+'_l'+str(multipole)+'.dat'))
    if verb: print 'reading from '+filename

    # Get the column names first
    with open(filename,'r') as f:
        cols = f.readline().split()
    
    # Check that all is as expected
    if '#' not in cols: 
        raise IOError('A header beginning with # is expected in file: '+filename)
    elif cols.index('k')!=1: raise IOError('Expecting k as first column in file: '+filename)
    cols.pop(0); cols.pop(0)

    # Get the data
    data = np.loadtxt(filename)
    
    return data[:,0], data[:,1:], cols

def into_1arr(ls,ks,mlps,mlps_lin,nms,nms_lin):

    # Precalculate lengths
    nks = [len(k) for k in ks]
    nktot = sum(nks)
    nncol = None; nlcol=None
    for ml,m in zip(mlps_lin,mlps):
        if nncol is None or nlcol is None:
            nlcol = len(ml[0,:])
            nncol = len(m[0,:])
        elif nlcol+nncol != len(ml[0,:])+len(m[0,:]):
            raise IOError("Multipoles have different number of columns. That won't work!")
        elif nlcol!= len(nms_lin) or nncol!=len(nms):
            raise IOError("Term array sizes differ from lists of columns. That won't work!")

    # Make the empty arrays in advance
    mpls = np.zeros(nktot)
    kvals = np.zeros(nktot)
    terms = np.zeros((nktot,nlcol+nncol))
    
    for nk,l,k,mpl,mpl_lin in zip(nks,ls,ks,mlps,mlps_lin):
        mpls[nk*l/2:nk*(l/2+1)] = l
        kvals[nk*l/2:nk*(l/2+1)] = k
        terms[nk*l/2:nk*(l/2+1),:nlcol] = mpl_lin
        terms[nk*l/2:nk*(l/2+1),nlcol:] = mpl

    return mpls, kvals, terms, nms_lin+nms

# Check the multipoles - not foolproof, but at least something
def check_mlps(mlps):
    for l in mlps:
        if l not in MULTIPOLES.keys(): raise Exception('We can only take multipoles of ' + str(MULTIPOLES.keys()) + '!')
    return True

# See how many times the kvals reverse -> this + 1 should tell you how many multipoles are stored in a vector
def count_mlps(kvals):
    return np.sum(np.diff(kvals)<0)+1

# Extract only one of the repeated k-arrays (return the index - same assumption as above)
def get_k_once(kvals, index=False):
    ends = [len(kvals)]
    ends[:0] = [i+1 for i in np.nonzero(np.diff(kvals)<0)[0]]
    index_once = np.arange(ends[0])
    if index:
        return index_once
    else:
        return kvals[index_once]

def read_zbEFT(basepath, mlps=[0,2,4], eft=FEFT, lin=FLIN):

    terms = [None]*len(mlps)
    terms_lin = [None]*len(mlps)
    kvals_nl = [None]*len(mlps)
    kvals_lin = [None]*len(mlps)

    names = None; nms_l = None
    for i,mlp in enumerate(mlps):

        kvals_lin[i], terms_lin[i], names_lin = read_file(basepath+lin, mlp)
        kvals_nl[i], terms[i], names_nl = read_file(basepath+eft, mlp)

        if names is None or nms_l is None:
            names = names_nl; nms_l = names_lin
        elif np.any(kvals_nl[i]!=kvals_lin[i]) or np.any(nms_l!=names_lin) or np.any(names!=names_nl):
            raise IOError('Mismatch in the k or column names in the files at'+basepath+'!')        

    return into_1arr(mlps,kvals_lin,terms,terms_lin,names_nl,nms_l)  
    