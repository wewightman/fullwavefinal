import numpy as np
import genfibers as gf
import muscleobjs as mo
import bmfrmtrig as bft
import json

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def prerout(seed: int, probefile="l74.json", ppw: int=16, rfat = 1.2E-3, sigfat = 300E-6, rfib = 300E-6, sigfib = 100E-6, dz = 1000E-6, dzsig = 100E-6, nfib=120, nfat=1200):
    """Generate parameter files needed for fullwave simulation launched by matlab"""
    logger.info("Loading probe file")
    # Load parameter file
    with open(probefile, 'r') as f:
        probe = json.load(f)

    c0 = 1540
    f0 = probe['impulseResponse']['f0']
    lam = c0/f0
    Nx = int(np.ceil(ppw*(probe['pitch']*probe['noElements']-probe['kerf'])/lam))
    Nz = int(np.ceil(ppw*40E-3/lam))


    logger.info("Generating field grid")
    # Generate field grid
    xmin = 0
    xmax = Nx*lam/ppw
    x = np.linspace(xmin, xmax, Nx)

    zmin = 0
    zmax = Nz*lam/ppw
    z = np.linspace(zmin, zmax, Nz)

    y = 0

    X, Y, Z = np.meshgrid(x, y, z)
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    logger.info("Select fibers")
    # Generate fibers
    N = nfib
    rng = np.random.default_rng(seed)
    dz_c = rng.normal(dz, dzsig, N)
    dz_c[dz_c<dzsig] = dzsig
    z_c = np.cumsum(dz_c) + 5E-3
    x_c = rng.uniform(xmax/2-5E-3, xmax/2+5E-3, N)
    y_c = rng.uniform(-dz/3, dz/3, N)
    r0 = rng.normal(rfib, sigfib, N)
    r0[r0<sigfib] = sigfib
    r1 = rng.uniform(50E-6, 150E-6, N)

    fmap = 2*gf.gengeneral(points, gf.sphere, [30E-3, r0[0]+2*r1[0], r0[0]+2*r1[0]], [0, 0, 0], [x_c[0], 0, z_c[0]])
    for ind in range(N-1):
        fmap -= gf.gengeneral(points, gf.sphere, [30E-3, r0[ind+1], r0[ind+1]], [0, 0, 0], [x_c[ind+1], y_c[ind+1], z_c[ind+1]])
        fmap += 2*gf.gengeneral(points, gf.sphere, [30E-3, r0[ind+1]+2*r1[ind+1], r0[ind+1]+2*r1[ind+1]], [0, 0, 0], [x_c[ind+1], y_c[ind+1], z_c[ind+1]])

    fmap = fmap

    fmap[(fmap==3)] = 2
    fmap[fmap>3] = 0

    logger.info("Gennerate cellulite map")
    # generate cellulite map
    N = nfat
    x = rng.uniform(xmin-5E-3, xmax+5E-3, N)
    z = rng.uniform(zmin-5E-3, zmax+5E-3, N)
    r = rng.normal(rfat, sigfat, 3*N).reshape((N, 3))
    r[r <= sigfat] = sigfat
    theta = rng.uniform(0, 2*np.pi, 3*N).reshape((N, 3))

    fatmask = np.zeros(points.shape[1])
    for ind in range(N):
        fatmask += gf.gengeneral(points, gf.gaussian, r[ind,:], theta[ind,:], [x[ind], 0, z[ind]])

    fatmask = fatmask

    x = rng.uniform(xmin-5E-3, xmax+5E-3, N)
    z = rng.uniform(zmin-5E-3, zmax+5E-3, N)
    r = rng.normal(rfat, sigfat, 3*N).reshape((N, 3))
    r[r <= sigfat] = sigfat
    theta = rng.uniform(0, 2*np.pi, 3*N).reshape((N, 3))

    fatmask1 = np.zeros(points.shape[1])
    for ind in range(N):
        fatmask1 += gf.gengeneral(points, gf.gaussian, r[ind,:], theta[ind,:], [x[ind], 0, z[ind]])

    fatmask1 = fatmask1


    fascia = (((fatmask>np.percentile(fatmask, 47.5)) & (fatmask < np.percentile(fatmask, 52.5))) 
              | ((fatmask1>np.percentile(fatmask1, 47.5)) & (fatmask1 < np.percentile(fatmask1, 52.5))))

    fmap[fmap==0] = 2*fascia[fmap==0]
    
    logger.info(f"Save Map nx:{Nx}, nz:{Nz}")
    fmap.astype(np.int8).tofile("map.bin")

    params = {
        'c0':c0,
        'f0':f0,
        'lam':lam,
        'Nx':Nx,
        'Nz':Nz,
        'ppw':ppw,
        'probe':probe
    }

    with open('params.json', 'w') as f:
        json.dump(params, f)

def pre_genmatmap(sampparams='model_params.json', simparams="field_params.json", probefile='probe.json', matkey="matkey.json"):
    logger.info("Loading probe file...")
    with open(probefile, 'r') as f:
        probe = json.load(f)
    
    logger.info("Loading field file...")
    with open(simparams, 'r') as f:
        fieldparams = json.load(f)

    logger.info("Loading sample parameter file...")
    with open(sampparams, 'r') as f:
        modelparams = json.load(f)

    logger.info("Loading meterial key file...")
    with open(matkey, 'r') as f:
        matkeys = json.load(f)

    logger.info("Generating sample set")
    set = mo.SampSet(**modelparams)

    c0 = fieldparams['c0']
    f0 = probe["centerFrequency"]
    ppw = fieldparams['ppw']
    rot = fieldparams['rot']
    lam = c0/f0
    Nx_half = int(np.ceil(ppw*(probe['pitch']*(probe['noElements']-1)/2)/lam))
    Nx = 2*Nx_half+1
    Nz = int(np.ceil(ppw*40E-3/lam))

    logger.info("Generating field grid")
    xmin = -Nx_half*lam/ppw
    xmax = Nx_half*lam/ppw
    x = np.linspace(xmin, xmax, Nx)

    zmin = 0
    zmax = Nz*lam/ppw
    z = np.linspace(zmin, zmax, Nz)

    y = 0

    X, Y, Z = np.meshgrid(x, y, z)
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    logger.info("Rotating field grid")
    T_rot = np.array(
        [[np.cos(rot), -np.sin(rot), 0],
         [np.sin(rot), np.cos(rot), 0],
         [0, 0, 1]]
    )
    points = T_rot @ points

    logger.info("Sampling material model")
    model = set(points)

    model[model == 0] = matkeys['inter']

    logger.info(f"Save Map nx:{Nx}, nz:{Nz}")
    model.astype(np.int8).tofile("map.bin")

    params = {
        'c0':c0,
        'f0':f0,
        'lam':lam,
        'Nx':Nx,
        'Nz':Nz,
        'ppw':ppw,
        'probe':probe
    }

    with open('matparams.json', 'w') as f:
        json.dump(params, f)

def postbeamform(filename='channels', savename = None, filepath='', c = 1540, extent = [-20E-3, 20E-3, 1.5E-3, 40E-3], dpx = [150E-6, 40E-6], fnum=1.5):
    with open(filepath + filename + ".json", 'r') as f:
        params = json.load(f)
    
    if savename is None:
        savename = filename + "_bf"
    
    # load binary file
    import ctypes as ct
    data = np.fromfile(filepath + filename + ".bin", dtype=np.float64).astype(ct.c_float).ctypes.data_as(ct.POINTER(ct.c_float))

    nele = params['nele']
    dele = params['dele']
    Ts = params['dT']
    theta = params['theta']
    fs = 1/Ts

    # Generate element grid
    xele = dele*(np.arange(nele) - (nele-1)/2)
    yele = 0
    zele = 0
    Xele, Yele, Zele = np.meshgrid(xele, yele, zele)
    eles = np.array([Xele.flatten(), Yele.flatten(), Zele.flatten()]).T

    # Recon parameters
    dz = dpx[1]
    zmin = extent[2]
    zmax = extent[3]
    dx = dpx[0]
    xmin = extent[0]
    xmax = extent[1]

    # Generate recon grid
    xgrid = np.arange(xmin, xmax, dx)
    ygrid = 0 #np.arange(ymin, ymax, dy)
    zgrid = np.arange(zmin, zmax, dz)
    Xgrid, Ygrid, Zgrid = np.meshgrid(xgrid, ygrid, zgrid)
    field = np.array([Xgrid.flatten(), Ygrid.flatten(), Zgrid.flatten()]).T

    tref = np.sin(theta)*xele/c
    tref = tref - np.min(tref)

    alphas = np.repeat(theta, nele)

    # initialize beamformer
    funcparams = {}
    funcparams['c'] = c
    funcparams['fnum'] = fnum
    funcparams['points'] = field
    funcparams['trefs'] = tref
    funcparams['refs'] = eles
    funcparams['alphas'] = alphas
    funcparams['nsamp'] = params['nT']
    funcparams['fs'] = fs
    funcparams['tstart'] = 0
    funcparams['parallel'] = True
    funcparams['ncores'] = 8

    bmfrm = bft.PWBeamformer(**funcparams)
    
    formed = bmfrm(data)

    formed.tofile(filepath + savename + ".bin")
