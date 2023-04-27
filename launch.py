import numpy as np
from string import Template
import json
import os

workroot = "/work/wew12/fullwave/"
repodir = "/hpc/group/ultrasound/wew12/repos/fullwavefinal/"
Nsteer = 11
dtheta = 1.5*np.pi/180
steers = dtheta*(np.arange(Nsteer) - (Nsteer-1)/2)
tilts = np.pi*np.linspace(60, 80, 3)/180
Nrot = 2
Nseed = 1
rots = np.pi*np.arange(Nrot)/Nrot

rbuns = [250E-6, 500E-6, 1E-3, 2E-3]
fracs = [0.45, 0.55]

# default material properties
matlegend = {
    'skin'  :1,
    'fat'   :2,
    'con'   :3,
    'musc'  :4,
    'inter' :5,
    'params': {
                # c,    rho,  alpha, b/a, sigma
        'skin'  :[1615, 1090, 0.35,  7.9, 0.030],
        'fat'   :[1465,  985, 0.40,  8.5, 0.010],
        'musc'  :[1580, 1050, 0.74,  6.6, 0.001],
        'con'   :[1613, 1120, 1.57,  0.0, 0.015],
        'inter' :[1543, 1027, 0.069, 0.0, 0.001]
    }
}

# default field selection parameters
default_field = {
    'rot':0, 
    'seed':0, 
    'ppw':16,
    'c0':1540
}

# default model parameters
default_model = {
        'width':50E-3,
        'tskin':1.5E-3,
        'mskin':matlegend['skin'],
        'tfat':5E-3,
        'mfat':matlegend['fat'],
        'tcon':1.5E-3,
        'mcon':matlegend['con'],
        'tmusc':32E-3,
        'mmusc':matlegend['musc'],
        'minter':matlegend['inter'],
        'theta':tilts[0],
        'rbun':0.5E-3,
        'sigbun':0.1E-3,
        'rfib':200E-6,
        'rfrac':0.45E-3,
        'sigfrac':0.05E-3,
        'seed':0
    }

def launch_matmaps():
    # iterate through seed, fiber tilt angle, bundle radius, level of overlap, and rotation angle
    for seed in range(Nseed):
        seeddir = workroot + f"seed_{seed:d}/"
        default_model['seed'] = seed
        default_field['seed'] = seed
        for tilt in tilts:
            tiltdir = seeddir + f"tilt_{180*tilt/np.pi:3.03f}/"
            default_model['theta'] = tilt
            for rbun in rbuns:
                bundir = tiltdir + f"rfibmm_{1E3*rbun}/"
                default_model['rbun'] = rbun
                for frac in fracs:
                    fracdir = bundir + f"frac_{100*frac:0.01f}/"
                    default_model['rfrac'] = frac
                    for rot in rots:
                        default_field['rot'] = rot
                        rotdir = fracdir + f"rot_{180*rot/np.pi:3.03f}/"
                        if not os.path.exists(rotdir): os.makedirs(rotdir)
                        os.chdir(rotdir)                        

                        # save image plane selections parameters and recon parameters
                        with open(rotdir+'field_params.json', 'w') as f:
                            json.dump(default_field, f)

                        # save material properties
                        with open(rotdir+'matkey.json', 'w') as f:
                            json.dump(matlegend, f)

                        # save parameter json
                        with open(rotdir + "model_params.json", 'w') as f:
                            json.dump(default_model, f)
                        
                        # copy scripts to directory
                        command = f"cp {repodir}prepostroutines.py {rotdir}\n"
                        command +=  f"cp {repodir}genfibers.py {rotdir}\n"
                        command +=  f"cp {repodir}muscleobjs.py {rotdir}\n"
                        command +=  f"cp {repodir}l74.json {rotdir}probe.json\n"
                        command +=  f"cp {repodir}shellscripts/gen_matmap.sh {rotdir}\n"
                        command += f"sbatch {rotdir}gen_matmap.sh"
                        os.system(command)
                        
def launch_simmaps():
    # iterate through seed, fiber tilt angle, bundle radius, level of overlap, and rotation angle
    for seed in range(Nseed):
        seeddir = workroot + f"seed_{seed:d}/"
        for tilt in tilts:
            tiltdir = seeddir + f"tilt_{180*tilt/np.pi:3.03f}/"
            for rbun in rbuns:
                bundir = tiltdir + f"rfibmm_{1E3*rbun}/"
                for frac in fracs:
                    fracdir = bundir + f"frac_{100*frac:0.01f}/"
                    for rot in rots:
                        rotdir = fracdir + f"rot_{180*rot/np.pi:3.03f}/"
                        for steer in steers:
                            steerdir = rotdir + f"steer_{180*steer/np.pi:3.03f}/"
                            if not os.path.exists(steerdir): os.makedirs(steerdir)
                            os.chdir(steerdir)
                            
                            # Save imaging parameters
                            with open("impar.json", 'w') as f:
                                json.dump({"theta":steer}, f)
                            
                            # copy scripts to directory
                            command = f"cp {rotdir}matparams.json {steerdir}\n"
                            command +=  f"cp {rotdir}matparams.json {steerdir}\n"
                            command +=  f"cp {rotdir}map.bin {steerdir}\n"
                            command +=  f"cp {rotdir}matkey.json {steerdir}\n"
                            command +=  f"cp {repodir}readmap.m {steerdir}\n"
                            command +=  f"cp {repodir}readjson.m {steerdir}\n"
                            command +=  f"cp {repodir}runfullwave.m {steerdir}\n"
                            command +=  f"cp {repodir}shellscripts/sim_matmap.sh {steerdir}\n"
                            command += f"sbatch {steerdir}sim_matmap.sh"
                            os.system(command)

def launch_beamform():
    # iterate through seed, fiber tilt angle, bundle radius, level of overlap, and rotation angle
    for seed in range(Nseed):
        seeddir = workroot + f"seed_{seed:d}/"
        for tilt in tilts:
            tiltdir = seeddir + f"tilt_{180*tilt/np.pi:3.03f}/"
            for rbun in rbuns:
                bundir = tiltdir + f"rfibmm_{1E3*rbun}/"
                for frac in fracs:
                    fracdir = bundir + f"frac_{100*frac:0.01f}/"
                    for rot in rots:
                        rotdir = fracdir + f"rot_{180*rot/np.pi:3.03f}/"
                        for steer in steers:
                            steerdir = rotdir + f"steer_{180*steer/np.pi:3.03f}/"
                            if not os.path.exists(steerdir):
                                print(f"{steerdir} doesnt exist... continuing...")
                                continue
                            os.chdir(steerdir)
                            
                            command =  f"cp {repodir}shellscripts/sim_beamform.sh {steerdir}\n"
                            command += f"sbatch {steerdir}sim_beamform.sh"
                            os.system(command)


