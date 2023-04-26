import numpy as np
import json
import os

workroot = "oop/"
Nsteer = 11
dtheta = 1.5*np.pi/180
steer = dtheta*(np.arange(Nsteer) - (Nsteer-1)/2)
tilts = np.pi*np.linspace(60, 80, 3)/180
Nrot = 2
Nseed = 1
rots = np.pi*np.arange(Nrot)/Nrot

rbuns = [250E-6, 500E-6, 1E-3, 2E-3]
fracs = [0.45, 0.55]

matlegend = {
    'mskin':1,
    'mfat':2,
    'mcon':3,
    'mmusc':4,
    'minter':5
}

default_model = {
        'width':40E-3,
        'tskin':1.5E-3,
        'mskin':matlegend['mskin'],
        'tfat':5E-3,
        'mfat':matlegend['mfat'],
        'tcon':1.5E-3,
        'mcon':matlegend['mcon'],
        'tmusc':32E-3,
        'mmusc':matlegend['mmusc'],
        'minter':matlegend['minter'],
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
                        rotdir = fracdir + f"rot_{180*rot/np.pi:3.03f}/"
                        if not os.path.exists(rotdir): os.makedirs(rotdir)

                        # save image plane selections parameters
                        with open(rotdir+'sel_params.json', 'w') as f:
                            json.dump({'rot':rot, 'seed':seed, 'ppw':16}, f)

                        # save parameter json
                        with open(rotdir + "model_params.json", 'w') as f:
                            json.dump(default_model, f)
                        
                        # copy scripts to directory
                        command = f"cp prepostroutines.py {rotdir}\n"
                        command +=  f"cp genfibers.py {rotdir}\n"
                        command +=  f"cp muscleobjs.py {rotdir}\n"
                        command +=  f"cp ./launch_matmap.sh {rotdir}\n"
                        command +=  f"cp l74.json {rotdir}\n"
                        command += f"sbatch {rotdir}launch_matmap.sh"
                        os.system(command)


