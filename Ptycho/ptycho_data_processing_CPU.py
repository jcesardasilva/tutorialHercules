#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reconstruction script for ptypy version 0.4.0

"""
# standard packages import
import argparse
import os
import re
import subprocess
import sys

import warnings
# to ignore some warnings in ptypy
warnings.filterwarnings('ignore')

# Third packages imports
import numpy as np
import ptypy
from ptypy import io
from ptypy.core import Ptycho
from ptypy import utils as u

"""
Reconstruction script for already prepared far-field ptychography data taken at ID16A beamline - ESRF
"""
# ---------------- Do not edit below this line --------------------------
def eff_distance(sx_pos, focus_pos, z_12):
    """
    Calculates the effective distance using Fresnel Scaling Theorem
    """
    z1 = np.abs(sx_pos - focus_pos)  # the focus to sample distance
    z2 = z_12 - z1  # the sample to detector distance
    eff_dist = z1 * z2 / z_12
    print("The effective sample-to-detector distance is {:0.05f}".format(eff_dist))
    return eff_dist

def record_originalpos():
    coords = []
    print("Saving original positions")
    for pname, pod in P.pods.items():
        # Save real position
        coords.append(np.copy(pod.ob_view.coord))

    np.savetxt("positions_theory.txt", coords)

def str2bool(v):
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

# parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "-f",
    "--files",
    #help="Path to folder containing the files for phase retrieval (Required)",
    help="Prepared .ptyd file for phase retrieval (Required)",
    type=str,
    required=True,
)


parser.add_argument(
    "-m",
    "--probemodes",
    help="Number of probe modes (default = 1)",
    default=1,
    type=int,
)

parser.add_argument(
    "-lo",
    "--loadpreviousobj",
    help="Activate the loading of a previous reconstruction as intial guess for the object (default=False)",
    default=False,
    type=str2bool,
)

parser.add_argument(
    "-lp",
    "--loadpreviousprobe",
    help="Activate the loading of a previous reconstruction as intial guess for the probe ss(default=False)",
    default=False,
    type=str2bool,
)

parser.add_argument(
    "-r",
    "--recons",
    help="File .ptyr containing previous reconstruction for initial guess \
            (only active if --loadprevious is True)",
    default=None,
    type=str,
)

parser.add_argument(
    "-DM",
    "--numiterDM",
    help="Number of iteration of Difference Map algorithm (default=1000)",
    default=1000,  # num_iter_DM,
    type=int,
)

parser.add_argument(
    "-ML",
    "--numiterML",
    help="Number of iteration of Maximum Likelihood optimization (default=200)",
    default=200,  # num_iter_ML,
    type=int,
)

parser.add_argument(
    "-pr",
    "--posref",
    help="Activate the position refinement (default=False)",
    default=False,
    type=str2bool,
)

args = parser.parse_args()

# Change of variables
filename = args.files
posref = args.posref
num_iter_DM = args.numiterDM
num_iter_ML = args.numiterML
num_pmodes = args.probemodes
load_previousobj = ["recon" if args.loadpreviousobj == True else None][0]
load_previousprobe = ["recon" if args.loadpreviousprobe == True else None][0]
resultname0 = args.recons
if load_previousobj or load_previousprobe:
    if resultname0 == None:
        raise ValueError("No file with previous results")
else:
    if resultname0 is not None:
        raise ValueError("Load previous recon file is not activated")
    # resultname0 = None

# check arguments
parser.parse_args()

## file specifics
scanname = os.path.splitext(filename)[0] #re.sub('_\d{4}$','',scanname) #scan name
dfilename = filename #ptyd filename

### calculation of the effective distance (Fresnel Scaling Theorem)
# quick fix to calculate z1
ptydfile=io.h5read(filename,'info')['info']
sx = ptydfile['sx']
focus_position = ptydfile['sx0']
z12 = ptydfile['distance']
dist_parallel = eff_distance(sx, focus_position, z12)
aperture_size = 120e-6 # in meters
focal_length = 0.1 # in meters
print('The effective sample-to-detector distance is {:0.05f}'.format(dist_parallel))

##### Ptypy parameter tree #####

## General parameter container
p = u.Param()
p.verbose_level = 3						## Verbosity level
p.data_type = "single"					## Reconstruction floating number precision
p.run = scanname						## Reconstruction identifier
p.ipython_kernel = False				## New feature?
#p.frames_per_block = 200                ## GPU only: It should be the same as the number of frames


## Global parameters for I/O
p.io = u.Param()
#p.io.home = "./" 						## Base directory for all I/O
p.io.rfile = "recons/%(run)s/%(run)s_%(engine)s.ptyr" ## Reconstruction file name (or format string)
#p.io.rfile = "recons/%(run)s/%(run)s_%(engine)s_%(iterations)04d.ptyr" ## Reconstruction file name (or format string)

## Auto-save options
p.io.autosave = u.Param()
p.io.autosave.active = True 			## Activation switch to save dump files
p.io.autosave.interval = -1			## Auto-save interval for dump files
#p.io.autosave.rfile = "dumps/%(run)s/%(run)s_%(engine)s_%(iterations)04d.ptyr" ## (10) Auto-save file name (or format string)

## Server / Client parameters
p.io.interaction = u.Param()
p.io.interaction.server = u.Param()
p.io.interaction.server.active = True                  ## Activation switch
p.io.interaction.server.address = "tcp://127.0.0.1"    ## The address the server is listening to.
p.io.interaction.server.port = 5560                    ## The port the server is listening to.
p.io.interaction.server.connections = 10               ## Number of concurrent connections on the server

## Plotting client parameters
p.io.autoplot = u.Param()
#p.io.autoplot.imfile = "plots/%(run)s/%(run)s_%(engine)s_%(iterations)04d.png" ## (17) Plot images file name (or format string)
p.io.autoplot.interval = 1                      ## Number of iterations between plot updates
p.io.autoplot.threaded = False                  ## Live plotting switch
p.io.autoplot.layout = 'nearfield'              ## Options for default plotter or template name
p.io.autoplot.dump = False                      ## Switch to dump plots as image files
p.io.autoplot.make_movie = False                ## Produce reconstruction movie after the reconstruction.

## Scan parameters
p.scans = u.Param() ## param container for instances of scan parameters
p.scans.ID16A = u.Param()
p.scans.ID16A.name = "Full" #"BlockFull"  or "Full"
p.scans.ID16A.propagation = "farfield"  ## Propagation type: "nearfield" or "farfield"

# Data parameters
#p.scans.ID16A.if_conflict_use_meta = True ## Give priority to metadata relative to input parameters
p.scans.ID16A.data= u.Param() ## Data preparation parameters
p.scans.ID16A.data.name = "PtydScan" ## Name of the PtyScan subclass to use
p.scans.ID16A.data.source = 'file' ## Describes where to get the data from.
p.scans.ID16A.data.dfile = dfilename ## Prepared data file path

# Illumination parameters
p.scans.ID16A.illumination = u.Param() ## illumination model parameters (probe)
p.scans.ID16A.illumination.aperture = u.Param() ## Beam aperture parameters
p.scans.ID16A.illumination.aperture.form = "rect" ## One of None, 'rect' or 'circ'
p.scans.ID16A.illumination.aperture.size = aperture_size ## Aperture width or diameter
p.scans.ID16A.illumination.propagation=u.Param() ## Parameters for propagation after aperture plane
p.scans.ID16A.illumination.propagation.parallel = dist_parallel ## Parallel propagation distance
p.scans.ID16A.illumination.propagation.focussed = focal_length ## propagation distance from aperture to focus

## mixed states - orthogonal modes
# this diversity will rescale the probe modes in order to taken into account the difference in photons
p.scans.ID16A.illumination.diversity = u.Param()
p.scans.ID16A.illumination.diversity.power = 0.1
p.scans.ID16A.illumination.diversity.noise = (1,2)
## coherence parameters
p.scans.ID16A.coherence = u.Param()
p.scans.ID16A.coherence.num_probe_modes = num_pmodes ## Number of probe modes
p.scans.ID16A.coherence.num_object_modes = 1 ## Number of object modes
#p.scans.ID16A.coherence.enegies = [1.0] # main energies
#p.scans.ID16A.coherence.spectrum = [1.0] # List of energies / wavelength relative to mean energy / wavelength
#p.scans.ID16A.coherence.probe_dispersion = None # If True, the same probe is used for all energies
#p.scans.ID16A.coherence.object_dispersion = None # If True, the same object is used for all energies


### Parameters to load from previous reconstruction
# if load model from previous reconstruction, use recon parameter
p.scans.ID16A.illumination.model = None ## None, 'recon', 'stxm'
if p.scans.ID16A.illumination.model == 'recon':
    p.scans.ID16A.illumination.recon = u.Param() ## parameters to load from previous reconstruction
    p.scans.ID16A.illumination.recon.rfile = resultname0 #"\*.ptyr" ## path to a .ptyr compatible file
    p.scans.ID16A.illumination.diversity = None
    p.scans.ID16A.illumination.propagation = u.Param()
    p.scans.ID16A.illumination.propagation.parallel = None
    #p.scans.ID16A.illumination.photons = None # Total number of photons
else:
    p.scans.ID16A.illumination.recon = None

## Object initiation parameters
p.scans.ID16A.sample= u.Param() ## initial object modelization parameters

### Parameters to load from previous reconstruction
# if load model from previous reconstruction, use recon parameter
p.scans.ID16A.sample.model = None ## None, 'recon', 'stxm'
if p.scans.ID16A.sample.model == 'recon':
    p.scans.ID16A.sample.recon = u.Param() ## parameters to load from previous reconstruction
    p.scans.ID16A.sample.recon.rfile = resultname0 #"\*.ptyr" ## path to a .ptyr compatible file
    p.scans.ID16A.sample.diversity = None
    p.scans.ID16A.sample.process = None # Model processing parameter. It must be None when loading previous data
else:
    p.scans.ID16A.sample.recon = None
    p.scans.ID16A.sample.fill = 1.0 + 1j * 0.0 # Default fill value

## Container for instances of "engines"  parameters"
p.engines = u.Param()

## Parameters for Difference map engine
## First engine = DM
p.engines.engine00 = u.Param()
p.engines.engine00.name = "DM"             ## Engine name
p.engines.engine00.numiter = num_iter_DM          ## number of iterations of DM
p.engines.engine00.numiter_contiguous = 1        ## Number of iteractions without interruption
p.engines.engine00.probe_support = None           ## Fraction of valid probe area (circular in probe frame)
p.engines.engine00.alpha = 1                      ## Difference map parameter
p.engines.engine00.probe_update_start = 2         ## Number of iteractions before probe update starts
p.engines.engine00.update_object_first = True     ## if False update object before probe
p.engines.engine00.overlap_converge_factor = 0.02 ## Threshold for interruption of the inner overlap loop
p.engines.engine00.overlap_max_iterations = 100   ## Maximum of iterations for the overlap constraint inner loop
p.engines.engine00.probe_inertia = 0.01           ## Weight of the current probe estimate in the update
p.engines.engine00.object_inertia = 0.1           ## Weight of the current object in the update
p.engines.engine00.fourier_relax_factor = 0.02    ## if rms error of model vs diffraction data is smaller than this fraction, Fourier constraint is met (set this value higher for noisy data)
p.engines.engine00.obj_smooth_std = 5            ## Gaussian smoothing (in pixel ) of the currect object prior to update
p.engines.engine00.clip_object = (0.,1.2)    ## Clip object amplitude into this interval

## position refinement
if posref:
    p.engines.engine00.position_refinement = u.Param()
    p.engines.engine00.position_refinement.start = 100  # 50
    p.engines.engine00.position_refinement.stop = 300  # 300
    p.engines.engine00.position_refinement.interval = 2
    p.engines.engine00.position_refinement.nshifts = 8
    p.engines.engine00.position_refinement.amplitude = 50e-9  # 6e-7
    p.engines.engine00.position_refinement.max_shift = 200e-9  # 6e-7
    p.engines.engine00.position_refinement.record = True
else:
    p.engines.engine00.position_refinement = False

## Second engine = ML
p.engines.engine01 = u.Param()                   ## Default second engine
p.engines.engine01.name = 'ML'            ## parameters for Maximum Likehood engine
p.engines.engine01.numiter = num_iter_ML         ## total number of iteractions
p.engines.engine01.numiter_contiguous = 1        ## number of iteration without interruption
p.engines.engine01.ML_type = "gaussian"             ## Likelihood model. One of 'gaussian','poison' or 'euclid'
p.engines.engine01.probe_support = None          ## Fraction of valid probe are (circular) in probe frame
p.engines.engine01.floating_intensities = True   ## if True, allow for adaptative rescalling of one diffraction pattern intensities (to correct for incident beam intensity fluctuations).
p.engines.engine01.intensity_renormalization = 1 ## A rescaling of the intensity so they can be interpreted as Poisson counts
p.engines.engine01.reg_del2 = True			   ## (142) Whether to use a Gaussian prior (smoothing) regularizer.
p.engines.engine01.reg_del2_amplitude = 1e-2#0.0001 		## (143) Amplitude of the Gaussian prior if used.
p.engines.engine01.smooth_gradient = 1          ## Smoothing preconditioner. If 0, not used, if >0 gaussian filter, if <0 Hann window.
p.engines.engine01.scale_precond = True          ## Whether to use the object/probe scaling preconditioner
p.engines.engine01.probe_update_start = 2        ## Number of iterations before probe update starts
#p.engine = u.Param()
#p.engine.ML = u.Param()
#p.engine.ML.floating_intensities = True

P = Ptycho(p,level=5)
if posref: record_originalpos()
P.print_stats()
P.finalize()
