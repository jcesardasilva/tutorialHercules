# coding: utf-8

# standard packages
import math

# third party packages
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np
import numpy.matlib as npm

def sphere(array_size,radius=1, showFig = False):
    theta, phi = np.linspace(0, 2 * np.pi, array_size+1), np.linspace(0, np.pi, array_size+1)
    THETA, PHI = np.meshgrid(theta, phi)
    X = radius * np.sin(PHI) * np.cos(THETA)
    Y = radius * np.sin(PHI) * np.sin(THETA)
    Z = radius * np.cos(PHI)
    if showFig:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='3d')
        plot = ax.plot_surface(X, Y, Z, color='b')
        ax.set_axis_off()
        ax.set_facecolor((0,0,0))
        plt.show(block=False)
    return X,Y,Z

def sphere2(array_size,radius=1, showFig = False):
    u, v = np.linspace(0, 2 * np.pi, array_size+1), np.linspace(0, np.pi, array_size+1)
    X = radius * np.outer(np.cos(u),np.sin(v))
    Y = radius * np.outer(np.sin(u),np.sin(v))
    Z = radius * np.outer(np.ones(np.size(u)),np.cos(v))
    if showFig:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='3d')
        plot = ax.plot_surface(X, Y, Z, color='b')
        ax.set_facecolor((0,0,0))
        plt.show(block=False)
    return X,Y,Z

def drawmol(mol, backgroundColor = 'k'):
    """
    # DRAWMOL(MOL)
    # Render simple molecules. Large molecules will require
    # considerable computing power.
    # Colors:
    # White: hydrogen
    # Black: carbon
    # Red: oxygen
    # Blue: nitrogen
    # Yellow: sulphur
    # Grey: other.

    #Based on DRAWPDB by Joe Hicklin, September 2001.
    """
    # creating size of the spheres
    if len(mol['x'])>50:
        #X,Y,Z = sphere2(10)
        sphere_size = 10
    else:
        #X,Y,Z = sphere2(20)
        sphere_size = 20

    # Set view parameters for all subplots.
    azimuth = 90
    altitude = 90

    # preparing the canvas
    plt.close('all')
    plt.ion()
    fig = plt.figure(1,figsize =(8,8))
    ax = fig.add_subplot(1,1,1,projection='3d')
    ax.view_init(altitude,azimuth)

    # loop over the mols
    for ii in range(len(mol['x'])):
        if mol['element'][ii] == 'H':
            #colorBead='w'
            colorBead = (.9,1,.9)
            r = 0.5
        elif mol['element'][ii] == 'C':
            #colorBead='k'
            colorBead = (0.25,0.25,0.25)
            r = 0.85
        elif mol['element'][ii] == 'O':
            #colorBead='r'
            colorBead = (1,0,0)
            r = 0.95
        elif mol['element'][ii] == 'N':
            #colorBead='b'
            colorBead = (0,0,1)
            r = 0.9
        elif mol['element'][ii] == 'S':
            #colorBead='y'
            colorBead = (1,1,0)
            r = 1.0
        else:
            colorBead=(0.6,0.6,0.6)
            r = 0.9
        # creating the sphere with radius value
        X,Y,Z = sphere2(sphere_size,radius=r)

        # preparing the light
        rgb = np.ones((Z.shape[0],Z.shape[1],3))
        light = LightSource(90, 45)

        # give the color to the shaded surface
        illuminatedSurface = light.shade_rgb(rgb*colorBead,Z)

        # adding beads with different color and radius to plots
        ax.plot_surface(mol['x'][ii]+X,mol['y'][ii]+Y,mol['z'][ii]+Z,
            rstride=1, cstride=1, linewidth=0,
            color=colorBead,
            #facecolors =illuminatedSurface,
            antialiased=False)
    if backgroundColor == 'black':
        ax.set_facecolor((0,0,0))
    elif backgroundColor == 'white':
        ax.set_facecolor((1,1,1,))
    elif backgroundColor == 'grey':
        ax.set_facecolor((0.5,0.5,0.5))
    else:
        print('Background color not implement yet. Using black')
        ax.set_facecolor((0,0,0))
    xmin, xmax, ymin, ymax = ax.axis()
    zmin = 0.8*np.min((xmin, ymin))
    zmax = 0.8*np.max((xmax,ymax))
    ax.set_xlabel('x [Angstrom]',fontsize = 14)
    ax.set_ylabel('y [Angstrom]',fontsize = 14)
    plt.tight_layout()
    ax.auto_scale_xyz([xmin,xmax],[ymin,ymax],[zmin,zmax])
    #ax.set_axis_off()
    plt.show(block=False)

def pdbreadatom(filename):
    """
    #Reads a PDB-file and stores selected information from the ATOM field.

    # COLUMNS        DATA TYPE       FIELD         DEFINITION
    # ---------------------------------------------------------------------------------
    #  1 -  6        Record name     "ATOM  "
    #  7 - 11        Integer         serial        Atom serial number.
    # 13 - 16        Atom            name          Atom name.
    # 17             Character       altLoc        Alternate location indicator.
    # 18 - 20        Residue name    resName       Residue name.
    # 22             Character       chainID       Chain identifier.
    # 23 - 26        Integer         resSeq        Residue sequence number.
    # 27             AChar           iCode         Code for insertion of residues.
    # 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
    #                                              Angstroms.
    # 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
    #                                              Angstroms.
    # 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
    #                                              Angstroms.
    # 55 - 60        Real(6.2)       occupancy     Occupancy.
    # 61 - 66        Real(6.2)       tempFactor    Temperature factor.
    # 73 - 76        LString(4)      segID         Segment identifier, left-justified.
    # 77 - 78        LString(2)      element       Element symbol*, right-justified.
    # 79 - 80        LString(2)      charge        Charge on the atom.

    # Details for the atom name:
    # 13 - 14    Chemical symbol* - right justified, except for hydrogen atoms
    # 15         Remoteness indicator (alphabetic)
    # 16         Branch designator (numeric)

    #Element and chemical symbol both refer to the corresponding entry in the
    #periodic table.
    """
    atom = {}
    with open(filename,'r') as fid:
        # ~ lines = [line.rstrip('\n') for line in fid if line.startswith('ATOM ') or line.startswith('HETATM ')]
        lines = [line.rstrip('\n') for line in fid if line.startswith('ATOM ')]
    for ii in range(len(lines)):
        atom[ii] = {}
        atom[ii]['name'] = lines[ii][13:16].rstrip()#lines[ii].split()[1]
        atom[ii]['x'] = eval(lines[ii][31:38].rstrip())#eval(lines[ii].split()[2])
        atom[ii]['y'] = eval(lines[ii][39:46].rstrip())#eval(lines[ii].split()[3])
        atom[ii]['z'] = eval(lines[ii][47:54].rstrip())#eval(lines[ii].split()[4])
        atom[ii]['tempFactor'] = eval(lines[ii][61:66].rstrip())#eval(lines[ii].split()[5])
        atom[ii]['element'] = lines[ii][77:78].rstrip()#lines[ii].split()[6]

    #If the file doesn't have an element entry (which it should), use the first characters of
    #the atom name (which should be the same symbol, right aligned). This is not
    #a very failsafe check.
    if atom[0]['element']=='' or atom[0]['element']!='0':
        print('Creating dictionary of elements from the name of the atoms')
        atom['element'] = np.array([atom[ii]['name'][0] for ii in range(len(lines))])
    else:
        atom['element'] = np.array([atom[ii]['element'] for ii in range(len(lines))])

    #Duplicate some fields for compatibility with some functions...
    atom['IDP'] = np.array([atom[ii]['tempFactor'] for ii in range(len(lines))])
    atom['x'] = np.array([atom[ii]['x'] for ii in range(len(lines))])
    atom['y'] = np.array([atom[ii]['y'] for ii in range(len(lines))])
    atom['z'] = np.array([atom[ii]['z'] for ii in range(len(lines))])
    atom['crd'] = np.vstack((atom['x'],atom['y'],atom['z']))
    return atom

def moltrans(mol, H, K, L):
    sizeH = H.shape
    nrCrds = np.prod(sizeH)
    nrAtoms = len(mol['x'])
    maxSzMb = 32

    # don't create matrices larger than a given size
    lim = np.ceil(maxSzMb*2**20/(8*nrAtoms))

    # initializing arrays
    F = np.zeros((1,nrCrds)).astype(np.complex)

    # raveling
    H = H.ravel()
    K = K.ravel()
    L = L.ravel()

    ii=0
    while ii<nrCrds:
        idx = np.arange(ii,np.min((ii+lim-1,nrCrds))).astype(np.int)
        ii += lim
        f = scatteringfactor(mol['element'],H[idx], K[idx],L[idx])
        f = f*debyewallerfactor(mol['IDP'],H[idx], K[idx],L[idx])
        F[:,idx] = structurefactor(f,mol['crd'],H[idx], K[idx],L[idx])

    return np.flipud(np.fliplr(F.T.reshape(sizeH)))

def readsf(sffile):
    """
    Coefficients for analytical approximation to scattering factors

     f(\sin(\theta)/\lambda) = \sum_{i=1}^4 a_i
                \exp(-b_i*(\sin(\theta)/\lambda)^2) + c
    order of the parameters in the files:
        name a1 b1 a2 b2 a3 b3 a4 b4 c
    """
    with open(sffile,'r') as fid:
        lines = [line.rstrip('\n') for line in fid if not line.startswith('#')]
    sf = {}
    for ii in range(len(lines)):
        sf[ii] = {}
        sf[ii]['a'] = np.zeros(4)
        sf[ii]['b'] = np.zeros(4)
        sf[ii]['label'] = lines[ii].split()[0]
        sf[ii]['a'][0] = eval(lines[ii].split()[1])
        sf[ii]['b'][0] = eval(lines[ii].split()[2])
        sf[ii]['a'][1] = eval(lines[ii].split()[3])
        sf[ii]['b'][1] = eval(lines[ii].split()[4])
        sf[ii]['a'][2] = eval(lines[ii].split()[5])
        sf[ii]['b'][2] = eval(lines[ii].split()[6])
        sf[ii]['a'][3] = eval(lines[ii].split()[7])
        sf[ii]['b'][3] = eval(lines[ii].split()[8])
        sf[ii]['c'] = eval(lines[ii].split()[9])
    sf['label'] = np.array([sf[ii]['label'] for ii in range(len(sf))])
    return sf

def scatteringfactor(molid,H,K,L):
    """
     f = scatteringfactor(id,HKL)
     f = scatteringfactor(id,H,K,L)
    id is a string with the chemical symbols of N atoms. (Nx1)
    HKL are n scattering vectors. (nx3)
    Output is a matrix with scattering factors. (Nxn)

    TODO:
    The sf structure should perhaps be input to the function.
    """
    # read scattering factors from file
    sf = readsf('atomsf.lib')
    HKL = np.vstack((H,K,L))
    atomTypes = np.unique(molid)#['element'])

    # Pre-allocate matrix
    f = np.zeros((molid.shape[0], HKL.shape[1])) #Pre-allocate matrix.
    stols = 1/4*np.sum(HKL**2,axis=0) #Pre-calculate (sin(theta)/lambda)^2.

    for ii in range(len(atomTypes)):
        sfnr = np.where(sf['label'] == atomTypes[ii])[0][0]
        idx = np.where(molid == atomTypes[ii])[0]
        a = sf[sfnr]['a']
        b = sf[sfnr]['b']
        c = sf[sfnr]['c']
        fType = c + np.exp(-np.matmul(stols[:,np.newaxis],b[np.newaxis,:])).dot(a[:,np.newaxis])
        f[idx,:] = npm.repmat(fType.T, idx[:,np.newaxis].shape[0],idx[:,np.newaxis].shape[1])

    #As a "sort of" correction for high resolution, which may yield negative
    #values for heavier elements, all negative values are set to eps.

    f[np.where(f<0)] = np.spacing(1)

    return f#.reshape(HKL.shape)

def debyewallerfactor(IDP,H,K,L):
    """
    function  T = debyewallerfactor(B,H)
    IDP are isotropic displacement parameters for N atoms. (Nx1)
    HKL are n scattering vectors. (Nx3).
    Output is a matrix with the atomic isotropic Gaussian Debye-Waller factors. (Nx1)
    """
    HKL = np.vstack((H,K,L))
    stols = np.sum(HKL**2,axis=0)/4 #Precompute (sin(theta)/lambda)^2. (Nx1)
    return np.exp(-IDP[:,np.newaxis].dot(stols[:,np.newaxis].T))

def structurefactor(f,R,H,K,L):
    """
    function F = structureFactor(f,R,HKL)
    f contains scattering factors for N atoms and n scattering vectors. (Nxn)
    R are coordinates for N atoms. (Nx3)
    HKL are the coordinates of n scattering vectors. (nx3)
    Output is a matrix with structure factors. (1xN)
    """
    HKL = np.vstack((H,K,L))
    phase = (-2*np.pi*R.T).dot(HKL)

    return np.sum(f*np.exp(1j*phase),axis=0)

def getViolations(rho, support):
    """
    :brief:             get the indices of the pixels that have values that violate the constraints
    :param rho:         the (unfinished) reconstruction image
    :param support:     the support area of the electron density
    :return:            the indices of the pixels that violate constraints
    """
    assert(rho.shape == support.shape)

    rhoInside = rho*support
    rhoOutside = rho*(1-support)
    inside = list(np.where(rhoInside < 0))      # constraint violation in area inside the support
    outside = list(np.where(rhoOutside != 0))   # and area outside the support (constraint points for phase)

    result = np.hstack((inside, outside))
    return tuple(result)                        # to be able to index with the result it must be a tuple
