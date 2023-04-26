import genfibers as gf
import numpy as np

class Fiber():
    """Elipsoid fibers with shells"""
    def __init__(self, rshell, mshell:int, rinner, linner, minner:int, thetas, origin):
        self.rshell = rshell
        self.mshell = mshell
        self.rinner = rinner
        self.linner = linner
        self.minner = minner
        self.thetas = thetas
        self.origin = origin

    def __call__(self, points):
        sigmas = [self.rshell + self.linner, self.rshell + self.rinner, self.rshell + self.rinner]
        matmap = 0
        matmap += self.mshell * gf.gengeneral(points, gf.sphere, sigmas, self.thetas, self.origin)

        sigmas = [self.linner, self.rinner, self.rinner]
        matmap += (self.minner - self.mshell) * gf.gengeneral(points, gf.sphere, sigmas, self.thetas, self.origin)

        return matmap
    
class Slab():
    """Slab function for a material"""
    def __init__(self, extent, mslab:int, origin):
        """defines a patch of skin materialed ovoids"""
        self.extent = extent
        self.mslab = mslab
        self.origin = origin

    def __call__(self, points):
        sigmas = [ex/2 for ex in self.extent]
        return self.mslab * gf.gengeneral(points, gf.cube, sigmas, [0,0,0], self.origin)
    
class SampSet():
    def __init__(self, width, tskin, mskin, tfat, mfat, tcon, mcon, tmusc, mmusc, minter, theta, rbun, sigbun, rfib, rfrac, sigfrac, seed=0):
        """
        width: width of domain in x/y
        tskin: skin thickness
        tfat: fat thickness
        tmusc: muscle thickness
        theta: angle of muscle fibers from vertical (no horizontal)
        rbun: radius of fiber bundle
        sigbun: standard deviation of budle radius
        rfib: thisckness of fiber sheath
        rfrac: the fraction of spacing between fibers that the radius takes up
        sigfrac: variation from grid spacing
        """

        # generate skin slab
        skinpar = {
            'extent':[width, width, tskin],
            'origin':[0, 0, tskin/2],
            'mslab':mskin
        }
        self.skin = Slab(**skinpar)
        
        # generate fat slab
        fatpar = {
            'extent':[width, width, tfat],
            'origin':[0, 0, tskin + tfat/2],
            'mslab':mfat
        }
        self.fat = Slab(**fatpar)

        # generate fat slab
        fatpar = {
            'extent':[width, width, tfat],
            'origin':[0, 0, tskin + tfat/2],
            'mslab':mfat
        }
        self.fat = Slab(**fatpar)

        # generate connective tissue slab
        conpar = {
            'extent':[width, width, tcon],
            'origin':[0, 0, tskin + tfat + tcon/2],
            'mslab':mcon
        }
        self.con = Slab(**conpar)

        # generate muscle interstitial tissue slab
        intpar = {
            'extent':[width, width, tmusc],
            'origin':[0, 0, tskin + tfat + tcon + tmusc/2],
            'mslab':minter
        }
        self.int = Slab(**intpar)

        # generate muscel fibers within the xy grid (hex packed)
        rng = np.random.default_rng(seed)
        dx = (rbun+rfib)/(rfrac*np.cos(theta))
        dy = np.cos(np.pi/6)*(rbun+rfib)/rfrac
        l = (2*tcon + tmusc)/np.cos(theta)
        xcoords = np.arange(-(width+l)/2, (width+l)/2, dx)
        xcoords += rng.normal(0, dx*sigfrac, len(xcoords))
        ycoords = np.arange(-(width+l)/2, (width+l)/2, dy)
        ycoords += rng.normal(0, dy*sigfrac, len(ycoords))
        rbuns = rng.normal(rbun, sigbun, len(xcoords)*len(ycoords))


        self.fibers = []
        for ix in range(len(xcoords)):
            for iy in range(len(ycoords)):
                if iy % 2 == 0:
                    offset = 0
                else:
                    offset = dx/2
                self.fibers.append(Fiber(rfib, mcon, rbuns[(ix+1)*iy], l/2, mmusc, [0, np.pi/2-theta, 0], [xcoords[ix]+offset, ycoords[iy], tskin+tfat+tcon+tmusc/2]))

    def __call__(self, points):
        matmap = np.zeros(points.shape[1], dtype=np.int8)
        
        # set the skin
        skinmap = self.skin(points)
        skinmask = skinmap != 0
        matmap[skinmask] = skinmap[skinmask]

        # set the fat
        fatmap = self.fat(points)
        fatmask = fatmap != 0
        matmap[fatmask] = fatmap[fatmask]

        # set the interstitial tissue
        intmap = self.int(points)
        intmask = intmap != 0
        matmap[intmask] = intmap[intmask]

        # Overlay each fiber on the model
        for fiber in self.fibers:
            fibmap = fiber(points)
            fibmask = fibmap != 0
            matmap[fibmask] = fibmap[fibmask]

        # set the connective tissue layer
        conmap = self.con(points)
        conmask = conmap != 0
        matmap[conmask] = conmap[conmask]

        return matmap

    
def test_slab():
    import matplotlib.pyplot as plt

    slabparams = {
        'extent':[10, 8, 3],
        'mslab':7,
        'origin':[0, 0, 0]
    }

    slab = Slab(**slabparams)

    x = np.linspace(-15, 15, 256)
    y = 0
    z = np.linspace(-15, 15, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = slab(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-15, 15, 15, -15])
    plt.xlabel("X")
    plt.ylabel("Z")
    plt.colorbar()

    y = np.linspace(-15, 15, 256)
    x = 0
    z = np.linspace(-15, 15, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = slab(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-15, 15, 15, -15])
    plt.xlabel("Y")
    plt.ylabel("Z")
    plt.colorbar()

    x = np.linspace(-15, 15, 256)
    z = 0
    y = np.linspace(-15, 15, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = slab(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-15, 15, 15, -15])
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.colorbar()
    
def test_fibers():
    import matplotlib.pyplot as plt
    fibparams = {
        'rshell': 0.5,
        'mshell': 4,
        'rinner': 1,
        'linner': 15,
        'minner': 3,
        'thetas': np.pi*np.array([45, 30, 0])/180,
        'origin': [0, 0, 1]
    }

    fib = Fiber(**fibparams)
    
    x = np.linspace(-15, 15, 256)
    y = 0
    z = np.linspace(-15, 15, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = fib(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-15, 15, 15, -15])
    plt.xlabel("X")
    plt.ylabel("Z")
    plt.colorbar()

    y = np.linspace(-15, 15, 256)
    x = 0
    z = np.linspace(-15, 15, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = fib(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-15, 15, 15, -15])
    plt.xlabel("Y")
    plt.ylabel("Z")
    plt.colorbar()

    x = np.linspace(-15, 15, 256)
    z = 0
    y = np.linspace(-15, 15, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = fib(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-15, 15, 15, -15])
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.colorbar()

def test_samp():
    sampparams = {
        'width':40,
        'tskin':3,
        'mskin':1,
        'tfat':5,
        'mfat':2,
        'tcon':1.5,
        'mcon':3,
        'tmusc':30.5,
        'mmusc':4,
        'minter':5,
        'theta':60*np.pi/180,
        'rbun':0.5,
        'sigbun':0.1,
        'rfib':0.1,
        'rfrac':0.45,
        'sigfrac':0.05
    }
    set = SampSet(**sampparams)
    
    x = np.linspace(-20, 20, 256)
    y = 0
    z = np.linspace(0, 40, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = set(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-20, 20, 40, 0])
    plt.xlabel("X")
    plt.ylabel("Z")
    plt.colorbar()

    y = np.linspace(-20, 20, 256)
    x = -6
    z = np.linspace(0, 40, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = set(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-20, 20, 40, 0])
    plt.xlabel("Y")
    plt.ylabel("Z")
    plt.colorbar()

    x = np.linspace(-20, 20, 256)
    z = 20
    y = np.linspace(-20, 20, 256)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.array([X.flatten(), Y.flatten(), Z.flatten()])

    mapped = set(points)
    mapped = mapped.reshape((256,256)).T

    plt.figure()
    plt.imshow(mapped, extent=[-20, 20, 20, -20])
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.colorbar()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # test_fibers()
    # test_slab()
    test_samp()
    plt.show()