import genfibers as gf

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
        sigmas = [self.rshell + self.linner/2, self.rshell + self.rinner, self.rshell + self.rinner]
        matmap = 0
        matmap += self.mshell * gf.gengeneral(points, gf.sphere, )