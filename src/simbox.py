import numpy as np

class Box:
    def __init__(self,vectors=([1.0,0.0,0.0],
                               [0.0,1.0,0.0],
                               [0.0,0.0,1.0])):
        self.boxshape = 'Triclinic'
        self.units = 'Angstroms'
        dim=(3,1)
        self.i = np.reshape(vectors[0],dim)
        self.j = np.reshape(vectors[1],dim)
        self.k = np.reshape(vectors[2],dim)

    def transform(self,matrix):
        self.i = np.matrix(matrix)*self.i
        self.j = np.matrix(matrix)*self.j
        self.k = np.matrix(matrix)*self.k

        
class OrthoBox(Box):
    def __init__(self,vectors=([1.0,0.0,0.0],
                               [0.0,1.0,0.0],
                               [0.0,0.0,1.0])):
        super().__init__(vectors)
        self.boxshape = 'Orthorhombic'
        self.lx = self.i[0,0]
        self.ly = self.j[1,0]
        self.lz = self.k[2,0]
        self.sidelens = [self.lx,self.ly,self.lz] 
        self.a1 = 90.00e0
        self.a2 = 90.00e0
        self.a3 = 90.00e0
<<<<<<< HEAD
        self.angles = [self.a1,self.a2,self.a3]
=======
>>>>>>> b1b2ab01c2bb2b0ccfa586e91701178bc6a8a7e8
        
    def _transform(self,size):
        matrix = np.matrix([[size[0],0,0],
                            [0,size[1],0],
                            [0,0,size[2]]])
        self.transform(matrix)
        self.lx *= size[0]
        self.ly *= size[1]
        self.lz *= size[2]
        self.sidelens = [self.lx,self.ly,self.lz]
