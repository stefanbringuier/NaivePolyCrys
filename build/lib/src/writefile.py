from crystal import Crystal

class WriteData:

    def __init__(self,filename="out.config"):

        self.filename = filename
        self.formt = ''
        self.extension = ''
        self.lammpsstyle = ''
        
    def writelammpsdata(self,polycrystal,**kwargs):
        if "filename" in kwargs:
            self.filename = kwargs["filename"]
            
        self.extension = 'data'
        self.formt = 'LAMMPS data'
        self.lammpsstyle = 'molecule'
        
        fullfilename = self.filename+"."+self.extension
        f = open(fullfilename,'w')
        natoms = polycrystal.natoms
        nspecies = polycrystal.nspecies
        
        #Species from crystal for grain ID=1, one may want to modify so
        #that info of each grain crystal type is stored.
        typecrys = polycrystal.graintypes[0]
        f.write('# %s %s Polycrystalline %s with %i grains\n' %(self.formt,
                                                             self.lammpsstyle,
                                                             typecrys,
                                                            polycrystal.numcrystallites)) 
        f.write('\n')
        f.write('%i atoms \n' %(natoms))
        f.write('%i atom types \n' %(nspecies))

        #Turn into generic function
        if polycrystal.boxshape == 'Orthorhombic':
            f.write('0.0 %f xlo xhi \n' %(polycrystal.lx))
            f.write('0.0 %f ylo yhi \n' %(polycrystal.ly))
            f.write('0.0 %f zlo zhi \n' %(polycrystal.lz))

        f.write('\n')
        f.write('Atoms \n')
        f.write('\n')
        for i in range(natoms):
            gid = polycrystal.polycrystal.get_tags()[i]+1
            ityp = polycrystal.polycrystal.get_chemical_symbols()[i]
            tid = polycrystal.chemmap[ityp]
            x,y,z = polycrystal.polycrystal.positions[i,:]
            f.write('%i %i %i %f %f %f \n' %(i+1,gid,tid,x,y,z))
        f.close()
