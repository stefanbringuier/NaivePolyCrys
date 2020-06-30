import sys
sys.path.append("../../build/lib/src")
from grains import Grains
from writefile import WriteData
ngrains = 4
ax = ay = az = 5*5.41
simbox = ([ax,0.0,0.0],
          [0.0,ay,0.0],
          [0.0,0.0,az])
grains = Grains(simbox,ngrains)
polycrystal= grains.makepolycrystal(['Silicon','Silicon','Germanium','Germanium'])
fileobject = WriteData('SiGex')
fileobject.writelammpsdata(polycrystal)
