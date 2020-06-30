import sys
sys.path.append("../../build/lib/src")
from grains import Grains
from writefile import WriteData
ngrains = 6
ax = ay = az = 12*3.15
simbox = ([ax,0.0,0.0],
          [0.0,ay,0.0],
          [0.0,0.0,az])
grains = Grains(simbox,ngrains)
polycrystal= grains.makepolycrystal(['Tungsten' for i in range(ngrains)])
fileobject = WriteData('Tungsten')
fileobject.writelammpsdata(polycrystal)
