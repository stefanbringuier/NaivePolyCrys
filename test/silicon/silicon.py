import sys
<<<<<<< HEAD
sys.path.append("../../build/lib/src")
from grains import Grains
from writefile import WriteData
ngrains = 6
ax = ay = az = 8*5.41
simbox = ([ax,0.0,0.0],
          [0.0,ay,0.0],
          [0.0,0.0,az])
grains = Grains(simbox,3)
polycrystal= grains.makepolycrystal(['Silicon' for i in range(ngrains)])
=======
sys.path.append("../build/lib/src")
from grains import Grains
from writefile import WriteData
ngrains = 3
ax = ay = az = 10*5.41
simbox = ([ax,0.0,0.0],
          [0.0,ay,0.0],
          [0.0,0.0,az])
polycrystal = Grains(simbox,3)
polycrystal.makepolycrystal(['Silicon' for i in range(ngrains)])
>>>>>>> b1b2ab01c2bb2b0ccfa586e91701178bc6a8a7e8
fileobject = WriteData('Silicon')
fileobject.writelammpsdata(polycrystal)
