import sys
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
fileobject = WriteData('Silicon')
fileobject.writelammpsdata(polycrystal)
