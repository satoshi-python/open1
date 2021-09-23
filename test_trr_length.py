
import MDAnalysis as MDA

for i in range(10):
    f = MDA.Universe("000" + str(i) + ".trr")
    print("000" + str(i) + ".trr:", len(f.trajectory))
for i in range(10,100):
    f = MDA.Universe("00" + str(i) + ".trr")
    print("00" + str(i) + ".trr:", len(f.trajectory))
for i in range(100,200):
    f = MDA.Universe("0" + str(i) + ".trr")
    print("0" + str(i) + ".trr:", len(f.trajectory))
