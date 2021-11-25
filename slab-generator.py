from ase.build import add_vacuum
from numpy.linalg import norm
import ase.io.vasp
import numpy as np


def arrX(angle):
    Rx = np.array([[1,     0,    0],
                  [ 0,  angle,  1],
                  [ 0,   -1,  angle]])
    return(Rx)

def arrY(angle):
    Ry = np.array([[angle,  0,  1],
                   [  0,  1,  0],
                   [  -1,  0,  angle]])
    return(Ry)

def arrZ(angle):
    Rz = np.array([[ angle,  1,  0],
                  [ -1,  angle,  0],
                  [  0,  0,  1]])
    return(Rz)


family = {

 '100'  : arrY(0),
 '010'  : arrX(0),
 '001'  : np.eye(3,3),
 '110'  : arrZ(1)@arrX(0),
 '101'  : arrY(1),
 '011'  : arrX(1),
 '111'  : [-arrX(-2)@arrY(-1),-arrZ(-1),arrX(1)@arrY(1)]

}

def generate_slab(face,structure,bulkdims,dims,vac):
    dir = "structures/"
    slab = ase.io.vasp.read_vasp(dir+"POSCAR_{}_{}.vasp".format(structure,face))

    bx,by,bz = bulkdims[0],bulkdims[1],bulkdims[2]

    if face == '111':
        R1,R2,R3 = family[face][0][:,0],family[face][1][:,1],family[face][2][:,2]
        R = np.array([R1,R2,R3]).T
    else:
        R = family[face]

    bulk = ase.io.vasp.read_vasp("POSCAR")
    a,b,c = bulk.cell

    primed = np.array([a,b,c])@R
    a,b,c = norm(primed[:,0]),norm(primed[:,1]),norm(primed[:,2])

    slab.cell = np.array([a/bx,0,0]),np.array([0,b/by,0]),np.array([0,0,c/bz])
    print(" 1. bulk dimensions applied")
    slab = slab*dims
    print(" 2. slab supercell generated")
    add_vacuum(slab,vac)
    print(" 3. vacuum layer of {} A added".format(vac))
    ase.io.vasp.write_vasp("POSCAR1",slab,direct=True,sort=True)

    slab = ase.io.vasp.read_vasp("POSCAR1")
    print(" 4. updating atomic data")

    # read in bulk structure
    with open("POSCAR") as f:
        lines = f.readlines()
        f.close()

    symbols = lines[5].split()
    typecount = np.asarray(lines[6].split(),dtype=int)
    Nbulk = sum(typecount)
    fracs = np.array([x/Nbulk for x in typecount])
    Nslab = len(slab)

    slabcount = fracs*Nslab
    slabcount = np.array([int(x) for x in slabcount])

    if sum(slabcount) != Nslab:
        add = Nslab - sum(slabcount)
        select = [i for i in range(len(slabcount)) if slabcount[i] == max(slabcount)][0]
        slabcount[select] += add

    with open("POSCAR1") as f:
        lines = f.readlines()
        f.close()

    lines[5] = " "+(' '.join(symbols)+'\n')
    newline = " "+'  '.join(str(num) for num in slabcount)+"\n"
    lines[6] = newline
    with open("POSCAR1",'w') as f:
        f.writelines(lines)
        f.close()

    print(" 5. {} {} surface slab of dims {} completed".format(lattice,face,slabdims))
    print(" ---------------------- procedure complete ------------------------")


print("\n |--------------------------- Welcome ---------------------------|")
print(" | Here is a list of supported surfaces, crystal lattices, and   |"+\
"\n | atom/unit cell                                                |")
print(" |---------------------------------------------------------------|")
#print("|                                                               |")
print(" | Lat.  | 1 0 0 | 0 1 0 | 0 0 1 | 1 1 0 | 1 0 1 | 0 1 1 | 1 1 1 |")
print(" |---------------------------------------------------------------|")
print(" |  sc   |   1   |   1   |   1   |   2   |   2   |   2   |   6   |")
print(" | bcc   |   2   |   2   |   2   |   4   |   4   |   4   |   12  |")
print(" | fcc   |   4   |   4   |   4   |   8   |   8   |   8   |   24  |")
print(" | hcp   |   4   |   4   |   4   |   8   |   8   |   8   |   24  |")
print(" |---------------------------------------------------------------|")
print(" ----------------------------- inputs ----------------------------")
print(" *if bulk is not supercell provide: 1 1 1")
face = input(" 1. choose a face: ")
lattice = input(" 2. lattice geometry of bulk: ")
bulkdims = np.asarray(input(" 3. if bulk structure is a supercell\n    provide the supercell dimensions* (x y z): ").split(),dtype=int)
slabdims = np.asarray(input(" 4. surface slab dimensions (x y z): ").split(),dtype=int)
vacuum = int(input(" 5. vacuum thickness: "))
print(" ---------------------- procedure initiating ----------------------")
bulkdims = (bulkdims[0],bulkdims[1],bulkdims[2])
slabdims = (slabdims[0],slabdims[1],slabdims[2])


generate_slab(face,lattice,bulkdims,slabdims,vacuum)
print(" ---------------------- please see POSCAR1 ------------------------")
