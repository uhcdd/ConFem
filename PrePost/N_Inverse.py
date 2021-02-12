import sys
sys.path.insert(1, '../src')
from numpy import array, zeros
from numpy.linalg import norm, det, inv
from ConFemBasics import SamplePoints

# shape function for C3D8 brick element
def N3D8( r, s, t):
    X = array([(1.-r)*(1.-s)*(1.-t)*0.125, (1.+r)*(1.-s)*(1.-t)*0.125, (1+r)*(1+s)*(1.-t)*0.125, (1.-r)*(1.+s)*(1.-t)*0.125, (1.-r)*(1.-s)*(1.+t)*0.125, (1.+r)*(1.-s)*(1.+t)*0.125, (1+r)*(1+s)*(1.+t)*0.125, (1.-r)*(1.+s)*(1.+t)*0.125])
    return X

if __name__ == '__main__':
    # shape function inverse for C3D8 brick
    IntType  = 2
    IntOrder = 1
    IntLen   = 8
    NMat = zeros((IntLen,IntLen), dtype=float)
    for i in range(IntLen):
        r = SamplePoints[IntType,IntOrder,i][0]
        s = SamplePoints[IntType,IntOrder,i][1]
        t = SamplePoints[IntType,IntOrder,i][2]
        XXX = N3D8(r, s, t)
        print(i,r,s,t,'__',XXX, sum(XXX))
        for j in range(IntLen):
            NMat[i,j] = XXX[j]
    print(det(NMat))
    NMatI = inv(NMat)
    print(NMatI)
    for i in range(IntLen):
        print(sum(NMatI[i]))