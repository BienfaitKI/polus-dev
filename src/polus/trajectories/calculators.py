import os


def RotateGeometry(Ref_geom,Test_geom,rotMethod="K"):
    import numpy as np
    import sys
    from scipy.spatial.transform import Rotation as R
    if rotMethod == "R":
        rot, _, _ = R.align_vectors(Ref_geom, Test_geom, return_sensitivity=True)
        ROT_geom = np.array(rot.apply(Test_geom))
    else:
        Ref_geom, Test_geom = np.array(Ref_geom), np.array(Test_geom)
        R, c, t  = kabsch_umeyama(Ref_geom,Test_geom)
        ROT_geom = np.array([t + c * R @ b for b in Test_geom])

    return ROT_geom



def ComputeRMSD(geom1,geom2,rotMethod,rotated,weights_):
    from polus.trajectories.globals import weightsVect
    import numpy as np
    import math
    import sys
    from scipy.spatial.transform import Rotation as R
    if not isinstance(geom1,np.ndarray):
        g1 = np.array(geom1)
    else:
        g1 = geom1
    if not isinstance(geom2,np.ndarray):
        g2 = np.array(geom2)
    else:
        g2 = geom2
    # Apply weight correction
    assert g1.shape == g2.shape
    # Rotate structures if necessary (not yet rotated)
    #if not rotated:
    if not rotated:
        if rotMethod == "R":
            rot, _, _ = R.align_vectors(g1, g2, return_sensitivity=True)
            diff = np.array(g1) - np.array(rot.apply(g2))
        else:
            R, c, t  = kabsch_umeyama(g1,g2)
            g2 = np.array([t + c * R @ b for b in g2])
            diff = g1 - g2
    else:
        diff = g1 - g2
    # Apply weights correction
    for i in range(diff.shape[0]):
        diff[i,:] = math.sqrt(weights_[i])*diff[i,:]

    rmsd = np.sqrt((diff * diff).sum() / np.array(g2).shape[0]) 

    return rmsd

#rigid_transform_3D function taken from:
# https://gist.github.com/oshea00/dfb7d657feca009bf4d095d4cb8ea4be
# Needs fixing
def kabsch(A, B, scale=False):
    import numpy as np
    import sys
    assert len(A) == len(B)
    N = A.shape[0];  # total points
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))
    # dot is matrix multiplication for array
    if scale:
        #H = np.transpose(BB) * AA / N
        H = np.matmul(np.tranpose(BB) * AA) / N
    else:
        #H = np.transpose(BB) * AA
        H = np.matmul(np.transpose(BB), AA) 
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T * U.T
    # special reflection case
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T * U.T
    if scale:
        varA = np.var(A, axis=0).sum()
        c = 1 / (1 / varA * np.sum(S))  # scale factor
        t = -R * (centroid_B.T * c) + centroid_A.T
    else:
        c = 1
        t = -R * centroid_B.T + centroid_A.T
        #print(t)
        #sys.exit()
    return c, R, t

# Needs fixing
def kabsch_(A, B, scale=False):
    import numpy as np
    assert len(A) == len(B)
    N = A.shape[0];  # total points
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))
    # dot is matrix multiplication for array
    if scale:
        H = np.transpose(BB) * AA / N
    else:
        H = np.transpose(BB) * AA
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T * U.T
    # special reflection case
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T * U.T
    if scale:
        varA = np.var(A, axis=0).sum()
        c = 1 / (1 / varA * np.sum(S))  # scale factor
        t = -R * (centroid_B.T * c) + centroid_A.T
    else:
        c = 1
        t = -R * centroid_B.T + centroid_A.T
    return c, R, t

def kabsch_umeyama(A, B,scale=False):
    import numpy as np
    assert A.shape == B.shape
    n, m = A.shape

    EA = np.mean(A, axis=0)
    EB = np.mean(B, axis=0)
    VarA = np.mean(np.linalg.norm(A - EA, axis=1) ** 2)

    H = ((A - EA).T @ (B - EB)) / n
    U, D, VT = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U) * np.linalg.det(VT))
    S = np.diag([1] * (m - 1) + [d])

    R = U @ S @ VT
    c = VarA / np.trace(np.diag(D) @ S)
    if not scale:
        c = 1
    t = EA - c * R @ EB

    return R, c, t

