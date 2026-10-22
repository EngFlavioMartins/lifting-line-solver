import numpy as np


def initializeVars(ncp):
    """ 
    Small function to initialize the variables. It just saves space, nothing else.
    """

    GammaNew = np.zeros(ncp)
    a        = np.zeros(ncp)
    aline    = np.zeros(ncp) 
    r_R      = np.zeros(ncp) 
    Fnorm    = np.zeros(ncp) 
    Ftan     = np.zeros(ncp) 
    Gamma    = np.zeros(ncp)
    alpha    = np.zeros(ncp)
    inflowAngle = np.zeros(ncp)

    return GammaNew, a, aline, r_R, Fnorm, Ftan, Gamma, alpha, inflowAngle



def update_Gamma_single_ring(Gamma, GammaNew, ConvWeight):
    """
    Update the circulation associated to each filament of each ring
    """
    return (1 - ConvWeight)*GammaNew + ConvWeight*Gamma



def velocity_induced_single_ring(ring, Xcp):
    """
    Calculate induced velocity at a control point Xcp from the four filaments of a vortex ring
    """
    
    CORE = 1e-3 # Vortex core size

    # All (xyz)-coordinates on the intradorse trailing edge vortices:
    xyz1 = np.vstack((ring.p1[0], ring.p1[1], ring.p1[2]))

    # All (xyz)-coordinates on the extradorse trailing edge vortices:
    xyz2 = np.vstack((ring.p2[0], ring.p2[1], ring.p2[2]))

    ntimes = len(ring.p1[0])-1

    vel = np.zeros(3) #, velocity_3D_from_vortex_filament(ring.Gamma, X1[:,0], X2[:,0], Xcp, CORE)

    for it in range(ntimes):

        # Determine the four points forming a vortex ring:
        P1 = xyz1[:,it]
        P2 = xyz2[:,it]
        P3 = xyz2[:,it+1] 
        P4 = xyz1[:,it+1] 

        vel += velocity_3D_from_vortex_filament(ring.Gamma, P1, P2, Xcp, CORE)
        vel += velocity_3D_from_vortex_filament(ring.Gamma, P2, P3, Xcp, CORE)
        vel += velocity_3D_from_vortex_filament(ring.Gamma, P3, P4, Xcp, CORE)
        vel += velocity_3D_from_vortex_filament(ring.Gamma, P4, P1, Xcp, CORE)

    return vel



## Determine velocity induced by a vortex filament at a point p
def velocity_3D_from_vortex_filament(GAMMA, X1, X2, P, CORE):

    R1 = np.sqrt((P[0] - X1[0])**2 + (P[1] - X1[1])**2 + (P[2] - X1[2])**2)
    R2 = np.sqrt((P[0] - X2[0])**2 + (P[1] - X2[1])**2 + (P[2] - X2[2])**2)

    R1XR2_X =  (P[1]-X1[1]) * (P[2]-X2[2]) - (P[2]-X1[2]) * (P[1]-X2[1])
    R1XR2_Y = -(P[0]-X1[0]) * (P[2]-X2[2]) + (P[2]-X1[2]) * (P[0]-X2[0])
    R1XR2_Z =  (P[0]-X1[0]) * (P[1]-X2[1]) - (P[1]-X1[1]) * (P[0]-X2[0])
    
    R1XR_SQR = R1XR2_X**2 + R1XR2_Y**2 + R1XR2_Z**2

    R0R1 = (X2[0]-X1[0]) * (P[0]-X1[0]) + (X2[1]-X1[1]) * (P[1]-X1[1]) + (X2[2]-X1[2]) * (P[2]-X1[2])
    R0R2 = (X2[0]-X1[0]) * (P[0]-X2[0]) + (X2[1]-X1[1]) * (P[1]-X2[1]) + (X2[2]-X1[2]) * (P[2]-X2[2])

    if R1XR_SQR < CORE ** 2:
        R1XR_SQR = CORE ** 2

    elif R1 < CORE:
        R1 = CORE # Gamma = 0

    elif R2 < CORE:
        R2 = CORE # Gamma = 0
        
    K = (GAMMA / (4 * np.pi * R1XR_SQR)) * (R0R1/R1 - R0R2/R2)

    U = K * R1XR2_X
    V = K * R1XR2_Y
    W = K * R1XR2_Z

    Vel = [U, V, W]

    return Vel






# ===================================================================
# Vectorised induction (same physics and core treatment as above)
# ===================================================================

def ring_segments(ring):
      """Start and end points (n_seg, 3) of every filament of a horseshoe ring, in the order
      used by velocity_induced_single_ring: for each wake station, P1-P2, P2-P3, P3-P4, P4-P1."""
      xyz1, xyz2 = np.asarray(ring.p1, float).T, np.asarray(ring.p2, float).T     # (n_pts, 3)
      P1, P2, P3, P4 = xyz1[:-1], xyz2[:-1], xyz2[1:], xyz1[1:]
      starts = np.stack([P1, P2, P3, P4], axis=1).reshape(-1, 3)
      ends   = np.stack([P2, P3, P4, P1], axis=1).reshape(-1, 3)
      return starts, ends


def filament_velocities(P, X1, X2, CORE=1e-3):
      """Velocity per unit circulation induced at points P (n, 3) by straight filaments X1 -> X2
      (m, 3); returns (n, m, 3). Same formula and core clipping as velocity_3D_from_vortex_filament."""
      r1 = P[:, None, :] - X1[None, :, :]
      r2 = P[:, None, :] - X2[None, :, :]
      R1 = np.sqrt(np.einsum("nmk,nmk->nm", r1, r1))
      R2 = np.sqrt(np.einsum("nmk,nmk->nm", r2, r2))
      cross = np.cross(r1, r2)
      sq = np.einsum("nmk,nmk->nm", cross, cross)
      r0 = (X2 - X1)[None, :, :]
      R0R1 = np.einsum("nmk,nmk->nm", np.broadcast_to(r0, r1.shape), r1)
      R0R2 = np.einsum("nmk,nmk->nm", np.broadcast_to(r0, r2.shape), r2)
      # the scalar version clips with if / elif / elif: only the first matching condition applies
      m1 = sq < CORE ** 2
      m2 = ~m1 & (R1 < CORE)
      m3 = ~m1 & ~m2 & (R2 < CORE)
      sq = np.where(m1, CORE ** 2, sq)
      R1 = np.where(m2, CORE, R1)
      R2 = np.where(m3, CORE, R2)
      K = (R0R1 / R1 - R0R2 / R2) / (4 * np.pi * sq)
      return K[..., None] * cross


def induction_matrices(bladepanels, trailingVortexes, chunk=16):
      """Matrices (ncp, nrings) of the velocity components induced at each control point by
      each ring at unit circulation. Equivalent to calling velocity_induced_single_ring for
      every pair, but evaluated with NumPy in blocks of control points."""
      cps = np.array([p.cp for p in bladepanels], float)
      segs = [ring_segments(ring) for ring in trailingVortexes]
      nseg = segs[0][0].shape[0]
      X1 = np.concatenate([s for s, _ in segs])
      X2 = np.concatenate([e for _, e in segs])
      M = np.empty((len(cps), len(segs), 3))
      for i in range(0, len(cps), chunk):
            v = filament_velocities(cps[i:i + chunk], X1, X2)          # (c, nrings * nseg, 3)
            M[i:i + chunk] = v.reshape(v.shape[0], len(segs), nseg, 3).sum(axis=2)
      return M[..., 0], M[..., 1], M[..., 2]
