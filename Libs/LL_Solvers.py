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



