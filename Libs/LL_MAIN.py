import numpy as np

from Libs.LL_Solvers  import initializeVars
from Libs.LL_Solvers  import induction_matrices
from Libs.LL_ConstructObjs import createAirfoils
from Libs.LL_ConstructObjs import loadBladeElement

def liftingLineSolver(mySolutionVector, bladepanels, trailingVortexes, Uwind, Omega, R, airfoil):
      
      # This code solves a lifting line model of a horizontal axis rotor
      # Inputs:
      #       rotor_wake_system: data structure that contains the geometry of the horseshoe vortex rings
      # and the control points at the blade
      #       wind: unperturbed wind velocity, also known as U_infinity
      #       Omega: rotational velocity of the rotor
      #       rotorradius: the radius of the rotor


      polar_alpha, polar_cl, polar_cd = createAirfoils(airfoil)
      polar_alpha, polar_cl, polar_cd = (np.asarray(x, float) for x in (polar_alpha, polar_cl, polar_cd))

      ncp  = len(bladepanels)      # Get number of control points

      # Variables for maximum number of iterations and convergence criteria
      Niterations = 100
      ConvWeight  = 0.1
      GammaCrit =   0.005

      # Induction matrices: velocity at each control point due to each ring at unit circulation
      MatrixU, MatrixV, MatrixW = induction_matrices(bladepanels, trailingVortexes)

      Gamma, a, aline, r_R, Fnorm, Ftan, Gamma, alpha, inflowAngle = initializeVars(ncp)

      cps   = np.array([p.cp for p in bladepanels], float)                # (ncp, 3)
      chord = np.array([p.chord for p in bladepanels], float)
      twist = np.array([p.twist for p in bladepanels], float)
      Vrot  = -np.cross([Omega, 0, 0], cps)                               # blade motion, seen from the blade
      r     = np.sqrt(np.einsum("ij,ij->i", cps, cps))
      azimdir = np.cross(np.stack([-1 / r, np.zeros_like(r), np.zeros_like(r)], axis=1), cps)

      # Every ring starts at unit circulation, as left by the induction loop of the original code
      ringGamma = np.ones(ncp)

      ## Calculate solution through an iterative process
      for kiter in range(Niterations):

            u, v, w = MatrixU @ ringGamma, MatrixV @ ringGamma, MatrixW @ ringGamma
            U = np.stack([Uwind + u + Vrot[:, 0], v + Vrot[:, 1], w + Vrot[:, 2]], axis=1)
            Vazim  = np.einsum("ij,ij->i", azimdir, U)
            Vaxial = U[:, 0]

            Fnorm, Ftan, newGamma, alpha, inflowAngle = loadBladeElement(Vaxial, Vazim, polar_alpha, polar_cl,
                                                                         polar_cd, chord, twist)

            ## Update variables:
            a = -(u + Vrot[:, 0]) / Uwind
            aline = np.ones(ncp) if Omega == 0 else U[:, 1] / (r * Omega) - 1   # ignore aline when not rotating
            r_R = r / R
            dGamma = newGamma - Gamma

            ## Check the convergence:
            top = np.max(np.abs(dGamma))
            bot = np.max(np.abs(Gamma)) + 1
            if top/bot < GammaCrit and kiter > 5:
                  print("Done.")
                  break

            ## Set new estimate of circulation on vortex rings:
            ringGamma = Gamma + ConvWeight * dGamma
            Gamma = ringGamma.copy()

      for j, jring in enumerate(trailingVortexes):
            jring.Gamma = ringGamma[j]

      mySolutionVector.a = a
      mySolutionVector.aline = aline
      mySolutionVector.r_R   = r_R
      mySolutionVector.Fnorm = Fnorm
      mySolutionVector.Ftan  = Ftan
      mySolutionVector.Gamma = Gamma
      mySolutionVector.alpha = alpha
      mySolutionVector.inflowAngle = inflowAngle

      ## Output results of converged solution
      return mySolutionVector
