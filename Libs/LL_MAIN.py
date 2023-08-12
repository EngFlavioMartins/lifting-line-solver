import numpy as np

from Libs.LL_Solvers  import initializeVars
from Libs.LL_Solvers  import velocity_induced_single_ring
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

      ncp  = len(bladepanels)      # Get number of control points

      # Variables for maximum number of iterations and convergence criteria
      Niterations = 100
      ConvWeight  = 0.1
      GammaCrit =   0.005

      # Initialize and calculate matrices for velocity induced by horseshoe vortex rings
      MatrixU = np.zeros((ncp,ncp)) # Matrix of induction, for velocity component in x-direction
      MatrixV = np.zeros((ncp,ncp)) # Matrix of induction, for velocity component in y-direction
      MatrixW = np.zeros((ncp,ncp)) # Matrix of induction, for velocity component in z-direction

      Gamma, a, aline, r_R, Fnorm, Ftan, Gamma, alpha, inflowAngle = initializeVars(ncp)

      for i, ipanel in enumerate(bladepanels):

            induced_velocity = np.zeros(3)

            for j, jring in enumerate(trailingVortexes):

                  jring.Gamma  = 1
                  induced_velocity = velocity_induced_single_ring(jring, ipanel.cp)

                  MatrixU[i,j] = induced_velocity[0]
                  MatrixV[i,j] = induced_velocity[1]
                  MatrixW[i,j] = induced_velocity[2]

      ## Calculate solution through an iterative process
      dGamma = np.zeros_like(Gamma)

      for kiter in range(Niterations):

            for i, ipanel in enumerate(bladepanels):

                  u, v, w = 0, 0, 0

                  for j, jring in enumerate(trailingVortexes):

                        u += MatrixU[i,j] * jring.Gamma
                        v += MatrixV[i,j] * jring.Gamma
                        w += MatrixW[i,j] * jring.Gamma

                  Vrot = - np.cross([Omega, 0, 0], ipanel.cp) # Wind, relatively to the blade

                  # Wind speed, relatively to the blade:

                  U = [Uwind + u + Vrot[0], 
                               v + Vrot[1], 
                               w + Vrot[2]]
 
                  r = np.sqrt( np.dot(ipanel.cp, ipanel.cp) )
                  azimdir = np.cross([-1/r, 0, 0], ipanel.cp)
                  Vazim   = np.dot(azimdir, U)
                  Vaxial  = np.dot([1, 0, 0], U)

                  #print('Vaxial:' + str(Vaxial) )
                  #print('Vazim:' + str(Vazim) )

                  temploads = loadBladeElement(Vaxial, Vazim, polar_alpha, polar_cl, polar_cd, ipanel.chord, ipanel.twist)


                  ## Update variables:
                  a[i] = -(u + Vrot[0]) / Uwind
                  if Omega == 0: # avoid dividing by zero. Ignore aline, in this case
                        aline[i] = 1
                  else:
                        aline[i] = U[1] / (r * Omega) - 1
                  r_R[i]   = r / R
                  Fnorm[i] = temploads[0]
                  Ftan[i]  = temploads[1]
                  alpha[i] = temploads[3]
                  inflowAngle[i] = temploads[4]

                  dGamma[i] = temploads[2] - Gamma[i] # difference in Gamma


            ## Check the convergence:     
            top = np.max(np.abs(dGamma))
            bot = np.max(np.abs(Gamma)) + 1
            if top/bot < GammaCrit and kiter > 5:
            #      print('Number of iterations: ' + str(kiter) )
            #      print('Error: ' + str(top/bot) )
                  print("Done.")
                  break

            ## Set new estimate of circulation on vortex rings:
            for j, jring in enumerate(trailingVortexes):
                  jring.Gamma = Gamma[j] + ConvWeight * dGamma[j]
                  Gamma[j] = jring.Gamma

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