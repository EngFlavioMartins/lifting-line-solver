import numpy as np
import pandas as pd

from Libs.LL_Classes import bladePanel
from Libs.LL_Classes import trailingVortex
from Libs.LL_Classes import solutionVector


# =====================================================


def construstGeometry(ChordCoefs, TwistCoefs, r, dr, dt, Uwind, Omega, nb=3, ncp=40, ndt=10):
      """
      Discretizes the geometry into panels using 'cosine' method.
      
      Parameters
      ----------
      N: integer, optional
            Number of blades;
            default: 3.
      n: integer, optional
            Number of panels;
            default: 40.
      r:  1D array of floats
            radial coordinates of the control points;
      
      Returns
      -------
      panels: 1D Numpy array of Panel objects.
            The list of panels
      """
      
      # Determine twist angle, theta [rad]:
      theta = (2 * np.pi / nb) * np.arange(0, nb, 1, int)

      # Create panels
      bladepanels = np.empty(nb*ncp, dtype=object)

      trailingVortexes = np.empty(nb*ncp, dtype=object)

      Radius = r[-1] + dr[-1]/2 # radius of the blade
      
      for blade in range(nb):
            for panel in range(ncp):
                  panelID = int(panel+blade*ncp)
                  bladepanels[panelID] = bladePanel(ChordCoefs, TwistCoefs, r[panel], dr[panel], Radius, theta[blade])
                  trailingVortexes[panelID] = trailingVortex(ChordCoefs, TwistCoefs, r[panel], dr[panel], Radius, theta[blade], Omega, Uwind, dt, ndt, panel)

      # Pre-allocate the solution vector:
      mySolutionVector = solutionVector(ncp)  

      return bladepanels, trailingVortexes, mySolutionVector



def createRadius(Rroot, R, ncp, method='linear'):

      if method == 'linear':
            r, dr = np.linspace(Rroot, R, num=ncp+1, retstep=True)

            rm = (r[1:] + r[:-1])/2 # location of the control points
            dr = dr*np.ones_like(rm)


      ## Apply cosine distribuition to points:
      elif method == 'cosine':
            theta = np.linspace(0, np.pi, num=ncp+1) # array of angles
            r = (R+Rroot)/2 - (R-Rroot)*0.5*np.cos(theta)

            rm = (r[1:] + r[:-1])/2 # location of the control points
            dr = r[1:] - r[:-1]

      return rm, dr
      
      
      
def loadBladeElement(Vnorm, Vtan, polar_alpha, polar_cl, polar_cd, chord, twist):
      """
      Calculates the aerodynamic loading at each blade element

      INPUTS:
            * Ua          -- axial speed of wind at blade section             --> (Nx1) array, [m/s]
            * Ut          -- tangential speed of the blade section            --> (Nx1) array, [m/s]
            * chord       -- chord length of each blade sectio                --> (Nx1) array, [m]
            * twist       -- twist angle of the blade section                 --> (Nx1) array, [rad]
            * polar_alpha -- polar AoA data from the airfoil database         --> (Nx1) array, [rad]
            * polar_cl    -- polar lift coeff. data from the airfoil database --> (Nx1) array
            * polar_cd    -- polar drag coeff. data from the airfoil database --> (Nx1) array
      OUTPUTS:
            * Fa          -- axial force acting on the rotor                  --> (Nx1) array, [N/m]
            * Ft          -- tangential force acting on the rotor             --> (Nx1) array, [N/m]
            * Circulation -- circulation value around each blade section ele. --> (Nx1) array, [m2/s]
            * alpha       -- AoA of each blade section ele.                   --> (Nx1) array, [rad]
            * inflowangle_deg -- inflow angle felt by blade section           --> (Nx1) array, [degree]
            * L           -- lift force acting on the rotor                   --> (Nx1) array, [N/m]
            * D           -- drag force acting on the rotor                   --> (Nx1) array, [N/m]
      """

      Umag2 = Vnorm**2 + Vtan**2

      inflowAngle = np.arctan2(Vnorm,Vtan)

      alpha = (twist + inflowAngle)*180/np.pi

      cl = np.interp(alpha, polar_alpha, polar_cl) # 2*np.pi*(alpha*np.pi/180) 
      cd = np.interp(alpha, polar_alpha, polar_cd)

      Lift = 0.5*Umag2*cl*chord
      Drag = 0 #0.5*Umag2*cd*chord

      Fnorm = Lift*np.cos(inflowAngle) + Drag*np.sin(inflowAngle)
      Ftan  = Lift*np.sin(inflowAngle) - Drag*np.cos(inflowAngle)

      Gamma = 0.5*np.sqrt(Umag2)*cl*chord

      return Fnorm , Ftan, Gamma, alpha, inflowAngle * 180/np.pi



def createAirfoils(airfoil_file):
    airfoil_data=pd.read_csv(airfoil_file, header=0, 
                            names = ["alfa", "cl", "cd", "cm"],  sep=',') # sep='\s+'

    polar_alpha = airfoil_data['alfa'][:]
    polar_cl = airfoil_data['cl'][:]
    polar_cd = airfoil_data['cd'][:]

    return polar_alpha, polar_cl, polar_cd



