import numpy as np


def blade(ChordCoefs, TwistCoefs, r, R):
      ca, cb, cc = ChordCoefs[0], ChordCoefs[1], ChordCoefs[2]
      ta, tb, tc = TwistCoefs[0], TwistCoefs[1], TwistCoefs[2]
      
      Chord = ca*(1 - cb*(r/R) ) + cc  # chord sizing distribuition, m
      Twist = ta*(1 - tb*(r/R) ) + tc  # blade twist angle (0=vertical, 90=horizontal), degress

      return Chord, Twist * np.pi/180



class bladePanel:
      """
      Contains information related to a rotor blade panel.
      """
      def __init__(self, ChordCoefs, TwistCoefs, r, dr, R, theta):
            """
            Initializes the panel.
            
            For a given r-coordinate of the panel's center, a chord-length c, and a panel thickness s,
            determines the 4 coordinates of the rotor's panel edges.

            Inputs
            ----------
            r:     float - radial position of the panel's center
            c:     float - chord length of the panel
            ds:    float - spanwise-thickness of the panel
            twist: float - twist angle of the panel
            theta: float - angle of the blade relatively to the x-plane
            
            Parameters
            ----------
            p1: [x,y,z]-coordinate of the first end-point.
            p2: [x,y,z]-coordinate of the second end-point.
            p3: [x,y,z]-coordinate of the third end-point.
            p4: [x,y,z]-coordinate of the fourth end-point.
            cp: [x,y,z]-coordinate of the control point
            normalDir: normal direction of the panel
            tanDir:    tangent direction of the panel
            """

            c1, twist1 = blade(ChordCoefs, TwistCoefs, r-dr/2, R)
            c2, twist2 = blade(ChordCoefs, TwistCoefs, r+dr/2, R)

            p1 = [-0.25 * c1 * np.sin(-1 * twist1), r-dr/2,  0.25 * c1 * np.cos(twist1)]
            p2 = [-0.25 * c2 * np.sin(-1 * twist2), r+dr/2,  0.25 * c2 * np.cos(twist2)]
            p3 = [ 0.75 * c2 * np.sin(-1 * twist2), r+dr/2, -0.75 * c2 * np.cos(twist2)] 
            p4 = [ 0.75 * c1 * np.sin(-1 * twist1), r-dr/2, -0.75 * c1 * np.cos(twist1)]

            ## Set panel's 4-points

            self.p1 = [p1[0], p1[1] * np.cos(theta) - p1[2] * np.sin(theta), p1[1] * np.sin(theta) + p1[2] * np.cos(theta)]
            self.p2 = [p2[0], p2[1] * np.cos(theta) - p2[2] * np.sin(theta), p2[1] * np.sin(theta) + p2[2] * np.cos(theta)]
            self.p3 = [p3[0], p3[1] * np.cos(theta) - p3[2] * np.sin(theta), p3[1] * np.sin(theta) + p3[2] * np.cos(theta)]
            self.p4 = [p4[0], p4[1] * np.cos(theta) - p4[2] * np.sin(theta), p4[1] * np.sin(theta) + p4[2] * np.cos(theta)]


            ## These are the coordinates, normal and tangent directions of the the panel's control point
            coord   = [0, r, 0]
            normal  = [     np.cos(theta), 0, -1 * np.sin(theta)]
            tangent = [-1 * np.sin(theta), 0, -1 * np.cos(theta)]

            ## Set control point, normal and tangential directions of panel:

            self.cp = [0, coord[1] * np.cos(theta) - coord[2] * np.sin(theta), 
                          coord[1] * np.sin(theta) + coord[2] * np.cos(theta)]
            
            self.normDir = [normal[0], normal[1] * np.cos(theta) - normal[2] * np.sin(theta), 
                                       normal[1] * np.sin(theta) + normal[2] * np.cos(theta)]
            
            self.tanDir  = [tangent[0], tangent[1] * np.cos(theta) - tangent[2] * np.sin(theta), 
                                        tangent[1] * np.sin(theta) + tangent[2] * np.cos(theta)]
            
            ## These are the chord length and blade twist at the control point
            ## It is is average value of the blade's section
            chord = (c1+c2)/2
            twist = (twist1+twist2)/2

            self.chord = chord
            self.twist = twist
            self.Gamma = 1.0    # Circulation value
            self.vt = 0.0       # tangential velocity
            #self.cp = 0.0       # pressure coefficient
            
      def show(self):
            print("cp, normal-dir, tang-dir, tang-vel., Circulation")
            print(self.cp, self.normDir, self.tanDir, self.vt, self.Gamma)






class trailingVortex:
      """
      Contains information related to a wake panels.
      """
      def __init__(self, ChordCoefs, TwistCoefs, r, dr, R, theta, Omega, Uwind, dt, ndt, panel):
        
            """
            Wake panel design:

            (r=R)
            ^
            |
            |     Gamma
            P2 ------------+
            |              |
            | Gamma        | Gamma
            |              | 
            P1 ------------+
            |     Gamma
            |
            +----------------------------> x
            (r=0)

            Initializes the panel.
        
            For a given r-coordinate of the panel's center, a chord-length c, and a panel thickness s,
            determines the 4 coordinates of the rotor's panel edges.

            Inputs
            ----------
            n:     float     - number of panels for the discretization
            N:     float     - number of blade elements
            TSR:   float     - tip speed ratio
            r:     Nx1 array - radial positions of the panels' center
            c:     Nx1 array - chord lengths of the panel
            dr:    Nx1 array - spanwise-thickness of the panel
            twist: Nx1 array - twist angle of the panel
          
            Parameters
            ----------
            xxxx
            """

            ndt += 2 # account for leading edge and control point's 2 extra points

            c1, twist1 = blade(ChordCoefs, TwistCoefs, r-dr/2, R)
            c2, twist2 = blade(ChordCoefs, TwistCoefs, r+dr/2, R)

            P1 = np.zeros((3,ndt))
            P2 = np.zeros((3,ndt))

            for n in range(ndt):
                 
                  ## Connecting to the bounding vortex:
                  if n == 0:
                        p1 = [0, r-dr/2,  0]
                        p2 = [0, r+dr/2,  0]

                        P1[:,n] = [p1[0], p1[1] * np.cos(theta) - p1[2] * np.sin(theta), 
                                          p1[1] * np.sin(theta) + p1[2] * np.cos(theta)]
                        P2[:,n] = [p2[0], p2[1] * np.cos(theta) - p2[2] * np.sin(theta), 
                                          p2[1] * np.sin(theta) + p2[2] * np.cos(theta)]
            
                  
                  ## Trailing edge vortex:
                  elif n == 1:
                        p1 = [ c1 * np.sin(-1 * twist1), r-dr/2, -1 * c1 * np.cos(twist1)]
                        p2 = [ c2 * np.sin(-1 * twist2), r+dr/2, -1 * c2 * np.cos(twist2)] 

                        P1[:,n] = [p1[0], p1[1] * np.cos(theta) - p1[2] * np.sin(theta), 
                                          p1[1] * np.sin(theta) + p1[2] * np.cos(theta)]
                        P2[:,n] = [p2[0], p2[1] * np.cos(theta) - p2[2] * np.sin(theta), 
                                          p2[1] * np.sin(theta) + p2[2] * np.cos(theta)]
                  
                  ## Remaining points:
                  elif n>1:
                        angle1 = theta + (n-2)*dt*Omega
                        angle2 = theta + (n-2)*dt*Omega
                        
                        P1[:,n] = [p1[0]+(n-2)*dt*Uwind, p1[1] * np.cos(angle1) - p1[2] * np.sin(angle1), p1[1] * np.sin(angle1) + p1[2] * np.cos(angle1)]
                        P2[:,n] = [p2[0]+(n-2)*dt*Uwind, p2[1] * np.cos(angle2) - p2[2] * np.sin(angle2), p2[1] * np.sin(angle2) + p2[2] * np.cos(angle2)]

            self.p1 = P1
            self.p2 = P2
            self.Gamma = 1

            # attribute an idea (panel=0,1,...) to associate the current ring with the panel-th control point:
            self.panel = panel 

      def show(self):
            print("p1, p2, Circulation")
            print(self.p1, self.p2, self.Gamma)





class solutionVector:
      """
      Contains information related to the output solution vector
      """
      def __init__(self, ncp):

            self.a     = np.zeros(ncp)
            self.aline = np.zeros(ncp)
            self.r_R   = np.zeros(ncp)
            self.Fnorm = np.zeros(ncp)
            self.Ftan  = np.zeros(ncp)
            self.Gamma = np.zeros(ncp)
            self.alpha = np.zeros(ncp)

      def updateSolutionVector(self, a, aline, r_R, Fnorm, Ftan, Gamma):

            self.a = a
            self.aline = aline
            self.r_R   = r_R
            self.Fnorm = Fnorm
            self.Ftan  = Ftan
            self.Gamma = Gamma
            self.alpha = alpha