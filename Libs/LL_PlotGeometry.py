import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

def createFigure():
      """
      Function to create an axed object
      """
      fig = plt.figure(figsize=(6,6), dpi=120)
      ax = fig.add_subplot( projection='3d')
      ax.set_proj_type('persp', focal_length=0.2)

      #plt.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)

      # make the panes transparent
      ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
      ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
      ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
      
      ## make the grid lines transparent
      ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
      ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
      ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

      return ax



def plotRotor(ax, bladepanels):
      """
      Function to plot a beaultiful presentation of the rotor's blade
      """
      for panel in bladepanels:
            xs = [ [panel.p1[0], panel.p2[0]], [panel.p4[0], panel.p3[0]] ]
            ys = [ [panel.p1[1], panel.p2[1]], [panel.p4[1], panel.p3[1]] ]
            zs = [ [panel.p1[2], panel.p2[2]], [panel.p4[2], panel.p3[2]] ]

            xs, ys, zs = np.array(xs), np.array(ys), np.array(zs)
            ax.plot_surface(xs, ys, zs, alpha=0.9, color = 'grey', zorder=99, shade=True, rasterized=True)
            

def plotTrailingVortexes(ax, trailingVortexes, npanel):
      """
      Function to plot the wake. Plot first line of p1 and all other lines from p2.
      """
      
      for i, tv in enumerate(trailingVortexes):
            p1 = tv.p1
            p2 = tv.p2

            ax.plot3D(p1[0], p1[1], p1[2], '.-', markersize=0.3, linewidth=0.2, color=cm.plasma(i/npanel), rasterized=True )
            ax.plot3D(p2[0], p2[1], p2[2], '.-', markersize=0.3, linewidth=0.2, color=cm.plasma(i/npanel), rasterized=True )