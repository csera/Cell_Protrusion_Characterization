# csera
# v3.0.0

import numpy as np
import sympy as sp
import plotly.graph_objs as go

# ---------------------------------------------------------------------
# Global system properties 
# ---------------------------------------------------------------------

# System properties: protrusion & bilayer
R0 = np.array([200,400,600]) # Radius of the protrusion's cylindrical body; in nm
L = 5000 # Protrusion length; in nm
KbPore = 15 # Bending constant of a generic bilayer; in kT
KbProtr = KbPore*2  # Bending constant of the double-membrane protrusion; in kT
    # Crude estimation that treats the double-bilayer structure as single unit
H0 = 0 # Membrane spontaneous curvature; in 1/nm
h = 7.5 # Bilayer thickness; in nm

# ---------------------------------------------------------------------
# Global sim parameters 
# ---------------------------------------------------------------------
rtMax = 40 # Max rt to compute for; 40 based on Chizmadzhev (1995); nm
rtSteps = rtMax*4 # Num steps to take from [1,rtMax] nm for val of rt; nm
rt = np.linspace(1,rtMax,num=rtSteps) 
    # Radius of "fold" to form pore rim as meas from inner monolayer surface; in nm
    # Range starts at ~upper limit of intermembr distance for spontaneous fusion and 
    # extends to rtMax
m = h+rt # Distance b/w outer leaflet & intermembrane midplane; nm
print("rtSteps = ",rtSteps)

stepSize = 0.1 # Step size of pore radius dilation; nm
max_r = 50 # Maximum pore radius to do computations for; nm
rEffSteps = int(np.ceil(max_r/stepSize + 1)) # Num steps to take over expanding rEffCap
    # +1 so that I preserve num of steps after adding initial values
print("rEffSteps = ",rEffSteps)

# ---------------------------------------------------------------------
# Key equations/fxns for G calculation
# ---------------------------------------------------------------------
# Equations come from Christoph's paper
def helfrichStress(k1,k2,Kb):
    '''
    Helfrich bending *stress*
    
    Args:
        k1 (float): Principal curvature 1; 1/nm
        k2 (float): Principal curvature 2; 1/nm
        Kb (float): Helfrich bending constant for membrane; kT
    Returns:
        float: The Helfrich bending *stress* (ie energy for an infinitesimal dA) 
    '''
    return 0.5*Kb*(k1+k2)**2 # Note units = energy/unit area = kT/nm^2

def xi(r0,rtVal,mVal,theta):
    '''
    Core function as defined in Christoph's paper
    
    Args:
        r0 (float): Fusion pore radius meas from intermembrane midplane; nm
        rtVal (float): Radius of "fold" to form pore rim; nm
        mVal (float): h+rtVAl = dist b/w outer leaflet & intermemb midplane; nm
        theta (float): Pore incline angle; radians
    Returns:
        float: Val of fxn xi from Christoph's paper
    '''
    return (r0 + mVal*np.cos(theta))/rtVal;
            
def W(wXi,theta):
    '''
    Core function as defined in Christoph's paper
    
    Args:
        wXi float: Relevant value of xi
        theta float: Pore incline angle; radians
    Returns:
        float: Val of fxn W from Christoph's paper
    '''
    xiSqr = wXi**2
    sqrtXiSqr = np.sqrt(xiSqr - 1)
    
    return ( (2*xiSqr)/sqrtXiSqr * (np.arctan((np.tan(theta/2 + np.pi/4) * sqrtXiSqr)/(wXi-1))
            - np.arctan((np.tan(theta/2) * sqrtXiSqr)/(wXi-1))) )

def T(tXi,rtVal,theta):
    '''
    Core function as defined in Christoph's paper
    
    Args:
        tXi (float): Relevant value of xi
        rtVAl (float): Radius of "fold" to form pore rim; nm
        theta (float): Pore incline angle; radians
    Returns:
        float: Val of fxn T from Christoph's paper
    '''
    cosTheta = np.cos(theta)
    
    return ( -4*cosTheta - H0*rtVal*(np.pi*tXi - 4*cosTheta)
            + H0**2 * rtVal**2 * (np.pi*tXi/2 - cosTheta) )

def GPore(r0,rtVal,mVal,theta,Kb):
    '''
    Calculates energy of an inclined toroidal fusion pore from Christoph's 2011 paper
    
    Args:
        r0 (float): Fusion pore radius meas from intermembrane midplane; nm
        rtVal (float): Radius of "fold" to form pore rim; nm
        mVal (float): h+rtVAl = dist b/w outer leaflet & intermemb midplane; nm
        theta (float): The angle of pore incline; radians
        Kb (float): Helfrich bending constant for membrane; kT
    Returns:
        float: The bending energy for an inclined toroidal pore; units set by Kb
    '''
    thisXi = xi(r0,rtVal,mVal,theta)
    return np.pi * Kb * (W(thisXi,theta) + W(thisXi,-theta) + 2*T(thisXi,rtVal,theta))

# ---------------------------------------------------------------------
# Custom functions
# ---------------------------------------------------------------------
def surfMinPath(rEff,rt,G):
    '''
    Finds the least energetically costly path for pore expansion on a given surface. Minimizes with respect to 
    rEffCap. Ie finds the rt value with the location and value of the minimum GTotCap for each rEffCap slice
    
    Args:
        rEff (Numpy array): Array of length rEffSteps containing the effective pore radii
        rt (Numpy array): Array of length rtSteps containing the intermembrane distances
        GTot (Numpy array): rtSteps x rEffSteps array containing the energy for each [rt,rEff] coordinate
        
    Returns:
        Numpy array: rEffStepsx3 Numpy array of floats giving the coordinates for the minium path 
            [:,0] = rEff. [:,1] =  rt. [:,2] = G
    '''
    # Note: An alternative aproach would be "water flow" pathfinding in which the program starts at a certain 
    # point/initial condition (eg rEffCap=4, rt=1) then just goes to wherever the lowest neighboring G is on 
    # the surface w/ a greater rEffCap. I did not use this approach since it would be v sensitive to init 
    # conditions and would be computationally taxing.
    print("Finding minimum path")
    
    GMinPathIndex = np.argmin(G,axis=0) #Arr of indices of min GTotCap in each row
    print("rt indices of min path =",GMinPathIndex)
    
    GMinPath = np.zeros((rEffSteps,3))
    
    # Get data for min path
    for i in range(0,rEffSteps):
        GMinPath[i,0] = rEff[i]
        GMinPath[i,1] = rt[GMinPathIndex[i]]
        GMinPath[i,2] = G[GMinPathIndex[i],i]
    
    print('GMinPath shape: ',GMinPath.shape)
    return GMinPath

def surfMaxPath(rEff,rt,G):
    '''
    Finds the most energetically costly path for pore expansion on the surface with respect to rt. Ie it finds 
    the max GTot value and location for each value of rt.
    
    Args:
        rEff (Numpy array): Array of length rEffSteps containing the effective pore radii
        rt (Numpy array): Array of length rtSteps containing the intermembrane distances
        G (Numpy array): rtSteps x rEffSteps array containing the energy for each [rt,rEff] coordinate
            
    Returns:
        Numpy array: rtStepsx3 Numpy array of floats giving the coordinates for the maximum path 
            [:,0] = rEff. [:,1] =  rt. [:,2] = G        
    '''
    print("Finding maximum path")

    # Find max path
    GMaxPathIndex = np.argmax(G,axis=1) #Arr of indices of min GTotCap in each row
    print("rt indices of max path = ",GMaxPathIndex)

    GMaxPath = np.zeros((rtSteps,3))
    
    # Get data for max path
    for i in range(0,rtSteps):
        GMaxPath[i,0] = rEff[GMaxPathIndex[i]]
        GMaxPath[i,1] = rt[i]
        GMaxPath[i,2] = G[i,GMaxPathIndex[i]] # the [rEffCap,rt] indices
    
    print('GMaxPath shape: ',GMaxPath.shape)
    return GMaxPath
    
print("Set-up completed")

# ---------------------------------------------------------------------
# Set-up: protrusion cap variables
# ---------------------------------------------------------------------
rEffCap = np.linspace(0,max_r,num=rEffSteps) # Effective radius of pore; nm
    # This is the pore radius at the narrowest point which < r0 if theta != 0
    # Note: linspace is start and end point inclusive by default. 
    # Thus there are rEffSteps-2 evenly spaced points in the open interval (0,max_r)
    
# Set up the arrays of size rtSteps x rEffSteps to hold my computations
GCap = np.zeros((R0.size,rtSteps,rEffSteps)) # Bending *energy* of all non-pore membrane in the cap; kT
    # Dim 1 R0 value done for
    # Dim 2 rt value done for (coordinate 1)
    # Dim 3 rEffCap value done for (coordinate 2)
GCapPore = np.zeros((R0.size,rtSteps,rEffSteps)) # Energy of the fusion pore (Helfrich bending only); kT
    # Dim 1 R0 value done for
    # Dim 2 rt value done for (coordinate 1)
    # Dim 3 rEffCap value done for (coordinate 2)
GTotCap = np.zeros((R0.size,rtSteps,rEffSteps)) # Total Gb of cap = GCap + GCapPore; kT
    # Dim 1 R0 value done for
    # Dim 2 rt value done for (coordinate 1)
    # Dim 3 rEffCap value done for (coordinate 2)
groundGTotCap = np.zeros((R0.size,rtSteps,rEffSteps)) # GTotCap w/ (r,0,0) translated to Gb = 0; kT

# Paths to trace on surfaces
GMinPathCap = np.zeros((R0.size,rEffSteps,3)) # Coordinates for the min G path
    # Dim 1 R0 value done for
    # Dim 2 (length rEffSteps) = index moving down rEff
    # Dim 3 (length 3) = spatial coordinates for the min in that slice
        # 0 = x (ie rEffCap), 1 = y (ie rt), 2 = z (GTotCap) 
groundGMinPathCap = np.zeros((R0.size,rEffSteps,3)) # GMinPathCap w/ start Gb set = 0
GMaxPathCap = np.zeros((R0.size,rtSteps,3)) # Coordinates for the max G path
    # Dim 1 R0 value done for
    # Dim 2 (length rtSteps) = index moving down rt
    # Dim 3 (length 3) = spatial coordinates for the min in that slice
        # 0 = x (ie rEffCap), 1 = y (ie rt), 2 = z (GTotCap) 

# Set generic view options for 3D surfaces
view = dict(
    center=dict(x=0,y=-0.1,z=-0.1),
    eye=dict(x=1,y=2,z=0.75)
)

# ---------------------------------------------------------------------
# Iterate calculations and surface plots over R0's
# ---------------------------------------------------------------------
for r in range(0,R0.size):
    print('------------------------------------------------------------')
    
    # ---------------------------------------------------------------------
    # Set/Compute initial conditions for referencing in computations
    # ---------------------------------------------------------------------
    # Note that array[:,n] refs calumn n in all rows in array
    stressCap = helfrichStress(1/R0[r],1/R0[r],KbProtr) # The static, uniform stress of the cap membrane
    GCap[r,:,0] = 2*np.pi*(R0[r]**2) * stressCap # Init Gb of cap; kT
        # is ACap * CapStress where ACap == SA of a hemisphere of radius R
    GCapPore[r,:,0] = 0
    GTotCap[r,:,0] = GCap[r,0,0]+GCapPore[r,0,0]

    # ---------------------------------------------------------------------
    # Find GCap, GCapPore, and thus GTotCap for each value of rt and rCapPore
    # ---------------------------------------------------------------------
    # First iterate over rCapPore then rt; nesting order doesn't matter since everything uses both
    print("Computing energies for R0 = ",R0[r])
    for i in range(1,rEffSteps):
        for rtIndex in range(0,rtSteps):
            #Vars for each iteration of loop
            theta = np.arcsin(rEffCap[i]/R0[r]) # Angle of pore incline where 0 is norm to z axis; radians
                # For theta, neglect difference b/w r0 and rEffCap[i]
            r0 = rEffCap[i] + m[rtIndex]*(1-np.cos(theta)) # r0 as defined in the figures in the set-up
                # What we actually plug into Christoph & Rob's solution
            hBand = R0[r]*np.cos(theta)
            #print("hBand(",i,") = ",hBand)

            # Find Gb of remaining cap
            A = 2*np.pi*R0[r]*hBand # Area of remaining cap
            GCap[r,:,i] = A*stressCap

            # Find Gb of the inclined pore
            GCapPore[r,rtIndex,i] = GPore(r0,rt[rtIndex],m[rtIndex],theta,KbPore)
            GTotCap[r,rtIndex,i] = GCap[r,rtIndex,i]+GCapPore[r,rtIndex,i]

    # ---------------------------------------------------------------------
    # Find MINimal G expansion path for cap pore
    # ---------------------------------------------------------------------
    GMinPathCap[r,:,:] = surfMinPath(rEffCap,rt,GTotCap[r,:,:])

    # ---------------------------------------------------------------------
    # Find MAXimal G expansion path
    # ---------------------------------------------------------------------
    # This is so that I can...
    #     a) Check that there is a peak past which expansion becomes spontaneous as in 
    #         the Chizmadzhev paper
    #     b) Get a sense of the pore size at which expansion becomes spontaneous
    #     c) Get a sense of the energetic barrier that must be overcome. This won't 
    #           necessarily be the actual value of the activation energy, but it will 
    #         probably be in the ballpark since the surface isn't super curvey
    #     Since there's that huge upsweep at low rt, I'm going to maximize with respect to rt rather
    #     than with respect to rEffCap, the approach used for finding the min path
    print('GTotCap shape: ',GTotCap.shape)
    print('GTotCap[r,:,:] shape: ',GTotCap[r,:,:].shape)
    print('GMaxPathCap shape: ',GMaxPathCap.shape)
    print('GMaxPathCap[r,:,:] shape:',GMaxPathCap[r,:,:].shape)
    GMaxPathCap[r,:,:] = surfMaxPath(rEffCap,rt,GTotCap[r,:,:])

    # ---------------------------------------------------------------------
    # Graph the results for the cap pore
    # ---------------------------------------------------------------------
    print("Graphing results for cap w/ R0 = ",R0[r])
    
    totCapSurf = go.Figure()
    # Graph GTotCap (GCap + GPore) surface
    print("Graphing total surface")
    totCapSurf.add_trace(go.Surface(
        #coloraxis=go.layout.Coloraxis(dict(
        #    title="Neg G", "coloraxis1")),
        colorbar=go.surface.ColorBar(title='Net G'),
        colorscale='matter',
        contours={"z":{"show":True,"start":np.max(GTotCap[r,:,:])/20,
                       "end":np.max(GTotCap[r,:,:])*0.95,
                       "size":np.max(GTotCap[r,:,:])/20}},
        name='Net G',
        x=rEffCap,y=rt,z=GTotCap[r,:,:],opacity=0.9))
    # Graph curve of min GTotCap for expanding rEffCap
    print("Graphing minimal G path")
    totCapSurf.add_trace(go.Scatter3d(
        marker=dict(
            color='blue',
            size=2),
        mode='lines',
        name='Min G path',
        x=GMinPathCap[r,:,0],y=GMinPathCap[r,:,1],z=GMinPathCap[r,:,2]))
    # Graph curve of max GTotCap for expanding rEffCap
    print("Graphing maximal G path")
    totCapSurf.add_trace(go.Scatter3d(
        marker=dict(
            color='red',
            size=2),
        mode='lines',
        name='Max G path',
        x=GMaxPathCap[r,:,0],y=GMaxPathCap[r,:,1],z=GMaxPathCap[r,:,2]))
    totCapSurf.update_layout(coloraxis_showscale=False,
        scene = dict(
                        xaxis_title='Pore radius (nm)',
                        yaxis_title='Intermembrane dist (nm)',
                        zaxis_title='G (kT)'),
        legend_orientation="h", # Make legend horizontal
        margin=dict(r=10, b=10, l=10, t=10),
        height=600,width=700)
    totCapSurf.update_layout(scene_camera=view)
    totCapSurf.show()
    totCapSurf.write_image(f'output data/Cap Gtot R0 {R0[r]}.pdf')
    totCapSurf.write_image(f'output data/Cap Gtot R0 {R0[r]}.svg')
    totCapSurf.write_image(f'output data/Cap Gtot R0 {R0[r]}.png')

    # Graph GCap component of GTotCap surface
    print("Graphing residual cap bending energy alone")
    bendCapSurf = go.Figure()
    bendCapSurf.add_trace(go.Surface(
        colorbar=go.surface.ColorBar(title='Residual cap'),
        colorscale='dense',
        contours={"z":{"show":True,"start":np.max(GCap[r,:,:])/2,
                       "end":np.max(GCap[r,:,:])*0.9,
                       "size":np.max(GCap[r,:,:])/5}},
        name="Residual cap",
        x=rEffCap,y=rt,z=GCap[r,:,:],opacity=0.9))
    bendCapSurf.update_layout(coloraxis_showscale=False,
        scene = dict(
                        xaxis_title='Pore radius (nm)',
                        yaxis_title='Intermembrane dist (nm)',
                        zaxis_title='G (kT)'),
        margin=dict(r=10, b=10, l=10, t=10),
        height=600,width=700)
    bendCapSurf.update_layout(scene_camera=view)
    bendCapSurf.show()
    bendCapSurf.write_image(f'output data/Cap hemisphere Gbend R0 {R0[r]}.pdf')
    bendCapSurf.write_image(f'output data/Cap hemisphere Gbend R0 {R0[r]}.svg')
    bendCapSurf.write_image(f'output data/Cap hemisphere Gbend R0 {R0[r]}.png')
    
    
    # ---------------------------------------------------------------------
    # "Ground" G values for the cap pore
    # ---------------------------------------------------------------------
    # First find the minimum value on GTotCap for rEffCap>0
    # To reduce the search space, we'll take advantage of the previous computation
    # of the minimum G path for pore expansion; we'll only search for the min w/i
    # GMinPathCap

    # First find the minimum energy for the cap
    #capMinG = np.amin(GMinPathCap[1:,2])
    #print('capMinG = ',capMinG)
    capGround = GMinPathCap[r,1,2]
    print('capGround = ',capGround)

    # Next "ground" the min path and the surface by that val
    groundGMinPathCap = GMinPathCap
    groundGMinPathCap[r,:,2] = groundGMinPathCap[r,:,2]-capGround
    
    groundGTotCap[r,:,:] = GTotCap[r,:,:]-capGround
#     print('GTotCap = ',GTotCap)
#     print('groundGTotCap = ',groundGTotCap)


print('------------------------------------------------------------')
# Use same sim params as for the cap calculation
# ---------------------------------------------------------------------
# Set-up: flat membrane vars
# ---------------------------------------------------------------------
# Set up the array of size rtSteps x rEffSteps to hold my computations
# Note: Since assuming H0 = 0, the flat membrane has no deformation energy and will thus
# be neglected in these calculations
GPoreFlat = np.zeros((rtSteps,rEffSteps)) # Gb of pore in flat membrane; kT
rEffFlat = np.linspace(0,max_r,num=rEffSteps) # Effective radius of pore; nm
    # This is the pore radius at the narrowest point. Since theta == 0, rEffFlat = r0
    # Note: linspace is start and end point inclusive by default. 
    # Thus there are rEffSteps-2 evenly spaced points in the open interval (0,max_r)
GMinPathFlat = np.zeros((rEffSteps,3)) # Coordinates for the min G path
    # Dim 1 (length rEffSteps) = index moving down rEffFlat
    # Dim 2 (length 3) = spatial coordinates for the min in that slice
        # 0 = x (ie rEffFlat), 1 = y (ie rt), 2 = z (GPoreFlat) 
GMaxPathFlat = np.zeros((rEffSteps,3)) # Coordinates for the max G path
    # Dim 1 (length rEffSteps) = index moving down rt
    # Dim 2 (length 3) = spatial coordinates for the min in that slice
        # 0 = x (ie rEffFlat), 1 = y (ie rt), 2 = z (GPoreFlat) 

# ---------------------------------------------------------------------
# Set initial conditions for referencing in computations
# ---------------------------------------------------------------------
# Note that array[:,n] refs calumn n in all rows in array
GPoreFlat[:,0] = 0

# ---------------------------------------------------------------------
# Find GPoreFlat for each value of rt and rPoreFlat
# ---------------------------------------------------------------------
# First iterate over rPoreFlat then rt; nesting order doesn't matter since everything uses both
print("Computing energies for FLAT pore")
for i in range(1,rEffSteps):
    for rtIndex in range(0,rtSteps):
        #Vars for each iteration of loop
        theta = 0 # Angle of pore incline; radians
            # This is constant and 0 since flat membrane

        # Find Gb of the flat pore
        GPoreFlat[rtIndex,i] = GPore(rEffFlat[i],rt[rtIndex],m[rtIndex],theta,KbPore)

print("Max GPoreFlat = ",np.max(GPoreFlat))

# ---------------------------------------------------------------------
# Find MINimal G expansion path for flat pore
# ---------------------------------------------------------------------
GMinPathFlat = surfMinPath(rEffFlat,rt,GPoreFlat)

# ---------------------------------------------------------------------
# Find MAXimal G expansion path
# ---------------------------------------------------------------------
# This is so that I can...
#     a) Check that there is a peak past which expansion becomes spontaneous as in 
#         the Chizmadzhev paper
#     b) Get a sense of the pore size at which expansion becomes spontaneous
#     c) Get a sense of the energetic barrier that must be overcome. This won't 
#           necessarily be the actual value of the activation energy, but it will 
#         probably be in the ballpark since the surface isn't super curvey
#     Since there's that huge upsweep at low rt, I'm going to maximize with respect to rt rather
#     than with respect to rEffFlat, the approach used for finding the min path
GMaxPathFlat = surfMaxPath(rEffFlat,rt,GPoreFlat)

# ---------------------------------------------------------------------
# Graph the results for the flat pore
# ---------------------------------------------------------------------
print("Graphing results")

# Graph GPoreFlat surface
flatSurf = go.Figure()
print("Graphing total surface")
flatSurf.add_trace(go.Surface(
    #coloraxis=go.layout.Coloraxis(dict(
    #    title="Neg G", "coloraxis1")),
    colorbar=go.surface.ColorBar(title='Net G'),
    colorscale='matter',
    contours={"z":{"show":True,"start":np.max(GPoreFlat)/20,"end":np.max(GPoreFlat)*0.95,
                   "size":np.max(GPoreFlat)/20}},
    name='Net G',
    x=rEffFlat,y=rt,z=GPoreFlat,opacity=0.9))
# Graph curve of min GPoreFlat for expanding rEffFlat
print("Graphing minimal G path")
flatSurf.add_trace(go.Scatter3d(
    marker=dict(
        color='blue',
        size=2),
    mode='lines',
    name='Min G path',
    x=GMinPathFlat[:,0],y=GMinPathFlat[:,1],z=GMinPathFlat[:,2]))
# Graph curve of max GPoreFlat for expanding rEffFlat
print("Graphing maximal G path")
flatSurf.add_trace(go.Scatter3d(
    marker=dict(
        color='red',
        size=2),
    mode='lines',
    name='Max G path',
    x=GMaxPathFlat[:,0],y=GMaxPathFlat[:,1],z=GMaxPathFlat[:,2]))
flatSurf.update_layout(coloraxis_showscale=False,
    scene = dict(
                    xaxis_title='Pore radius (nm)',
                    yaxis_title='Intermembrane dist (nm)',
                    zaxis_title='G (kT)'),
    legend_orientation="h", # Make legend horizontal
    margin=dict(r=10, b=10, l=10, t=10),
    height=600,width=700)
flatSurf.update_layout(scene_camera=view)
flatSurf.show()    
flatSurf.write_image(f'output data/Flat Gtot.pdf')
flatSurf.write_image(f'output data/Flat Gtot.svg')
flatSurf.write_image(f'output data/Flat Gtot.png')

# ---------------------------------------------------------------------
# "Ground" G values for the flat pore
# ---------------------------------------------------------------------
# Follow same procedure as for the cap pore

# restrGMinPathFlat = np.zeros(rEffSteps)
# restrGMinPathFlat = GMinPathFlat[:,2]

# restrGMinPathFlat = np.delete(restrGMinPathFlat,0)

# First find the minimum energy for the flat pore
# flatMinG = np.amin(restrGMinPathFlat)
#flatMinG = np.amin(GMinPathFlat[1:,2])
#print('flatMinG = ',flatMinG)
flatGround = GMinPathFlat[1,2]
print('flatGround = ',flatGround)

# Next "ground" the min path and the surface by that val
#groundGMinPathFlat = GMinPathFlat[1:,2]-flatMinG
groundGMinPathFlat = GMinPathFlat
groundGMinPathFlat[:,2] = groundGMinPathFlat[:,2]-flatGround
# print('GMinPathCap[:,2] = ',GMinPathCap[:,2])
# print('groundGMinPathCap[:,2] = ',groundGMinPathCap[:,2])

groundGPoreFlat = GPoreFlat[:,:]-flatGround
print('GPoreFlat = ',GPoreFlat)
print('groundGPoreFlat = ',groundGPoreFlat)


# ---------------------------------------------------------------------
# COMPARISONS
# ---------------------------------------------------------------------
# Graph "grounded" energy surfaces for cap pore vs flat pore
# ---------------------------------------------------------------------
compSurf = go.Figure()
for r in range(0,R0.size):
    print("Graphing cap pore grounded tot G surface for R0 = ",R0[r])
    print('groundGTotCap[r,:,:] shape: ',groundGTotCap[r,:,:].shape)
    compSurf.add_trace(go.Surface(
        colorbar=go.surface.ColorBar(title='Net G'),
        colorscale='matter',
        contours={"z":{"show":True,"start":np.max(groundGTotCap[r,:,:])/20,
                       "end":np.max(groundGTotCap[r,:,:])*0.95,
                       "size":np.max(groundGTotCap[r,:,:])/20}},
        name=f'Net G, cap, R0 = {R0[r]}', #use f-string for form style insertion of data
        x=rEffCap,y=rt,z=groundGTotCap[r,:,:],opacity=0.95))
    print("Graphing min G path for cap pore")
    compSurf.add_trace(go.Scatter3d(
        marker=dict(
            color='red',
            size=2),
        mode='lines',
        name=f'Min G path, cap, R0 = {R0[r]}',
        x=groundGMinPathCap[r,1:,0],y=groundGMinPathCap[r,1:,1],z=groundGMinPathCap[r,1:,2]))
    
print("Graphing flat pore energy surface")
compSurf.add_trace(go.Surface(
    colorbar=go.surface.ColorBar(title='Net G'),
    colorscale='dense',
    contours={"z":{"show":True,"start":np.max(groundGPoreFlat)/20,"end":np.max(groundGPoreFlat)*0.95,
                   "size":np.max(groundGPoreFlat)/20}},
    name='Net G, flat',
    x=rEffFlat,y=rt,z=groundGPoreFlat,opacity=0.95))
print("Graphing min G path for flat pore")
compSurf.add_trace(go.Scatter3d(
    marker=dict(
        color='blue',
        size=2),
    mode='lines',
    name='Min G path, flat',
    x=GMinPathFlat[1:,0],y=GMinPathFlat[1:,1],z=GMinPathFlat[1:,2]))

compSurf.update_layout(coloraxis_showscale=False,
    scene = dict(
                    xaxis_title='Pore radius (nm)',
                    yaxis_title='Intermembrane dist (nm)',
                    zaxis_title='G (kT)'),
    legend_orientation="h", # Make legend horizontal
    margin=dict(r=10, b=10, l=10, t=10),
    height=600,width=700)
compSurf.update_layout(scene_camera=view)
compSurf.show()


# ---------------------------------------------------------------------
# Graph "grounded" energy min G paths cap pore vs flat pore
# ---------------------------------------------------------------------
print("Plotting minimum paths: G vs rEff")
compMinCurves = go.Figure()
# Skip rEff = 0 so that range stays nice for looking at the rest of the graph
compMinCurves.add_trace(go.Scatter(x=rEffCap[1:], y=groundGMinPathFlat[1:,2], 
                        name="Flat pore", mode='lines'))

for r in range(0,R0.size):
    compMinCurves.add_trace(go.Scatter(x=rEffCap[1:], y=groundGMinPathCap[r,1:,2], 
                        name=f'Cap pore, R0 ={R0[r]}', mode='lines'))
# compMinCurves.add_trace(go.Scatter(x=rEffCap[1:], y=GMinPathCap[1:,2], name="Cap pore",
#                         mode='lines'))

compMinCurves.update_layout(title="Grounded minimum G paths",
                        xaxis_title="Pore radius (nm)", yaxis_title="G (kT)")
compMinCurves.show()
compMinCurves.write_image(f'output data/Cap v flat pore - min G_tot path.pdf')
compMinCurves.write_image(f'output data/Cap v flat pore - min G_tot path.svg')
compMinCurves.write_image(f'output data/Cap v flat pore - min G_tot path.png')

# ---------------------------------------------------------------------
# Find indices of rt and rEff boundaries
# ---------------------------------------------------------------------
# numpy.argmax() and numpy.argmin() can take arrays w/ a comparison operator as their 
# input. Eg numpy.argmax(arr>=3) will give the first index in a flattened arr where
# the value >= 3.  See note in markdown above explaining comparison operators in python
minRtIndex = np.argmax(rt>14)
print('rt>14 @ index of',minRtIndex)

minREffIndex = np.argmax(rEffCap<0.5)
print('rEffCap>0.5 @ index of ',minREffIndex)
maxREffIndex = np.argmin(rEffCap<16.5)
print('rEffCap<16.5 up to index ',maxREffIndex)

# ---------------------------------------------------------------------
# Graph "grounded" energy surfaces for cap vs flat - restr domains
# ---------------------------------------------------------------------
compSurfRestr = go.Figure()
for r in range(0,R0.size):
    print("Graphing cap pore grounded total energy surface, R0=",R0[r])
    compSurfRestr.add_trace(go.Surface(
        colorbar=go.surface.ColorBar(title='Net G'),
        colorscale='matter',
        contours={"z":{"show":True,"start":np.max(groundGTotCap)/20,"end":np.max(groundGTotCap)*0.95,
                       "size":np.max(groundGTotCap)/20}},
        name=f'Cap, R0 = {R0[r]}',
        x=rEffCap[1:maxREffIndex],y=rt[minRtIndex:],z=groundGTotCap[r,minRtIndex:,1:maxREffIndex]+5,
        opacity=0.9))
    print("Graphing min G path for cap pore, R0=",R0[r])
    compSurfRestr.add_trace(go.Scatter3d(
        marker=dict(
            color='red',
            size=2),
        mode='lines',
        name=f'Min G path, cap, R0={R0[r]}',
        x=groundGMinPathCap[r,1:maxREffIndex,0],
        y=groundGMinPathCap[r,1:maxREffIndex,1],
        z=groundGMinPathCap[r,1:maxREffIndex,2]+5))

print("Graphing flat pore energy surface")
compSurfRestr.add_trace(go.Surface(
    colorbar=go.surface.ColorBar(title='Net G'),
    colorscale='dense',
    contours={"z":{"show":True,"start":np.max(groundGPoreFlat)/20,"end":np.max(groundGPoreFlat)*0.95,
                   "size":np.max(groundGPoreFlat)/20}},
    name='Flat',
    x=rEffFlat[1:maxREffIndex],y=rt[minRtIndex:],z=groundGPoreFlat[minRtIndex:,1:maxREffIndex],
    opacity=0.9))
print("Graphing min G path for flat pore")
compSurfRestr.add_trace(go.Scatter3d(
    marker=dict(
        color='blue',
        size=2),
    mode='lines',
    name='Min G path, flat',
    x=groundGMinPathFlat[1:maxREffIndex,0],
    y=groundGMinPathFlat[1:maxREffIndex,1],
    z=groundGMinPathFlat[1:maxREffIndex,2]))

compSurfRestr.update_layout(coloraxis_showscale=False,
    scene = dict(
                    xaxis_title='Pore radius (nm)',
                    yaxis_title='Intermembrane dist (nm)',
                    zaxis_title='G (kT)'),
    legend_orientation="h", # Make legend horizontal
    margin=dict(r=10, b=10, l=10, t=10),
    height=600,width=700)
compSurfRestr.update_layout(scene_camera=view)
compSurfRestr.show()
compSurfRestr.write_image(f'output data/Cap v flat pore - G_tot restr domain.pdf')
compSurfRestr.write_image(f'output data/Cap v flat pore - G_tot restr domain.svg')
compSurfRestr.write_image(f'output data/Cap v flat pore - G_tot restr domain.png')

# ---------------------------------------------------------------------
# Graph "grounded" energy min G paths cap vs flat - restr rEff
# ---------------------------------------------------------------------
print("Plotting minimum paths: G vs rEff")
compMinCurvesRestr = go.Figure()
for r in range(0,R0.size):
    compMinCurvesRestr.add_trace(go.Scatter(x=rEffCap[1:maxREffIndex], y=groundGMinPathCap[r,1:maxREffIndex,2], 
                            name=f'Cap pore, R0 = {R0[r]}', mode='lines'))
compMinCurvesRestr.add_trace(go.Scatter(x=rEffCap[1:maxREffIndex], y=groundGMinPathFlat[1:maxREffIndex,2], 
                        name="Flat pore", mode='lines'))
compMinCurvesRestr.update_layout(title="Grounded minimum G paths",
                        xaxis_title="Pore radius (nm)", yaxis_title="G (kT)")
compMinCurvesRestr.show()
compMinCurvesRestr.write_image(f'output data/Cap vs flat pore - min G_tot path restr domain.pdf')
compMinCurvesRestr.write_image(f'output data/Cap vs flat pore - min G_tot path restr domain.svg')
compMinCurvesRestr.write_image(f'output data/Cap vs flat pore - min G_tot path restr domain.png')