# Curtis Sera
# 2020-09-14, v1.0
# 
# Adaptation of my Mathematica "Membrane approach 3.nb" to something that is 
# open source and has easier graphics editing options

import numpy as np
import plotly.graph_objs as go

# First, define functions
def hydrationP(P0,c,d):
    '''
    Hydration repulsion pressure
    
    Args:
        P0 (float): Hydration P at zero separation; Pa
        c (float): hydration decay constant; m
        d (float): intermembrane separation distance; m
    Returns:
        float: The hydration pressure in Pa at distance d
    '''
    return P0 * np.exp(-d/c)

def randUndulationP(kT,P0,Kb,a,d):
    '''
    Repulsive undulation pressure according to Rand's equation
    
    Args:
        kT (float): k_b * T in Kelvin; J/kT
        P0 (float): Hydration P at zero separation; Pa
        Kb (float): Helfrich bending constant; kT
        a (float): Undulatory decay constant; m
        d (float): intermembrane separation distance; m
    Returns:
        float: The undulation pressure in Pa at distance d
    '''
    return ((np.pi*kT)/(32*a))*np.sqrt(P0/(Kb*kT*a))*np.exp(-d/(2*a))

def vdwP(kT,H,t,d):
    '''
    Van der Waals interaction pressure
    
    Args:
        kT (float): k_b * T in Kelvin; J/kT
        H (float): Haymaker constant; kT
        t (float): lipid bilayer thickness; m
        d (float): intermembrane separation distance; m
    Returns:
        float: The van der Waals interaction pressure at distance d
    '''
    return (-(H*kT)/(6*np.pi))*((1/d**3)-(2/(d+t)**3)+(1/d+2*t)**3)

# Define system constants
kT = 4.28e-21 # Avg thermal E @ 37C; J
P0 = 1e12 # Hydration pressure; Pa
c = 0.2e-9 # Hydration decay constant; m
a = 0.2e-9 # Undulation decay constant; m
Kb = 15 # Helfrich bending constant; kT
H = 1 # Haymaker constant; kT
t = 5e-9 # lipid bilayer thickness; m

# Define computational params
dMax = 15e-9 # max intermembrane separation (d) to check; m
n = 1000 # num of points to evaluate the functions for

dSpan = np.linspace(0,dMax,num=n+1)
dSpan = np.delete(dSpan,0)

Ph = np.zeros(n) # Hydration pressure; Pa
Pu = np.zeros(n) # Undulation pressure; Pa
Pv = np.zeros(n) # Van der Waals pressure; Pa
Ptot = np.zeros(n) # Ph+Pu+Pv; Pa

Pabm = np.zeros(n) # P fom ABM; Pa

for x in range(0,n):
    d = dSpan[x]
    
    Ph[x] = hydrationP(P0,c,d)
    Pu[x] = randUndulationP(kT,P0,Kb,a,d)
    Pv[x] = vdwP(kT,H,t,d)
    
    Ptot[x] = Ph[x]+Pu[x]+Pv[x]
    
    Pabm[x] = 3.565e5 # Upper bound est done by hand on 2019-06-11; Pa

# Plot over whole domain (0,dMax]
plotDflt = go.Figure()
plotDflt.add_trace(go.Scatter(x=dSpan, y=Ph, name='Hydration', 
                                  line=dict(color='skyblue', dash='dash'),
                                  mode='lines'))
plotDflt.add_trace(go.Scatter(x=dSpan, y=Pu, name='Undulation', 
                                  line=dict(color='salmon', dash='dash'),
                                  mode='lines'))
plotDflt.add_trace(go.Scatter(x=dSpan, y=Pv, name='Van der Waals',
                                  line=dict(color='darkorchid', dash='dash'),
                                  mode='lines'))
plotDflt.add_trace(go.Scatter(x=dSpan, y=Ptot, name='Net interaction',
                                  line=dict(color='blue'),
                                  mode='lines'))
plotDflt.update_layout(title="Interactions between pure lipid bilayers",
                       # Set xaxis ticks manually to plot in nm rather than m
                       xaxis = dict(
                           tickmode='array',
                           tickvals=[0,2e-9,4e-9,6e-9,8e-9,10e-9,12e-9,14e-9],
                           ticktext=[0,2,4,6,8,10,12,14]
                       ),
                       xaxis_title="d (nm)", yaxis_title="Pressure (Pa)")
plotDflt.show()

# Plot over a restricted range to better see behavior
plotCeil = go.Figure()
plotCeil.add_trace(go.Scatter(x=dSpan, y=Ph, name='Hydration', 
                                  line=dict(color='skyblue', dash='dash'),
                                  mode='lines'))
plotCeil.add_trace(go.Scatter(x=dSpan, y=Pu, name='Undulation', 
                                  line=dict(color='salmon', dash='dash'),
                                  mode='lines'))
plotCeil.add_trace(go.Scatter(x=dSpan, y=Pv, name='Van der Waals',
                                  line=dict(color='orchid', dash='dash'),
                                  mode='lines'))
plotCeil.add_trace(go.Scatter(x=dSpan, y=Ptot, name='Net interaction',
                                  line=dict(color='blue'),
                                  mode='lines'))
plotCeil.update_layout(title="Interactions between pure lipid bilayers",
                       xaxis = dict(
                           tickmode='array',
                           tickvals=[0,2e-9,4e-9,6e-9,8e-9,10e-9,12e-9,14e-9],
                           ticktext=[0,2,4,6,8,10,12,14]
                       ),
                       xaxis_title="d (m)", yaxis_title="Pressure (Pa)")
plotCeil.update_yaxes(range=[-3000,4000])
plotCeil.show()

# Plot over a restricted range vs ABM pressure
plotCeilVAbm = go.Figure()
plotCeilVAbm.add_trace(go.Scatter(x=dSpan, y=Ph, name='Hydration', 
                                  line=dict(color='skyblue', dash='dash'),
                                  mode='lines'))
plotCeilVAbm.add_trace(go.Scatter(x=dSpan, y=Pu, name='Undulation', 
                                  line=dict(color='salmon', dash='dash'),
                                  mode='lines'))
plotCeilVAbm.add_trace(go.Scatter(x=dSpan, y=Pv, name='Van der Waals',
                                  line=dict(color='orchid', dash='dash'),
                                  mode='lines'))
plotCeilVAbm.add_trace(go.Scatter(x=dSpan, y=Ptot, name='Net interaction',
                                  line=dict(color='blue'),
                                  mode='lines'))
plotCeilVAbm.add_trace(go.Scatter(x=dSpan, y=Pabm, name='ABM',
                                  line=dict(color='gold'),
                                  mode='lines'))
plotCeilVAbm.update_layout(title="Bilayer Interactions vs ABM",
                       xaxis = dict(
                           tickmode='array',
                           tickvals=[0,2e-9,4e-9,6e-9,8e-9,10e-9,12e-9,14e-9],
                           ticktext=[0,2,4,6,8,10,12,14]
                       ),
                       xaxis_title="d (m)", yaxis_title="Pressure (Pa)")
plotCeilVAbm.update_yaxes(range=[-3000,4e5])
plotCeilVAbm.show()

# Plot over a restricted range vs ABM pressure isolating the net interaction and ABM
plotCeilVAbmIso = go.Figure()
plotCeilVAbmIso.add_trace(go.Scatter(x=dSpan, y=Ptot, name='Net interaction',
                                  line=dict(color='blue'),
                                  mode='lines'))
plotCeilVAbmIso.add_trace(go.Scatter(x=dSpan, y=Pabm, name='ABM',
                                  line=dict(color='gold'),
                                  mode='lines'))
plotCeilVAbmIso.update_layout(title="Bilayer Interactions vs ABM",
                       xaxis = dict(
                           tickmode='array',
                           tickvals=[0,2e-9,4e-9,6e-9,8e-9,10e-9,12e-9,14e-9],
                           ticktext=[0,2,4,6,8,10,12,14]
                       ),
                       xaxis_title="d (m)", yaxis_title="Pressure (Pa)")
plotCeilVAbmIso.update_yaxes(range=[-3000,4e5])
plotCeilVAbmIso.show()