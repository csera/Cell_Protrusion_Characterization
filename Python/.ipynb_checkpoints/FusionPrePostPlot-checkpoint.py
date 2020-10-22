# Curtis Sera
# v1.0 (for v2.0 of data generating scripts)
# 2020-10-21

import numpy as np
import plotly.graph_objs as go

capDGb = np.genfromtxt("output data/Cap DGb 2-0.csv",delimiter=',')
cap_rt = np.genfromtxt("output data/Cap rt 2-0.csv", delimiter=',')

bodyDGb = np.genfromtxt("output data/Body DGb 2-0.csv",delimiter=',')
body_rt = np.genfromtxt("output data/Body rt 2-0.csv", delimiter=',')

baseDGb = np.genfromtxt("output data/Base DGb 2-0.csv",delimiter=',')
base_rt = np.genfromtxt("output data/Base rt 2-0.csv", delimiter=',')

print('base_rt = ',base_rt)
print('baseDGb = ',baseDGb)
print('bodyDGb = ',bodyDGb)

print('Plotting energies')

fig = go.Figure()
fig.add_trace(go.Scatter(x=cap_rt, y=capDGb, mode='markers + lines', name='Cap'))
fig.add_trace(go.Scatter(x=body_rt, y=bodyDGb, mode='markers + lines', name='Body'))
fig.add_trace(go.Scatter(x=base_rt, y=baseDGb, mode='markers + lines', name='Base'))
fig.update_layout(title="$Compiled: \\quad \Delta G_b = G_b^{post} - G_b^{pre}$",
                        xaxis_title="$r_{t,rim} \\, (nm)$", yaxis_title="Energy (kT)")
fig.show()

logFig = go.Figure()
logFig.add_trace(go.Scatter(x=cap_rt, y=capDGb, mode='markers + lines', name='Cap'))
logFig.add_trace(go.Scatter(x=body_rt, y=bodyDGb, mode='markers + lines', name='Body'))
logFig.add_trace(go.Scatter(x=base_rt, y=baseDGb, mode='markers + lines', name='Base'))
logFig.update_layout(title="$Compiled: \\quad \Delta G_b = G_b^{post} - G_b^{pre}$",
                        xaxis_title="$r_{t,rim} \\, (nm)$", yaxis_title="Energy (kT)")
logFig.update_yaxes(type='log')
logFig.show()

ratlFig = go.Figure()
ratlFig.add_trace(go.Scatter(x=cap_rt, y=1/capDGb, mode='markers + lines', name='Cap'))
ratlFig.add_trace(go.Scatter(x=body_rt, y=1/bodyDGb, mode='markers + lines', name='Body'))
ratlFig.add_trace(go.Scatter(x=base_rt, y=1/baseDGb, mode='markers + lines', name='Base'))
ratlFig.update_layout(title="$Compiled: \\quad \\frac{1}{\Delta G_b} = \\frac{1}{G_b^{post} - G_b^{pre}}$",
                        xaxis_title="$r_{t,rim} \\, (nm)$", yaxis_title="$Energy^{-1} \\, (kT^{-1})$")
ratlFig.show()

ratlFigSqrt2 = go.Figure()
ratlFigSqrt2.add_trace(go.Scatter(x=cap_rt, y=1/(capDGb**2), mode='markers + lines', name='Cap'))
ratlFigSqrt2.add_trace(go.Scatter(x=body_rt, y=1/(bodyDGb**2), mode='markers + lines', name='Body'))
ratlFigSqrt2.add_trace(go.Scatter(x=base_rt, y=1/(baseDGb**2), mode='markers + lines', name='Base'))
ratlFigSqrt2.update_layout(title="$Compiled: \\quad \\frac{1}{(\Delta G_b)^{2}} = \\frac{1}{(G_b^{post} - G_b^{pre})^{2}}$",
                        xaxis_title="$r_{t,rim} \\, (nm)$", yaxis_title="$Energy^{-2} \\, (kT^{-2})$")
ratlFigSqrt2.show()

ratlFigSqrt = go.Figure()
ratlFigSqrt.add_trace(go.Scatter(x=cap_rt, y=1/(capDGb**0.5), mode='markers + lines', name='Cap'))
ratlFigSqrt.add_trace(go.Scatter(x=body_rt, y=1/(bodyDGb**0.5), mode='markers + lines', name='Body'))
ratlFigSqrt.add_trace(go.Scatter(x=base_rt, y=1/(baseDGb**0.5), mode='markers + lines', name='Base'))
ratlFigSqrt.update_layout(title="$Compiled: \\quad \\frac{1}{(\Delta G_b)^{0.5}} = \\frac{1}{(G_b^{post} - G_b^{pre})^{0.5}}$",
                        xaxis_title="$r_{t,rim} \\, (nm)$", yaxis_title="$Energy^{-0.5} \\, (kT^{-0.5})$")
ratlFigSqrt.show()