# Curtis Sera
# v2.2 for 
#    v2.2.x of CapFusionPrePost, BodyFusionPrePost, and BaseFusionPrePost
# Using *MATPLOTLIB*, not Plotly since req LaTeX
# 2020-10-21

import numpy as np
import matplotlib.pyplot as pt

capDGb = np.genfromtxt("output data/Cap DGb 2-2.csv",delimiter=',')
cap_rt = np.genfromtxt("output data/Cap rt 2-2.csv", delimiter=',')

bodyDGb = np.genfromtxt("output data/Body DGb 2-2.csv",delimiter=',')
body_rt = np.genfromtxt("output data/Body rt 2-2.csv", delimiter=',')

baseDGb = np.genfromtxt("output data/Base min path 2-2.csv",delimiter=',')
base_rt = np.genfromtxt("output data/Base rt 2-2.csv", delimiter=',')

print('base_rt = ',base_rt)
print('baseDGb = ',baseDGb)
print('bodyDGb = ',bodyDGb)

print('Plotting energies')

print('making linear plot')
linFig = pt.figure()
linSP = linFig.add_subplot(111)
linSP.plot(cap_rt, capDGb, label='Cap')
linSP.plot(body_rt, bodyDGb, label='Body')
linSP.plot(base_rt, baseDGb, label='Base')
linSP.set_xlabel('$r_{t,protr} \\, (nm)$')
linSP.set_ylabel('$Energy \\, (kT)$')
linSP.set_title('$Compiled: \\quad \\Delta G_b = G_b^{post} - G_b^{pre}$')
linSP.legend()

print('making log-linear (semilog) plot')
logFig = pt.figure()
logSP = logFig.add_subplot(111)
logSP.plot(cap_rt, capDGb, label='Cap')
logSP.plot(body_rt, bodyDGb, label='Body')
logSP.plot(base_rt, baseDGb, label='Base')
logSP.set_xlabel('$r_{t,protr} \\, (nm)$')
logSP.set_ylabel('$Energy \\, (kT)$')
logSP.set_yscale('log')
logSP.set_title('$Compiled \\, (semilog\, scale): \\quad \\Delta G_b = G_b^{post} - G_b^{pre}$')
logSP.legend()

print("making reciprocal (ie rational) plot")
ratlFig = pt.figure()
ratlSP = ratlFig.add_subplot(111)
ratlSP.plot(cap_rt, 1/capDGb, label='Cap')
ratlSP.plot(body_rt, 1/bodyDGb, label='Body')
ratlSP.plot(base_rt, 1/baseDGb, label='Base')
ratlSP.set_xlabel('$r_{t,protr} \\, (nm)$')
ratlSP.set_ylabel('$Energy^{-1} \\, (kT^{-1})$')
ratlSP.set_title('$Compiled: \\quad \\frac{1}{\Delta G_b} = \\frac{1}{G_b^{post} - G_b^{pre}}$')
ratlSP.legend()

print('making reciprocal sqr')
ratlSqrFig = pt.figure()
ratlSqrSP = ratlSqrFig.add_subplot(111)
ratlSqrSP.plot(cap_rt, 1/(capDGb**2), label='Cap')
ratlSqrSP.plot(body_rt, 1/(bodyDGb**2), label='Body')
ratlSqrSP.plot(base_rt, 1/(baseDGb**2), label='Base')
ratlSqrSP.set_xlabel('$r_{t,protr} \\, (nm)$')
ratlSqrSP.set_ylabel('$Energy^{-2} \\, (kT^{-2})$')
ratlSqrSP.set_title('$Compiled: \\quad \\frac{1}{\Delta G_b^2} = \\frac{1}{(G_b^{post} - G_b^{pre})^2}$')
ratlSqrSP.legend()

print('making reciprocal sqrt')
ratlSqrtFig = pt.figure()
ratlSqrtSP = ratlSqrtFig.add_subplot(111)
ratlSqrtSP.plot(cap_rt, 1/(capDGb**0.5), label='Cap')
ratlSqrtSP.plot(body_rt, 1/(bodyDGb**0.5), label='Body')
ratlSqrtSP.plot(base_rt, 1/(baseDGb**0.5), label='Base')
ratlSqrtSP.set_xlabel('$r_{t,protr} \\, (nm)$')
ratlSqrtSP.set_ylabel('$Energy^{-1/2} \\, (kT^{-1/2})$')
ratlSqrtSP.set_title('$Compiled: \\quad \\frac{1}{\Delta G_b^{1/2}} = \\frac{1}{(G_b^{post} - G_b^{pre})^{1/2}}$')
ratlSqrtSP.legend()

pt.show()

# Save plots
print('saving plots as PDFs')
linFig.savefig("output data/Compiled DGb linear 2-2.pdf")
logFig.savefig("output data/Compiled DGb semilog 2-2.pdf")
ratlFig.savefig("output data/Compiled DGb rational 2-2.pdf")
ratlSqrFig.savefig("output data/Compiled DGb ratl-sqr 2-2.pdf")
ratlSqrtFig.savefig("output data/Compiled DGb ratlsqrt 2-2.pdf")
