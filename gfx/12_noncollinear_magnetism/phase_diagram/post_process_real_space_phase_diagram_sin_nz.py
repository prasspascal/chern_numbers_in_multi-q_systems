import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot parameter
mpl.pyplot.rcdefaults()
#plt.rcParams['figure.figsize'] = [20, 10]
plt.rcParams['savefig.facecolor'] = "white"
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['text.usetex'] = True
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'CMU serif'

fs = 16

phasediagram = np.load("./data/real_space_phase_diagram_sin_nz.npy")

plt.imshow(phasediagram.transpose(), origin='lower', aspect='auto',
           extent=(-2/np.sqrt(3),2/np.sqrt(3), 0, 2*np.pi), interpolation='nearest', cmap='RdBu')

cbar = plt.colorbar()
cbar.set_label(label=r'$\langle{\hat{n}_z}\rangle$', size=fs)
cbar.ax.tick_params(labelsize=14)

eps = 1e-10

n_m = 300

for delta in [-4*np.pi,-2*np.pi, 0, 2*np.pi,4*np.pi]:
    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( 3*(np.pi + np.arccos(m))) +delta)
    plt.plot(ms,phase,'--', color='lightgray')


    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( 3*(np.pi - np.arccos(m)))  +delta)
    plt.plot(ms,phase,'--', color='lightgray')

    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( (3*np.pi + np.arccos(m)))  +delta)
    plt.plot(ms,phase,'--', color='lightgray')

    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( (3*np.pi - np.arccos(m)))  +delta)
    plt.plot(ms,phase,'--', color='lightgray')

plt.contour(phasediagram.transpose(), levels=[-2,0,2],origin='lower',
           extent=(-2, 2, 0, 2*np.pi), colors=['black', 'black', 'black'])

plt.text(0.05, np.pi - 0.1, r'$\langle{\hat{n}_z}\rangle=0$',{'fontsize': 12})

xticks = [-1.0,-0.5,0.0,0.5,1.0]
xticks_labels = [r'$-1.0$',r'$-0.5$',r'$0.0$',r'$0.5$',r'$1.0$']
plt.xticks(xticks, xticks_labels, fontsize=14)

yticks = [0, np.pi, 2*np.pi]
yticks_labels = [r'$0$',r'$\pi$',r'$2\pi$']
plt.yticks(yticks, yticks_labels, fontsize=16)

ax = plt.gca()
ax.set_xlabel(r"$m$", fontsize=fs)
ax.set_ylabel(r"$\varphi$", fontsize=fs)
ax.set_xlim((-2/np.sqrt(3),2/np.sqrt(3)))
ax.set_ylim((0,2*np.pi))

plt.tight_layout()
plt.savefig("./plots/real_space_phase_diagram_sin_nz.png", dpi=300)
plt.clf()
