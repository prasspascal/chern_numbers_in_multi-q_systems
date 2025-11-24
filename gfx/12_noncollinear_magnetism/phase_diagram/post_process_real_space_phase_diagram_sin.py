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
mpl.rcParams['font.family'] = 'cmu serif'

fs = 16

phasediagram = np.load("./data/real_space_phase_diagram_sin.npy")

plt.imshow(phasediagram.transpose(), origin='lower', aspect='auto',
           extent=(-2/np.sqrt(3),2/np.sqrt(3), 0, 2*np.pi), interpolation='bicubic', cmap='RdBu')

cbar = plt.colorbar()
cbar.set_label(label=r'$\mathrm{deg}~\hat{\mathbf{n}}$', size=fs)
cbar.ax.tick_params(labelsize=14)

# plt.contour(phasediagram.transpose(), levels=[-2,-1,0,1,2],origin='lower',
#            extent=(-2, 2, 0, 2*np.pi), colors=['black', 'black', 'black','black','black'])

eps = 1e-10

n_m = 300
for delta in [-4*np.pi,-2*np.pi, 0, 2*np.pi,4*np.pi]:
    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( 3*(np.pi + np.arccos(m))) +delta)
    plt.plot(ms,phase,'-', color='black')


    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( 3*(np.pi - np.arccos(m)))  +delta)
    plt.plot(ms,phase,'-', color='black')

    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( (3*np.pi + np.arccos(m)))  +delta)
    plt.plot(ms,phase,'-', color='black')

    ms = np.linspace(-1+eps, 1-eps,n_m)
    phase = []
    for m in ms:
        phase.append( ( (3*np.pi - np.arccos(m)))  +delta)
    plt.plot(ms,phase,'-', color='black')

label_fs = 14

# -- deg = -1

plt.text(0.2, np.pi/2-0.2, r'$\mathrm{deg}\hat{\ \mathbf{n}}=\mp 1$',{'fontsize': label_fs})
plt.text(0.2, 3*np.pi/2 , r'$\mathrm{deg}\hat{\ \mathbf{n}}=\mp 1$',{'fontsize': label_fs})

# -- deg = 0

plt.text(-1.1, np.pi - 0.1, r'$\mathrm{deg}\hat{\ \mathbf{n}}=0$',{'fontsize': label_fs})
plt.text(0.65, 2*np.pi - 0.4, r'$\mathrm{deg}\hat{\ \mathbf{n}}=0$',{'fontsize': label_fs})
plt.text(0.65, 0.2, r'$\mathrm{deg}\hat{\ \mathbf{n}}=0$',{'fontsize': label_fs})

# -- deg = 1

plt.text(-0.75, np.pi/2+0.1 , r'$\mathrm{deg}\hat{\ \mathbf{n}}=\pm 1$',{'fontsize': label_fs})
plt.text(-0.75, 3*np.pi/2-0.3 , r'$\mathrm{deg}\hat{\ \mathbf{n}}=\pm 1$',{'fontsize': label_fs})

# -- deg = -2

plt.text(-0.4, 2*np.pi-0.6, r'$\mathrm{deg}\hat{\ \mathbf{n}}=\mp 2$',{'fontsize': label_fs, 'color': 'lightgray'})

plt.text(-0.4,  0.4, r'$\mathrm{deg}\hat{\ \mathbf{n}}=\mp 2$',{'fontsize': label_fs, 'color': 'lightgray'})

# -- deg = 2

plt.text(-0.1, np.pi-0.1, r'$\mathrm{deg}\hat{\ \mathbf{n}}=\pm 2$',{'fontsize': label_fs, 'color': 'lightgray'})

# ticks and labels

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
plt.savefig("./plots/real_space_phase_diagram_sin.png", dpi=300)
plt.clf()
