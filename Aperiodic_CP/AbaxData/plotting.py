import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace, log, mean
from os.path import isfile
from scipy.optimize import curve_fit

def func(x, a, b):
    return a * x + b

fig, axes = plt.subplots(1, figsize=(10,6))
axins = inset_axes(axes, width="50%", height="50%",loc=3, borderpad = 0.1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 20) 

k = 1
with open (f'./k={k}/last_params.txt', 'r') as f:
    for line in f.read().splitlines():
        param, value = line.split(',')
        if param == 'tmax':
            tmax = float(value)
        elif param == 'steps':
            steps = int(value)
        elif param == 'size':
            size = float(value)
        elif param == 'la':
            lamb_a = float(value)
        else:
            pass

param_string = '_t{:.1e}_la{:.2f}_s{}'.format(tmax, lamb_a, int(size))

if k == 1:
    colors = {
            'sup' : '#d81900',
            'inf' : '#af9e7c',
            'crit' : '#333333'
            }
elif k == 2:
    colors = {
            'sup' : '#1b59c4',
            'inf' : '#bf9c4a',
            'crit' : '#333333'
            }
elif k == 3:
    colors = {
            'sup' : '#307f1d',
            'inf' : '#9e6c4d',
            'crit' : '#333333'
            }
else:
    colors = {
            'sup' : '#f26419',
            'inf' : '#2f74a3',
            'crit' : '#333333'
            }

alphas = [ x for x in linspace(0.4, 1, steps) ]

# Main plot
with open(f'./k={k}/plot_targets{param_string}.txt', 'r') as f:
    for line in f.read().splitlines():
        step, file, stat = line.split(',')
        step = int(step)
        try:
            dfile = open(file, 'r')
            ts = []
            rs = []
            for line in dfile.read().splitlines():
                t, rho, dt, drho = line.split(',')
                ts.append(float(t))
                rs.append(float(rho))
            dfile.close()
            if stat != 'crit':
                if stat == 'active':
                    axes.loglog(ts, rs, c = colors['sup'], alpha = alphas[step])
                else:
                    axes.loglog(ts, rs, c = colors['inf'], alpha = alphas[step])
            else:
                axes.loglog(ts, rs, lw = 3, c = colors['crit'])
        except FileNotFoundError:
            pass

axes.set_xlim(1, tmax)
axes.set_ylim(1e-1,1)
axes.set_xlabel(r't', fontsize = 20)
axes.set_ylabel(r'$\lambda$', fontsize = 20)

# Inset plot
Infs = []
Sups = []

with open(f'./k={k}/lambdas_hist{param_string}.dat', 'r') as f:
    for line in f.read().splitlines():
        linf, lsup = line.split(',')
        Infs.append(float(linf))
        Sups.append(float(lsup))

Steps = [ s for s in range(steps) ]
axins.plot(
        Steps, Sups,
        lw = 2, ls = '--', c = colors['sup'],
        marker = 'o', label = r"$\lambda_{sup}$")

axins.plot(Steps, Infs,
        lw = 2,ls = '--', c = colors['inf'],
        marker = 'o', label = r"$\lambda_{inf}$")

Crits = [(s+i)/2 for s, i in list(zip(Sups,Infs)) ]
axins.plot(Steps, Crits,
        lw = 3, ls = '--', c = colors['crit'],
        marker = 'o', label = r"$\lambda_{crit}$")

lcrit = Crits[-1]
custom_lines = []
for clr in colors.values():
    custom_lines.append(Line2D([0], [0], color=clr, lw=4))
axes.legend(custom_lines, [r"$\lambda_{sup}$", r"$\lambda_{inf}$", r"$\lambda^* = $" + r"${:.4f}$".format(lcrit)], fontsize = 13)

#axins.hlines(y = 3.29785, xmin = 0, xmax = len(Steps) - 1,
        #ls='--',alpha=0.5, color = colors['crit'], lw = 1.)

axins.set_xlabel(r"$j$", fontsize = 13)
axins.xaxis.set_label_position("top")
axins.set_ylabel(r"$\lambda_j$", fontsize = 13)
axins.yaxis.set_label_position("right")
axins.legend(fontsize = 13)

axins.tick_params(
        labelsize = 10,
        left=False, right=True,
        bottom=False, top=True,
        labelleft=False, labelright=True,
        labeltop = True, labelbottom = False
        )

fig.tight_layout()

outfile = f'./Plots/acp_{param_string}.pdf'
if isfile(outfile):
    ans = input('> File already exists, do you with to overwrite it?(yes/no)')
    if ans == 'yes':
        plt.savefig(outfile)
    else:
        pass
else:
    plt.savefig(outfile)

plt.show()
