import subprocess
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from numpy import linspace, log, isnan
from random import randint
from time import sleep

def progress_bar(completed = 0, total = 1, text = '', color = 0, size = 10):
    '''
    INPUT: (Number of) completed tasks, total amount of tasks, text to display inside the bar, color and total size of the bar.
    OUTPUT: String of the progress bar
    '''
    offset = 2
    text = " " * offset + text
    perc = completed / total
    hilight = int( size * perc )
    bar = ''
    for index in range(hilight):
        if index < len(text):
            bar += "\033[48;5;{}m\033[38;5;0m{}\033[0m".format(color, text[index])
        else:
            bar += "\033[48;5;{}m\033[38;5;0m \033[0m".format(color)
    for index in range(hilight, size):
        if index < len(text):
            bar += "\033[48;5;0m\033[38;5;{}m{}\033[0m".format(color, text[index])
        else:
            bar += "\033[48;5;0m\033[38;5;0m \033[0m"
    bar += "  {}%".format(int(perc*100))
    return bar

def func(x, a, b, c):
    return a + (b * x) + c * (x ** 2) 
def Curvature(fname):
    try:
        dfile = open(fname, 'r')
        times = []
        rhos = []
        for line in dfile.read().splitlines():
            t, r, dt, dr = line.split(',')
            if float(t) == 0 or float(r) == 0:
                pass
            else:
                if isnan(float(t)) == False:
                    times.append(float(t))
                    rhos.append(float(r))
        fit_parans, fit_cov = curve_fit(func, log(times), log(rhos))
        #fit_parans = polyfit(log(times), log(rhos), 2)
        return fit_parans[2]
    except:
        print("file {} not found".format(fname))
        raise FileNotFoundError

def lin_func(x, a, b):
    return (a * x) + b
def Slope(fname, start = 0, end = -1):
    try:
        dfile = open(fname, 'r')
        times = []
        rhos = []
        for line in dfile.read().splitlines():
            t, r, dt, dr = line.split(',')
            if float(t) == 0 or float(r) == 0:
                pass
            else:
                if isnan(float(t)) == False:
                    times.append(float(t))
                    rhos.append(float(r))
        fit_parans, fit_cov = curve_fit(
                lin_func,
                log(times[start:end]),
                log(rhos[start:end])
                )
        return fit_parans[0]
    except:
        print("file {} not found".format(fname))
        raise FileNotFoundError

def IsActive(fname, s1 = 200, e1 = 300, s2 = -50, e2 = -1):
    if abs(Slope(fname, s1, e1)) > abs(Slope(fname, s2, e2)):
        return 'active'
    else:
        return 'inactive'

def RunSims(code, analysis, k, rho_0, lamb_a, lamb_b, tsup, size, sim_i, sim_f):
    process = subprocess.Popen(
            [
            code,
            str(k),
            str(rho_0),
            str(lamb_a), 
            str(lamb_b), 
            str(tsup),
            str(size),
            str(sim_i),
            str(sim_f),
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()
    exec_time = float(out)

    process = subprocess.Popen(
            [
            analysis,
            str(k),
            str(lamb_a), 
            str(lamb_b), 
            str(tsup),
            str(size),
            str(0),
            str(sim_f)
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()
    fname = out.decode("utf-8")  

    return (exec_time, fname)

fig, axes = plt.subplots(1, figsize=(10,6))
axins = inset_axes(axes, width="50%", height="50%",loc=3, borderpad = 0.1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 20) 

size = 1000
sims = 100
lamb_a = 2.4
rho_0 = 0.999
tinf = 0

k = 1

if k == 1:
    linf, lsup = (3, 8)
elif k == 2:
    linf, lsup = (3, 6)
elif k == 3:
    linf, lsup = (3, 5)
else:
    linf, lsup = (2, 8)

steps = 15
tmax = 1e4
sim_incr = 100

times = linspace(0, tmax, steps)
flicks = linspace(0.95, 0.99, steps+1)

Sups = []
Infs = []

code = './aperiodic_cp'
analysis = './aperiodic_analysis'
output = open('lambdas_hist.dat', 'w')

last_lsup = linf_reps = 0

lambdas = {}
lambdas[lsup] = 'active'
lambdas[linf] = 'inactive'

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

alphs = ( x for x in linspace(0.4, 1, steps))

for step in range(steps - 1):
    #print("\033[H\033[J")
    print("-" * 40)
    print(progress_bar(
        step,
        steps-1,
        f"Depth {step}",
        color = 2,
        size = 50)
        )
    alpha_level = alphs.__next__()

    if linf_reps == 2:
        linf_reps = 0
        linf *= flicks[step]

    Sups.append(lsup)
    Infs.append(linf)
    output.write(f"{lsup},{linf}\n")
    for lamb_b in [lsup, linf]:
        tsup = times[step + 1]

        if lamb_b == lsup:
            print("Testing superior value \033[38;5;7m{:4f}\033[0m... ".format(lamb_b), end = '')
            clr = colors['sup']
        else:
            print("Testing inferior value \033[38;5;7m{:4f}\033[0m... ".format(lamb_b), end = '')
            clr = colors['inf']

        sim_i = 0
        sim_f = sims + sim_i

        exec_time, fname = RunSims(code,
                analysis,
                k,
                rho_0,
                lamb_a,
                lamb_b,
                tsup,
                size,
                sim_i,
                sim_f
                )

        lambdas[lamb_b] = IsActive(fname)
        if lambdas[lamb_b] == 'active':
            print(f"\033[38;5;2;1m{lambdas[lamb_b].title()}.\033[0m   (Exec. Time = {float(exec_time)})")
        else:
            print(f"\033[38;5;1;1m{lambdas[lamb_b].title()}.\033[0m   (Exec. Time = {float(exec_time)})")

        try:
            dfile = open(fname, 'r')
            ts = []
            rs = []
            for line in dfile.read().splitlines():
                t, rho, dt, drho = line.split(',')
                ts.append(float(t))
                rs.append(float(rho))
            dfile.close()
            axes.loglog(ts, rs, c = clr, alpha = alpha_level)
        except FileNotFoundError:
            pass

    last_lsup, last_linf = lsup, linf

    if lambdas[lsup] == 'active':
        if lambdas[linf] == 'active':
            lsup = (lsup + linf) / 2.
        else:
            linf = (lsup + linf) / 2.
    else:
        lsup = linf = 0
        for lamb_b, status in sorted(lambdas.items(), key= lambda l: l[0]):
            if status == 'active':
                if linf == 0:
                    linf = lamb_b
                else:
                    lsup = lamb_b
                    break
    if linf == last_linf:
        linf_reps += 1
    else:
        linf_reps = 0

output.write(f"{lsup},{linf}\n")
output.close()

#print("\033[H\033[J")
print("-" * 40)
print(progress_bar(
    100,
    100,
    f"Done =D",
    color = 3,
    size = 50)
    )
lcrit = (lsup + linf) / 2
print(f"Critical lambda after {steps} steps: \033[38;5;3;1m{lcrit}\033[0m")

Steps = [ s for s in range(steps - 1) ]
axins.plot(
        Steps, Sups,
        lw = 2, ls = '--', c = colors['sup'],
        marker = 'o', label = r"$\lambda_{sup}$")

axins.plot(Steps, Infs,
        lw = 2,ls = '--', c = colors['inf'],
        marker = 'o', label = r"$\lambda_{inf}$")

axins.plot(Steps, [(s+i)/2 for s, i in list(zip(Sups,Infs)) ],
        lw = 3, ls = '--', c = colors['crit'],
        marker = 'o', label = r"$\lambda_{crit}$")

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

tsup = tmax
lamb_b = lcrit

sim_i = 0
sim_f = sims + sim_i

exec_time, fname = RunSims(code,
        analysis,
        k,
        rho_0,
        lamb_a,
        lamb_b,
        tsup,
        size,
        sim_i,
        sim_f
        )
try:
    dfile = open(fname, 'r')
    ts = []
    rs = []
    for line in dfile.read().splitlines():
        t, rho, dt, drho = line.split(',')
        ts.append(float(t))
        rs.append(float(rho))
    dfile.close()
    axes.plot(ts, rs,lw = 3, c = colors['crit'])
except FileNotFoundError:
    pass

custom_lines = []
for clr in colors.values():
    custom_lines.append(Line2D([0], [0], color=clr, lw=4))
axes.legend(custom_lines, [r"$\lambda_{sup}$", r"$\lambda_{inf}$", r"$\lambda^* = $" + r"${:.4f}$".format(lcrit)], fontsize = 13)

axes.set_xlim(1, tmax)
axes.set_ylim(1e-1,1)
axes.set_xlabel(r't', fontsize = 20)
axes.set_ylabel(r'$\lambda$', fontsize = 20)
fig.tight_layout()
plt.savefig(f'aperiodic_bissection_tmax_{tmax}_k{k}.pdf')
plt.show()
