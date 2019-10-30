import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from numpy import linspace, mean, log, isnan, sqrt
from scipy.optimize import curve_fit
from sys import argv
from os.path import isfile

fig, axes = plt.subplots(2, 1, figsize=(6,8))
fig2, axes2 = plt.subplots(1)

def func(x, a, b):
    return (a * x) + b

#Ncolors = 20
#colormap = plt.cm.tab20
#Ncolors = min(colormap.N,Ncolors)
#colors = ( colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors))

Ks = [ int(num) for num in argv[1:] ]

for k in Ks:
    if k == 1:
        colors = {
                'active' : "#f42c04",
                'inactive' : '#2f74a3',
                'crit' : '#333333'
                }
    elif k == 2:
        colors = ["#4a5d99", "#474338", "#ce8d66"] 
    elif k == 3:
        colors = {
                'active' : "#29c653",
                'inactive' : "#c47f4a",
                'crit' : "#333333"
                }
    elif k == 4:
        colors = {
                'active' : "#d64865",
                'inactive' : "#d8a15d",
                'crit' : "#333333"
                }
    else:
        colors = {
                'active' : '#f26419',
                'inactive' : '#2f74a3',
                'crit' : '#333333'
                }


    lambdas = set()
    tmaxs = {}
    files = {}
    regimes = {}

    with open(f'./k={k}/acp_bissec_results.dat', 'r') as f:
        next(f)
        for line in f.read().splitlines():
            step, dt, la, lb, tmax, size, sims, fname, regime = line.split(',')
            lambdas.add(lb)
            if lb in tmaxs.keys():
                if tmax > tmaxs[lb]:
                    tmaxs[lb] = tmax
                    files[lb] = fname
                    regimes[lb] = regime
                else:
                    pass
            else:
                tmaxs[lb] = tmax
                files[lb] = fname
                regimes[lb] = regime

    top_lbs = []
    btn_lbs = []

    for lb, regime in regimes.items():
        if regime == 'active':
            top_lbs.append(float(lb))
        elif regime == 'inactive':
            btn_lbs.append(float(lb))
        elif regime == 'crit':
            crit_lb = lb
        else:
            print("lambda without regime")
            exit(0)

    top_lbs.sort()
    btn_lbs.sort()

    d1 = dict(zip(top_lbs, linspace(0.1, 1, len(top_lbs))))
    d2 = dict(zip(btn_lbs[::-1], linspace(0.1, 1, len(btn_lbs))))
    alphas = {**d1, **d2}

    for lb in lambdas:
        #dfile = open(files[lb].replace('rho_av', 'surv_prob'), 'r')
        dfile = open(files[lb], 'r')
        ts = []
        ps = []
        for line in dfile.read().splitlines():
            #t, rho, dt, drho = line.split(',')
            t, prob = line.split(',')
            if True in (isnan(float(t)), isnan(float(prob))):
                pass
            else:
                ts.append(float(t))
                ps.append(float(prob))
        dfile.close()
        if regimes[lb] == 'crit':
            print(f"critical lambda: {lb}")
            fit_parans, fit_cov = curve_fit(
                        func,
                        log(ts[-50:-1]),
                        log(ps[-50:-1])
                        )
            print(f"delta: {fit_parans[0]}({sqrt(fit_cov[1][1])}) - delta limpo: 0.15946")
        if regimes[lb] != 'crit':
            axes[0].loglog(ts, ps, c = colors[regimes[lb]], alpha = alphas[float(lb)])
        else:
            lcrit = float(lb)
            axes[0].loglog(ts, ps, lw = 3, c = colors['crit'])
        #axes[0].loglog(ts, [t ** -.160 - 0.01 for t in ts])
        axes[0].set_xlabel(r"t", fontsize = 15)
        axes[0].set_ylabel(r"$P_s$", fontsize = 15)
        axes[0].set_xlim(6.3, 1.7e+6)
        axes[0].set_ylim(0.0215, 1)

        custom_lines = []
        for clr in colors.values():
            custom_lines.append(Line2D([0], [0], color=clr, lw=4))
        axes[0].legend(
                custom_lines, 
                    [
                        r"$\lambda_{sup} = $" + f"{float(top_lbs[-1]):.3f}",
                        r"$\lambda_{inf\ } = $" + f"{float(btn_lbs[0]):.3f}",
                        r"$\lambda_{crit} = $" + f"{float(crit_lb):.3f}"
                    ],
                fontsize = 12)


        dfile = open(files[lb].replace('surv_prob', 'rho_av'), 'r')
        ts = []
        rs = []
        for line in dfile.read().splitlines():
            t, rho, dt, drho = line.split(',')
            if True in (isnan(float(t)), isnan(float(rho))):
                pass
            else:
                if float(rho) > 0:
                    ts.append(float(t))
                    rs.append(float(rho))
        dfile.close()
        if regimes[lb] == 'crit':
            print(f"critical lambda: {lb}")
            fit_parans, fit_cov = curve_fit(
                        func,
                        log(ts[-50:-1]),
                        log(rs[-50:-1])
                        )
            print(f"theta: {fit_parans[0]}({sqrt(fit_cov[1][1])}) - theta limpo: 0.31369")
        if regimes[lb] != 'crit':
            axes[1].loglog(ts, rs, c = colors[regimes[lb]], alpha = alphas[float(lb)])
        else:
            lcrit = float(lb)
            axes[1].loglog(ts, rs, lw = 3, c = colors['crit'])
        axes[1].set_xlabel(r"t", fontsize = 15)
        axes[1].set_ylabel(r"$\rho$", fontsize = 15)
        axes[1].set_xlim(6.3, 1.7e+6)
        axes[1].set_ylim(5.0e-5, 0.04)

        if regimes[lb] != 'crit':
            #axes2.loglog(ts, [ r/p for r, p in list(zip(rs, ps))], color = colors[regimes[lb]], alpha = alphas[float(lb)])
            #axes2.loglog(ps, rs , color = colors[regimes[lb]], alpha = alphas[float(lb)])
            pass
        else:
            pass
            #axes2.loglog(ps, rs , color = colors[regimes[lb]], alpha = alphas[float(lb)])
        #axes2.loglog(linspace(.1,1,100), [1e-4 * p ** (-1.96712) for p in linspace(.1,1,100)])
        #axes2.loglog(ts, [ r/p for r, p in list(zip(rs, ps)) ])

#
#
#    rho_est = {} # lambda : rho_est
#    #rho_desv
#    corr_times = []
#    corr_rhos = []
#    epsilons = { lb : abs(lb - lcrit) for lb in top_lbs}
#    for lb in top_lbs:
#        dfile = open(files[lb].replace('surv_prob', 'rho_av'), 'r')
#        #dfile = open(files[str(lb)], 'r')
#        ts = []
#        rs = []
#        for line in dfile.read().splitlines():
#            t, rho, dt, drho = line.split(',')
#            ts.append(float(t))
#            rs.append(float(rho))
#        dfile.close()
#
#        rho_est[lb] = mean(rs[-100:-1])
#
#        found = False
#        hits = 0
#        for t, rho in list(zip(ts,rs)):
#            tol = 1.1
#            #if rho < tol * epsilons[lb] * rho_est[lb]:
#            if rho < tol * rho_est[lb]:
#                if hits > 1:
#                    corr_times.append(float(t))
#                    corr_rhos.append(float(rho))
#                    found = True
#                    break
#                else:
#                    hits += 1
#            else:
#                hits = 0 
#        if not found:
#            epsilons.pop(lb)
#
#        clr = colors['active']
##    axes.loglog(corr_times[-1], corr_rhos[-1], ls='',
##            marker='D', color = clr, markerfacecolor= clr,
##            markeredgewidth=1.5, markeredgecolor='black' )
##    axes.loglog([1e6], [rho_est[lb]], ls='',
##            marker='D', color = clr, markerfacecolor= clr,
##            markeredgewidth=1.5, markeredgecolor='black' )
##
##axes2.loglog(list(epsilons.values()), corr_times, ls = '--', marker = 'o')
##axes2.loglog(list(epsilons.values()), list(rho_est.values()), ls = '--', marker = 'o')
fig.tight_layout()
if isfile(f'bissec_result_k{Ks[0]}.pdf'):
    print("Overwrite existing plot?")
    ans = input("> ")
    if ans.lower() in ("y", "yes"):
        fig.savefig(f'bissec_result_k{Ks[0]}.pdf')
    else:
        pass
else:
    fig.savefig(f'bissec_result_k{Ks[0]}.pdf')

plt.show()

