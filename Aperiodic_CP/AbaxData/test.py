import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numpy import linspace, log, isnan
from scipy.optimize import curve_fit

#def func(x, a, b, c):
#    return a + (b * x) + c * (x ** 2) 
#def IsActive(fname):
#    try:
#        dfile = open(fname, 'r')
#        times = []
#        rhos = []
#        for line in dfile.read().splitlines():
#            t, r, dt, dr = line.split(',')
#            if float(t) == 0 or float(r) == 0:
#                pass
#            else:
#                if isnan(float(t)) == False:
#                    times.append(float(t))
#                    rhos.append(float(r))
#        fit_parans, fit_cov = curve_fit(func, log(times[500:]), log(rhos[500:]))
#        #fit_parans = polyfit(log(times), log(rhos), 2)
#        if fit_parans[2] >= 0:
#            return 'active'
#        else:
#            return 'inactive'
#    except:
#        print("file {} not found".format(fname))
#        raise FileNotFoundError

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

def IsActive(fname, s1 = 100, e1 = 200, s2 = -100, e2 = -1):
    if abs(Slope(fname, s1, e1)) > abs(Slope(fname, s2, e2)):
        return 'active'
    else:
        return 'inactive'

fig, axes = plt.subplots(1)

ts = []
rs = []
with open('rho_av.dat', 'r') as f:
    for line in f.read().splitlines():
        t, r, dt, dr = line.split(',')
        rs.append(float(r))
        ts.append(float(t))
axes.loglog(ts, rs)

print(IsActive('rho_av.dat'))

plt.show()
