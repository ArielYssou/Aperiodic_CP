import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from numpy import logspace, exp, log
from scipy.optimize import curve_fit
from random import uniform

def func(x, a, b, c):
    return a + (b * x) + c * (x ** 2) 
def Curvature(times, rhos, begin = 0, end= -1):
    fit_parans, fit_cov = curve_fit(
            func,
            log(times[begin:end]),
            log(rhos[begin:end]))
    return fit_parans[2]

def lin_func(x, a, b):
    return (a * x) + b
def Slope(times, rhos, start = 0, end = -1):
    fit_parans, fit_cov = curve_fit(
            lin_func,
            log(times[start:end]),
            log(rhos[start:end])
            )
    return fit_parans[0]

def IsActive(times, rhos, s1 = 200, e1 = 300, s2 = -100, e2 = -1, method = 'curvature'):
    if method == 'lin':
        reference_slope = abs(Slope(times, rhos, s1, e1))
        final_slope = abs(Slope(times, rhos, s2, e2))
        tolerance = 1.01
        if final_slope <= tolerance * reference_slope:
            return 'active'
        else:
            return 'inactive'
    elif method == 'curvature':
        if Curvature(times, rhos, s1, e2) > 0:
            return 'active'
        else:
            return 'inactive'

def rho_t(lb, t, mu = 1,rho_init = 0.15):
    return ((exp((mu-lb)*t)/rho_init) + (lb/(mu-lb))*(exp((mu-lb)*t)-1)) ** (-1)

def Noises(points):
    n = 0
    amp = 0.001
    while n < points:
        yield exp(uniform(1 - amp, 1 + amp))
        amp *= 1.007
        n += 1

fig, axes = plt.subplots(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 20)
colors = []

mu= 1.0
lc = 1.00000000001
epsilon = 0.00003
curves = 10
points = 1000
tmax  = 5
times = logspace(0 , tmax, points)

# Critical Curve
clr ="#474338"  
lb = lc
noises = Noises(points)
rhos = [ rho_t(lb, t) * noises.__next__() for t in times ]
print(f"Critical curve: {IsActive(times,rhos)}")
axes.loglog(times, rhos, color = clr, ls = '-')

# Active Curve
clr ="#f28123"  
lb = (1.0 + epsilon) * lc
noises = Noises(points)
rhos = [ rho_t(lb, t) * noises.__next__() for t in times ]
print(f"Active curve: {IsActive(times,rhos)}")
axes.loglog(times, rhos, color = clr, ls = '-')

# Inactive Curve
clr = "#f42c04" 
lb = (1.0 - epsilon) * lc
noises = Noises(points)
rhos = [ rho_t(lb, t) * noises.__next__() for t in times ]
print(f"Inactive curve: {IsActive(times,rhos)}")
axes.loglog(times, rhos, color = clr, ls = '-')

custom_lines = []
for clr in colors:
    custom_lines.append(Line2D([0], [0], color=clr, lw=4))
axes.legend(custom_lines, [r"$\lambda < \lambda_c$", r"$\lambda = \lambda_c$", r"$\lambda > \lambda_c$"], fontsize = 12)

axes.set_xlim(1, 10 ** tmax)
#axes.set_ylim(0.0005, 0.15)
axes.set_xlabel(r't', fontsize=20)
axes.set_ylabel(r"$\rho$", fontsize=20)
fig.tight_layout()
plt.savefig("clean.pdf")
plt.show()

exit(0)



