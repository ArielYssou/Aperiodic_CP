import subprocess
from scipy.optimize import curve_fit
from numpy import linspace, logspace, log, polyfit, isnan
from random import randint
from os import listdir
from os.path import isfile, isdir, join
from time import sleep

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

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
def Curvature(fname, begin = 0, end= -1):
    try:
        dfile = open(fname, 'r')
        times = []
        rhos = []
        for line in dfile.read().splitlines():
            t, r = line.split(',')
            if float(t) == 0 or float(r) == 0:
                pass
            else:
                if isnan(float(t)) == False:
                    times.append(float(t))
                    rhos.append(float(r))
        fit_parans, fit_cov = curve_fit(func, log(times[begin:end]), log(rhos[begin:end]))
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
            t, r = line.split(',')
            #t, r, dt, dr = line.split(',')
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

def IsActive(fname, s1 = 100, e1 = 150, s2 = -50, e2 = -1, method = 'lin'):
    if not isfile(fname):
        return 'active'
    elif file_len(fname) <= 200:
        return 'inactive'
    else:
        pass
    
    crit1 = 'active'
    crit2 = 'active'
     
    # Method 1: Compare the slop in the middle to the final slope
    reference_slope = abs(Slope(fname, s1, e1))
    final_slope = abs(Slope(fname, s2, e2))
    tolerance = 1.0001
    if final_slope <= tolerance * reference_slope:
        crit1 = 'active'
    else:
        crit1 = 'inactive'

    # Method 2: 
    if Curvature(fname, s1, e2) > 0:
        crit2 = 'active'
    else:
        crit2 = 'inactive'

    if 'inactive' in (crit1, crit2):
        if 'active' in (crit1, crit2):
            tolerance = 1.0001
            reference_slope = abs(Slope(fname, s1, e1))
            final_slope = abs(Slope(fname, -20, -1))
            tolerance = 1.01
            if final_slope <= tolerance * reference_slope:
                return 'active'
            else:
                return 'inactive'
        else:
            return 'inactive'
 
    else:
        return 'active'

def RunSims(code, k, rho_0, delta_t, lamb_a, lamb_b, tsup, size, sim_i, sim_f, jobs):
    process = subprocess.Popen(
            [
            "bash",
            "aurora.sh",
            "-e",
            code,
            str(k),
            str(rho_0),
            str(delta_t),
            str(lamb_a), 
            str(lamb_b), 
            str(tsup),
            str(size),
            str(sim_i),
            str(sim_f),
            str(jobs)
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()

def RunAnalysis(analysis, k, delta_t, la, lb, tsup, size, sim_f):
    process = subprocess.Popen(
            [
            analysis,
            str(k),
            str(delta_t),
            str(la), 
            str(lb), 
            str(tsup),
            str(size),
            str(0),
            str(sim_f)
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()
    fname = out.decode("utf-8")  
    return fname

def BissecStep(code, analysis, k, rho_0, delta_t, lamb_a, lamb_b, tsup, size, sim_i, sim_f, jobs):
    RunSims(
        code,
        k,
        rho_0,
        delta_t,
        lamb_a,
        lamb_b,
        tsup,
        size,
        sim_i,
        sim_f,
        jobs
        )

    while True:
        if Finished(code):
            break
        else:
            MoveQueue()

    fname = RunAnalysis(
            analysis,
            k,
            delta_t,
            lamb_a,
            lamb_b,
            tsup,
            size,
            sim_f
            )

    return fname

def Finished(code):
    process = subprocess.Popen(
            [
            "bash",
            "aurora.sh",
            "-s",
            code.split('/')[-1]
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()
    if float(out) == 1:
        return False
    else:
        return True

def MoveQueue():
    process = subprocess.Popen(
            [
            "bash",
            "aurora.sh",
            "-m"
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()


def CheckGateway(remote_host, remote_basedir):
    dir = remote_basedir
    process = subprocess.Popen(
            [
            "ssh",
            remote_host,
            "bash",
            "{}/gateway.sh".format(dir)
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()
    if out.decode('utf-8') == 'open':
        print("Gateway is open")
        return True
    else:
        print("Gateway is closed. Aborting")
        return False

def UploadFile(file, remote_host, remote_basedir):
    # Creates additional dirs
    additional_dir = ''
    ignore_dirs = [
            'Simulacoes',
            'data_aperiodic'
            ]
    dirs = file.replace('/home/ariel/','').rsplit('/')
    for dir in dirs[:-1]: # The last element is the file name
        if dir not in ignore_dirs:
            additional_dir += '/' + dir

    remote_fulldir = remote_basedir + additional_dir
    if additional_dir != '':
        process = subprocess.Popen(
                [
                "ssh",
                remote_host,
                "mkdir",
                "-p",
                remote_fulldir,
                ],
                stdout = subprocess.PIPE)
        out, err = process.communicate()

    # Upload file
    process = subprocess.Popen(
            [
            "scp",
            file,
            "{}:{}".format(remote_host, remote_fulldir),
            ],
            stdout = subprocess.PIPE)
    out, err = process.communicate()
    return 0

def ExistingSims(k, lamb_a, lamb_b, tmax, size):
    target = '/home/ariel/Simulacoes/data_aperiodic'
    target = target + '/k={}'.format(int(k))
    target = target + '/Lambda_a={:.3f}'.format(lamb_a)
    target = target + '/Lambda_b={:.8f}'.format(lamb_b)
    target = target + '/Tmax={:e}'.format(tmax)
    target = target + '/Size={}'.format(int(size))

    if not isdir(target):
        return 0
    else:
        pass

    files = [f for f in listdir(target) if isfile(join(target, f))]

    for f in ['rho_av.dat', 'surv_prob.dat']:
        try:
            files.remove(f)
        except ValueError:
            pass
    for f in files:
        if 'rho' not in f:
            files.remove(f)

    return len(files)

def Bissection(k = 2, lsup = 200, linf = 100.0, steps = 15, resume = False):
    size = 50000
    tmax = 1e6
    delta_t = 1
    sims = 50000
    sim_increment = 0
    lamb_a = 1.3
    rho_0 = -1
    tinf = 0
    jobs = 150

#    if not resume:
#        if k == 1:
#            linf, lsup = (3, 8)
#        elif k == 2:
#            linf, lsup = (3, 8)
#        elif k == 3:
#            linf, lsup = (3, 8)
#        else:
#            linf, lsup = (3, 8)
    lmid = (lsup + linf) / 2

    times = linspace(0, tmax, steps)
    #flicks = linspace(0.95, 0.96, steps+1)
    #flicks = linspace(0.95, 0.99, steps+1)

    Sups = []
    Infs = []

    root = '/home/ariel/Simulacoes'
    remote_host = "ariel@ariel-Inspiron-5421"
    remote_basedir = "/home/ariel/Desktop/Mestrado/Aperiodic_CP/Aperiodic_Bissection/AbaxData"

    code = '{}/acp_abax'.format(root)
    analysis = '{}/acp_analysis'.format(root)

    root += '/data_aperiodic/k={}'.format(k)

    out_fname = '{}/acp_bissec_results.dat'.format(root)

    if isfile(out_fname):
        if resume:
            out_file = open(out_fname, 'a')
        else:
            print("Warning, there is a set of existing data, overwrite? (yes/no)")
            ans = input("> ")
            if ans.lower() in ('y','yes'):
                out_file = open(out_fname, 'w+')
            else:
                return
    else:
        out_file = open(out_fname, 'w+')

    transfer_files = [
            out_fname
            ]

    out_file.write("step, delta_t, lamb_a, lamb_b, tmax, size, sims, fname, regime\n")

    regimes = {}
    regimes[linf] = 'inactive'
    regimes[lsup] = 'active'

    for step in range(steps - 1):
        print("-" * 40)
        print(progress_bar(
            step,
            steps-1,
            "Depth {}".format(step),
            color = 2,
            size = 50)
            )

        for lamb_b in [lmid]:
            tsup = times[step + 1]

            print("Testing value \033[38;5;7m{:4f}\033[0m... ".format(lamb_b), end = '')

            sim_i = 0
            sim_f = sims + sim_i

            if k == 0:
                lamb_a = lamb_b
            else:
                pass

            fname = BissecStep(
                code,
                analysis,
                k,
                rho_0,
                delta_t,
                lamb_a,
                lamb_b,
                tsup,
                size,
                sim_i,
                sim_f,
                jobs
                )

            regimes[lamb_b] = IsActive(fname)
            if regimes[lamb_b] == 'active':
                print("\033[38;5;2;1m{}.\033[0m".format(regimes[lamb_b].title()))
            else:
                print("\033[38;5;1;1m{}.\033[0m".format(regimes[lamb_b].title()))
            transfer_files.append(fname)
            fname_alt = fname 
            transfer_files.append(fname_alt.replace('surv_prob', 'rho_av'))
            transfer_files.append(fname_alt)
            fname = fname.replace('/home/ariel/','')
            fname = fname.replace('Simulacoes/','')
            fname = fname.replace('data_aperiodic/','')
            fname = remote_basedir + '/' + fname

            out_file.write("{},{},{},{},{},{},{},{},{}\n".format(step,delta_t, lamb_a, lamb_b, tsup, size, sim_f, fname, regimes[lamb_b]))

        if regimes[lmid] == 'active':
            if regimes[linf] == 'active':
                lsup = linf
                for lamb, status in sorted(regimes.items(), key= lambda l: l[0])[::-1]:
                    if status == 'inactive':
                        linf = lamb
            else:
                lsup = lmid
        else:
            if regimes[lsup] == 'active':
                linf = lmid
            else:
                lsup = linf
                for lamb, status in sorted(regimes.items(), key= lambda l: l[0]):
                    if status == 'inactive':
                        linf = lamb
        lmid = (lsup + linf) / 2

    print("-" * 40)
    print(progress_bar(
        100,
        100,
        "Done =D",
        color = 3,
        size = 50)
        )

    regimes[lmid] = 'crit'
    
    lamb_b = lmid

    fname = BissecStep(
        code,
        analysis,
        k,
        rho_0,
        delta_t,
        lamb_a,
        lamb_b,
        tsup,
        size,
        sim_i,
        sim_f,
        jobs
        )

    transfer_files.append(fname)
    fname_alt = fname 
    transfer_files.append(fname_alt.replace('surv_prob', 'rho_av'))
    transfer_files.append(fname_alt)
    fname = fname.replace('/home/ariel/','')
    fname = fname.replace('Simulacoes/','')
    fname = fname.replace('data_aperiodic/','')
    fname = remote_basedir + '/' + fname

    out_file.write("{},{},{},{},{},{},{},{},{}\n".format(step, delta_t, lamb_a, lamb_b, tsup, size, sim_f, fname, regimes[lamb_b]))

    out_file.close()

    CheckGateway(remote_host, remote_basedir)
    for file in transfer_files:
        UploadFile(file, remote_host, remote_basedir)

if __name__ == '__main__':
    from sys import argv

    if '-b' in argv:
        Bissection()
    elif '-i' in argv:
        Increase_Resolution()
    elif '-c' in argv:
        Continue_Bissec()
    else:
        Bissection()
