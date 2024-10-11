import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import cm
from concurrent.futures import ThreadPoolExecutor

# Constants
V_init = np.array([0.0, 0.0, 0.0])
freq = 2450000000.0
B = np.array([0.0, 0.0, 0.0875])
E = np.array([1.0, 0.0, 0.0])
E2 = np.array([0.0, -1.0, 0.0])
e = 1.6021766208e-19
me = 9.10938356e-31
runlength = 16000

# Range: 1000000000 - 4000000000, Step size: 50000000
res16000_var_freq = [5423777887.860481, 6054527909.020197, 6474045787.760949, 6197663167.538217, 5734210911.973459, 6111349919.325712, 6970928784.201629, 7477551546.594364, 7120915614.107124, 6610325638.7028055, 7231773432.785794, 8460691273.110773, 9146504930.323093, 8710246431.986973, 8160856647.57918, 9231723028.700577, 11169783830.639965, 12294212230.409092, 11851639490.995533, 11374237321.598768, 13555832039.606537, 17363237843.39732, 20061438617.07363, 20386929567.54864, 21189479434.14722, 28843442384.732937, 44273325032.09847, 66352206279.49256, 109257636323.60173, 262257170921.40237, 101022275222.27974, 62397214681.50242, 42164670339.63267, 26902153609.301514, 18422354092.338448, 17300696231.888443, 17565613892.950676, 15326626475.14969, 11376353091.739325, 8680504457.726845, 9027994651.785418, 9919011050.057579, 9083256186.981844, 6941197398.2908125, 5449056345.575133, 5939356810.211611, 6782056645.775941, 6324324015.092768, 4860615641.900285, 3853120310.2923026, 4336824337.690223, 5083492094.287296, 4778899837.479259, 3662844852.9859967, 2910588890.9077992, 3362421825.686735, 4023755043.9221845, 3796102170.4041314, 2890099753.074362, 2293575152.879548]

# Function to compute final velocity magnitude for a given frequency (frt)
def compute_velocity_f(frt):
    V = np.copy(V_init)
    for i in range(runlength):
        F = e * ((E * math.cos(2.0 * math.pi * frt * (i / 1000000000000.0))) + np.cross(V, B))
        A = F / me
        V += A * (1.0 / 1000000000000.0)  # Update velocity
    final_velocity = np.sqrt(V.dot(V)) / (runlength / 1000000000000.0)  # Magnitude of final velocity
    return frt, final_velocity

def compute_velocity_B(b):
    V = np.copy(V_init)
    for i in range(runlength):
        F = e * ((E * math.cos(2.0 * math.pi * freq * (i / 10000000000000.0))) + np.cross(V, b))
        A = F / me
        V += A * (1.0 / 10000000000000.0)  # Update velocity
    final_velocity = np.sqrt(V.dot(V)) / (runlength / 10000000000000.0)  # Magnitude of final velocity
    return b, final_velocity

def fapp(x, l):
    return 2.62e11 / (1.0 + np.square( (x - 2470000000)/(l * freq) ) )

def compute_fapp_err(l):
    errs = 0.0
    for x in range(0, len(res16000_var_freq)):
        errs += np.square( fapp((1000000000.0 + 50000000.0 * x), l) - res16000_var_freq[x])
    return errs / len(res16000_var_freq);

def minimize_fapp_err():
    rl = 0.0000001
    rr = 1.0
    for i in range(0, 100):
        rm = (rl + rr) / 2
        el = compute_fapp_err(rl)
        er = compute_fapp_err(rr)
        #em = compute_fapp_err(rm)
        if el > er:
            rl = rm
        else:
            rr = rm
    l = rl
    x = np.linspace(1000000000.0, 4000000000.0, 1000)
    xpre = np.linspace(1000000000.0, 4000000000.0, len(res16000_var_freq))

    fig, ax1 = plt.subplots()
    ax1.set_ylabel('avg a in m/s^2')
    ax1.set_xlabel('f in Hz')
    ax1.plot(x, np.log10(fapp(x, l)), label="approx., l = {}".format(l))
    ax1.plot(xpre, np.log10(res16000_var_freq), label="sim")
    plt.legend(loc="upper left")
    plt.show()


def par_sim():
    lAf = []
    lAb = []
    lF = []
    lB = []

    with ThreadPoolExecutor() as executor:
        # Submit tasks to thread pool
        #futures = [executor.submit(compute_velocity_f, 50000000.0 * f) for f in range(0, 200)]
        futures = [executor.submit(compute_velocity_f, 50000000.0 * f) for f in range(20, 80)] # 1GHz-4gHz
        #futures = [executor.submit(compute_velocity_f, 2250000000.0 + 10000000.0 * f) for f in range(0, 50)] # 2.25GHz-2.75GHz
        # Collect results as they complete
        for future in futures:
            frt, final_velocity = future.result()
            lF.append(frt)
            lAf.append(final_velocity)
            print("{0} : {1}".format(frt, final_velocity))

        # futuresb = [executor.submit(compute_velocity_B, np.array([0.0, 0.0, 0.0875 + (f * 0.0025)])) for f in range(0, 105)]
    
        # for future in futuresb:
        #     b, final_velocity = future.result()
        #     lB.append(b)
        #     lAb.append(final_velocity)
        #     print("{0} : {1}".format(b[2], final_velocity))
    fig, ax1 = plt.subplots()
    ax1.ylabel('avg a in m/s^2 after 1 s')
    ax1.xlabel('f in Hz')
    ax1.plot(lF, lAf)
    ax2 = ax1.twiny()
    ax2.xlabel('B in T')
    ax2.plot(lB, lAb)
    plt.show()

def vis_electron_path(b, f, steps):
    B = np.array([0.0, 0.0, b])
    Pos = np.array([0.0, 0.0, 0.0])
    V = np.array([0.0, 0.0, 0.0])
    Px = []
    Py = []
    Vh = []
    for i in range(steps):
        p = 2.0 * math.pi * f * (i / 1000000000000.0)
        F = e * (( E * math.cos(p) + E2 * math.sin(p) ) + np.cross(V, B))
        A = F / me
        V += A * (1.0 / 1000000000000.0)  # Update velocity
        Pos += V * (1.0 / 1000000000000.0)  # Update position
        Px.append(Pos[0])
        Py.append(Pos[1])
        Vh.append(np.sqrt(V.dot(V)))
    maxV = max(Vh)
    x = np.linspace(0, maxV, len(Vh))
    # plt.scatter(Px, Py, c=cm.inferno(Vh/maxV))
    plt.plot(Px, Py)
    plt.show()

# minimize_fapp_err()

vis_electron_path(0.0875, 2450000000, 20000)
