import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange 

def expcos(x):
    return np.exp(x)*np.cos(x)

def num_int(N:int,x_inf:float = 0.,x_sup:float = np.pi, fprint:bool = False, func:callable = expcos):
    samp_points = np.linspace(x_inf,x_sup,N)
    func_vals = np.array([func(x) for x in samp_points])
    
    res = np.trapezoid(func_vals,samp_points)
    if fprint:
        print(f"The resulting value of the integral for {N} samples is {res:.16f}")
    return res

if __name__ == "__main__":
    N = np.logspace(4,26, 46, base=2)
    x_inf = 0
    x_sup = np.pi/2
    I_true = (np.exp(np.pi/2)-1)/2

    epsrel = np.zeros(len(N))
    for its in trange(len(N)):
        I_res = num_int(round(N[its]),x_inf,x_sup)
        epsrel[its] = np.abs(I_res/I_true-1)

    # plotting
    fig, ax = plt.subplots()
    plt.plot(N, epsrel, marker='o')
    plt.title(r"$\int_0^{\pi/2}e^{x}\cos(x)dx$")
    plt.xlabel("N")
    plt.ylabel(r"$\epsilon_{rel}$")
    ax.set_xscale('log', base=2)
    plt.yscale("log")

    plt.savefig("epsvN_python.pdf", format="pdf", bbox_inches="tight")
    plt.savefig("epsvN_python.png", format="png", bbox_inches="tight")
