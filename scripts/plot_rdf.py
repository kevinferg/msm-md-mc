import numpy as np
import matplotlib.pyplot as plt

def read_data(filename):
    arr = np.loadtxt(filename)
    r = arr[:, 0]
    g = arr[:, 1]
    return r, g

def make_plot(r, g, filename = None):
    plt.figure(figsize=(4, 3),dpi=200)
    
    plt.plot(r, g, 'b-', label = "g(r)", linewidth = 1.5)

    plt.xlim(0,np.max(r))
    plt.ylim(0,np.max(g)*1.05)
    plt.xlabel("r")
    plt.ylabel("g(r)")

    plt.grid("on")
    if filename is not None:
        plt.savefig(filename, bbox_inches = "tight")
    plt.show()
    plt.close()
    

def main():
    in_file = "rdf.txt"
    out_file = "fig/rdf.png"

    r, g = read_data(in_file)
    make_plot(r, g, out_file)


if __name__ == "__main__":
    main()
