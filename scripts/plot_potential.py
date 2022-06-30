import numpy as np
import matplotlib.pyplot as plt

def read_data(filename):
    arr = np.loadtxt(filename)
    r = arr[:, 0]
    U = arr[:, 1]
    F = arr[:, 2]
    return r, U, F

def make_plot(r, U, F, filename = None):
    plt.figure(figsize=(4, 3))
    
    plt.plot(r, U, 'b-', label = "U, Potential", linewidth = 1.5)
    plt.plot(r, F, 'r--', zorder = 0, label = "F, Force", linewidth = 2)

    plt.xlim([0,4])
    plt.ylim([-2.5,3])
    plt.xlabel("r")
    plt.legend()

    plt.grid("on")
    if filename is not None:
        plt.savefig(filename, bbox_inches = "tight")
    plt.show()
    plt.close()
    

def main():
    in_file = "potential.log"
    out_file = "fig/morse_potential.png"

    r, U, F = read_data(in_file)
    make_plot(r, U, F, out_file)



if __name__ == "__main__":
    main()
