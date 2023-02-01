import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def zetaPlot(Array):
    Array = Array.swapaxes(0, 1)[:, 250:350+1]
    plt.figure(figsize=(10, 6))
    plt.tight_layout()
    # plt.title(r"$\zeta$ $[\frac{1}{s}]$" f",  t = {t * dt} s", fontsize=14)
    plt.xlabel("x [km]")
    plt.ylabel("z [km]")
    plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
    plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
    plt.contourf(Array, levels=np.arange(-0.2, 0.2+0.02, 0.02), extend='both', cmap=cm.twilight_shifted)
    cbar = plt.colorbar(pad=0.05)
    cbar.set_ticks(np.arange(-0.2, 0.2+0.02, 0.02))
    plt.savefig(f"../graphs/zeta/test.png", dpi=300)
    plt.show()
    plt.close()

