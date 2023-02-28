import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import xarray as xr


######## Should Be Tuned ########### 
left_range, right_range, split_range = 250, 350+1, 10

left_th, right_th, split_th = -16, 16, 11

left_zeta, right_zeta, split_zeta = -0.03, 0.03, 11

wind_skip = 2
wind_scale = 300
wind_speed = 10

left_qc, right_qc, split_qc = 0, 3, 4

DT = 0.1
xSize, zSize = 150000, 15000
DX = DZ = 250

LEAP = 250
######## Should Be Tuned ###########

NX, NZ = int(xSize / DX), int(zSize / DZ)

cmap = cm.twilight_shifted
DPI = 300

############################# qc+qr+th+u+w ##########################
def qc_qr_th_u_w(t):
    data = xr.open_dataset(f"../outputs/nc/{t}.nc")
    X, Z = np.meshgrid(np.arange(right_range-left_range), np.arange(NZ))

    # load data
    th = data['th'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]
    u = data['u'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]
    w = data['w'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]
    qc = data['qc'].to_numpy().swapaxes(0, 1)[:, left_range:right_range] * 1000
    qr = data['qr'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]

    if t == 0:
        rho = data['rho'].to_numpy()[1:-1]
        tb = data['tb'].to_numpy()[1:-1]
        qvb = data['qvb'].to_numpy()[1:-1]
        qvsb = data['qvsb'].to_numpy()[1:-1]

        fig, ax = plt.subplots(1, 4, figsize=(10, 6), gridspec_kw={'width_ratios': [1, 1, 1, 6]}, dpi=DPI)
        # fig.tight_layout()

        z = np.arange(125, zSize, DZ)[1:-1] / 1000

        ax[1].get_yaxis().set_visible(False)
        ax[2].get_yaxis().set_visible(False)
        ax[3].get_yaxis().set_visible(False)

        ax[0].plot(rho, z)
        ax[0].set_ylim(0, 15)
        ax[0].set_title("Init. " + r"$\overline{\rho}$")
        ax[0].set_ylabel("z [km]")
        ax[0].set_xlabel(r"$\overline{\rho}$  [$\frac{kg}{m^3}$]")

        ax[1].plot(tb, z)
        ax[1].set_ylim(0, 15)
        ax[1].set_title("Init. " + r"$\overline{\theta}$")
        ax[1].set_xlabel(r"$\overline{\theta}$  [K]")

        ax[2].plot(qvb, z)
        ax[2].plot(qvsb, z)
        ax[2].set_ylim(0, 15)
        ax[2].legend([r"$\overline{q_v}$", r"$\overline{q_{vs}}$"])
        ax[2].set_title("Init. " + r"$\overline{q_v}$&$\overline{q_{vs}}$")
        ax[2].set_xlabel(r"mixing ratio  [$\frac{kg}{kg}$]")

        # deal with qr
        qrx, qrz = [], []
        for i in range(qr.shape[0]):
            for k in range(qr.shape[1]):
                if qr[i][k] >= 0.001:
                    qrx.append(i)
                    qrz.append(k)

        # plot th
        ax[3].set_title(f"t = {t * DT} s", fontsize=14)
        ax[3].set_xlabel("x [km]")
        ax[3].set_xticks(np.linspace(0, right_range-left_range-1, split_range), np.linspace(DX*left_range/1000, DX*right_range/1000, split_range))
        CF = ax[3].contourf(th, levels=np.linspace(left_th, right_th, split_th), extend='both', cmap=cmap)
        cbar = plt.colorbar(CF, pad=0.05)
        cbar.set_ticks(np.linspace(left_th, right_th, split_th))
        cbar.set_label(r"$\theta^{'}$ [K]")

        # plot u,w 
        Q = ax[3].quiver(X[::wind_skip, ::wind_skip], Z[::wind_skip, ::wind_skip], u[::wind_skip, ::wind_skip], w[::wind_skip, ::wind_skip], angles='xy', units="width", scale=wind_scale)
        qk = plt.quiverkey(Q, 0.8, 0.9, wind_speed, f"{wind_speed}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')

        # plot qc
        CS = ax[3].contour(qc, levels=np.linspace(left_qc, right_qc, split_qc), colors='yellow', linewidths=1.5)
        plt.clabel(CS, inline=1, fontsize=14, fmt="%1.0f")

        # plot qr
        scatter = ax[3].scatter(qrz, qrx, c='deepskyblue', s=5, marker='|')

        plt.show()

    else:
        # deal with qr
        qrx, qrz = [], []
        for i in range(qr.shape[0]):
            for k in range(qr.shape[1]):
                if qr[i][k] >= 0.001:
                    qrx.append(i)
                    qrz.append(k)

        # plot th
        plt.figure(figsize=(10, 6), dpi=DPI)
        plt.title(f"t = {t * DT} s", fontsize=14)
        plt.xlabel("x [km]")
        plt.ylabel("z [km]")
        plt.xticks(np.linspace(0, right_range-left_range-1, split_range), np.linspace(DX*left_range/1000, DX*right_range/1000, split_range))

        plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
        CF = plt.contourf(th, levels=np.linspace(left_th, right_th, split_th), extend='both', cmap=cmap)
        cbar = plt.colorbar(pad=0.05)
        cbar.set_ticks(np.linspace(left_th, right_th, split_th))
        cbar.set_label(r"$\theta^{'}$ [K]")

        # plot u,w 
        Q = plt.quiver(X[::wind_skip, ::wind_skip], Z[::wind_skip, ::wind_skip], u[::wind_skip, ::wind_skip], w[::wind_skip, ::wind_skip], angles='xy', units="width", scale=wind_scale)
        qk = plt.quiverkey(Q, 0.7, 0.9, wind_speed, f"{wind_speed}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')

        # plot qc
        CS = plt.contour(qc, levels=np.linspace(left_qc, right_qc, split_qc), colors='yellow', linewidths=1.5)
        plt.clabel(CS, inline=1, fontsize=14, fmt="%1.0f")

        # plot qr
        scatter = plt.scatter(qrz, qrx, c='deepskyblue', s=5, marker='|')

        # plot legend
        h1, _ = CS.legend_elements()
        h2, _ = CF.legend_elements()
        plt.legend([h1[0], scatter], [r"$q_c\quad[\frac{g}{kg}]$", r"$q_r > 1 [\frac{g}{kg}]$"], loc='upper right')
    plt.savefig(f"../graphs/qc+qr+th+u+w/{int(t/LEAP)}.png", dpi=DPI)
    plt.close()
    
############################# zeta ##########################
def zeta(t):
    data = xr.open_dataset(f"../outputs/nc/{t}.nc")
    X, Z = np.meshgrid(np.arange(right_range-left_range), np.arange(NZ))

    # load data
    zeta = data['zeta'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]
    u = data['u'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]
    w = data['w'].to_numpy().swapaxes(0, 1)[:, left_range:right_range]

    # plot zeta
    plt.figure(figsize=(10, 6), dpi=DPI)
    plt.title(r"$\zeta$ $[\frac{1}{s}]$" + f",  t = {t * DT} s", fontsize=14)

    plt.xlabel("x [km]")
    plt.ylabel("z [km]")
    plt.xticks(np.linspace(0, right_range-left_range-1, 11), np.linspace(DX*left_range/1000, DX*right_range/1000, 11))
    plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
    
    plt.contourf(zeta, levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    cbar = plt.colorbar(pad=0.05)
    cbar.set_ticks(np.linspace(left_zeta, right_zeta, split_zeta))
    
    # plot u,w
    Q = plt.quiver(X[::wind_skip, ::wind_skip], Z[::wind_skip, ::wind_skip], u[::wind_skip, ::wind_skip], w[::wind_skip, ::wind_skip], angles='xy', units="width", scale=wind_scale)
    qk = plt.quiverkey(Q, 0.7, 0.9, wind_speed, f"{wind_speed}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')

    plt.savefig(f"../graphs/zeta/{int(t/LEAP)}.png", dpi=DPI)