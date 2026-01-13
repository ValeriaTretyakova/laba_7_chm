import matplotlib.pyplot as plt
import numpy as np

BASIS = "shannon"
N = 512


def get_data(fn):
    try:
        with open(fn, 'r') as f:
            return np.fromstring(f.read(), sep=' ')
    except:
        return None

for s in range(1, 5):
    low = get_data(f"{BASIS}_low_{s}.txt")
    high = get_data(f"{BASIS}_high_{s}.txt")
    p = get_data(f"{BASIS}_p_{s}.txt")
    orig = get_data("original.txt")

    if low is not None:
        fig, axs = plt.subplots(2, 2, figsize=(10, 6))
        fig.suptitle(f"Этап {s}")
        x = np.linspace(0, N, len(low))

        axs[0, 0].plot(orig, 'm.-', lw=0.5, ms=2)
        axs[0, 0].set_title("Исходный сигнал")

        axs[0, 1].plot(x, high, 'c.-', lw=0.5, ms=4)
        axs[0, 1].set_title(f"ВЧ коэффициенты <z, psi_-{s}, k>")

        axs[1, 0].plot(p, 'y.-', lw=0.5, ms=2)
        axs[1, 0].set_title(f"Восстановленный сигнал P_-{s}(z)")

        axs[1, 1].plot(x, low, 'b.-', lw=0.5, ms=4)
        axs[1, 1].set_title(f"НЧ коэффициенты <z, psi_-{s}, k>")

        for ax in axs.flat:
            ax.grid(True, linestyle='--', alpha=0.6)

        plt.tight_layout()

p2_filtered = get_data(f"{BASIS}_p_2.txt")
orig = get_data("original.txt")

if p2_filtered is not None:
    plt.figure(figsize=(12, 6))
    plt.title(f"Восстановленный сигнал P_-1(z) после обнуления <z, psi_-2, k>", fontsize=12)

    plt.plot(orig, '.-', color='black', alpha=0.3, ms=1, label="Исходный сигнал")

    plt.plot(p2_filtered, 'r.-', lw=2, ms=4, label=f"Отфильтрованный сигнал")

    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()

plt.show()