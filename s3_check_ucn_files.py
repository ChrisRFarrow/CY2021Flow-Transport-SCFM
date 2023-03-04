import numpy as np
import matplotlib.pyplot as plt


'''
- xxx


'''


def plot_conc_vs_time(file):
    ifile = f'output/{file}.dat'  # 4 cy2021 ecf
    arr = np.loadtxt(ifile)
    #arr[0, 0] = -999
    #arr[270:663, 640:1164] = 0

    # Save
    np.savetxt(f'output/{file}_revised.dat',
               arr, delimiter=" ", fmt='%.6f')

    arr2 = np.ma.masked_less_equal(arr, 0)

    #plt.imshow(arr2, interpolation='none')
    fig, ax = plt.subplots(nrows=1, ncols=1)
    # vmin=1e-9, vmax=1e5
    im1 = ax.imshow(arr2, cmap='jet', aspect='auto')
    # add space for colour bar
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    fig.colorbar(im1, cax=cbar_ax)
    # plt.show()
    plt.savefig(f'output/{file}.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    # Plumes in CY2020 for CY2021 ECF Runs
    list_files = ['C14_20201', 'No3_20201', 'TCE_20201',
                  'HexCr_20201', 'SR90_20201', 'Trit_20201']

    list_files = ['No3_20201']
    for file in list_files:
        plot_conc_vs_time(file)

    '''
    # Plumes in CY2019 for CY2021 Runs
    list_files = ['C14_2019', 'No3_2019', 'TCE_2019',
                  'HexCr_2019', 'SR90_2019', 'Trit_2019']
    for file in list_files:
        plot_ic(file)
    '''
