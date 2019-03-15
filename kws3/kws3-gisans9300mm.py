import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import bornagain as ba
from os.path import join
from ruamel.yaml import YAML


sensitivity_matrix = "/datadisk/data/KWS3/p15749/sens.sens"
datapath = '/home/mary/build/bornagain/users/uliana/kws3/gisans9m'

# critical angles
a_Si = 0.596   # degrees
a_Au = 0.878   # degrees
a_Cu = 1.059   # degrees


# data files
f65952 = '00065952_0000_p15749_first-sample-9m-friday_HRD_standard.det'        # 0 deg, large slit
f65953 = '00065953_0000_p15749_first-sample-9m-friday_HRD_standard.det'        # 0 deg, large slit, to sum with f65952
f66023 = '00066023_0000_p15749_first-sample-9m-sarturday_HRD_standard.det'     # 0 deg, small slit
f66024 = '00066024_0000_p15749_first-sample-9m-sarturday_HRD_standard.det'     # 0 deg, small slit, to sum with f66023
f66028 = '00066028_0000_p15749_sample-rot45dgr-9m-sunday_HRD_standard.det'     # 45 deg, small slit
f66029 = '00066029_0000_p15749_sample-rot90dgr-9m-sunday_HRD_standard.det'     # 90 deg, small slit
f66055 = '00066055_0000_p15749_final-sample_HRD_standard.det'                  # 0 deg, 1hour
f66057 = '00066057_0000_p15749_final-sample-90dgr_HRD_standard.det'            # 90 deg?, 1 hour
f66058 = '00066058_0000_p15749_final-sample-90dgr_HRD_standard.det'            # 90 deg?, to sum with f66057?


def calc_q(ai, af):
    return 2.0*np.pi*(np.sin(ai*ba.deg) + np.sin(af*ba.deg))/1.28


def load_data(fname):
    rawdata = np.loadtxt(fname)
    sens = np.loadtxt(sensitivity_matrix)
    data = rawdata*sens

    return np.rot90(data, 3)


def plot_2d_raw(fname, zmax=100.0):
    data = load_data(join(datapath, fname))
    plt.imshow(data, norm=matplotlib.colors.LogNorm(0.1, zmax), aspect='auto', cmap='jet')
    plt.title(fname)
    plt.show()


def plot_all_raw():
    for fname in [f65952, f65953, f66023, f66024, f66028, f66029, f66055, f66057, f66058]:
        plot_2d_raw(fname, 100)


def plot_large_mirror_raw():
    fname = '0006603{}_0000_p15749_large-mirror_HRD_standard.det'
    for i in range(3):
        plot_2d_raw(fname.format(i), 100)


def get_kws3_detector(u0=45.3, v0=29.9):
    """
    Creates and returns KWS-3 detector
    detector dimensions: 90x90 mm
    SDD: 9200 mm (GISANS)
    """
    detector = ba.RectangularDetector(256, 90.0, 256, 90.0)
    detector.setPerpendicularToDirectBeam(9200.0, u0, v0)
    return detector


def get_simulation(ai=0.56, u0=45.3, v0=29.9):
    """
    Returns a GISAXS simulation with beam and detector defined
    wavelength 12.8 angstrom
    incident angle 0.0 degrees (transmission)
    """
    simulation = ba.GISASSimulation()
    simulation.setBeamParameters(12.8*ba.angstrom, ai*ba.deg, 0.0*ba.deg)
    simulation.setDetector(get_kws3_detector(u0, v0))
    simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(5.0, 5.0))
    simulation.setBeamIntensity(1000)
    distr_1 = ba.DistributionGaussian(1.28*ba.nm, 0.1)
    simulation.addParameterDistribution("*/Beam/Wavelength", distr_1, 50, 2.0, ba.RealLimits.positive())
    simulation.getOptions().setIncludeSpecular(True)
    return simulation


def plot_all_2d():
    data65952 = load_data(join(datapath, f65952))
    data65953 = load_data(join(datapath, f65953))
    data66023 = load_data(join(datapath, f66023))
    data66024 = load_data(join(datapath, f66024))
    data66028 = load_data(join(datapath, f66028))
    data66029 = load_data(join(datapath, f66029))
    data66055 = load_data(join(datapath, f66055))
    data66057 = load_data(join(datapath, f66057))
    data66058 = load_data(join(datapath, f66058))

    datalist = [(65952, data65952), (65953, data65953), (66023, data66023), (66024, data66024), (66028, data66028),
                (66029, data66029), (66055, data66055), (66057, data66057), (66058, data66058)]

    simulation = get_simulation()
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)

    plt.style.use('seaborn-talk')
    plt.suptitle("Sample S3098 SDD=9200 mm", fontsize=16)
    i = 1
    for run, data in datalist:
        ax = plt.subplot(3, 3, i)
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        im = plt.imshow(data,
                        norm=matplotlib.colors.LogNorm(0.1, 50),
                        extent=axes_limits,
                        aspect='auto', cmap='jet')

        cb = plt.colorbar(im)
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.ylabel(r'$Q_z$ (a. u.)')
        plt.title(r"Run {}".format(run))
        i += 1

    plt.show()


def plot_all_1d():
    data65952 = load_data(join(datapath, f65952))
    data65953 = load_data(join(datapath, f65953))
    data66023 = load_data(join(datapath, f66023))
    data66024 = load_data(join(datapath, f66024))
    data66028 = load_data(join(datapath, f66028))
    data66029 = load_data(join(datapath, f66029))
    data66055 = load_data(join(datapath, f66055))
    data66057 = load_data(join(datapath, f66057))
    data66058 = load_data(join(datapath, f66058))

    datalist = [(65952, data65952), (65953, data65953), (66023, data66023), (66024, data66024), (66028, data66028),
                (66029, data66029), (66055, data66055), (66057, data66057), (66058, data66058)]

    simulation = get_simulation()
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)

    plt.style.use('seaborn-talk')
    plt.suptitle("Sample S3098 SDD=9200 mm", fontsize=16)
    i = 1
    for run, data in datalist:
        plt.subplot(3, 3, i)
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        qy = np.linspace(axes_limits[0], axes_limits[1], data.shape[0])
        counts = np.sum(data, axis=0)
        plt.semilogy(qy, counts, color='k', marker='o', markersize=5, linestyle='None')
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.ylabel(r'$I$ (a. u.)')
        plt.title(r"Run {}".format(run))
        i += 1

    plt.show()


def plot_summary():
    data00deg = load_data(join(datapath, f66023))   # 0 degrees
    data45deg = load_data(join(datapath, f66028))   # 45 degrees
    data90deg = load_data(join(datapath, f66029))   # 90 degrees

    simulation = get_simulation()
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)
    qy = np.linspace(axes_limits[0], axes_limits[1], data00deg.shape[0])

    sum00deg = np.sum(data00deg, axis=0)
    sum45deg = np.sum(data45deg, axis=0)
    sum90deg = np.sum(data90deg, axis=0)

    mon00deg = get_monitor_counts("{}.yaml".format(os.path.splitext(f66023)[0]))
    mon45deg = get_monitor_counts("{}.yaml".format(os.path.splitext(f66028)[0]))
    mon90deg = get_monitor_counts("{}.yaml".format(os.path.splitext(f66029)[0]))

    plt.semilogy(qy, sum00deg/mon00deg, color='k', marker='o', markersize=5, linestyle='None', label=r"$\omega=0^\circ$")
    plt.semilogy(qy, sum45deg/mon45deg, color='b', marker='s', markersize=5, linestyle='None', label=r"$\omega=45^\circ$")
    plt.semilogy(qy, sum90deg/mon90deg, color='g', marker='*', markersize=5, linestyle='None', label=r"$\omega=90^\circ$")

    ymin, ymax = 2e-07, 2.5e-06
    plt.axvline(x=0.0125, color='yellow', linestyle='--')
    plt.axvline(x=0.02, color='yellow', linestyle='--')
    plt.fill_between(qy, ymin, ymax, where=np.logical_and(qy >= 0.0125, qy <= 0.02), facecolor='yellow', alpha=0.3)
    plt.text(0.0135, 1.2e-06, "300-500 nm", fontsize=14)
    plt.axvline(x=0.003, color='0.7', linestyle=':')
    plt.axvline(x=0.006, color='0.7', linestyle=':')
    plt.fill_between(qy, ymin, ymax, where=np.logical_and(qy >= 0.003, qy <= 0.006), facecolor='0.7', alpha=0.3)
    plt.text(0.003, 6.0e-07, r"1-2 $\mu$m", fontsize=14)
    plt.axvline(x=-0.0125, color='yellow', linestyle='--')
    plt.axvline(x=-0.02, color='yellow', linestyle='--')
    plt.text(-0.019, 1.2e-06, "300-500 nm", fontsize=14)
    plt.fill_between(qy, ymin, ymax, where=np.logical_and(qy >= -0.02, qy <= -0.0125), facecolor='yellow', alpha=0.3)
    plt.axvline(x=-0.003, color='0.7', linestyle=':')
    plt.axvline(x=-0.006, color='0.7', linestyle=':')
    plt.fill_between(qy, ymin, ymax, where=np.logical_and(qy >= -0.006, qy <= -0.003), facecolor='0.7', alpha=0.3)
    plt.text(-0.006, 6.0e-07, r"1-2 $\mu$m", fontsize=14)

    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$I$ (a. u.)')
    plt.ylim([ymin, ymax])
    plt.title(r"Sample rotation summary")
    plt.legend(loc='upper right', fontsize=14)

    plt.show()


def get_monitor_counts(fname):
    yaml = YAML()
    with open(join(datapath, fname)) as f:
        metadata = yaml.load(f)
        for d in metadata['measurement']['devices']:
            if d['name'] == 'mon1':
                return d['value'][0]


def plot_0deg_summary():
    # load data
    data23 = load_data(join(datapath, f66023))   # 0 degrees
    data24 = load_data(join(datapath, f66024))   # 0 degrees
    data55 = load_data(join(datapath, f66055))   # 0 degrees
    data52 = load_data(join(datapath, f65952))   # 0 degrees
    data53 = load_data(join(datapath, f65953))   # 0 degrees
    data_large_slit = data52 + data53
    data_small_slit = data23 + data24

    # get monitor counts
    mon55 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66055)[0]))
    mon23 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66023)[0]))
    mon24 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66024)[0]))
    mon52 = get_monitor_counts("{}.yaml".format(os.path.splitext(f65952)[0]))
    mon53 = get_monitor_counts("{}.yaml".format(os.path.splitext(f65953)[0]))

    simulation = get_simulation()
    axes_limits = ba.get_axes_limits(simulation.result(), ba.AxesUnits.QSPACE)
    qy = np.linspace(axes_limits[0], axes_limits[1], data23.shape[0])

    simulation_1 = get_simulation(u0=45.0)
    axes_limits_1 = ba.get_axes_limits(simulation_1.result(), ba.AxesUnits.QSPACE)
    qy_1hour = np.linspace(axes_limits_1[0], axes_limits_1[1], data55.shape[0])

    sum_large_slit = np.sum(data_large_slit, axis=0)/(mon52 + mon53)
    sum_small_slit = np.sum(data_small_slit, axis=0)/(mon23 + mon24)
    sum_1hour = np.sum(data55, axis=0)/mon55

    plt.style.use('seaborn-talk')
    plt.subplot(1, 2, 1)
    plt.semilogy(qy, sum_small_slit, color='k', marker='o', markersize=5, linestyle='None', label=r"Small slit")
    plt.semilogy(qy, sum_large_slit, color='b', marker='s', markersize=5, linestyle='None', label=r"Large slit")
    plt.semilogy(qy_1hour, sum_1hour, color='g', marker='*', markersize=5, linestyle='None', label=r"1 hour")

    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$I$ (a. u.)')
    plt.ylim([2.0e-07, 4.0e-05])
    plt.title(r"Sample rotation $0^\circ$ summary", fontsize=16)
    plt.legend(loc='upper right', fontsize=14)

    plt.subplot(1, 2, 2)
    img = plt.imread(join(datapath, 'rot0deg.jpg'))
    plt.imshow(img)
    plt.axis('off')

    plt.show()


def plot_90deg_summary():
    # load data
    data29 = load_data(join(datapath, f66029))   # 90 degrees
    data57 = load_data(join(datapath, f66057))   # 90 degrees
    data58 = load_data(join(datapath, f66058))   # 90 degrees
    data_1hour = data57 + data58

    # get monitor counts
    mon57 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66057)[0]))
    mon58 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66058)[0]))
    mon29 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66029)[0]))

    simulation = get_simulation()
    axes_limits = ba.get_axes_limits(simulation.result(), ba.AxesUnits.QSPACE)
    qy = np.linspace(axes_limits[0], axes_limits[1], data29.shape[0])

    simulation_1 = get_simulation(u0=46.0)
    axes_limits_1 = ba.get_axes_limits(simulation_1.result(), ba.AxesUnits.QSPACE)
    qy_1hour = np.linspace(axes_limits_1[0], axes_limits_1[1], data57.shape[0])

    sum_29 = np.sum(data29, axis=0)/mon29
    sum_1hour = np.sum(data_1hour, axis=0)/(mon57 + mon58)

    plt.style.use('seaborn-talk')
    plt.subplot(1, 2, 1)
    plt.semilogy(qy, sum_29, color='k', marker='o', markersize=5, linestyle='None', label=r"Over the night")
    plt.semilogy(qy_1hour, sum_1hour, color='b', marker='s', markersize=5, linestyle='None', label=r"1 hour")

    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$I$ (a. u.)')
    plt.ylim([2.0e-07, 5.0e-06])
    plt.title(r"Sample rotation $90^\circ$ summary", fontsize=16)
    plt.legend(loc='upper right', fontsize=14)

    plt.subplot(1, 2, 2)
    img = plt.imread(join(datapath, 'rot90deg.jpg'))
    plt.imshow(img)
    plt.axis('off')

    plt.show()


def plot_45deg_summary():
    # load data
    data28 = load_data(join(datapath, f66028))   # 45 degrees

    # get monitor counts
    mon28 = get_monitor_counts("{}.yaml".format(os.path.splitext(f66028)[0]))

    simulation = get_simulation()
    axes_limits = ba.get_axes_limits(simulation.result(), ba.AxesUnits.QSPACE)
    qy = np.linspace(axes_limits[0], axes_limits[1], data28.shape[0])

    sum_28 = np.sum(data28, axis=0)/mon28

    plt.style.use('seaborn-talk')
    plt.subplot(1, 2, 1)
    plt.semilogy(qy, sum_28, color='k', marker='o', markersize=5, linestyle='None', label=r"Over the night")

    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$I$ (a. u.)')
    plt.ylim([2.0e-07, 2.0e-06])
    plt.title(r"Sample rotation $45^\circ$ summary", fontsize=16)
    plt.legend(loc='upper right', fontsize=14)

    plt.subplot(1, 2, 2)
    img = plt.imread(join(datapath, 'rot45deg.jpg'))
    plt.imshow(img)
    plt.axis('off')

    plt.show()


if __name__ == '__main__':
    # plot_large_mirror_raw()
    # plot_all_2d()
    # plot_all_1d()
    plot_0deg_summary()
    plot_45deg_summary()
    plot_90deg_summary()
    plot_summary()
