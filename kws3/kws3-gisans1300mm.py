import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import bornagain as ba
from os.path import join
from ruamel.yaml import YAML


sensitivity_matrix = "/datadisk/data/KWS3/p15749/sens.sens"
datapath = '/home/mary/build/bornagain/users/uliana/kws3/gisans1m'

# critical angles
a_Si = 0.596   # degrees
a_Au = 0.878   # degrees
a_Cu = 1.059   # degrees


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


def plot_2d_sum(fnames, zmax=100.0):
    result = np.zeros((256, 256))
    for fname in fnames:
        data = load_data(join(datapath, fname))
        result += data

    plt.imshow(result, norm=matplotlib.colors.LogNorm(0.1, zmax), aspect='auto', cmap='jet')
    plt.show()


def get_kws3_detector():
    """
    Creates and returns KWS-3 detector
    detector dimensions: 90x90 mm
    SDD: 9300 mm (for SANS)
    """
    u0 = 45.3  # in mm
    v0 = 30.2  # in mm
    detector = ba.RectangularDetector(256, 90.0, 256, 90.0)
    detector.setPerpendicularToDirectBeam(1300.0, u0, v0)
    return detector


def get_simulation(ai=0.56):
    """
    Returns a GISAXS simulation with beam and detector defined
    wavelength 12.8 angstrom
    incident angle 0.0 degrees (transmission)
    """
    simulation = ba.GISASSimulation()
    simulation.setBeamParameters(12.8*ba.angstrom, ai*ba.deg, 0.0*ba.deg)
    simulation.setDetector(get_kws3_detector())
    simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(5.0, 5.0))
    simulation.setBeamIntensity(1000)
    distr_1 = ba.DistributionGaussian(1.28*ba.nm, 0.1)
    simulation.addParameterDistribution("*/Beam/Wavelength", distr_1, 50, 2.0, ba.RealLimits.positive())
    simulation.getOptions().setIncludeSpecular(True)
    return simulation


def get_rotation_angle(fname):
    yaml = YAML()
    with open(join(datapath, fname.format('yaml'))) as f:
        metadata = yaml.load(f)
        for d in metadata['measurement']['devices']:
            if d['name'] == 'sam_rot':
                return d['value']


def plot_s3278_2d():
    fname = "0006543{r1}_000{r2}_p15749_S3278-GiSANS-0.6dgr-rot-scan-m15-15dgr_HRD_standard.{ext}"
    runs = [(0, 1), (1, 1), (2, 2), (3, 3), (4, 4), (5, 1), (6, 2), (7, 3), (8, 0)]
    simulation = get_simulation()
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)

    plt.style.use('seaborn-talk')
    plt.suptitle("Sample S3278 SDD=1300 mm", fontsize=16)
    for i in runs:
        data = load_data(join(datapath, fname.format(r1=i[0], r2=i[1], ext='det')))
        omega = get_rotation_angle(join(datapath, fname.format(r1=i[0], r2=i[1], ext='yaml')))
        ax = plt.subplot(3, 3, i[0]+1)
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        im = plt.imshow(data,
                        norm=matplotlib.colors.LogNorm(0.1, 1000),
                        extent=axes_limits,
                        aspect='auto', cmap='jet')

        cb = plt.colorbar(im)
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.ylabel(r'$Q_z$ (nm$^{-1}$)')
        plt.title(r"Run 0006543{n}, $\omega={om}^\circ$".format(n=i[0], om=omega))

    plt.show()


def plot_s3278_qz():
    fname = "0006543{r1}_000{r2}_p15749_S3278-GiSANS-0.6dgr-rot-scan-m15-15dgr_HRD_standard.{ext}"
    runs = [(0, 1), (1, 1), (2, 2), (3, 3), (4, 4), (5, 1), (6, 2), (7, 3), (8, 0)]
    simulation = get_simulation()
    result = simulation.result()

    plt.style.use('seaborn-talk')
    plt.suptitle(r"Sample S3278 SDD=1300 mm, slice along $\alpha_f$", fontsize=16)
    for i in runs:
        data = load_data(join(datapath, fname.format(r1=i[0], r2=i[1], ext='det')))
        omega = get_rotation_angle(join(datapath, fname.format(r1=i[0], r2=i[1], ext='yaml')))
        hist = result.histogram2d(ba.AxesUnits.DEGREES)
        hist.setContent(data)
        ax = plt.subplot(3, 3, i[0]+1)
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        zslice = hist.projectionY(0.0)
        plt.semilogy(zslice.getBinCenters(), zslice.getBinValues(), color='k', marker='.', markersize=5, linestyle='None')
        plt.axvline(x=0.56, color='r', linestyle='--')
        ycoord = 12000
        if i[0] == 8:
            ycoord=5000
        plt.annotate(r'$\alpha_i=0.56^{\circ}$', xy=(0.56, ycoord), xytext=(0.9, 2000),
                     fontsize=14, arrowprops=dict(facecolor='red', shrink=0.0, width=0.5, headwidth=4))
        plt.xlabel(r'$\alpha_f$ ($^{\circ}$)')
        plt.title(r"Run 0006543{n}, $\omega={om}^\circ$".format(n=i[0], om=omega))

    plt.show()


def plot_s3278_2d_roi():
    fname = "0006543{r1}_000{r2}_p15749_S3278-GiSANS-0.6dgr-rot-scan-m15-15dgr_HRD_standard.{ext}"
    runs = [(0, 1), (1, 1), (2, 2), (3, 3), (4, 4), (5, 1), (6, 2), (7, 3), (8, 0)]
    simulation = get_simulation()
    result = simulation.result()
    axes_limits = [-0.05, 0.05, 0.06, 0.205]

    plt.style.use('seaborn-talk')
    plt.suptitle("Sample S3278 SDD=1300 mm", fontsize=16)
    for i in runs:
        data = load_data(join(datapath, fname.format(r1=i[0], r2=i[1], ext='det')))
        omega = get_rotation_angle(join(datapath, fname.format(r1=i[0], r2=i[1], ext='yaml')))
        hist = result.histogram2d(ba.AxesUnits.QSPACE)
        hist.setContent(data)
        h = hist.crop(-0.05, 0.06, 0.05, 0.205)
        ax = plt.subplot(3, 3, i[0]+1)
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        im = plt.imshow(h.array(),
                        norm=matplotlib.colors.LogNorm(0.1, 500),
                        extent=axes_limits,
                        aspect='auto', cmap='jet')

        cb = plt.colorbar(im)
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.ylabel(r'$Q_z$ (nm$^{-1}$)')
        plt.title(r"Run 0006543{n}, $\omega={om}^\circ$".format(n=i[0], om=omega))

    plt.show()


def plot_s3098_1deg(ai=1.2):
    fname = '00065{r}_00{s:02d}_p15749_SiO2,_Si,_Au,_Cu_SiO2,_Si,_Au,_Cu_HRD_standard.{ext}'
    runs = range(1)  # put 22 to view the all 22 files
    simulation = get_simulation(ai=ai)
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)
    # ROI, nm^-1
    xmin, ymin, xmax, ymax = -0.05, 0.02, 0.05, 0.17
    qSi = calc_q(ai, a_Si)
    qAu = calc_q(ai, a_Au)
    qCu = calc_q(ai, a_Cu)

    plt.style.use('seaborn-talk')
    grid = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3)
    for i in runs:
        data = load_data(join(datapath, fname.format(r=483+i, s=i, ext='det')))
        omega = get_rotation_angle(join(datapath, fname.format(r=483+i, s=i, ext='yaml')))

        plt.suptitle(r"Run 00065{n}, $\omega={om}^\circ$".format(n=483+i, om=omega), fontsize=16)
        # ===================
        # raw data
        # ===================
        plt.subplot(grid[0, 0])
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        im = plt.imshow(data,
                        norm=matplotlib.colors.LogNorm(0.1, 500),
                        extent=axes_limits,
                        aspect='auto', cmap='jet')

        cb = plt.colorbar(im)
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.ylabel(r'$Q_z$ (nm$^{-1}$)')
        plt.title(r"Raw data")

        # ===================
        # raw data ROI
        # ===================
        plt.subplot(grid[0, 1])
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        hist = result.histogram2d(ba.AxesUnits.QSPACE)
        hist.setContent(data)
        h = hist.crop(xmin, ymin, xmax, ymax)
        im = plt.imshow(h.array(),
                        norm=matplotlib.colors.LogNorm(0.1, 500),
                        extent=[xmin, xmax, ymin, ymax],
                        aspect='auto', cmap='jet')

        cb = plt.colorbar(im)
        plt.axhline(y=qSi, color='0.3', linestyle='--', linewidth=1)
        plt.axhline(y=qAu, color='b', linestyle='--', linewidth=1)
        plt.axhline(y=qCu, color='g', linestyle='--', linewidth=1)
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.ylabel(r'$Q_z$ (nm$^{-1}$)')
        plt.title(r"Raw data (ROI)")

        # ===================
        # qz slice (degrees)
        # ===================
        plt.subplot(grid[1, 0])
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        hist = result.histogram2d(ba.AxesUnits.DEGREES)
        hist.setContent(data)
        zslice = hist.projectionY(0.0)
        plt.semilogy(zslice.getBinCenters(), zslice.getBinValues(), color='k', marker='o', markersize=5, linestyle='None')
        plt.axvline(x=ai, color='r', linestyle='--')
        plt.annotate(r'$\alpha_i={}^\circ$'.format(ai), xy=(ai, 300), xytext=(0.9, 500),
                     fontsize=14, arrowprops=dict(facecolor='red', shrink=0.0, width=0.5, headwidth=4))
        plt.xlabel(r'$\alpha_f$ ($^{\circ}$)')
        plt.title(r"Slice along $\alpha_f$")

        # ===================
        # qy slice (nm^-1)
        # ===================
        plt.subplot(grid[1, 1])
        plt.subplots_adjust(wspace=0.2, hspace=0.3)
        hist = result.histogram2d(ba.AxesUnits.QSPACE)
        hist.setContent(data)
        ysliceSi = hist.projectionX(qSi)
        ysliceCu = hist.projectionX(qCu)
        ysliceAu = hist.projectionX(qAu)
        plt.semilogy(ysliceSi.getBinCenters(), ysliceSi.getBinValues(), color='0.3', marker='o', markersize=7,
                     linestyle='None', label=r'$Q_c(Si)={:.3f} $'.format(qSi) + r'nm$^{-1}$')
        plt.semilogy(ysliceAu.getBinCenters(), ysliceAu.getBinValues(), color='b', marker='s', markersize=7,
                     linestyle='None', label=r'$Q_c$(Au)={:.3f} '.format(qAu) + r'nm$^{-1}$')
        plt.semilogy(ysliceCu.getBinCenters(), ysliceCu.getBinValues(), color='g', marker='*', markersize=9,
                     linestyle='None', label=r'$Q_c$(Cu)={:.3f} '.format(qCu) + r'nm$^{-1}$')
        plt.xlim([0.0, 0.2])
        plt.ylim([0.7, 10])
        plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
        plt.title(r"Slices along $Q_y$")
        plt.legend(loc='upper right', fontsize=12)

        plt.show()


if __name__ == '__main__':
    # plot_s3278_2d()
    # plot_s3278_qz()
    # plot_s3278_2d_roi()
    plot_s3098_1deg(0.57)
