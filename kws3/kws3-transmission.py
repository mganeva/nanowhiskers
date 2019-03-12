import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import bornagain as ba
from os.path import join
from ruamel.yaml import YAML


sensitivity_matrix = "/datadisk/data/KWS3/p15749/sens.sens"
datapath = 'kws3/transmission'

# transmission data files
ec_fname1 = '00064696_0000_p15749_EC_HRD_standard.{}'          # 5x5 mm, 9.5 m
ec_fname2 = '00064699_0000_p15749_EC_HRD_standard.{}'          # 5x5 mm, 9.5 m
ec_fname3 = '00064701_0000_p15749_EC_HRD_standard.{}'          # 5x5 mm, 9.5 m, 1 hour
s3278_fname1 = '00064697_0000_p15749_S-3278_HRD_standard.{}'   # 5x5 mm, 9.5 m
s3278_fname2 = '00064702_0000_p15749_S-3278_HRD_standard.{}'   # 5x5 mm, 9.5 m, 1 hour
s3098_fname1 = '00064698_0000_p15749_S-3098_HRD_standard.{}'   # 5x5 mm, 9.5 m
s3098_fname2 = '00064700_0000_p15749_S-3098_HRD_standard.{}'   # 5x5 mm, 9.5 m, 1 hour


def get_ec_factor():
    """
    given: measured EC of thickness d2=0.05 cm
    measured sample of thickness d1=0.038 cm
    attenuation law: I = I0*exp(-Sigma*d)
    thus, I1/I2 = I1*exp(-Sigma*(d2 - d1))
    Sigma_Si = 0.169 cm^-1 (https://www.ncnr.nist.gov/resources/activation/)
    :return: exp(-Sigma*(d2 - d1))
    """
    return np.exp(-0.169*(0.05 - 0.038))


def load_data(fname):
    rawdata = np.loadtxt(fname)
    sens = np.loadtxt(sensitivity_matrix)
    data = rawdata*sens

    return np.rot90(data, 3)


def plot_2d_raw(fname, zmax=100.0):
    data = load_data(join(datapath, fname.format('det')))
    plt.imshow(data, norm=matplotlib.colors.LogNorm(0.1, zmax), aspect='auto', cmap='jet')
    plt.show()


def plot_2d_sum(fnames, zmax=100.0):
    result = np.zeros((256, 256))
    for fname in fnames:
        data = load_data(join(datapath, fname.format('det')))
        result += data

    plt.imshow(result, norm=matplotlib.colors.LogNorm(0.1, zmax), aspect='auto', cmap='jet')
    plt.show()


def get_monitor_counts(fname):
    yaml = YAML()
    with open(join(datapath, fname.format('yaml'))) as f:
        metadata = yaml.load(f)
        for d in metadata['measurement']['devices']:
            if d['name'] == 'mon1':
                return d['value'][0]


def plot_2d_reduced(fname, zmax=100.0):
    # load data
    data = load_data(join(datapath, fname.format('det')))
    # load empty can
    ec_data = load_data(join(datapath, ec_fname3.format('det')))

    # read monitor counts
    ec_mon1 = get_monitor_counts(ec_fname3)
    data_mon1 = get_monitor_counts(fname)

    # normalize data and subtract empty can
    data_norm = data/data_mon1
    ec_norm = ec_data/ec_mon1
    factor = get_ec_factor()
    ec_scaled = ec_norm/factor
    result = data_norm - ec_scaled
    np.savetxt("qqq.txt", result)

    plt.imshow(result, norm=matplotlib.colors.LogNorm(1.0e-11, zmax), aspect='auto', cmap='jet')
    plt.show()


def get_kws3_detector():
    """
    Creates and returns KWS-3 detector
    detector dimensions: 90x90 mm
    SDD: 9500 mm
    """
    u0 = 45.7  # in mm
    v0 = 48.5  # in mm
    detector = ba.RectangularDetector(256, 90.0, 256, 90.0)
    detector.setPerpendicularToDirectBeam(9500.0, u0, v0)
    return detector


def get_simulation():
    """
    Returns a GISAXS simulation with beam and detector defined
    wavelength 12.8 angstrom
    incident angle 0.0 degrees (transmission)
    """
    simulation = ba.GISASSimulation()
    simulation.setBeamParameters(12.8*ba.angstrom, 90.0*ba.deg, 0.0*ba.deg)
    simulation.setDetector(get_kws3_detector())
    simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(5.0, 5.0))
    simulation.setBeamIntensity(1.0e-4)
    distr_1 = ba.DistributionGaussian(1.28*ba.nm, 0.1)
    simulation.addParameterDistribution("*/Beam/Wavelength", distr_1, 50, 2.0, ba.RealLimits.positive())
    simulation.getOptions().setIncludeSpecular(True)
    return simulation


def reduce(fname1, fname2):
    # load data
    raw1 = load_data(join(datapath, fname1.format('det')))
    raw2 = load_data(join(datapath, fname2.format('det')))
    raw = raw1 + raw2  # to get more statistics

    # load empty can
    ec_data = load_data(join(datapath, ec_fname3.format('det')))

    # read monitor counts
    ec_mon1 = get_monitor_counts(ec_fname3)
    mon1 = get_monitor_counts(fname1) + get_monitor_counts(fname2)

    # normalize data and subtract empty can
    normalized = raw/mon1
    ec_norm = ec_data/ec_mon1
    factor = get_ec_factor()
    ec_scaled = 0.95*ec_norm/factor                  # to avoid oversubtracting

    # reduced data arrays
    data = normalized  - ec_scaled
    return data


def get_i_of_q(data, nbins=256):
    """
    calculates I(Q) where Q = sqrt(Qx^2 + Qy^2)
    :param data:
    :return:
    """
    simulation = get_simulation()
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)
    shape = data.shape
    x = np.linspace(axes_limits[0], axes_limits[1], shape[0])
    y = np.linspace(axes_limits[2], axes_limits[3], shape[1])

    # xx, yy = np.meshgrid (x, y)
    # q = np.sqrt(xx**2 + yy**2)
    # result = np.array([q.flatten(), data.flatten()]).transpose()
    result = []
    for i in range(shape[0]):
        for j in range(shape[1]):
            q = np.sqrt(x[i]**2 + y[j]**2)
            result.append([q, data[i, j]])
    result = np.array(result)
    # sort data
    result = result[result[:, 0].argsort()]

    # bin data
    bins = np.linspace(0.0, np.max(q), nbins)
    indices = np.digitize(result[:, 0], bins)

    a = []
    for i in range(bins.size):
        idx = np.where(indices==i)
        if idx[0].size > 0:
            a.append([bins[i], np.mean(result[:, 1][idx])])
    return np.array(a)


def plot_2d_reduced_qspace(zmin=1.0e-09):
    # reduced data arrays
    s3098_data = reduce(s3098_fname1, s3098_fname2)
    s3278_data = reduce(s3278_fname1, s3278_fname2)

    # convert to qspace
    simulation = get_simulation()
    s3098_hist = simulation.result().histogram2d(ba.AxesUnits.QSPACE)
    s3278_hist = simulation.result().histogram2d(ba.AxesUnits.QSPACE)
    s3098_hist.setContent(s3098_data)
    s3278_hist.setContent(s3278_data)

    # plot the result
    plt.style.use('seaborn-talk')
    grid = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3)

    # ==============
    # S3098
    # ===============
    ax = plt.subplot(grid[0, 0])
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    im = plt.imshow(s3098_hist.array(),
                    norm=matplotlib.colors.LogNorm(zmin, 01.0e-06),
                    extent=[s3098_hist.getXmin(), s3098_hist.getXmax(), s3098_hist.getYmin(), s3098_hist.getYmax()],
                    aspect='auto', cmap='jet')

    cb = plt.colorbar(im)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$Q_x$ (nm$^{-1}$)')
    plt.title("S3098")

    # ==============
    # S3278
    # ===============
    ax = plt.subplot(grid[0, 1])
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    im = plt.imshow(s3278_hist.array(),
                    norm=matplotlib.colors.LogNorm(zmin, 01.0e-06),
                    extent=[s3278_hist.getXmin(), s3278_hist.getXmax(), s3278_hist.getYmin(), s3278_hist.getYmax()],
                    aspect='auto', cmap='jet')

    cb = plt.colorbar(im)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$Q_x$ (nm$^{-1}$)')
    plt.title("S3278")
    # =======
    # I(Q)
    # =======
    plt.subplot(grid[1, 0:])
    s3098_slice = get_i_of_q(s3098_data)
    s3278_slice = get_i_of_q(s3278_data)
    plt.semilogy(s3098_slice[:, 0], s3098_slice[:, 1], color='k', marker='.', markersize=7, linestyle='None', label='S3098')
    plt.semilogy(s3278_slice[:, 0], s3278_slice[:, 1], color='b', marker='o', markersize=5, linestyle='None', label='S3278')
    plt.axvline(x=0.003, color='0.7', linestyle='--', linewidth=1)
    plt.axvline(x=0.0055, color='r', linestyle='--', linewidth=1)
    plt.axvline(x=0.006, color='0.7', linestyle='--', linewidth=1)
    plt.annotate(r'{} $\mu$m'.format(1.14), xy=(0.0055, 1.3e-07), xytext=(0.007, 4.0e-07),
                 fontsize=14, arrowprops=dict(facecolor='red', shrink=0.0, width=0.5, headwidth=4))
    plt.axvline(x=0.020, color='g', linestyle='--', linewidth=1)
    plt.annotate(r'{} nm'.format(314), xy=(0.02, 0.17e-07), xytext=(0.017, 0.4e-07),
                 fontsize=14, arrowprops=dict(facecolor='green', shrink=0.0, width=0.5, headwidth=4))
    plt.axvline(x=0.0115, color='g', linestyle='--', linewidth=1)
    plt.annotate(r'{} nm'.format(546), xy=(0.0115, 0.17e-07), xytext=(0.013, 0.4e-07),
                 fontsize=14, arrowprops=dict(facecolor='green', shrink=0.0, width=0.5, headwidth=4))
    x = s3278_slice[:, 0]
    plt.fill_between(x, 0, 1, where=np.logical_and(x >= 0.003, x <= 0.0061), facecolor='0.7', alpha=0.5)
    plt.xlim(0.0, 0.024)
    plt.ylim(1.0e-10, 1.0e-06)
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.xlabel(r'$Q$ (nm$^{-1}$)')
    plt.ylabel(r'$I(Q)$, a.u.')
    plt.legend(loc='upper right', fontsize=18)

    plt.show()


def get_sample(radius=153.4, d=1000.0):
    # Defining Materials
    material_1 = ba.HomogeneousMaterial("Air", 0.0, 0.0)
    material_3 = ba.MaterialBySLD("Si", 2.0737e-06, -2.3758e-11)
    material_2 = ba.MaterialBySLD("Au", 4.6665e-06, -1.6205e-08)

    # Defining Layers
    layer_1 = ba.Layer(material_1)

    # Defining Form Factors
    formFactor = ba.FormFactorCone6(radius * ba.nm, 300.0 * ba.nm, 81.0 * ba.deg)

    # Defining Particles
    particle = ba.Particle(material_3, formFactor)

    interference = ba.InterferenceFunctionRadialParaCrystal(d*ba.nm, 2e3*ba.nm)
    pdf = ba.FTDistribution1DGauss(250.0 * ba.nm)
    interference.setProbabilityDistribution(pdf)

    # Defining Particle Layouts and adding Particles
    layout_1 = ba.ParticleLayout()
    layout_1.addParticle(particle, 1.0)
    layout_1.setTotalParticleSurfaceDensity(0.01)
    layout_1.setInterferenceFunction(interference)

    # Adding layouts to layers
    layer_1.addLayout(layout_1)

    # Defining Multilayers
    multiLayer_1 = ba.MultiLayer()
    multiLayer_1.addLayer(layer_1)
    return multiLayer_1


def plot_sim_qspace(expdata, title, radius=150.0, d=1000.0, zmin=1.0e-09):
    simulation = get_simulation()

    # experimental data
    s_hist = simulation.result().histogram2d(ba.AxesUnits.QSPACE)
    s_hist.setContent(expdata)

    # run simulation
    sample = get_sample(radius=radius, d=d)
    simulation.setSample(sample)
    simulation.setTerminalProgressMonitor()
    simulation.runSimulation()
    result = simulation.result()
    axes_limits = ba.get_axes_limits(result, ba.AxesUnits.QSPACE)

    # plot the result
    plt.style.use('seaborn-talk')
    grid = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3)
    # ==============
    # experimental data
    # ===============
    ax = plt.subplot(grid[0, 0])
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    im = plt.imshow(s_hist.array(),
                    norm=matplotlib.colors.LogNorm(zmin, 01.0e-06),
                    extent=axes_limits,
                    aspect='auto', cmap='jet')

    cb = plt.colorbar(im)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$Q_x$ (nm$^{-1}$)')
    plt.title("{} experiment".format(title))

    # ==============
    # Simulation
    # ===============
    ax = plt.subplot(grid[0, 1])
    plt.subplots_adjust(wspace=0.2, hspace=0.3)

    im = plt.imshow(result.array(),
                    norm=matplotlib.colors.LogNorm(zmin, 01.0e-06),
                    extent=axes_limits, aspect='auto', cmap='jet')

    cb = plt.colorbar(im)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$Q_x$ (nm$^{-1}$)')
    plt.title("Simulation")

    # =======
    # I(Q)
    # =======
    plt.subplot(grid[1, 0:])
    qy = np.array(result.axis(0, ba.AxesUnits.QSPACE))
    exp_slice = np.mean(s_hist.array(), axis=0)
    sim_slice = np.mean(result.array(), axis=0)
    plt.semilogy(qy, exp_slice, color='k', marker='.', markersize=7, linestyle='None', label='Experiment')
    plt.semilogy(qy, sim_slice, color='b', linestyle='-', label='Simulation')
    plt.axvline(x=0.003, color='0.7', linestyle='--', linewidth=1)
    plt.axvline(x=0.0055, color='r', linestyle='--', linewidth=1)
    plt.axvline(x=0.006, color='0.7', linestyle='--', linewidth=1)
    plt.annotate(r'{} $\mu$m'.format(1.14), xy=(0.0055, 1.3e-07), xytext=(0.007, 4.0e-07),
                 fontsize=14, arrowprops=dict(facecolor='red', shrink=0.0, width=0.5, headwidth=4))
    plt.axvline(x=0.020, color='g', linestyle='--', linewidth=1)
    plt.annotate(r'{} nm'.format(314), xy=(0.02, 0.17e-07), xytext=(0.017, 0.4e-07),
                 fontsize=14, arrowprops=dict(facecolor='green', shrink=0.0, width=0.5, headwidth=4))
    plt.axvline(x=0.0115, color='g', linestyle='--', linewidth=1)
    plt.annotate(r'{} nm'.format(546), xy=(0.0115, 0.17e-07), xytext=(0.013, 0.4e-07),
                 fontsize=14, arrowprops=dict(facecolor='green', shrink=0.0, width=0.5, headwidth=4))
    plt.fill_between(qy, 0, 1, where=np.logical_and(qy >= 0.0029, qy <= 0.0061), facecolor='0.7', alpha=0.5)
    plt.xlim(0.0, 0.024)
    plt.ylim(1.0e-10, 1.0e-06)
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.xlabel(r'$Q_y$ (nm$^{-1}$)')
    plt.ylabel(r'$I(Q_y)$, a.u.')
    plt.legend(loc='upper right', fontsize=18)
    plt.title(r"Integrated over $Q_x$")

    plt.show()


if __name__ == '__main__':
    # plot_2d_sum([ec_fname1, ec_fname2, ec_fname3], 100)
    # plot_2d_raw(ec_fname3, 100000)
    # plot_2d_reduced(s3278_fname1, 1.0e-06)
    plot_2d_reduced_qspace()
    s3098_data = reduce(s3098_fname1, s3098_fname2)
    plot_sim_qspace(s3098_data, 'S3098', 546.0, 1500.0)
