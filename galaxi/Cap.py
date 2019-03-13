import numpy
import bornagain as ba
from bornagain import deg, angstrom, nm, kvector_t
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['image.cmap'] = 'jet'


def get_sample():
    # Defining Materials
    material_1 = ba.HomogeneousMaterial("example01_Air", 0.0, 0.0)
    material_2 = ba.HomogeneousMaterial("Au", 3.53665637e-05, 2.9383311e-06)
    material_3 = ba.HomogeneousMaterial("Si", 5.73327e-06, 1.006366e-07)

    # Defining Layers
    layer_1 = ba.Layer(material_1)
    layer_2 = ba.Layer(material_3)

    # Defining Form Factors
    formFactor_1 = ba.FormFactorTruncatedSphere(159.0*nm, 244.0*nm, 106.0*nm)

    # Defining Particles
    particle_1 = ba.Particle(material_2, formFactor_1)
    particle_1_position = kvector_t(0.0*nm, 0.0*nm, 333.0*nm)
    particle_1.setPosition(particle_1_position)

    # Defining Particle Layouts and adding Particles
    layout_1 = ba.ParticleLayout()
    layout_1.addParticle(particle_1, 1.0)
    layout_1.setTotalParticleSurfaceDensity(0.01)

    # Adding layouts to layers
    layer_1.addLayout(layout_1)

    # Defining Multilayers
    multiLayer_1 = ba.MultiLayer()
    multiLayer_1.addLayer(layer_1)
    multiLayer_1.addLayer(layer_2)
    return multiLayer_1


def get_simulation():
    simulation = ba.GISASSimulation()
    
    detector = ba.RectangularDetector(981, 168.732, 1043, 179.396)
    detector.setPerpendicularToDirectBeam(3530.0, 104.558, 39.422)
    simulation.setDetector(detector)
    
    simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(0.043, 0.043))
    simulation.setBeamParameters(0.1341*nm, 0.2*deg, 0.0*deg)
    simulation.setBeamIntensity(5.0e+06)
    simulation.setTerminalProgressMonitor()
    simulation.setRegionOfInterest(40.0, 50.0, 160.0, 170.0)
    return simulation


def run_simulation():
    sample = get_sample()
    simulation = get_simulation()
    simulation.setSample(sample)
    simulation.runSimulation()
    return simulation.result()


def plot(result):
    plt.figure(figsize=(12.80, 10.24))
    plt.subplot()
    ba.plot_colormap(result, units=ba.AxesUnits.QSPACE, title="Q-space",
                     xlabel=r'$Q_{y} [1/nm]$', ylabel=r'$Q_{z} [1/nm]$', zmin=0.2, zmax=1e+5)
    plt.savefig('cup.jpg')


if __name__ == '__main__':
    result = run_simulation()
    plot(result)