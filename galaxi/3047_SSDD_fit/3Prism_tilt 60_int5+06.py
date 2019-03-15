import numpy
import bornagain as ba
from bornagain import deg, angstrom, nm, kvector_t
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['image.cmap'] = 'jet'


j=0
while j<361:
    def get_sample():
        # Defining Materials
        material_1 = ba.HomogeneousMaterial("Air", 0.0, 0.0)
        material_2 = ba.HomogeneousMaterial("Si", 5.73327e-06, 1.006366e-07)

        # Defining Layers
        layer_1 = ba.Layer(material_1)
        layer_2 = ba.Layer(material_2)

        # Defining Form Factors
        formFactor_1 = ba.FormFactorPrism6(159.0 * nm, 300.0 * nm)

        # Defining Particles
        particle_1 = ba.Particle(material_2, formFactor_1)
        particle_1_rotation = ba.RotationY(60.0 * deg)
        particle_1.setRotation(particle_1_rotation)
        particle_1_position = kvector_t(0.0 * nm, 0.0 * nm, 138 * nm)
        particle_1.setPosition(particle_1_position)

        # Defining composition of particles at specific positions
        z1 = j
        z2 = j+120
        z3 = j+240

        particleComposition_1 = ba.ParticleComposition()
        particleComposition_1.addParticle(particle_1)
        particleComposition_1_rotation = ba.RotationZ(z1 * deg)
        particleComposition_1.setRotation(particleComposition_1_rotation)
        particleComposition_2 = ba.ParticleComposition()
        particleComposition_2.addParticle(particle_1)
        particleComposition_2_rotation = ba.RotationZ(z2 * deg)
        particleComposition_2.setRotation(particleComposition_2_rotation)
        particleComposition_3 = ba.ParticleComposition()
        particleComposition_3.addParticle(particle_1)
        particleComposition_3_rotation = ba.RotationZ(z3 * deg)
        particleComposition_3.setRotation(particleComposition_3_rotation)


        # Defining Particle Layouts and adding Particles
        layout_1 = ba.ParticleLayout()
        layout_1.addParticle(particleComposition_1, 1.0)
        layout_1.addParticle(particleComposition_2, 1.0)
        layout_1.addParticle(particleComposition_3, 1.0)

        layout_1.setTotalParticleSurfaceDensity(0.001)

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
        detector.setPerpendicularToDirectBeam(830.0, 106.1068, 47.644)
        simulation.setDetector(detector)

        simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(0.043, 0.043))
        simulation.setBeamParameters(0.1341 * nm, 0.2 * deg, 0.0 * deg)
        simulation.setBeamIntensity(5e+06)
        simulation.setTerminalProgressMonitor()
        simulation.setRegionOfInterest(0, 48, 168.73, 179.396)
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
        plt.savefig('{id:03d}.jpg'.format(id=j))


    if __name__ == '__main__':
        result = run_simulation()
        plot(result)


    j = j + 5