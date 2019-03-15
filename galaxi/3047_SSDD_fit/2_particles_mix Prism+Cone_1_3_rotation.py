import numpy
import bornagain as ba
from bornagain import deg, angstrom, nm, kvector_t
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['image.cmap'] = 'jet'

nslices_1 = 5
nslices_2 = 15

j=0
while j<121:
    def get_sample():
        # Defining Materials
        material_1 = ba.HomogeneousMaterial("example01_Air", 0.0, 0.0)
        material_2 = ba.HomogeneousMaterial("Si", 5.73327e-06, 1.006366e-07)

        # Defining Layers
        layer_1 = ba.Layer(material_1)
        layer_2 = ba.Layer(material_2)

        particleComposition_1 = ba.ParticleComposition()

        for i in range(nslices_1):
            r = 159
            z = i * 15 * nm
            y = z + 15 * nm

            # Defining Form Factors
            formFactor_1 = ba.FormFactorCone6(r, 5.0 * nm, 68.0 * deg)
            formFactor_2 = ba.FormFactorCone6(r, 10.0 * nm, 78.0 * deg)

            # Defining Particles
            particle_1 = ba.Particle(material_2, formFactor_1)
            particle_1_rotation = ba.RotationY(180.0 * deg)
            particle_1.setRotation(particle_1_rotation)
            particle_1_position = kvector_t(0.0 * nm, 0.0 * nm, y * nm)
            particle_1.setPosition(particle_1_position)
            particle_2 = ba.Particle(material_2, formFactor_2)
            particle_2_position = kvector_t(0.0 * nm, 0.0 * nm, z * nm)
            particle_2.setPosition(particle_2_position)

            # Defining composition of particles at specific positions

            particleComposition_1.addParticle(particle_1)
            particleComposition_1.addParticle(particle_2)

        for i in range(nslices_2):
            r = 159 * nm - i * 2
            z = i * 15 * nm
            z2 = z + nslices_1 * 15*nm
            y = z + 15 * nm
            y2 = y + + nslices_1 * 15*nm

            # Defining Form Factors
            formFactor_1 = ba.FormFactorCone6(r, 5.0 * nm, 68.0 * deg)
            formFactor_2 = ba.FormFactorCone6(r, 10.0 * nm, 78.0 * deg)

            # Defining Particles
            particle_1 = ba.Particle(material_2, formFactor_1)
            particle_1_rotation = ba.RotationY(180.0 * deg)
            particle_1.setRotation(particle_1_rotation)
            particle_1_position = kvector_t(0.0 * nm, 0.0 * nm, y2 * nm)
            particle_1.setPosition(particle_1_position)
            particle_2 = ba.Particle(material_2, formFactor_2)
            particle_2_position = kvector_t(0.0 * nm, 0.0 * nm, z2 * nm)
            particle_2.setPosition(particle_2_position)

            # Defining composition of particles at specific positions

            particleComposition_1.addParticle(particle_1)
            particleComposition_1.addParticle(particle_2)

        particleComposition_1_rotation = ba.RotationZ(j * deg)
        particleComposition_1.setRotation(particleComposition_1_rotation)

        # Defining Particle Layouts and adding Particles
        layout_1 = ba.ParticleLayout()
        layout_1.addParticle(particleComposition_1, 1.0)
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

        simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(0.43, 0.43))
        simulation.setBeamParameters(0.134 * nm, 0.2 * deg, 0.0 * deg)
        simulation.setBeamIntensity(5.0e+06)
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



