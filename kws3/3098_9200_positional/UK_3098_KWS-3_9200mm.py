import bornagain as ba
from bornagain import deg, angstrom, nm, kvector_t, nm2
from matplotlib import pyplot as plt
import matplotlib
import  numpy as np

matplotlib.rcParams['image.cmap'] = 'jet'

i=0
while i < 61:

    def get_sample():
        # Defining Materials
        material_1 = ba.HomogeneousMaterial("Air", 0.0, 0.0)
        material_2 = ba.MaterialBySLD("Au", 4.6665e-06, -1.6205e-08)
        material_3 = ba.MaterialBySLD("Si", 2.0737e-06, -2.3758e-11)

        # Defining Layers

        layer_1 = ba.Layer(material_1)
        layer_2 = ba.Layer(material_3)

        formFactor_1 = ba.FormFactorCone6(159.0 * nm, 640.0 * nm, 86.0 * deg)
        formFactor_2 = ba.FormFactorCone6(157.0 * nm, 640.0 * nm, 86.0 * deg)
        formFactor_3 = ba.FormFactorTruncatedSphere(115.0 * nm, 160.0 * nm, 0.0 * nm)

        particle_1 = ba.Particle(material_2, formFactor_1)
        particle_2 = ba.Particle(material_3, formFactor_2)
        particle_3 = ba.Particle(material_2, formFactor_3)
        particle_3_position = kvector_t(0.0 * nm, 0.0 * nm, 640.0 * nm)
        particle_3.setPosition(particle_3_position)

        # Defining Core Shell Particles

        particleCoreShell_1 = ba.ParticleCoreShell(particle_2, particle_1)
        particleCoreShell_1_rotation = ba.RotationZ(0 * deg)
        particleCoreShell_1.setRotation(particleCoreShell_1_rotation)

        # Defining composition of particles at specific positions
        particleComposition_1 = ba.ParticleComposition()
        particleComposition_1.addParticle(particleCoreShell_1)
        particleComposition_1.addParticle(particle_3)
        particleComposition_1_rotation = ba.RotationX(0.0 * deg)
        particleComposition_1.setRotation(particleComposition_1_rotation)

        # Defining Interference Functions
        interference_1 = ba.InterferenceFunction2DLattice(1500.0 * nm, 1500.0 * nm, 120.0 * deg, i * deg)
        interference_1_pdf = ba.FTDecayFunction2DCauchy(1000.0 * nm, 1000.0 * nm, 0.0 * deg)
        interference_1.setDecayFunction(interference_1_pdf)
        interference_1.setPositionVariance(500.0 * nm2)

        # Defining Particle Layouts and adding Particles
        layout_1 = ba.ParticleLayout()
        layout_1.addParticle(particleComposition_1, 1.0)
        layout_1.setInterferenceFunction(interference_1)
        layout_1.setTotalParticleSurfaceDensity(0.000001)

        # Adding layouts to layers
        layer_1.addLayout(layout_1)

        # Defining Multilayers
        multiLayer_1 = ba.MultiLayer()
        multiLayer_1.addLayer(layer_1)
        multiLayer_1.addLayer(layer_2)
        return multiLayer_1

    def get_simulation():
        simulation = ba.GISASSimulation()
    
        detector = ba.RectangularDetector(256, 90.0, 256, 90.0)
        detector.setPerpendicularToDirectBeam(9200.0, 45.3, -140)
        simulation.setDetector(detector)
    
        simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(0.12, 0.12))
        simulation.setBeamParameters(1.28*nm, 0.56*deg, 0.0*deg)
        simulation.setBeamIntensity(1.0e+04)
        simulation.setTerminalProgressMonitor()
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
                         xlabel=r'$Q_{y} [1/nm]$', ylabel=r'$Q_{z} [1/nm]$', zmin=0.001, zmax=1e+5)
        plt.savefig('1_{id:03d}.png'.format(id=i))
        #plt.show()


    if __name__ == '__main__':
        result = run_simulation()
        plot(result)
        arr = result.array()
        np.savetxt("intensity1_{id:03d}.txt".format(id=i), arr)


    i=i+5