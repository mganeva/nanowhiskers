import bornagain as ba
from bornagain import deg, angstrom, nm, kvector_t
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

        # Defining Form Factors
        formFactor_1 = ba.FormFactorCone6(90.0*nm, 270.0*nm, 75.0*deg)
        formFactor_2 = ba.FormFactorCone6(88.0*nm, 270.0*nm, 75.0*deg)

        # Defining Particles
        particle_1 = ba.Particle(material_2, formFactor_1)
        particle_2 = ba.Particle(material_3, formFactor_2)

        # Defining Core Shell Particles

        particleCoreShell_1 = ba.ParticleCoreShell(particle_2, particle_1)
        particleCoreShell_1_rotation = ba.RotationZ(i*deg)
        particleCoreShell_1.setRotation(particleCoreShell_1_rotation)

        # Defining Particle Layouts and adding Particles
        layout_1 = ba.ParticleLayout()
        layout_1.addParticle(particleCoreShell_1, 1.0)
        layout_1.setTotalParticleSurfaceDensity(0.0001)

        # Defining Roughness Parameters
       # layerRoughness_1 = ba.LayerRoughness(1.0, 0.3, 5.0*nm)

        # Adding layouts to layers
        layer_1.addLayout(layout_1)

        # Defining Multilayers
        multiLayer_1 = ba.MultiLayer()
        multiLayer_1.addLayer(layer_1)
        multiLayer_1.addLayer(layer_2)
       # multiLayer_1.addLayerWithTopRoughness(layer_2, layerRoughness_1)
        return multiLayer_1


    def get_simulation():
        simulation = ba.GISASSimulation()
    
        detector = ba.RectangularDetector(256, 90.0, 256, 90.0)
        detector.setPerpendicularToDirectBeam(1300.0, 45.0, 29.0)
        simulation.setDetector(detector)
    
        #simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(0.34, 0.34))
        simulation.setBeamParameters(1.28*nm, 0.6*deg, 0.0*deg)
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
                         xlabel=r'$Q_{y} [1/nm]$', ylabel=r'$Q_{z} [1/nm]$', zlabel=None)
        plt.savefig('1_{id:03d}.png'.format(id=i))
        #plt.show()


    if __name__ == '__main__':
        result = run_simulation()
        plot(result)
        arr = result.array()
        np.savetxt("intensity_{id:03d}.txt".format(id=i), arr)


    i=i+5