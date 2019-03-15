import matplotlib
from matplotlib import pyplot as plt
import bornagain as ba

#matplotlib.rcParams['image.cmap'] = 'jet'

wavelength = 1.3414*ba.angstrom
alpha_i = 0.18*ba.degree

# detector setup as given from instrument responsible
pilatus_npx, pilatus_npy = 981, 1043
pilatus_pixel_size = 0.172  # in mm
detector_distance = 830.0  # in mm
beam_xpos, beam_ypos = 616.9, 277.3  # in pixels

crop_xmin, crop_ymin, crop_xmax, crop_ymax = -5.8, 0, 3.8, 7.3  # crop coordinates, nm^-1


def create_detector():
    """
    Creates and returns GALAXY detector
    """
    u0 = beam_xpos*pilatus_pixel_size  # in mm
    v0 = beam_ypos*pilatus_pixel_size  # in mm
    detector = ba.RectangularDetector(pilatus_npx, pilatus_npx*pilatus_pixel_size, pilatus_npy, pilatus_npy*pilatus_pixel_size)
    detector.setPerpendicularToDirectBeam(detector_distance, u0, v0)

    return detector

def create_simulation():
    """
    Creates and returns GISAS simulation with beam and detector defined
    """
    simulation = ba.GISASSimulation()
    simulation.setBeamParameters(wavelength, alpha_i, 0.0)
    simulation.setDetector(create_detector())
    return simulation



def load_real_data(hist, filename):
        """
        Fill histogram representing our detector with intensity data from tif file.
        Returns cropped version of it, which represent the area we are interested in.
        """
        hist.load(filename)
        return hist.crop(crop_xmin, crop_ymin, crop_xmax, crop_ymax)           # if you don't want to crop


def plot_data(filename, k):
        """
        Load data and plot results
        """
        simulation = create_simulation()
        hist = simulation.result().histogram2d(ba.AxesUnits.QSPACE)
        result = load_real_data(hist, filename)
        plt.figure(figsize=(12.80, 10.24))
        #plt.subplot()
        # showing the result
        im = plt.imshow(result.array()+1, norm=matplotlib.colors.LogNorm(1.0, result.getMaximum()), cmap='jet', aspect='auto')
        cb = plt.colorbar(im)
        cb.set_label(r'Intensity (arb. u.)', size=16)
        plt.xlabel(r'$Q_y (nm^{-1})$', fontsize=16)
        plt.ylabel(r'$Q_z (nm^{-1})$', fontsize=16)
        #plt.show()
        plt.savefig('{id:04d}.jpg'.format(id=k))
        print("{id:04d}.png".format(id=k))


if __name__ == '__main__':
    for k in range(2949, 3022):
        filename1 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Real_Date_3047_SSDD/Koneva_3047_SSDD_3{id}.tif".format(id=k)
        plot_data(filename=filename1, k=k)

