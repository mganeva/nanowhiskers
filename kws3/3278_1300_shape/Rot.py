import numpy as np
import matplotlib.pyplot as plt

b = np.loadtxt('00065431_0001_p15749_S3278-GiSANS-0.6dgr-rot-scan-m15-15dgr_HRD_standard.det', dtype=int)
c = np.rot90(b, k=3)
np.savetxt('00065952_3098.txt', c, fmt='%d')
