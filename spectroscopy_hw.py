import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
     
def Gauss(x, a, x0, sigma): return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def Emission_Fit(xpixel, intensity, mask):
    mean = sum(xpixel[mask]*intensity[mask])/sum(intensity[mask])
    sigma = np.sqrt(sum(intensity[mask] * (xpixel[mask] - mean)**2) / sum(intensity[mask]))
    popt, pcov = curve_fit(Gauss, xpixel[mask], intensity[mask], p0=[max(intensity[mask]), mean, sigma])

    # Plot the data
    '''
    plt.plot(xpixel[mask], intensity[mask], 'ro')
    plt.plot(xpixel[mask], Gauss(xpixel[mask], *popt), 'k-')
    plt.show()
    '''
    return popt

# Import the data
data_dir = '/media/fmendez/E637-B25B/Fall_2020/Obstech/HW/Spectroscopy/Data/'

# Perform the dark subtraction for each lamp
ne = fits.getdata(data_dir+'neon_2s.FIT')
hg = fits.getdata(data_dir+'mercury_0d1s.FIT')
h = fits.getdata(data_dir+'hydrogen_1s.FIT')

# Plot the reduced image
norm = ImageNormalize(hg, interval=ZScaleInterval())
'''
plt.title('Neon Lamp')
plt.imshow(h, cmap='gray', origin='lower', norm=norm)
plt.show()
'''
# Mask the spectrum and create the 1D version
ne_mask = ne[350:450, ::]
hg_mask = hg[350:450, ::]
h_mask = h[350:450, ::]
'''
plt.subplot(3,1,1)
plt.title('Neon Spectrum', size=15)
plt.imshow(ne_mask, cmap='gray', origin='lower', norm=norm)

plt.subplot(3,1,2)
plt.title('Mercury Spectrum', size=15)
plt.imshow(hg_mask, cmap='gray', origin='lower', norm=norm)

plt.subplot(3,1,3)
plt.title('Hydrogen Spectrum', size=15)
plt.imshow(h_mask, cmap='gray', origin='lower', norm=norm)
plt.tight_layout()
plt.show()
'''
ne_counts = np.sum(ne_mask, axis=0)
ne_xpix = np.arange(ne_counts.size)

hg_counts = np.sum(hg_mask, axis=0)
hg_xpix = np.arange(hg_counts.size)

h_counts = np.sum(h_mask, axis=0)
h_xpix = np.arange(h_counts.size)

# Plot the 1D spectrum
'''
plt.subplot(3,1,1)
plt.title('Neon 1D Spectrum', size=15)
plt.ylabel('Intensity', size=15)
plt.plot(ne_xpix, ne_counts, 'k-')

plt.subplot(3,1,2)
plt.title('Mercury 1D Spectrum', size=15)
plt.ylabel('Intensity', size=15)
plt.plot(hg_xpix, hg_counts, 'k-')

plt.subplot(3,1,3)
plt.title('Hydrogen 1D Spectrum', size=15)
plt.ylabel('IntensityX', size=15)
plt.xlabel(' Pixels', size=15)
plt.plot(h_xpix, h_counts, 'k-')
plt.tight_layout()
plt.show()
'''
# Isolate the emission lines to perform the Gaussian fit
ne1 = (ne_xpix > 12.5) & (ne_xpix < 29.5)
ne2 = (ne_xpix > 28.5) & (ne_xpix < 44.5)
hg1 = (hg_xpix > 63.) & (hg_xpix < 73.)
hg2 = (hg_xpix > 73.) & (hg_xpix < 83.)
hg3 = (hg_xpix > 217.) & (hg_xpix < 234.)
hg4 = (hg_xpix > 745.) & (hg_xpix < 765.)
hg5 = (hg_xpix > 898.) & (hg_xpix < 914.)
h1 = (h_xpix > 507.) & (h_xpix < 521.)
h2 = (h_xpix > 758.) & (h_xpix < 770.)

# Perform the Gaussian fit and calculate the pix position of the peak
x_ne1 = Emission_Fit(ne_xpix, ne_counts, ne1)[1]
x_ne2 = Emission_Fit(ne_xpix, ne_counts, ne2)[1]
x_hg1 = Emission_Fit(hg_xpix, hg_counts, hg1)[1]
x_hg2 = Emission_Fit(hg_xpix, hg_counts, hg2)[1]
x_hg3 = Emission_Fit(hg_xpix, hg_counts, hg3)[1]
x_hg4 = Emission_Fit(hg_xpix, hg_counts, hg4)[1]
x_hg5 = Emission_Fit(hg_xpix, hg_counts, hg5)[1]
x_h1 = Emission_Fit(h_xpix, h_counts, h1)[1]
x_h2 = Emission_Fit(h_xpix, h_counts, h2)[1]

# Create the arrays for the xpix positions and corresponding wavelengths
x_pos = np.array([x_ne2, x_hg1, x_hg2, x_hg3, x_hg4, x_hg5, x_h1, x_h2])
waves = np.array([585.249, 579.065, 576.959, 546.074, 435.835, 404.656, 486.1, 434.04])

# Perform the linear fit to the data and plot it
line_fit = np.poly1d(np.polyfit(x_pos, waves, 1.))
'''
plt.title('X Pix Positions vs Wavelength Linear Fit', size=15)
plt.xlabel('X Pix Positions', size=15)
plt.ylabel('Wavelengths (nm)', size=15)
plt.plot(x_pos, waves, 'ro')
plt.plot(x_pos, line_fit(x_pos), 'k--')
plt.tight_layout()
plt.show()
'''
# Plot each lamp spectrum with the wavelength solution
'''
plt.subplot(3,1,1)
plt.title('Neon Spectrum (Wave Calibrated)', size=15)
plt.ylabel('Flux', size=15)
plt.plot(line_fit(ne_xpix), ne_counts, 'k-')

plt.subplot(3,1,2)
plt.title('Mercury Spectrum (Wave Calibrated)', size=15)
plt.ylabel('Flux', size=15)
plt.plot(line_fit(hg_xpix), hg_counts, 'k-')

plt.subplot(3,1,3)
plt.title('Hydrogen Spectrum (Wave Calibrated)', size=15)
plt.ylabel('Flux', size=15)
plt.xlabel('Wavelength (nm)', size=15)
plt.plot(line_fit(h_xpix), h_counts, 'k-')
plt.tight_layout()
plt.show()
'''
# Compute the measured wavelengths and measure the residuals and the RMS
measured_waves = np.array([585.448, 579.05, 576.92, 545.900, 435.97, 404.555, 486.076, 434.19])
residuals = (waves - measured_waves)/waves

rms_res = np.std(residuals)

print(min(line_fit(h_xpix)), max(line_fit(h_xpix)))
print('The measured centroids are', x_pos)
print('The Known wavelengths are', waves)
print('The measured wavelengths are', measured_waves)
print('The RMS of the Residuals is', rms_res)

plt.title('Wavelength Residuals', size=15)
plt.ylabel('Residuals', size=15)
plt.xlabel('X Pix Positions', size=15)
plt.plot(x_pos, residuals, 'ko')
plt.tight_layout()
plt.show()

# Resources for line list
# Neon: https://www.researchgate.net/figure/Neon-emission-spectra-solid-line-from-dust-free-dc-positive-column-thin-line-from_fig1_323359177
# Mercury: http://hyperphysics.phy-astr.gsu.edu/hbase/quantum/atspect2.html#:~:text=The%20prominent%20mercury%20lines%20are,weak%20line%20at%20491.604%20nm.
# Hydrogen: http://chemed.chem.purdue.edu/genchem/topicreview/bp/ch6/bohr.html










