import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.io import fits
from scipy import stats
from astropy.modeling import models, fitting
import os

######################################################### READ NOISE MEASUREMENTS ############################################################################
print('These are the measurements for Read Noise:')

# Create the path to data and file list
data_dir = '/media/fmendez/E637-B25B/Fall_2020/Obstech/CCD_Lab/Data/Group3/'
file_list = os.listdir(data_dir)

# Create the master arrays for the exposures
exp_0p1s = []
exp_1s = []
exp_10s = []

for file in file_list:
    # Extract the exposure times from the header info and the image data
    exp_time = fits.getheader(data_dir+file)['EXPTIME']
    image = fits.getdata(data_dir+file).astype(float)

    # Append the images to the corresponding master array
    if exp_time == 0.1:
        exp_0p1s.append(image)

    if exp_time == 1.0:
        exp_1s.append(image)

    if exp_time == 10.:
        exp_10s.append(image)

exp_0p1s = np.array(exp_0p1s)
exp_1s = np.array(exp_1s)
exp_10s = np.array(exp_10s)

# Create the STD arrays
std_array_0p1s = np.std(exp_0p1s, axis=0)
std_array_1s = np.std(exp_1s, axis=0)
std_array_10s = np.std(exp_10s, axis=0)

# Compute the mean, median, and uncertainty
mean_0p1s, median_0p1s, uncertainty_0p1s = np.mean(std_array_0p1s), np.median(std_array_0p1s), np.std(std_array_0p1s)/np.sqrt(len(std_array_0p1s))
mean_1s, median_1s, uncertainty_1s = np.mean(std_array_1s), np.median(std_array_1s), np.std(std_array_1s)/np.sqrt(len(std_array_1s))
mean_10s, median_10s, uncertainty_10s = np.mean(std_array_10s), np.median(std_array_10s), np.std(std_array_10s)/np.sqrt(len(std_array_10s))

print('Measurements From Statistics:','\nMean 0.1s:', round(mean_0p1s, 4), 'Median 0.1s:', round(median_0p1s, 4), 'Uncertainty 0.1s:', round(uncertainty_0p1s, 2))
print('Mean 1s:', round(mean_1s, 4), 'Median 1s:', round(median_1s, 4), 'Uncertainty 1s:', round(uncertainty_1s, 2))
print('Mean 10s:', round(mean_10s, 4), 'Median 10s:', round(median_10s, 4), 'Uncertainty 10s:', round(uncertainty_10s, 2))

# Create the different histograms of the RMS arrays
heights_0p1s, borders_0p1s = np.histogram(np.ravel(std_array_0p1s), bins=10000)
widths_0p1s = np.diff(borders_0p1s)
centers_0p1s = (borders_0p1s[:-1] + widths_0p1s)

heights_1s, borders_1s = np.histogram(np.ravel(std_array_1s), bins=25000)
widths_1s = np.diff(borders_1s)
centers_1s = (borders_1s[:-1] + widths_1s)

heights_10s, borders_10s = np.histogram(np.ravel(std_array_10s), bins=2000)
widths_10s = np.diff(borders_10s)
centers_10s = (borders_10s[:-1] + widths_10s)

# Perform the Gaussian Fit
t_init_0p1s = models.Gaussian1D()
fit_t_0p1s = fitting.LevMarLSQFitter()
t_0p1s = fit_t_0p1s(t_init_0p1s, centers_0p1s, heights_0p1s)

t_init_1s = models.Gaussian1D()
fit_t_1s = fitting.LevMarLSQFitter()
t_1s = fit_t_1s(t_init_1s, centers_1s, heights_1s)

t_init_10s = models.Gaussian1D()
fit_t_10s = fitting.LevMarLSQFitter()
t_10s = fit_t_10s(t_init_10s, centers_10s, heights_10s)

x_fit_0p1s = np.linspace(borders_0p1s[0], borders_0p1s[-1], 10000)
x_fit_1s = np.linspace(borders_1s[0], borders_1s[-1], 10000)
x_fit_10s = np.linspace(borders_10s[0], borders_10s[-1], 10000)

# Plot the pixel distribution and the Gaussian Fit
'''
plt.subplot(1,3,1)
plt.title(r'$t=0.1 s$', size=15)
plt.ylabel('Frequency', size=15)
plt.xlim(0, 60)
plt.bar(centers_0p1s, heights_0p1s, width=widths_0p1s, color='red', label='Pixel Distribution')
plt.plot(x_fit_0p1s, t_0p1s(x_fit_0p1s), color='black', label='Gaussian Fit')

plt.subplot(1,3,2)
plt.title(r'$t=1 s$', size=15)
plt.xlabel('ADU/pixel', size=15)
plt.xlim(0, 60)
plt.bar(centers_1s, heights_1s, width=widths_1s, color='red', label='Pixel Distribution')
plt.plot(x_fit_1s, t_1s(x_fit_1s), color='black', label='Gaussian Fit')

plt.subplot(1,3,3)
plt.title(r'$t=10 s$',size=15)
plt.xlim(0, 60)
plt.bar(centers_10s, heights_10s, width=widths_10s, color='red', label='Pixel Distribution')
plt.plot(x_fit_10s, t_10s(x_fit_10s), color='black', label='Gaussian Fit')

plt.tight_layout()
plt.legend()
plt.show()
'''
# Measure the mean and sigma values from the Gaussian fit
gauss_mean_0p1s, gauss_uncertainty_0p1s = round(t_0p1s.mean[0], 4), round(t_0p1s.stddev[0]/np.sqrt(len(centers_0p1s)), 2)
gauss_mean_1s, gauss_uncertainty_1s = round(t_1s.mean[0], 4), round(t_1s.stddev[0]/np.sqrt(len(centers_1s)), 2)
gauss_mean_10s, gauss_uncertainty_10s = round(t_10s.mean[0], 4), round(t_10s.stddev[0]/np.sqrt(len(centers_10s)), 2)

print('\nMeasurements From Gaussian Fit:', '\nMean 0.1s:', gauss_mean_0p1s, 'Uncertainty 0.1s:', gauss_uncertainty_0p1s)
print('Mean 1s:', gauss_mean_1s, 'Uncertainty 1s:', gauss_uncertainty_1s)
print('Mean 10s:', gauss_mean_10s, 'Uncertainty 10s:', gauss_uncertainty_10s)

# Plot the read noise values
'''
plt.title('Read Noise as a Function of Time', size=15)
plt.ylabel('Read Noise', size=15)
plt.xlabel('Exposure Time (secs)', size=15)
plt.errorbar([0.1,1.,10.], [mean_0p1s, mean_1s, mean_10s], yerr=[uncertainty_0p1s, uncertainty_1s, uncertainty_10s], \
             fmt='ro', capsize=3, label='Statistical Measurements')
plt.errorbar([0.1,1.,10.], [gauss_mean_0p1s, gauss_mean_1s, gauss_mean_10s], \
             yerr=[gauss_uncertainty_0p1s, gauss_uncertainty_1s, gauss_uncertainty_10s], fmt='ko', capsize=3, label='Gaussian Fit Measurements')
plt.legend()
plt.tight_layout()
plt.show()
'''
##################################################################### DARK CURRENT MEASUREMENTS #############################################################
print('\nThese are the Measurements for Dark Current')

# Import the files fpr the bias level and crate the master Bias
biases = []

for file in file_list:
    # Extract the exposure times from the header info and the image data
    exp_time = fits.getheader(data_dir+file)['EXPTIME']
    image = fits.getdata(data_dir+file).astype(float)

    # Append the images to the corresponding master array
    if exp_time == 0.001:
        biases.append(image)

biases = np.asarray(biases)
master_bias = np.median(biases, axis=0)

# Calculate the Dark Current
times = []
means = []
uncertainties = []

for file in file_list:
    # Extract the exposure times from the header info and the image data
    exp_time = fits.getheader(data_dir+file)['EXPTIME']
    image = fits.getdata(data_dir+file).astype(float)

    # Ignore the exposure times < 0.2 and = 1:
    if exp_time >= 0.2 and exp_time != 1. and exp_time != 10.:
        # Subtract the master bias
        debiased = image - master_bias

        # Create the Gaussian fit to the ADU/pixel distribution
        heights, borders = np.histogram(np.ravel(debiased), bins=10000)
        widths = np.diff(borders)
        centers = (borders[:-1] + widths)

        t_init = models.Gaussian1D()
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, centers, heights)

        x_fit = np.linspace(borders[0], borders[-1], 10000)

        # Measure the mean value and uncertainty and append it to the lists
        mean = round(t.mean[0], 4)
        uncertainty = round(t.stddev[0]/np.sqrt(len(centers)), 2)

        times.append(exp_time)
        means.append(mean)
        uncertainties.append(uncertainty)

        print('Exposure Time:', exp_time, 's', 'Mean:', mean, 'Uncertainty:', uncertainty)

        # Plot the pixel distribution and the gaussian fit
        '''
        plt.title(r'$t=$'+str(exp_time)+r'$s$', size=15)
        plt.xlabel('ADU/pixels', size=15)
        plt.ylabel('Frequency', size=15)
        plt.xlim(-300, 300)
        plt.bar(centers, heights, width=widths, color='red', label='Pixel Distribution')
        plt.plot(x_fit, t(x_fit), color='black', label='Gaussian Fit')
        plt.legend()
        plt.tight_layout()
        plt.show()
        '''

# Create the linear fit of the mean values and the exposure times and plot it
polyfit_curve = np.poly1d(np.polyfit(times, means, 1.))
x_fit = np.arange(min(times), max(times))
'''
plt.title(r'Measurements of $dDN/dt$', size=15)
plt.xlabel(r'$t$ (s)', size=15)
plt.ylabel('Mean ADU/pixel', size=15)
plt.plot(x_fit, polyfit_curve(x_fit),  'k--')  
plt.errorbar(times, means, yerr=uncertainties, fmt='ro', capsize=3)
plt.tight_layout()
plt.show()
'''
# Measure dDN/dt and its uncertainty
dDN = np.polyfit(times, means, 1., cov=True)[0][0]
sigma_DN = np.polyfit(times, means, 1., cov=True)[1][0][0]

print('Dark Current = ',round(dDN, 4), '+/-', round(sigma_DN, 6))

################################################################# GAIN FACTOR MEASUREMENTS ##################################################################
print('\nThese are the measurements for the Gain Factor:')

# Import the data for the gain
data_dir = '/media/fmendez/E637-B25B/Fall_2020/Obstech/CCD_Lab/Data/Group_3/'
file_list = [g_file for g_file in os.listdir(data_dir) if g_file.startswith('gain') and g_file.endswith('.fits')]

# Etract the different exposure times from the flats
times = []

for file in file_list:
    # Extract exposure times
    exp_time = fits.getheader(data_dir+file)['EXPTIME']

    if exp_time not in times:
        times.append(exp_time)

# Measure the signal and the variance
signals = []
variances = []

for time in np.sort(np.array(times)):
    master = []
    
    for file in file_list:
        # Extract the exposure times (again) and the image data
        exp_time = fits.getheader(data_dir+file)['EXPTIME']
        image = fits.getdata(data_dir+file).astype(float)

        if exp_time == time:
            # Create the master image
            master.append(image)

    # Measure the variance and the signal and append it to the list
    rms_array = np.std(np.asarray(master), axis=0)
    mean_array = np.mean(np.asarray(master), axis=0)

    signal = np.mean(mean_array)
    variance = np.mean(rms_array**2.)

    print('Exposure Time:', time, 'Signal:', round(signal, 4), 'Variance:', round(variance, 4))

    signals.append(signal)
    variances.append(variance)

    # Clear the master list for the next exposure time
    master.clear()

# create the linear fot for the signal and the variance and plot it
signals = np.asarray(signals)
variances = np.asarray(variances)

polyfit_curve = np.poly1d(np.polyfit(signals, variances, 1.))
x_fit = np.arange(min(signals), max(signals))
'''
plt.title('Measurements of G', size=15)
plt.xlabel('Signal', size=15)
plt.ylabel('Variance', size=15)
plt.plot(x_fit, polyfit_curve(x_fit),  'k--')  
plt.plot(signals, variances, 'ro')
plt.tight_layout()
plt.show()
'''         
# Measure 1/G and its uncertainty
gain = 1./np.polyfit(signals, variances, 1., cov=True)[0][0]

slp, itc, rv, pv, stderr = stats.linregress(signals, variances)
gain_uncertainty = stderr/(np.sqrt(sum((signals-np.mean(signals))**2.)))

print('Gain:', round(gain, 4), '+/-', round(gain_uncertainty, 6))
        
###############################################################################3 LINEARITY ##################################################################
print('\nThese are the measurements of linearty:')

# Import the files for bright and fain liniarity
bright_file_list = [b_lin_file for b_lin_file in os.listdir(data_dir) if b_lin_file.startswith('lin_bright') and b_lin_file.endswith('.fits')]
faint_file_list = [f_lin_file for f_lin_file in os.listdir(data_dir) if f_lin_file.startswith('lin_faint') and f_lin_file.endswith('.fits')]

# Measure the signal as a function of exposure time for the bright exposures
print('Bright Linearity:')

bright_signals = []
bright_times = []

for file in bright_file_list:
    # Extract the exposure times and the image data
    exp_time = fits.getheader(data_dir+file)['EXPTIME']
    image = fits.getdata(data_dir+file).astype(float)

    bright_times.append(exp_time)

    # Measure the signal from the image
    signal = np.mean(image)
    bright_signals.append(signal)

    print('Exposure Time:', exp_time, 'Signal:', round(signal, 4))

# Measure the signal as a function of exposure time for the faint exposures
print('Faint Linearity:')

faint_signals = []
faint_times = []

for file in faint_file_list:
    # Extract the exposure times and the image data
    exp_time = fits.getheader(data_dir+file)['EXPTIME']
    image = fits.getdata(data_dir+file).astype(float)

    faint_times.append(exp_time)

    # Measure the signal from the image
    signal = np.mean(image)
    faint_signals.append(signal)

    print('Exposure Time:', exp_time, 'Signal:', round(signal, 4))
    
# Fit a line to the linearities, measure the residuals,  and plot them
bright_polyfit_curve = np.poly1d(np.polyfit(bright_times, bright_signals, 1.))
bright_x_fit = np.arange(min(bright_times), max(bright_times))

faint_polyfit_curve = np.poly1d(np.polyfit(faint_times, faint_signals, 1.))
faint_x_fit = np.arange(min(faint_times), max(faint_times))

bright_residuals = (np.array(bright_signals) - bright_polyfit_curve(bright_times))/np.array(bright_signals)
faint_residuals = (np.array(faint_signals) - faint_polyfit_curve(faint_times))/np.array(faint_signals)
'''
plt.title('Linearity Measurements', size=15)
plt.xlabel(r'$t$ (s)', size=15)
plt.ylabel('Signal', size=15)
plt.plot(faint_times, faint_signals, 'bo')
plt.plot(bright_times, bright_signals, 'ro')
plt.plot(bright_x_fit, bright_polyfit_curve(bright_x_fit),  'r--', label='Bright Linearity')  
plt.plot(faint_x_fit, faint_polyfit_curve(faint_x_fit),  'b--', label='Faint Linearity')
plt.legend()
plt.tight_layout()
plt.show()
'''
plt.subplot(2,1,1)
plt.title('Bright Linearity Residuals', size=15)
plt.xlabel('Signal', size=15)
plt.ylabel('Residuals', size=15)
plt.plot(bright_signals, bright_residuals, 'ro', label='Bright Residuals')

plt.subplot(2,1,2)
plt.title('Faint Linearity Residuals', size=15)
plt.xlabel('Signal', size=15)
plt.ylabel('Residuals', size=15)
plt.plot(faint_signals, faint_residuals, 'bo', label='Faint Residuals')
plt.tight_layout()
plt.show()





    





