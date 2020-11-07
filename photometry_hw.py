from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats, mad_std
from astropy.visualization import ZScaleInterval, ImageNormalize
import matplotlib.pyplot as plt
import numpy as np
import os
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
from photutils import datasets
from photutils import DAOStarFinder, find_peaks
from skimage.feature import *
from scipy.ndimage import interpolation as interp
import os
import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore', np.RankWarning)

def Master_Dark_Undebiased(data_list):
    # Convert the data list to an array and median combine it
    master_dark = np.median(np.asarray(data_list), axis=0)

    return master_dark

def Master_Flat_Undarked(data_list):
    # Convert the data list to an array and median combine it
    master_flat = np.median(np.asarray(data_list), axis=0)

    # Normalize the master flat frame
    norm_master_flat = master_flat/np.median(master_flat)

    return norm_master_flat

def Dark_Subtraction(file, calibration_list, et):
    
    for cal_file in calibration_list:
        # Extract the frame information
        mfr = fits.getheader(master_files_dir+cal_file)['FRAME']

        if mfr == 'Dark':
            md_et = fits.getheader(master_files_dir+cal_file)['EXPTIME']

            if et == md_et:
                print('Subtracting', cal_file, 'from', et)
                dark_subtracted = fits.getdata(data_dir+file) - fits.getdata(master_files_dir+cal_file)

                return dark_subtracted

def Flat_Division(file, calibration_list, fl):
    
    for cal_file in calibration_list:
        # Extract the frame information
        mfr = fits.getheader(master_files_dir+cal_file)['FRAME']

        if mfr == 'Flat':
            mf_fl = fits.getheader(master_files_dir+cal_file)['FILTER']

            if fl == mf_fl:
                print('Dividing', cal_file, 'from dark subtracted', fl, '\n')
                flat_divided = file/fits.getdata(master_files_dir+cal_file)

                return flat_divided

def Align_Image(images):
    # select the template image and align the target images
    template_image = images[0]
    aligned_images = []

    for i in range(len(images)):

        if i != 0:
            print('Aligning image', i+1, 'out of', len(images))
            # Select the target images
            target_image = images[i]

            # Calculate the shift of the images
            result, error, diffphase = register_translation(template_image, target_image, 5000)

            # Align the images based on the shift
            align_image = interp.shift(target_image, result)

            aligned_images.append(align_image)

    # Add the template image back to the aligned images list
    aligned_images.append(template_image)

    return aligned_images

def Sources_Positions(x_peaks, y_peaks, values):
    p = np.zeros((len(x_peaks), 3))

    for i in range(len(x_peaks)):
        p[i,0] = x_peaks[i]
        p[i,1] = y_peaks[i]
        p[i,2] = values[i]

    return p

def Peak_Positions(data, upper_lim, lower_lim):
    mean, median, std = sigma_clipped_stats(data, sigma=3.)
    treshold = median + (3.*std)
    peaks_table = find_peaks(data, treshold, box_size=11.)
    
    positions = Sources_Positions(peaks_table['x_peak'], peaks_table['y_peak'], peaks_table['peak_value'])

    x = []
    y = []
    val = []

    for i in range(len(positions[:,0])):

        if positions[:,2][i] > lower_lim and positions[:,2][i] < upper_lim:
            x.append(positions[:,0][i])
            y.append(positions[:,1][i])
            val.append(positions[:,2][i])
          
    return np.array(x), np.array(y), np.array(val)

def Gauss_Fit(data):
    # Convert the image to 1D
    data_1d = np.sum(data, axis=0)
    xpix = np.arange(data_1d.size)
        
    # Perform the Gaussian Fit
    t_init = models.Gaussian1D()
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, xpix, data_1d)

    x_fit = np.linspace(xpix[0], xpix[-1], 1000)

    # Plot the Gauss fit
    '''
    plt.title('Gaussian Fit')
    plt.plot(xpix, data_1d, 'ro', label='Star')
    plt.plot(x_fit, t(x_fit), color='black', label='Gaussian Fit')
    plt.legend()
    plt.show()
    '''
    return t

def Photometry(data, star_an, sky_an, et):
    # Measure the brightness of the star and the sky
    star_phot = aperture_photometry(data, star_an)
    sky_phot = aperture_photometry(data, sky_an)

    # Calculate the mean sky brightness and subtract it from the star
    sky_mean = sky_phot['aperture_sum']/sky_an.area
    sky_sum = sky_mean*star_an.area
    star_tot_brightness = star_phot['aperture_sum'] - sky_sum

    # Measure the flux and its error (Use the gain form the CCD Lab)
    gain = 0.0778
    flux = (star_tot_brightness*gain)/et
    flux_err = np.sqrt(flux)

    # Measure the Instrumental magnitude and its error
    inst_mag = -2.5*np.log10(flux)
    inst_mag_err = (2.5/np.log(10))*(flux_err/flux)

    return star_tot_brightness, flux, flux_err, inst_mag, inst_mag_err


date = '2020-10-17'
target = 'H PERSEI'
main_dir = '/media/fmendez/E637-B25B/Fall_2020/Obstech/HW/Photometry/'+date+'/'
data_dir = main_dir+'data/'

file_list = os.listdir(data_dir)
'''
# Extract the information from the header files
objects = []
exp_times = []
filters = []

for file in file_list:
    obj = fits.getheader(data_dir+file)['OBJECT']
    et = fits.getheader(data_dir+file)['EXPTIME']
    filt = fits.getheader(data_dir+file)['FILTER']
    fr = fits.getheader(data_dir+file)['FRAME']

    if obj not in objects and et not in exp_times and filt not in filters:
        objects.append(obj)
        exp_times.append(et)
        filters.append(filt)

    else:
        pass
    
    if obj == target and fr == 'Light':
        data = fits.getdata(data_dir+file)
        norm = ImageNormalize(data, interval=ZScaleInterval())
        
        plt.title(obj+' '+filt+' '+str(et)+' '+fr+' '+file)
        plt.imshow(data, cmap='gray', origin='lower', norm=norm)
        plt.show()

        keep = input('Keep File? ')
        
        if keep == 'n':
            print('Deleting', file)
            os.remove(data_dir+file)

        else:
            print('Keeping', file)
            pass
    
    
# Create the master DARK Frames frame
print('MASTER DARK FRAME CREATION')

for time in exp_times:
    print('Creating Master dark for', time, 'secs:')
    master_dark_array = []

    for file in file_list:
        et = fits.getheader(data_dir+file)['EXPTIME']
        fr = fits.getheader(data_dir+file)['FRAME']

        if time == et and fr == 'Dark':
            # Extract the image data and append it to the list
            darks_data = fits.getdata(data_dir+file)
            master_dark_array.append(darks_data)
            
    # Create the master Dark frame
    print(len(master_dark_array), 'frames found')
    master_dark = Master_Dark_Undebiased(master_dark_array)

    # Create and save the master dark frame
    save_dir = main_dir+'/reduced_files/'
    
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    
    hdu = fits.PrimaryHDU(master_dark)
    hdu.header['EXPTIME'] = time
    hdu.header['FRAME'] = 'Dark'

    hdul = fits.HDUList([hdu])

    if os.path.isfile(save_dir+'master_dark_'+str(time)+'s.fits'):
        print('Master Dark already saved')
        pass

    else:
        hdul.writeto(save_dir+'master_dark_'+str(time)+'s.fits')
        print('Master Dark Saved!')

    master_dark_array.clear()

# Create Master FLAT Frame
print('\nMASTER FLAT FRAME CREATION')

for filt in filters:
    print('Creating Master flat for filter', filt, ':')
    master_flat_array = []

    for file in file_list:
        f = fits.getheader(data_dir+file)['FILTER']
        fr = fits.getheader(data_dir+file)['FRAME']

        if filt == f and fr == 'Flat':
            # Extract the image data and append it to the list
            flat_data = fits.getdata(data_dir+file)
            master_flat_array.append(flat_data)
            
    # Create the master Flat frame
    print(len(master_flat_array), 'frames found')
    master_flat = Master_Flat_Undarked(master_flat_array)

    # Create and save the master flat frame
    save_dir = main_dir+'/reduced_files/'
    
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    
    hdu = fits.PrimaryHDU(master_flat)
    hdu.header['FILTER'] = filt
    hdu.header['FRAME'] = 'Flat'

    hdul = fits.HDUList([hdu])

    if os.path.isfile(save_dir+'master_flat_'+str(filt)+'.fits'):
        print('Master Flat already saved')
        pass

    else:
        hdul.writeto(save_dir+'master_flat_'+str(filt)+'.fits')
        print('Master Flat Saved!')

    master_flat_array.clear()

# Perform the data reduction                          
print('\nDATA REDUCTION')

master_files_dir = main_dir+'reduced_files/'
reduced_files = []
reduced_data = []

for file in file_list:
    # Extract the information
    obj = fits.getheader(data_dir+file)['OBJECT']
    et = fits.getheader(data_dir+file)['EXPTIME']
    filt = fits.getheader(data_dir+file)['FILTER']
    fr = fits.getheader(data_dir+file)['FRAME']

    if obj == target and fr == 'Light':
        # Perform the Dark subtraction
        mfile_list = os.listdir(master_files_dir)
        d_sub = Dark_Subtraction(file, mfile_list, et)
        
        # Perform Flat Division
        f_div = Flat_Division(d_sub, mfile_list, filt)

        # Plot the reduced file
        
        norm = ImageNormalize(f_div, interval=ZScaleInterval())
        
        plt.title(obj+' '+filt+' '+str(et)+' '+fr+' '+file)
        plt.imshow(f_div, cmap='gray', origin='lower', norm=norm)
        plt.show()
        
        reduced_files.append(file)
        reduced_data.append(f_div)

print('Data Reduction Complete! Total reduced files:', len(reduced_files))

# Perform Reduced Files Alignment
print('\nIMAGE ALIGNMENT')

aligned_data = Align_Image(reduced_data)

# Save the aligned frames
save_aligned_dir = main_dir+'aligned_files/'

if not os.path.isdir(save_aligned_dir):
    os.mkdir(save_aligned_dir)
    
for i in range(len(reduced_files)):
    # Extract the information
    obj = fits.getheader(data_dir+reduced_files[i])['OBJECT']
    et = fits.getheader(data_dir+reduced_files[i])['EXPTIME']
    filt = fits.getheader(data_dir+reduced_files[i])['FILTER']
    fr = fits.getheader(data_dir+reduced_files[i])['FRAME']

    data = aligned_data[i]

    # Create fits files
    hdu = fits.PrimaryHDU(data)
    
    hdu.header['OBJECT'] = obj
    hdu.header['EXPTIME'] = et
    hdu.header['FILTER'] = filt
    hdu.header['FRAME'] = fr

    hdul = fits.HDUList([hdu])

    if os.path.isfile(save_aligned_dir+'aligned_'+reduced_files[i]):
        print('Aligned file already saved')
        pass

    else:
        hdul.writeto(save_aligned_dir+'aligned_'+reduced_files[i])
        print('Aligned file Saved!')

# Create the Master science frames
print('\nMEDIAN COMBINE THE SCIENCE FRAMES')
save_aligned_list = os.listdir(save_aligned_dir)

for filt in filters:
    print('Creating Master science for filter', filt, ':')
    master_science_array = []

    for file in save_aligned_list:
        obj = fits.getheader(save_aligned_dir+file)['OBJECT']
        et = fits.getheader(save_aligned_dir+file)['EXPTIME']
        f = fits.getheader(save_aligned_dir+file)['FILTER']
        fr = fits.getheader(save_aligned_dir+file)['FRAME']

        if filt == f:
            # Extract the image data and append it to the list
            science_data = fits.getdata(save_aligned_dir+file)
            master_science_array.append(science_data)

        else:
            pass

    # Create the master Flat frame and plot it
    if len(master_science_array) > 0:
        print(len(master_science_array), 'frames found for filter', filt)
        master_science = np.median(np.asarray(master_science_array), axis=0)

        norm = ImageNormalize(master_science, interval=ZScaleInterval())
            
        plt.title(filt+' '+file)
        plt.imshow(master_science, cmap='gray', origin='lower', norm=norm)
        plt.show()

        # Create and save the master flat frame
        save_dir = main_dir+'/reduced_files/'
        
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)
        
        hdu = fits.PrimaryHDU(master_science)
        hdu.header['OBJECT'] = obj
        hdu.header['EXPTIME'] = et
        hdu.header['FILTER'] = filt
        hdu.header['FRAME'] = fr

        hdul = fits.HDUList([hdu])

        if os.path.isfile(save_dir+'master_science_'+str(filt)+'.fits'):
            print('Master science already saved')
            pass

        else:
            hdul.writeto(save_dir+'master_science_'+str(filt)+'.fits')
            print('Master science Saved!')

        master_science_array.clear()

    else:
        pass
'''
# Measure the FWHM, PSF photometry, and appreture photometry
data_dir = main_dir+'/reduced_files/'
data_files = os.listdir(data_dir)

for file in data_files:
    # Analyze only the science frames
    fr = fits.getheader(data_dir+file)['FRAME']

    # Extract the information of the header
    if fr == 'Light':
        obj = fits.getheader(data_dir+file)['OBJECT']
        et = fits.getheader(data_dir+file)['EXPTIME']
        f = fits.getheader(data_dir+file)['FILTER']

        # Extract the image data and find the ADU peaks
        data = fits.getdata(data_dir+file).astype(float)

        upper_lim = 70000.
        lower_lim = 20000.
        
        x, y, val = Peak_Positions(data, upper_lim, lower_lim)
        print('\n', obj, f, '- Peaks found:', len(x) )
        
        for i in range(len(x)):
            # Decide if the identified star should be measured or not
            norm = ImageNormalize(data, interval=ZScaleInterval())
            '''
            plt.title(obj+' '+f+' Filter')
            plt.imshow(data, cmap='gray', origin='lower', norm=norm)
            plt.plot(x[i], y[i], '+', c='darkblue', ms=20)
            plt.text(x[i]+20, y[i]-50, str(i+1), c='darkblue') 
            plt.show()
            '''
            measure = input('\nDo you wish to measure this star(y/n):')

            if measure == 'n':
                print('Peak', i+1, 'ignored')
                pass

            else:
                print('Measuring peak', i+1)

                # Crop the image to see the star only
                croped_data = data[int(y[i])-15:int(y[i])+15, int(x[i])-15:int(x[i])+15] 
                '''
                norm = ImageNormalize(croped_data, interval=ZScaleInterval())
                
                plt.title(obj+' '+f+' Peak Value: '+str(val[i]))
                plt.imshow(croped_data, cmap='gray', origin='lower', norm=norm)
                plt.show()
                '''
                # Measure the FWHM of the star
                fwhm = Gauss_Fit(croped_data).stddev[0]
                print('- The FWHM is', round(fwhm, 4), 'pixels')

                # Define the appertures for the annuli
                # Main Annulus
                star_an = CircularAperture((x[i], y[i]), r=1.3*fwhm)

                # Sky positions - Slect an area of the image where the are no stars
                sky_pos = (1400., 500.)
                sky_an = CircularAnnulus(sky_pos, r_in=2.*fwhm, r_out=3.*fwhm)
                
                plt.title(obj+' '+f)
                plt.imshow(data, cmap='gray', origin='lower', norm=norm)
                star_an.plot(color='yellow')
                sky_an.plot(color='orange')
                plt.show()
                
                # Measure star's total brightnes, flux, and instrumental magnitude
                star_tot_brightness, flux, flux_err, inst_mag, inst_mag_err = Photometry(data, star_an, sky_an, et)

                print('- The total star brightnes is', round(float(star_tot_brightness), 4,),
                      '\n- The star flux is', round(float(flux), 4), '+/-', round(float(flux_err), 1),
                      '\n- The star instrumental magnitude is', round(float(inst_mag), 4), '+/-', round(float(inst_mag_err), 6))

                    
                

                    
                    

                
              
                
        

        



    
        

        
            

    
    

    

    
    


        






















            
            

