import numpy as np
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize
import matplotlib.pyplot as plt
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
            result, error, diffphase = register_translation(template_image, target_image, 1000)

            # Align the images based on the shift
            align_image = interp.shift(target_image, result)

            aligned_images.append(align_image)

    # Add the template image back to the aligned images list
    aligned_images.append(template_image)

    return aligned_images


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
    

    
    
    
        

        
            

    
    

    

    
    


        






















            
            

