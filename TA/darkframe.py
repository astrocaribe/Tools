# Generate a dark frame using input fields.
# Inputs: 
#         infiles - The files from which darks are needed
#         a_size - Array size for the dark frame. Used in initial creation
# Output: 
#      darkframe - A darkframe (same dimensions as input images)
#                  containing the mean of all the input images
def darkframe(infiles, a_size):
    
    # Number of files
    num_files = len(infiles)
    
    # Initial empty darkgrid in dimensions of the incoming files
    darkgrid = np.zeros(a_size)
    
    # Iterate through number of input frames and extract
    # mean pixel counts
    for ii in range(num_files):
        data = pf.getdata(infiles[ii])
        darkgrid = darkgrid + data
        
    darkframe = darkgrid/num_files
    
    # Store the resulting dark frame in a FITS file
    hdu = pf.PrimaryHDU(darkframe)
    print "Saving generated dark frame to FITS file:"
    hdu.writeto('DarkFrame.fits', clobber=True)
    
    return darkframe