# Setup global path
path='201212-Alternate-TA-simulations/db_SKY_F140X_MIRROR/exposuresA1600/'

# Frame of interest
frame = 0

# Flag for generating / subtracting a dark
gen_dark = True

# Window half width for source zoom/analysis, A1600 aperture pixel center, 
# and N x N pixel window limits
whw = 16                            # in pixels
a1600_x, a1600_y = 988, 1416
xll, xul = a1600_x-whw, a1600_x+whw
yll, yul = a1600_y-whw, a1600_y+whw

# Setup filenames for the random offsets and peakups
# There are 5 random offsets and 10 peakups

# Build the filenames for random offsets
seq_ro = np.arange(5)
file_nodark = []
file_dark = []
title = []
nodark = 'F140X_MIRROR_A1600_centered_random_offsets_{0:02d}_000_491_nodark.fits'
dark = 'F140X_MIRROR_A1600_centered_random_offsets_{0:02d}_000_491_{1:02d}.fits'

for ii in seq_ro: file_nodark.append(nodark.format(seq_ro[ii]+1))
for ii in seq_ro: file_dark.append(dark.format(seq_ro[ii]+1, seq_ro[ii]))
for ii in seq_ro: title.append("Random_Offset_#{0:02d}".format(seq_ro[ii]+1))
    
# Build the filenames for peakups, and append to file list
seq_ro = np.arange(10)
nodark = 'F140X_MIRROR_A1600_peakup_{0:02d}_000_491_nodark.fits'
dark = 'F140X_MIRROR_A1600_peakup_{0:02d}_000_491_{1:02d}.fits'

for ii in seq_ro: file_nodark.append(nodark.format(seq_ro[ii]+1))
for ii in seq_ro: file_dark.append(dark.format(seq_ro[ii]+1, seq_ro[ii]))
for ii in seq_ro: title.append("Peak-up_#{0:02d}".format(seq_ro[ii]+1))


if gen_dark:     
    # Generate a dark frame for use with centroiding the "Real" source
    inputfiles = []
    print "Using the following files for dark frame generation:"
    for ii in range(5): 
        inputfiles.append(path+file_dark[ii])
        print file_dark[ii]
    
    # Generate the dark fram using the ae compiled list
    dark_frame = darkframe(inputfiles, (2048,2048))
    imshow(dark_frame[frame, :, :], cmap='gray')
else:
    print "No dark frame will be generated!"
    
# Generate a 2x2 diagnostic image for each file     
for ii in range(len(file_nodark)-14):

    # Create a 2x2 diagnostic image
    fig, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2, 3, figsize=(9., 6.))
    fig.suptitle(title[ii]+" (Frame #"+str(frame+1)+" of 3)")
    fig.text(.5, 0.05, file_nodark[ii],  va="center", ha="center")

    
    # Extract data from the "No dark" image...
    nodark_data = pf.getdata(path+file_nodark[ii])
    nodark_field = nodark_data[frame, :, :]
    print
    print file_nodark[ii] #, nodark_field.min(), nodark_field.max()
    
    # ... and plot the field and zoomed source
    ax1.imshow(nodark_field, cmap='gray', interpolation='nearest', vmax=5.0)
    ax1.set_title('"Naked-eye" Field', size='small')
    ax2.imshow(nodark_field, cmap='gray', interpolation='nearest', vmax=50.0)
    ax2.set_xlim(xll, xul)
    ax2.set_ylim(yll, yul)
    ax2.plot(a1600_x, a1600_y, marker='s', mec='red', mfc='none', ms=96, alpha=0.7)
    ax2.plot(a1600_x, a1600_y, marker='s', mec='blue', mfc='none', ms=6, alpha=0.7)
    
    # Include an X to mark brightest pixel...
    x_mark = np.where(nodark_field == nodark_field.max())
    print "Max pixel ('Naked-eye' image):", nodark_field.max(), "at ["+str(x_mark[1][0])+", "+str(x_mark[0][0])+"]"
    ax2.plot(x_mark[1][0], x_mark[0][0], marker='x', mec='green', ms=4, alpha=0.7)
    
    # ...and calculate and mark the flux-weighted centroid of the PSF (red filled circle)
    center = centroid(nodark_field[yll:yul, xll:xul])
    print "Naked-eye Centroid:", center, center[0]+xll, center[0]+yll
    ax2.plot(center[0]+xll, center[1]+yll, marker='o', mec='red', mfc='red', ms=6, alpha=0.3)
    ax2.set_title('"Naked-eye" Source (Zoomed)', size='small')
    
    # Plot the ramp-ups of the "Naked-eye" image brightest, centroid and aperture center pixel levels
    nodark_ramp = nodark_data[:, x_mark[0][0], x_mark[1][0]]
    centroid_ramp = nodark_data[:, center[0]+yll, center[1]+xll]
    ap_center_ramp = nodark_data[:, a1600_y, a1600_x]
    ax5.set_title('Normalized Pixel Ramp', size='small')
    ax5.plot(nodark_ramp/nodark_ramp.min(), 'gx-', alpha=0.5, label = "Brightest")
    ax5.plot(centroid_ramp/centroid_ramp.min(), 'ro-', alpha=0.5, label = "Centroid")
    ax5.plot(ap_center_ramp/ap_center_ramp.min(), 'bs-', alpha=0.5, label = "Ap. Center")
    ax5.legend(numpoints=1, frameon=False, prop={'size':'small'}, loc=2)
    
    # Extract data from the "Real", dark included image...
    dark_data = pf.getdata(path+file_dark[ii])
    dark_field = dark_data[frame, :, :]
    
    # Subtract the generated dark frame if flag set
    if gen_dark: dark_field = dark_field - dark_frame[frame, :, :]
    
    # ... and plot the field and zoomed source
    ax3.imshow(dark_field, cmap='gray', interpolation='nearest')
    ax3.set_title('"Real" Field', size='small')
    ax4.imshow(dark_field, cmap='gray', interpolation='nearest')
    ax4.set_xlim(xll, xul)
    ax4.set_ylim(yll, yul)
    ax4.plot(a1600_x, a1600_y, marker='s', mec='red', mfc='none', ms=96, alpha=0.7)
    ax4.plot(a1600_x, a1600_y, marker='s', mec='blue', mfc='none', ms=6, alpha=0.7)

    center = centroid(dark_field[yll:yul, xll:xul])
    print "Real Centroid:", center, center[0]+xll, center[0]+yll
    
    ax4.plot(center[0]+xll, center[1]+yll, marker='o', mec='red', mfc='red', ms=6, alpha=0.3)
    ax4.set_title('"Real" Source (Zoomed)', size='small')
  
    
    # Include an X to mark brightest pixel... (Need to edit!)
    x_mark = np.where(dark_field == dark_field[yll:yul, xll:xul].max())
    print yll, yul, xll, xul, dark_field[yll:yul, xll:xul].max()
    print "Max pixel ('Real' image):", dark_field.max(), "at", x_mark[1][0], x_mark[0][0]
    ax4.plot(x_mark[1][0], x_mark[0][0], marker='x', mec='green', ms=4, alpha=0.7)
    
    # Plot the ramp-ups of the "Real" image brightest, centroid and aperture center pixel levels
    # Begin by subtracting the generated dark exposure
    dark_data = dark_data - dark_frame
    dark_ramp = dark_data[:, x_mark[0][0], x_mark[1][0]]
    centroid_ramp = dark_data[:, center[0]+yll, center[1]+xll]
    ap_center_ramp = dark_data[:, a1600_y, a1600_x]
    ax6.set_title('Normalized Pixel Ramp', size='small')
    ax6.plot(dark_ramp/dark_ramp.min(), 'gx-', alpha=0.5)
    ax6.plot(centroid_ramp/centroid_ramp.min(), 'ro-', alpha=0.5)
    ax6.plot(ap_center_ramp/ap_center_ramp.min(), 'bs-', alpha=0.5)
 
    # Save figure tn external file    
    out_file = title[ii]+'.pdf'
    fig.savefig(out_file)
    
    print out_file, "is ready!"
    print