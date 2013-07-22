# *************************** Seup ***************************
# Setup global path and dark images path
path = '201212-Alternate-TA-simulations/db_SKY_F140X_MIRROR/exposuresA1600/'
dpath = '201212-Alternate-TA-simulations/db_SKY_F140X_MIRROR/exposuresA200/'

# Flag for generating / subtracting a dark
gen_dark = True

# Window half width for source zoom/analysis, A1600 aperture pixel center, 
# and N x N pixel window limits
whw = 16                            # in pixels
a1600_x, a1600_y = 988, 1416
xll, xul = a1600_x-whw, a1600_x+whw
yll, yul = a1600_y-whw, a1600_y+whw

# *************************** Seup ***************************

# *************************** Filenames ***************************
# Setup filenames for the random offsets and peakups
# There are 5 random offsets and 10 peakups

# Build the filenames for random offsets
seq = np.arange(5)
nodark = path+'F140X_MIRROR_A1600_centered_random_offsets_{0:02d}_000_491_nodark.fits'
dark = path+'F140X_MIRROR_A1600_centered_random_offsets_{0:02d}_000_491_{1:02d}.fits'
file_nodark = [nodark.format(seq[ii]+1) for ii in seq]
file_dark = [dark.format(seq[ii]+1, seq[ii]) for ii in seq]
title = ["Random_Offset_#{0:02d}".format(seq[ii]+1) for ii in seq]
    
# Build the filenames for peakups, and append to file list
seq = np.arange(10)
nodark = path+'F140X_MIRROR_A1600_peakup_{0:02d}_000_491_nodark.fits'
dark = path+'F140X_MIRROR_A1600_peakup_{0:02d}_000_491_{1:02d}.fits'

for ii in seq: file_nodark.append(nodark.format(seq[ii]+1))
for ii in seq: file_dark.append(dark.format(seq[ii]+1, seq[ii]))
for ii in seq: title.append("Peak-up_#{0:02d}".format(seq[ii]+1))

# *************************** Filenames ***************************

# *************************** Dark Frames ***************************
# If flag is set (gen_dark), generate a dark frame. Build the filenames
# in the same manner as previous section (5 random offsets, 10 peakups),
# then generate dark with darkframe_gen function...
if gen_dark:
    print "Generating a dark file..."
    seq = np.arange(5)
    dark = dpath+'F140X_MIRROR_A200_centered_random_offsets_{0:02d}_000_491_{1:02d}.fits'
    f_darks = [dark.format(seq[ii]+1, seq[ii]) for ii in seq]
    
    seq = np.arange(10)
    dark = dpath+'F140X_MIRROR_A200_peakup_{0:02d}_000_491_{1:02d}.fits'
    for ii in seq: f_darks.append(dark.format(seq[ii]+1, seq[ii]))

    # Generate the dark frame using the as compiled list
    dark_frame = darkframe_gen(f_darks)
    fig, ax = plt.subplots()
    cax = ax.imshow(dark_frame, cmap='gray')
    fig.colorbar(cax, ax=ax)
    fig.savefig("DarkFrame.pdf")      # Save as PDF for inspection
else:
    print "No dark frame will be generated!"
    
# *************************** Dark Frames ***************************

# *************************** Analysis ***************************

# Generate diagnostic images
for ii in xrange(len(file_nodark)-14):
    
    # Set up plot
    fig = plt.figure()
    axNaked = plt.subplot2grid((2, 2), (0, 0))
    axReal = plt.subplot2grid((2, 2), (0, 1))
    axResidual = plt.subplot2grid((2, 2), (1, 0), colspan = 2)

    # ********************* "Naked-eye" image *********************
    # Read in image...
    print "\nReading file "+file_nodark[ii][67:-1]+" ..."
    naked_image = readimage(file_nodark[ii])
    
    # Set limits and plot "Naked-eye" source (32 x 32 pixels)
    axNaked.set_xlim(xll, xul)
    axNaked.set_ylim(yll, yul)
    for tt in axNaked.get_xticklabels(): tt.set_fontsize(6)
    for tt in axNaked.get_yticklabels(): tt.set_fontsize(6)

    # Plot image and colorbar to have an idea of pixel count spread        
    nax = axNaked.imshow(naked_image[yll:yul, xll:xul], cmap='bone', interpolation='nearest')
    cb = fig.colorbar(nax, ax=axNaked)
    for tt in cb.ax.get_yticklabels(): tt.set_fontsize(6)
        
    # Draw the A1600 aperture (red square) and the aperture ceter (blue square)
    # and brightest pixel (green cross)
    axNaked.plot(a1600_x, a1600_y, marker='s', mec='red', mfc='none', ms=48)
    axNaked.plot(a1600_x, a1600_y, marker='s', mec='blue', mfc='none', ms=3)
    
    # Find and plot the centroid (red filled circle)
    center = centroid(naked_image[yll:yul, xll:xul])
    true_center = "[" + np.str(np.floor(center[0]+xll)) + ", " + np.str(np.floor(center[1]+yll)) + "]"
    print "'Naked-eye' centroid: ", center, true_center
    axNaked.plot(center[0]+xll, center[1]+yll, marker='o', mec='red', mfc='red', ms=3, alpha=0.3)
    
    # Find and mark the brightest pixel
    x_mark = np.where(naked_image == naked_image[yll:yul, xll:xul].max())
    print "Max pixel ('Naked-eye' image):", naked_image[yll:yul, xll:xul].max(), "at ["+str(x_mark[1][0])+", "+str(x_mark[0][0])+"]"
    axNaked.plot(x_mark[1][0], x_mark[0][0], marker='x', mec='green', ms=3, alpha=0.7)
    # ********************* "Naked-eye" image *********************
    
    
    # ********************* "Real" image *********************
    # Read in image...
    print "Reading file "+file_dark[ii][67:-1]+" ..."
    real_image = readimage(file_dark[ii])
    
    # Perform dark subtraction if flag set
    if gen_dark: real_image = real_image - dark_frame
    
    # Set limits and plot "Naked-eye" source (32 x 32 pixels)
    axReal.set_xlim(xll, xul)
    axReal.set_ylim(yll, yul)
    for tt in axReal.get_xticklabels(): tt.set_fontsize(6)
    for tt in axReal.get_yticklabels(): tt.set_fontsize(6)

    # Plot image and colorbar to have an idea of pixel count spread        
    nax = axReal.imshow(real_image, cmap='bone', interpolation='nearest')
    cb = fig.colorbar(nax, ax=axReal)
    for tt in cb.ax.get_yticklabels(): tt.set_fontsize(6)
        
    # Draw the A1600 aperture (red square) and the aperture ceter (blue square)
    # and brightest pixel (green cross)
    axReal.plot(a1600_x, a1600_y, marker='s', mec='red', mfc='none', ms=48)
    axReal.plot(a1600_x, a1600_y, marker='s', mec='blue', mfc='none', ms=3)
    
    # Find and plot the centroid (red filled circle)
    center = centroid(real_image[yll:yul, xll:xul])
    true_center = "[" + np.str(np.floor(center[0]+xll)) + ", " + np.str(np.floor(center[1]+yll)) + "]"
    print "'Real' centroid: ", center, true_center
    axReal.plot(center[0]+xll, center[1]+yll, marker='o', mec='red', mfc='red', ms=3, alpha=0.3)
    
    # Find and mark the brightest pixel (within 32 x 32 sub image)
    x_mark = np.where(real_image == real_image[yll:yul, xll:xul].max())
    print "Max pixel ('Naked-eye' image):", real_image[yll:yul, xll:xul].max(), "at ["+str(x_mark[1][0])+", "+str(x_mark[0][0])+"]"
    axReal.plot(x_mark[1][0], x_mark[0][0], marker='x', mec='green', ms=3, alpha=0.7)
    # ********************* "Real" image *********************
    
    outfile = "Test"+np.str(ii+1)+".pdf"
    fig.savefig(outfile)

# *************************** Analysis ***************************