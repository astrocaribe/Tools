# encoding: utf-8
# Import necessary modules
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfits as pf
from mpl_toolkits.axes_grid1 import ImageGrid

import logging
import pywcs

# Header
__author__ = "Tommy Le Blanc"
__version__ = "1.0.4"

# HISTORY
#    1. Nov 2013 - Vr. 1.0: Added initial PSF tools
#                          - Added mask_psf function
#                          - Added gen_central_thrumap function
#                          - Added gen_3x3_thrumap function
#                  Vr. 1.0.1: Minor updates to gen_central)thrumap and
#                             gen_3x3_thrumap functions
#    2. Feb 2014 - Vr. 1.0.2: Updates to gen_central_thrumap
#                          - Include either an input filename or a 
#                            numpy datacube in input call
#                          - Changed s_path to outfile to include path 
#                            info as well as filename
#    3. Mar 2014 - Vr. 1.0.3: Minor updates to gen_central_thrumap
#                          - Included an option to display a grid 
#                            pattern on map
#    4. Mar 2014 - Vr. 1.0.4: Minor updates to mask_image()
#                          - Updated shutter dimensions as per latest
#                            NIRSpec OCD
#                          - Updated calculation of 5x5 shutter mask 
#                            dimensions for subsequent computations.

# Utility definitions
# *********************** mask_psf ***********************
# Create and apply an MSA shutter mask given the x and y MSA center parameters,
# a centroid location, and an array size for the mask image
def mask_image(x, y, c, image, verbose=False):
    """
    Create and apply an MSA shutter mask pattern on an input PSF

    This functions generates and applied a NIRSpec MSA shutter mask on an input
    PSF, given the x and y MSA center parameters, a centroid location, and an 
    array size for the mask image. The resulting mask is centered on the input 
    c location, effectively zeroing out the PSF where the mask applies.

    Keyword arguments:
    x, y    -- Central MSA pointing coordinates
    c       -- Input PSF centroid coordinates
    image   -- Input PSF image
    verbose -- Toggle verbose mode. (default False).

    Output(s):
    masked_image -- The masked PSF
    shutter_flux -- The flux values, summed per 5 x 5 shutter

    Example use:

        >>    masked_psf, psf_shutter_flux = psfutils.mask_image(0.0, 0.0, center, psf)

        Generates a masked psf image along with a 5 x 5 summed shutter_flux array
        given an x, y MSA position of (0, 0) and a PSF centered at center, and an 
        initial image psf.
    """

    # MSA shutter width, height, and grid width
    # Dimensions are for 10X oversampled pixels, each = 0.01"
    #dx, dy, gridwidth = 20, 40, 6
    dx, dy, gridwidth = 20, 45, 6
    print('Shutter dimensions: ', dx, dy, gridwidth)


    # Locate pointing center based on input x and y centering parameters
    # (for MSA grid mask)
    centerx = c[0]-(x*(dx/2.))
    centery = c[1]-(y*(dy/2.))

    if verbose: print("Pointing center: [{},{}]".format(centerx, centery))
    
    # Model the MSA grid by zeroing out flux in grid pixels
    # Start by copying the scaled psf image
    msa_mask = np.ones(image.shape)
            
    # Calculate dimensions (in pixels) of a 5 x 5 MSA shutter array
    x_shutter_dim = (gridwidth*6) + (dx*5)
    y_shutter_dim = (gridwidth*6) + (dy*5)

    xstart = centerx - (x_shutter_dim/2.)   #x dim of MSA
    xu = centerx + (x_shutter_dim/2.) 
    ystart = centery - (y_shutter_dim/2.)   #y dim of MSA
    yu = centery + (y_shutter_dim/2.)

    # Preserve the x and y MSA start coordinates (Find a better way!)
    # Also store the x and y start/end positions for each cell in the
    # MSA
    xl = xstart
    yl = ystart

    s_xpos = []
    e_xpos = []
    s_ypos = []
    e_ypos = []
    
    #print("X/Y start/end limits: [{}, {}], [{}, {}]".format(xl, yl, xu, yu))
    #print()
    
    # Generate MSA vertical shutter mask
    for bb in np.arange(6):
        xv1 = xstart
        xv2 = xv1 + gridwidth
    
        s_xpos = np.append(s_xpos, xv2)
        e_xpos = np.append(e_xpos, xv1)
        msa_mask[yl:yu+1, xv1:xv2+1] = 0.
        xstart = xv2 + dx
        #print("(mask check): [{}, {}]".format(xv1, xv2))

    # Delete the last and first elements, respectively
    # of the start/end x position arrays
    s_xpos = np.delete(s_xpos, -1)    
    e_xpos = np.delete(e_xpos, 0)

    # Generate MSA horizontal shutter mask
    for cc in np.arange(6):
        yv1 = ystart
        yv2 = yv1 + gridwidth

        s_ypos = np.append(s_ypos, yv2)
        e_ypos = np.append(e_ypos, yv1)
        msa_mask[yv1:yv2+1, xl:xu+1] = 0.
        ystart = yv2 + dy

    # Delete the last and first elements, respectively
    # of the start/end y position arrays
    s_ypos = np.delete(s_ypos, -1)    
    e_ypos = np.delete(e_ypos, 0)
    
    # Apply the mask to the input image
    masked_image = msa_mask * image
    
    # Store and print the flux recovered in each if the 5 x 5 MSA shutters
    shutter_flux = np.empty((5, 5), dtype=np.float32)

    for tt in xrange(len(s_ypos)):
        for uu in xrange(len(s_xpos)):
            shutter_flux[tt, uu] = masked_image[s_ypos[tt]:e_ypos[tt]+1, s_xpos[uu]:e_xpos[uu]+1].sum()
            
    if verbose: print("Total flux per 5x5 shutter:")
    if verbose: print("===========================")
    if verbose: 
        for tt in xrange(len(s_xpos)): 
            print(shutter_flux[4-tt, :])    
            
    return masked_image, shutter_flux
# *********************** mask_psf ***********************


# *********************** gen_central_thrumap ***********************
def gen_central_thrumap(d_cube, wavelength, lvls, outfile='./ThruMap.pdf', infile=False, grid=False):
    """
    Generate a throughput map of the central shutter (from a 5x5 config).

    This function generate a throughput map (percentage of flux) in a 
    central shutter (of a 5x5 shutter configuration), given an input 
    datacube. This version assumes a shape (441, 5, 5), corresponding to 
    21 x and y dithers of a psf within the central shutter.

    Keyword arguments:
    d_cube     -- Thruput datacube
    wavelength -- Wavelength over which the map is being generated
    lvls       -- The contour levels to plot
    s_path     -- The save path for the generated map
    infile     -- File containing the thruput datacube (Optional). If set, 
                  d_cube is ignored
    grid       -- Include/exclude a grid over the map (Optional). If set, grid is
                  drawn.

    Output(s):
    fig        -- The generated throughput map figure

    Example use:

    1.    >>    psfutils.gen_central_thrumap(thruput_cube, 2.1, lvls=[0, 25, 50, 100],
                                     s_path='./Maps')

        Generates a throughput contour map of the central shutter using the input 
        file datacube.fits of a psf generated at 2.1 µm, with contour levels at 
        0, 25, 50, and 100. The map is stored in './Maps'.

    2.    >>    psfutils.gen_central_thrumap(thruput_cube, 2.1, lvls=[0, 25, 50, 100],
                                     s_path='./Maps', infile='datacube.fits')

        Generates a throughput contour map of the central shutter using the input 
        file datacube.fits of a psf generated at 2.1 µm, with contour levels at 
        0, 25, 50, and 100. The map is stored in './Maps'. Even though a thruput_cube 
        has been passed, it is ignored because an input file was included.

    """

    if infile:
        # Read in throughput cube
        hdu = pf.open(infile)
        steps = hdu[0].header['NSTEP']
        cube = hdu[1].data.astype(np.float64)
        hdu.close()

        # Reshape the central shutter values into 2 dimensions
        x = cube[:, 2, 2].reshape(steps,steps)
    else:
        steps = np.sqrt(d_cube.shape[0]).astype(np.int)
        x = d_cube[:, 2, 2].reshape(steps,steps)

    # Central shutter representation
    fig, ax = plt.subplots(figsize=(5,8))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('MSA Shutter Throughput Map ({} $\mu m$)'.format(wavelength))

    # If grid keyword set, include a grid layout over plot
    if grid: ax.grid(axis='both', which='both')

    cmap = cm.get_cmap('jet', 20)    # 20 discrete colors

    cax = plt.contourf(x, cmap=cmap, levels=lvls, antialiased=False, vmin=0.15, vmax=0.7, alpha=0.8)
    plt.axes().set_aspect(2.)

    # Extract and plot throughput max
    t_max = np.where(x == x.max().astype(np.float64))
    max_x, max_y = t_max[1][0], t_max[0][0]
    ax = plt.plot(max_x, max_y, marker='+', ms=10, mfc='black', mec='black')

    cnt = plt.contour(x, colors='black', linewidths=.5, levels=lvls)
    plt.clabel(cnt, fontsize=8)

    fig.savefig(outfile)

    print('MSA Central Shutter Throughput Map ({} µm) done ...'.format(wavelength))

    return fig
# *********************** gen_central_thrumap ***********************


# *********************** gen_3x3_thrumap ***********************
def gen_3x3_thrumap(infile, wavelength, lvls, s_path):
    """
    Generate a throughput map of the 3x3 central shutter region 
    (from a 5x5 config).

    This function generate a throughput map (percentage of flux) in an 
    inner 3x3 shutter (of a 5x5 shutter configuration), given an input 
    datacube. This version assumes a shape (441, 5, 5), corresponding to 
    21 x and y dithers of a psf within the central shutter.

    Keyword arguments:
    infile     -- Central MSA pointing coordinates
    wavelength -- Wavelength over which the map is being generated
    lvls       -- The contour levels to plot
    s_path     -- The save path for the generated map

    Output(s):
    fig        -- The generated throughput map figure

    Example use:

        >>    psfutils.gen_3x3_thrumap('datacube.fits', 2.1, lvls=[0, 25, 50, 100],
                                     s_path='./Maps')

        Generates a throughput contour map of the inner 3x3 shutter region 
        (from 5x5 shutter configuration) using the input file datacube.fits 
        of a psf generated at 2.1 µm, with contour levels at 0, 25, 50, and 
        100. The map is stored in './Maps'.
    """

    # Read in throughput cube
    hdu = pf.open(infile)
    steps = hdu[0].header['NSTEP']
    cube = hdu[1].data.astype(np.float64)

    # Setup figure and 3x3 grid
    fig = plt.figure(figsize=(6,12))
    fig.suptitle('3x3 Throughput Map ({} $\mu m$)'.format(wavelength))

    grid = ImageGrid(fig, 111, # similar to subplot(111)
                     nrows_ncols = (3, 3), # creates 3x3 grid of axes
                     axes_pad=0.25,  # pad between axes in inch.
                     aspect=False,
                     share_all=True,
                     label_mode='1',
                     cbar_mode='single',
                     cbar_set_cax=True)

    axisNum = 0

    cmap = cm.get_cmap('jet', 20)    # 20 discrete colors

    # Contour fill grids    
    for row in xrange(3):
        for column in xrange(3):
        
            x = cube[:, 3-row, 3-column].reshape(21,21)
            im = grid[axisNum].contourf(x, cmap=cmap, antialiased=False, levels=lvls)
            grid[axisNum].text(10,10,axisNum, color='white')
            axisNum += 1

    # Insert a colorbar
    cax = grid.cbar_axes[0].colorbar(im)

    # Save figure
    fig.savefig(s_path+'/Images3D/ThruMap_{}.jpg'.format(wavelength))
    
    hdu.close()
    print('MSA 3x3 Shutter Throughput Map ({} µm) done ...'.format(wavelength))
# *********************** gen_3x3_thrumap ***********************




# Print diagnostic load message
print("(psfutils): PSF Utilities Version {} loaded!".format(__version__))