# Import necessary modules
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfits as pf

import logging
import pywcs

# Header
__author__ = "Tommy Le Blanc"
__version__ = "1.0"

# HISTORY
#    1. Nov 2013 - Vr. 1.0: Added initial PSF tools
#                          - Added mask_psf function
#                          - Added gen_2d_thrumap function
#                          - Added gen_3d_thrumap function


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

        >>    masked_psf, psf_shutter_flux = mask_image(0.0, 0.0, center, psf)

        Generates a masked psf image along with a 5 x 5 summed shutter_flux array
        given an x, y MSA position of (0, 0) and a PSF centered at center, and an 
        initial image psf.
    """

    # Locate pointing center based on input x and y centering parameters
    # (for MSA grid mask)
    centerx = c[0]-(x*(dx/2.))
    centery = c[1]-(y*(dy/2.))

    if verbose: print("Pointing center: [{},{}]".format(centerx, centery))
    
    # Model the MSA grid by zeroing out flux in grid pixels
    # Start by copying the scaled psf image
    msa_mask = np.ones(image.shape)
            
    # Calculate dimensions (in pixels) of a 5 x 5 MSA shutter array
    xstart = centerx - (136/2.)   #x dim of MSA = 136px
    xu = centerx + (136/2.) 
    ystart = centery - (236/2.)   # ydim of MSA = 236px
    yu = centery + (236/2.)

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




# Print diagnostic load message
print("(psfutils): PSF Utilities Version {} loaded!".format(__version__))