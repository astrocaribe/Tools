# Import necessary modules
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf

# Header
__author__ = "Tommy Le Blanc"
__version__ = "1.2.1"

# HISTORY
#    1. Jan 2013 - Vr. 1.0: Added initial versions of centroid and bytescl
#                  Vr. 1.1: Changed/added features to centroid function
#    2. Sep 2013 - Vr. 1.2: 
#                          - Added readimage function
#                          - Added diff_image function
#                          - Added a1600_pixtosky function
#                          - Added Rotate2D function
#                          - Added ta_transform function
#                          - Changed the print function to the Python 3.0 format
#                            for future compatibility.
#    3. Sep 2013 - Vr. 1.2.1:
#                          - Updated readimage to read from ext=1

# Utility definitions
# *********************** centroid ***********************
def centroid(grid, gw, initial=None, debug=False, verbose=False):
    '''
    ABOUT: 
        This function takes the input image grid, as well as a grid width 
        (in pixel units) over which to do a final centroiding. An initial guess 
        is executed over the entire input image, and then another centroiding is 
        iterated given the initial computed guess. The final output is a rounded 
        centroid coordinate (rounded to account for any loss of accuracy due to 
        using floors throughout computation).

    INPUT(S):
        grid - Image for which centrod is to be computed. For now, an NxN image
            is required.

        gw - The grid width (in pixels) over which to compute the final centroid.
        initial - The optional initial guess coordinates. If none is supplied, 
            then the approximate center (+/- 1 px) is used. Note that the co-
            ordinates are recognized in normal order (i.e, initial=[x,y]).

        initial - An initial guess for targeted centroiding, in the form (x, y).
            Defaults to None if the keyword is omiited.

        debug - Toggle debugging mode. Defaults to False, in which case no 
            diagnostic image is produced or passed.

    OUTPUT(S):
        centroid - The computed centroid coordinates. The 
            final coordinates are relative to grid dimensions.

        fig - diagnostic figure generated to troubleshoot the centroiding 
            calculation. The figure contains 4 panels; i) the image (grid) for 
            which centroiding is being calculated, ii), iii) the weights and 
            weighted totals for x and y in grid, and iv) a scaled view of the 
            image with final centroid (if updated from an initial calculation).

    USAGE:
        Can be used in one of two ways:
        1. Debug mode:
             
             center, debugfig = centroid(img, cgw, debug=True)
             
             In this mode, the centroiding algorithm will estimate an initial centroid
             of img, given a centroiding window cgw. Setting debug=True produces a 
             diagnostic figure that is passed (along with the centroid) to the calling 
             script; this figure can be captured and printed for analysis.
        2. Centroid calculation mode:
             
             center = centroid(img, cgw, initial=init)

             In this mode, you can calculate a targeted centroid of img, given a 
             centroiding window cgw (of a few pixels), and an initial centroid guess. 
             The result is the calculated centroid for img in the form (x, y).

    HISTORY:
        Ver. 1.0: January 2013
            - Original function; calculated weights and weighted averages by 
            manually summing over row/columns.
        Ver. 2.0: February 2013
            - Optimizes computation by summing over row/column in one go.
        Ver. 2.1: 12th February, 2013
            - Modified to include a smaller centroiding computation. An 
            initial guess is computed over the entire NxN input grid, and then a 
            tighter computation over a smaller nxn pixel grid is conducted.
            - Included centroiding window edge detection when determining 
            new computation window
        Ver. 2.2: 18th February, 2013
            - Modified to include built in centroiding initial guessing by way
            of an optional keyword; initial. Instead of including an initial 
            guess mechanism, an "initial" keyword is added so that the function 
            can be used to perform an initial guess, and re-running this function
            with the initial guess will perform a more targeted cetroiding.
        Ver. 2.3: 7th March, 2013
            - Added a debug keyword for returning degub figures back to the main
            program for saving/viewing.
        Ver. 2.31: 7th August, 2013
            - Updated the debug portion of scripts by adding plot and figure titles.
        Ver. 2.4: 8th August, 2013
            - Included a "verbose" flag to turn off diagnostic output to screen.    
    '''
    __author__ = "Tommy Le Blanc"
    __version__ = "2.4"
    
    import logging
    
    # Set the initial guess coordinates
    if initial == None: 
        if verbose: print("Initial Guess (x,y)...")
        init_guess = ((grid.shape[1]/2) - 1, (grid.shape[0]/2) - 1)
    else:
        if verbose: print("Centroid...")
        init_guess = initial #(initial[1], initial[0])
        
    if verbose: print("(centroid) Centroid Search Position: ({:.2f}, {:.2f})".format(init_guess[0], init_guess[1]))    
    

    # ************* Centroiding *************
    
    # Define centroiding limits based on input grid width
    # specifications (grid width x/y lower and upper limits)
    gw_xl = init_guess[0] - (gw/2)
    gw_xu = gw_xl + gw - 1
    gw_yl = init_guess[1] - (gw/2)
    gw_yu = gw_yl + gw - 1

    
    # Check wether the new centroid window passes the edges, and if so, 
    # reset window to edges
    if gw_xl < 0:
        gw_xl = 0.
        gw_xu = gw_xl + gw - 1. 
    elif gw_xu >= grid.shape[1]:
        gw_xu = grid.shape[1] - 1.
        gw_xl = grid.shape[1] - gw 

    if gw_yl < 0:
        gw_yl = 0.
        gw_yu = gw_yl + gw - 1.
    elif gw_yu >= grid.shape[0]:
        gw_yu = grid.shape[0] - 1.
        gw_yl = grid.shape[0] - gw
        
    # Define optimized centroiding grid
    newGrid = grid[gw_yl:gw_yu+1, gw_xl:gw_xu+1]

    
    # ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ******
    if debug:
        fig = plt.figure(figsize=(8,8))
        axWeightX = plt.subplot(221)
    
        axImgScaled = plt.subplot(222)
        axImg = plt.subplot(223)
        axWeightY = plt.subplot(224)
    
        axWeightX.set_xlim(0, newGrid.shape[1]-1)
        axWeightY.set_ylim(0, newGrid.shape[0]-1)

        # Set display min/max
        if grid.shape[0] > 16:
            va, vb = 0, 1.e-4
        else:
            va, vb = 0, 0.15

        axImg.imshow(newGrid, cmap='gray', interpolation='nearest', vmin=va, vmax=vb)
        axImg.autoscale(axis='both', enable=False)
    
        axImgScaled.set_title('Centroid window: {}'.format(newGrid.shape))
        cax=axImgScaled.imshow(newGrid, cmap='jet', interpolation='nearest', vmin=va, vmax=vb)
        axImgScaled.autoscale(axis='both', enable=False)
    # ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ******
  
    
    #Diagnostic output
    if verbose: print("(centroid) Centroiding window: ({:.2f}, {:.2f}), ({:.2f}, {:.2f})".format(gw_xl, gw_yl, gw_xu, gw_yu))
    if verbose: print("(centroid) Total Pixel Count:", np.sum(newGrid))
    if verbose: print("(centroid) Average Pixel Count:", np.average(newGrid))
    if verbose: print("(centroid) Grid Shape:", newGrid.shape)

    # Weights and weighted averages
    weight_i = np.sum(newGrid, 0)
    weight_j = np.sum(newGrid, 1)
    weight_i = np.where(weight_i < 0, 0, weight_i)
    weight_j = np.where(weight_j < 0, 0, weight_j)
    
    weight_avg_i = weight_i * np.arange(len(newGrid))
    weight_avg_j = weight_j * np.arange(len(newGrid))
    
        
    # X and Y position of the centroid    
    x = np.sum(weight_avg_i)/np.sum(weight_i)
    y = np.sum(weight_avg_j)/np.sum(weight_j)
    centroid = (x+gw_xl, y+gw_yl)
    if verbose: print("(centroid) Raw X/Y and centroid: ({:.2f}, {:.2f})".format(x, y))
    # ************* Centroiding *************
  
    
    
    # ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ******
    if debug:
        temp_i = list(enumerate(weight_i))
        px, py = zip(*temp_i)
        axWeightX.plot(px,py, marker='o', mfc='blue', mec='blue', ms=4, alpha=0.5)

        temp_j = list(enumerate(weight_j))
        px, py = zip(*temp_j)
        axWeightY.plot(py,px, marker='o', mfc='blue', mec='blue', ms=4, alpha=0.5)

        temp_ai = list(enumerate(weight_avg_i))
        px, py = zip(*temp_ai)
        axWeightX.plot(px,py, marker='o', mfc='green', mec='green', ms=4, alpha=0.5)

        temp_aj = list(enumerate(weight_avg_j))
        px, py = zip(*temp_aj)
        axWeightY.plot(py,px, marker='o', mfc='green', mec='green', ms=4, alpha=0.5)
    
        axImg.plot(x, y, marker='o', mfc='red', mec='red', alpha=0.5)
        axImgScaled.plot(x, y, marker='o', mfc='red', mec='white', alpha=0.5)
    
        fig.tight_layout()
        fig.colorbar(cax, ax=axImgScaled, shrink=0.75, ticks=[va, 0, vb])
    # ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ****** DEBUG! ******
 
    
    
    if verbose: print("(centroid) New centroid: ({:.2f}, {:.2f})".format(centroid[0], centroid[1]))
    
    # If in debugging mode, return the figure as well, else just return the centroid
    if debug: return centroid, fig
    else: return centroid
# *********************** centroid ***********************


# *********************** bytescl ***********************
def bytescl(img, bottom, top):
    '''
    ABOUT:
        This function scales a pixel image with arbitrary (bottom, top) limits,  and 
        changes them to (0, 1), (i.e. - Changes the limit from bottom - top to 0 - 1).

    INPUT:
        bottom, top - Lower and upper limits of the original scale.

    OUTPUT:
        scl_img - Scaled image with new limits 0(min) - 1(max).

    HISTORY:
        Ver. 1: January 2013
            - Original function.    
    '''
    scl_img = (((top - bottom) * (img - img.min())) / (img.max() - img.min())) + bottom
    
    return scl_img
# *********************** bytescl ***********************


# *********************** readimage ***********************
# Extract an image from a multi-ramp integration FITS file
def readimage(infile):
    '''
    ABOUT:
        Extract am image from a multi-ramp integration FITS file. Currently, JWST NIRSpec FITS
        images consists of a 3-ramp integration, with each succesive image containing more photon
        counts than the next. This meathod of reading/combining the images ensures the elimination
        of cosmic rays that may occur in a few of the frames, as well (more importantly) reduce the
        impact of hot pixels by chossing the minimum pixel count on a pixel by pixel basis.

    INPUT(S):
        infile - A multi-frame FITS image. In this case, a 3-frame image is expected as per test specs.
                 **Note** The input image must be in ext=1 for this version to work.

    OUTPUT(S):
        omega - A combined FITS image that compares the differenced 3-frame image (generating alpha, beta
            image frames), and selecting the lowest count pixel between the two.

     HISTORY:
        Ver. 1: May 2013
            - Original function            
    '''
    
    # Read in input file, and generate the alpha and beta images
    # (for subtraction)
    master = pf.getdata(infile, 1)
    alpha = master[1, :, :] - master[0, :, :]
    beta = master[2, :, :] - master[1, :, :]
    
    # Generate a final image by doing a pixel to pixel check 
    # between alpha and beta images, storing lower value
    omega = np.where(alpha < beta, alpha, beta)
     
    print('(readimage): Image {} read ...'.format(infile))
        
    # Return the extracted image
    return omega
# *********************** readimage ***********************


# *********************** diff_image ***********************
# New method for conducting image differencing
def diff_image(im1, im2):
    '''
    ABOUT: 
           This function takes as input the images to be subtracted,
           and output the final subtracted image. The order in which
           the input images are input is important; the 2nd will be 
           subtracted from the first.
    INPUT:
           im1, im2 - The images for which differencing is to be 
           perfomred.
    OUTPUT: 
           output_im - Final difference image. Dimensions match that
           of im1 and im2.

    HISTORY:
        Ver. 1: May 2013
           - Original function
    '''
    
    # Subtract the min of im1 and im2 for im1
    output_im = im1 - np.where(im1 < im2, im1, im2)
    
    return output_im
# *********************** diff_image ***********************



# *********************** a1600_pixtosky ***********************
def a1600_pixtosky(angle_d, pixcrd, origin):
    '''
    ABOUT: 
        This function converts input pixels coordinates into sky coordinates for 
        the A1600 aperture on the NIRSpec MSA.

    INPUT:
        angle - The rotation angle between the aperture grid and the sky coordinates, 
                in degrees.

        (x_pix, y_pix) - The pixel coordinate that needs to be converted. This can
            be a single coordinate, or a numpy array of coordinates.

    OUTPUT:
        (x_sky, y_sky) - The computed centroid coordinates, rounded to nearest 1. The 
            final coordinates are relative to grid dimensions. Same length as 
            (x_pix, y_pix) above.

    HISTORY:
        Ver. 1.0: January 2013
            - Original function
    '''
    __author__ = "Tommy Le Blanc"
    __version__ = "1.0"
    
    # Import necessary modules
    import pywcs
    
    # Change angle from degrees to radians
    angle_r = np.radians(angle_d)
    
    # Create a new WCS object of two axes
    wcs = pywcs.WCS(naxis=2)

    # Setup a pixel projection
    wcs.wcs.crpix = [7.5, 7.5]
    wcs.wcs.crval = [321.536442, 5.689362]
    wcs.wcs.cdelt = np.array([0.106841547, 0.108564798])
    wcs.wcs.pc = ([np.sin(angle_r),np.cos(angle_r)],
                  [np.cos(angle_r),-np.sin(angle_r)])


    # Calculate new coordinates for each pair of input coordinates
    world = wcs.wcs_pix2sky(pixcrd, origin)
    
    return world
# *********************** a1600_pixtosky ***********************


# *********************** ta_transform ***********************
def ta_transform(pnts, mirror=False):
    '''
    ABOUT: 
           This function transforms extracted pixel positions to sky
           in x-y space. For this application, a mirror along the 
           aperture "horizaontal" is included (using the 'Mirror'
           keyword) to correctly represent the converted coordinates. 
    INPUT:
           pnts - A list of pixel coordinates. Can be only one pair, 
                  or several.
    OUTPUT: 
           transformed_pnts - Transformed pixel coordinates, in sky.

     HISTORY:
         Ver.1: May 2013
             - Original function
    '''
    # Import necessary modules
    import numpy as np

    # Function constants
    
    # Transformation derivatives (as per pixel)
    dxdx, dydx = .076839, -0.068522
    dxdy, dydy = .069291, .079071
    
    # Reference point in pixels and sky coordinates, respectively
    Xo_pix, Yo_pix = 989.204158, 1415.929442
    Xo_sky, Yo_sky = 321.536442, 5.689362
    
    
    # Split the entered tuples into seperate variables
    #print('pnts LENGTH!!! = {}'.format(len(pnts)))
    if len(pnts) > 1:
        X_pix, Y_pix = np.array(zip(*pnts))
    else:
        X_pix, Y_pix = pnts[0][0], pnts[0][1]

        
    print('(ta_tranform): X_pix {}'.format(X_pix))
    print('(ta_tranform): Y_pix {}'.format(Y_pix))
    print('(ta_tranform): type {}'.format(type(X_pix)))
    
    ## Convert angle to radians
    #ang = np.deg2rad(angle)
    
    # Calculate the coordinate transforms
    X_temp = Xo_sky + (dxdx*(X_pix-Xo_pix)) + (dxdy*(Y_pix-Yo_pix))
    Y_temp = Yo_sky + (dydx*(X_pix-Xo_pix)) + (dydy*(Y_pix-Yo_pix))
    
    
    if mirror:
        # Rotate/Mirror the tranformed coordinated
        if len(pnts) > 1:
            temp = np.array(zip(X_temp, Y_temp))
        else:
            temp = np.array((X_temp, Y_temp))
        
        cnt = np.array((321.536442, 5.689362))
        
        a, b = 0.992682, 0.120757
        transformed_pnts = np.dot(temp-cnt,np.array([[a,b],[b,-a]]))+cnt
        
        print('(ta_tranform): Mirroring at y=x included!\n')
    
    else:
        print('(ta_tranform): No mirroring included')
        if len(pnts) > 1:
            transformed_pnts = np.array(zip(X_temp, Y_temp))
        else:
            transformed_pnts = np.array((X_temp, Y_temp))
    
    return transformed_pnts
# *********************** ta_transform ***********************


# *********************** Rotate2D ***********************
def Rotate2D(pts,cnt,angle):
    '''
    ABOUT:
        Rotates a set of points in 2D space around a center point.
    
    INPUT(S):
        pts - A set/list of points that need to be rotated.
        cnt - The center of rotation. Coordinate in the form of (x, y) expected.
        angle - The angle of rotation that needs to be computed, in radians.

    OUTPUT(S):
        rot_pts - The resulting points after rotation has been calculated.

    HISTORY:
        Ver. 1: Jan 2013
            - Original function
    '''

    import numpy as np
    
    ang = np.deg2rad(angle)
    
    rot_pts = np.dot(pts-cnt,np.array([[np.sin(ang),np.cos(ang)],[np.cos(ang),np.sin(-ang)]]))+cnt
    return rot_pts
# *********************** Rotate2D ***********************




# Print diagnostic message
print("(tautils): TA Utilities Version {} loaded!".format(__version__))
