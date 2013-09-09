# Import necessary modules
import numpy as np
import matplotlib.pyplot as plt

# Header
__author__ = "Tommy Le Blanc"
__version__ = "1.1"

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
        centroid - The computed centroid coordinates, rounded to nearest 1. The 
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
        if verbose: print "Initial Guess (x,y)..."
        init_guess = ((grid.shape[1]/2) - 1, (grid.shape[0]/2) - 1)
    else:
        if verbose: print "Centroid..."
        init_guess = initial #(initial[1], initial[0])
        
    if verbose: print "(centroid) Centroid Search Position: ({:.2f}, {:.2f})".format(init_guess[0], init_guess[1])    
    

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
    if verbose: print "(centroid) Centroiding window: ({:.2f}, {:.2f}), ({:.2f}, {:.2f})".format(gw_xl, gw_yl, gw_xu, gw_yu)
    if verbose: print "(centroid) Total Pixel Count:", np.sum(newGrid)
    if verbose: print "(centroid) Average Pixel Count:", np.average(newGrid)
    if verbose: print "(centroid) Grid Shape:", newGrid.shape

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
    if verbose: print "(centroid) Raw X/Y and centroid: ({:.2f}, {:.2f})".format(x, y)
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
 
    
    
    if verbose: print "(centroid) New centroid: ({:.2f}, {:.2f})".format(centroid[0], centroid[1])
    
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



# Print diagnostic message
print "(tautils): TA Utilities Version {} loaded!".format(__version__)
