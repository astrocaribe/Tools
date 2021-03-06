from __future__ import print_function
import numpy as np

# Header
__author__ = "Tommy Le Blanc"
__version__ = "1.0.1"

# HISTORY
#    1. Feb 2014 - Vr. 1.0: Added initial versions the JWST target locate
#                           algorithm. Added the following:
#                           i) 


# Utility definitions
# *************************** checkbox_2D ***************************
def checkbox_2D(image, checkbox, debug=False):
    """
    Find the course location of an input psf by finding the 
    brightest checkbox.
    
    This function uses a 2 dimensional image as input, and
    finds the the brightest checkbox of given size in the
    image.
    
    Keyword arguments:
    image    -- 2 dimensional psf image
    checkbox -- A sliding partial filter that equal the sum
                of values in an n x n region centered on the
                current pixel, where n is an odd integer.
    
    Output(s):
    checkbox_ctr -- A tuple containing the brightest checkbox 
                    location.
    checkbox_hfw -- A tuple containing the checkbox halfwidth.
    
    Example usage:
    
        >> cb_cen, cb_hw = checkbox_2D(psf, 5)
        
        Find the location of the brightest checkbox, given a
        checkbox size of 5. Returns the brightest checkbox 
        center and halfwidths.
    """
    
    # Calculate the checkbox half-width
    chw = (checkbox - 1) / 2
    
    # Calculate the image size
    xsize, ysize = image.shape[1], image.shape[0]
    
    # Calculate the x and y widths of checkbox region
    xwidth, ywidth = xsize - checkbox + 1, ysize - checkbox + 1
    
    # If the checkbox size is not equal to both the X and Y sizes, 
    # find the pixel with the brightest checkbox
    if checkbox != xsize and checkbox != ysize:
        xpeak = 0
        ypeak = 0
        sumpeak = 0
        for ii in xrange(xsize - checkbox):
            for jj in xrange(ysize - checkbox):
                t = np.sum(image[jj:jj+checkbox, ii:ii+checkbox])
                if t > sumpeak:
                    xpeak = ii + chw + 1
                    ypeak = jj + chw + 1
                    sumpeak = t
        
        print('(checkbox_2D): Checkbox not equal to both x/ysize.')

    
    # If the checkbox size is equal to both the X and Y sizes
    if checkbox == xsize and checkbox == ysize:
        xpeak = xsize / 2
        ypeak = ysize / 2
        sumpeak = np.sum(image, axis=None)
        
        print('(checkbox_2D): Checkbox equal to x/ysize.')
        
    # Print calculated checkbox center, and sum within checkbox centroid

    # Find the checkbox region half-width in x and y
    xhw = (xwidth / 2) - 1
    yhw = (ywidth / 2) - 1
        
    if xpeak < xhw or xpeak > xsize - xhw or ypeak < yhw or ypeak > ysize - yhw:
        print('(checkbox_2D): WARNING - Peak too close to edge of image.')
        
#    NOTE: Use this section of the input image is a subset of a larger image
#          Not currently needed for this analysis
#    # Determine the center of the brightest checkbox, in extracted
#    # image coordinates
#    xpeak = xpeak + xhw
#    ypeak = ypeak + yhw

    # Debug messages
    if debug:
        print('(checkbox_2D): chw = ', chw)
        print('(checkbox_2D): xsize, ysize = {}, {}'.format(xsize, ysize))
        print('(checkbox_2D): xwidth, ywidth = {}, {}'.format(xwidth, ywidth))
        print('(checkbox_2D): xpeak, ypeak = {}, {}'.format(xpeak, ypeak))
        print('(checkbox_2D): sumpeak = ', sumpeak)
        print('(checkbox_2D): xhw, yhw = {}, {}'.format(xhw, yhw))
        
    checkbox_ctr = np.array((xpeak, ypeak))
    checkbox_hfw = np.array((xhw, yhw))

    return checkbox_ctr, checkbox_hfw
# *************************** checkbox_2D ***************************

# *************************** checkbox_1D ***************************
def checkbox_1D(image, checkbox, debug=False):
    """
    Find the course location of an flattened input psf by 
    finding the brightest checkbox.
    
    This function uses an image as input, flattens it into a 
    vector and finds the the brightest checkbox of given size 
    in the image.
    
    Keyword arguments:
    image    -- 2 dimensional psf image
    checkbox -- A sliding partial filter that equal the sum
                of values in an n x n region centered on the
                current pixel, where n is an odd integer.
    
    Output(s):
    xpeak -- The brightest checkbox location.
    xhw   -- Checkbox halfwidth.

    
    Example usage:
    
        >> cb_center, cb_hw = checkbox_1D(vector, 5)
        
        Find the location of the brightest checkbox in a vector, 
        given a checkbox size of 5. Returns the brightest checkbox 
        center and halfwidth.
    """
    
    # Collapse input image, currently onto X axis
    # Reshape to reflect collapse onto x axis
    vector = np.sum(image, axis=0)
    print('(checkbox_1D): Image collapsed into 1D vector.')
    
    # Calculate the checkbox half-width
    chw = (checkbox - 1) / 2

    
    # Calculate the image size
    xsize, ysize = image.shape[1], image.shape[0]
    
    # Calculate the x and y widths of checkbox region
    xwidth = xsize - checkbox + 1

    # If the checkbox size is not equal to both the X and Y sizes, 
    # find the pixel with the brightest checkbox
    if checkbox != xsize and checkbox != ysize:
        xpeak = 0
        ypeak = 1
        sumpeak = 0
        for ii in xrange(xsize - checkbox):
            t = np.sum(vector[ii:ii+checkbox])
            if t > sumpeak:
                xpeak = ii + 1
                sumpeak = t

        print('(checkbox_1D): Checkbox not equal to xsize.')
                
                
    # If the checkbox size is equal to both the X and Y sizes
    if checkbox == xsize:
        xpeak = xsize / 2
        sumpeak = np.sum(vector, axis=None)
        
        print('(checkbox_1D): Checkbox equal to xsize.')
        
    # Print checkbox center and peak around centroid region

    # Find the checkbox region half-width in x and y
    xhw = (xwidth / 2) - 1
        
    if xpeak < xhw or xpeak > xsize - xhw:
        print('(checkbox_1D): WARNING - Peak too close to edge of image.')
        
    
    # Debug messages
    if debug:
        print('(checkbox_1D): chw = ', chw)
        print('(checkbox_1D): xhw = ', xhw)
        print('(checkbox_1D): xsize = ', xsize)
        print('(checkbox_1D): xwidth = ', xwidth)
        print('(checkbox_1D): xpeak = ', xpeak)
        print('(checkbox_1D): sumpeak = ', sumpeak)
        
#    NOTE: Use this section of the input image is a subset of a larger image
#          Not currently needed for this analysis
#    # Determine the center of the brightest checkbox, in extracted
#    # image coordinates
#    xpeak = xpeak + xhw
    
    return xpeak, xhw
# *************************** checkbox_1D ***************************

# *************************** centroid_2D ***************************
def centroid_2D(image, checkbox_center, checkbox_halfwidth, max_iter=0, threshold=0, 
	            debug=False):
    """
    Fine location of the target by calculating the centroid for 
    the region centered on the brightest checkbox.
    
    Performs the centroid calculation on the checkbox region 
    calculated using the function checkbox_2D().
    
    Keyword arguments:
    image              -- 2 dimensional psf image
    checkbox_center    -- The location of the brightest checkox, 
                          in x and y, in the input image. (pix)
    checkbox_halfwidth -- The halfwidths in both x and y of the 
                          checkbox from centroid_2D(). (pix)
    max_iter           -- Max number of iterations for fine 
                          centroiding. Ignored if convergence is
                          reached beforehand.
    threshold          -- Threshold for successful convergence.
    
    
    Output(s):
    centroid  -- Tuple containing the location of the target, 
                 in the format [x, y].
    c_cum     -- The calculated flux sum within the checkbox 
                 region.
    
    Example usage:
    
        >> cb_centroid, cb_sum = centroid_2D(psf, cb_cen, cb_hw, 
                        max_iter=5, threshold=0.0001, debug=True)
        
        Find the fine centroid location of the target in psf, given
        both the checkbox center and halfwidths. Do a fine centroid 
        iteration of a maximum of 5, unless an threshold of 0.0001 
        is reached first
    """
    
    # First calculate centroid to use for the first iteration
    c_sum = 0
    xsum = 0
    ysum = 0
    
    convergence_flag = 'N/A'
    
    # Unpack the checkbox_center and checkbox_halfwidth into 
    # their appropriate variables
    xpeak, ypeak = checkbox_center
    xhw, yhw = checkbox_halfwidth    
    
    for ii in xrange(xpeak - xhw - 1, xpeak + xhw - 1):
        for jj in xrange(ypeak - yhw - 1, ypeak + yhw - 1):
            xloc = ii + 1
            yloc = jj + 1
            c_sum = c_sum + image[jj, ii]
            xsum = xsum + xloc * image[jj, ii]
            ysum = ysum + yloc * image[jj, ii]
            
    if debug:
        # Initial sum calculation (before iterations)
        print('(centroid_2D): Init. Sum (before iterations) = ', c_sum)

    if c_sum == 0:
        print('(centroid_2D): ERROR - divide by zero.')
        return
    else:
        xcen = xsum / c_sum
        ycen = ysum / c_sum
            
    # Iteratively calculate centroid until solution converges, using 
    # neighboring pixels to apply weighting...
    old_xcen = xcen
    old_ycen = ycen
    num_iter = 0
    
    for kk in xrange(max_iter):
        num_iter += 1
        c_sum = 0
        xsum = 0
        ysum = 0
        
        # Set up x and y centroid scanning ranges
        x_range = np.array((np.floor(old_xcen - xhw) - 1, np.ceil(old_xcen + xhw) - 1))
        y_range = np.array((np.floor(old_ycen - yhw) - 1, np.ceil(old_ycen + yhw) - 1))
                
        for ii in xrange(np.int(x_range[0]), np.int(x_range[1])):
            for jj in xrange(np.int(y_range[0]), np.int(y_range[1])):
                
                # Initalize weights to zero
                xweight = 0
                yweight = 0
                
                # Adjust weights given distance from current centroid
                xoff = np.abs((ii + 1) - old_xcen)
                yoff = np.abs((jj + 1) - old_ycen)
                
                # If within the original centroid box, set weight to 1
                # for both x and y.
                # If on the border, the scale weight
                if xoff <= xhw:
                    xweight = 1
                elif xoff > xhw and xoff < (xhw + 1):
                    xweight = xhw + 1 - xoff
                    
                if yoff <= yhw:
                    yweight = 1
                elif yoff > yhw and yoff < (yhw + 1):
                    yweight = yhw + 1 - yoff
                    
                # Compute cummulative weight
                weight = xweight * yweight
                
                # Calculate centroid
                xloc = ii + 1
                yloc = jj + 1
                
                c_sum = c_sum + image[jj, ii] * weight
                xsum = xsum + xloc * image[jj, ii] * weight
                ysum = ysum + yloc * image[jj, ii] * weight
                
        if c_sum == 0:
            print('(centroid_2D): ERROR - Divide by zero.')
        else:
            xcen = xsum / c_sum
            ycen = ysum / c_sum
                    
            # Check for convergence
            if np.abs(xcen - old_xcen) <= threshold and np.abs(ycen - old_ycen) <= threshold:
                convergence_flag = 'Success'
                break
            elif kk == max_iter:
                convergence_flag = 'Fail'
                break
            else:
                old_xcen = xcen
                old_ycen = ycen
    
    # Debug messages
    if debug:
        print('(centroid_2D): xpeak, ypeak = {}, {}'.format(xpeak, ypeak))
        print('(centroid_2D): xhw, yhw = {}, {}'.format(xhw, yhw))
        print('(centroid_2D): xcen, ycen = {}, {} '.format(xcen, ycen))        
        
                        
    print('(centroid_2D): Centroid = [{}, {}] for num_iter = {}.'.format(xcen-1, ycen-1, num_iter))
    print('(centroid_2D): Converged? ', convergence_flag)
          
    # -1 on both axes, as Python is 0 major    
    centroid = np.array((xcen-1, ycen-1))
    return centroid, c_sum
# *************************** centroid_2D ***************************



# *************************** centroid_1D ***************************
def centroid_1D(image, xpeak, xhw, debug=False):
    """
    Fine location of the target by calculating the centroid for 
    the region centered on the brightest checkbox.
    
    Performs the centroid calculation on the checkbox region 
    calculated using the function checkbox_1D().
    
    Keyword arguments:
    image -- 2 dimensional psf image
    xpeak -- The location of the brightest checkox in the 
             flattened vector.
    xhw   -- The halfwidth of the checkbox region calculated in
             checkbox_1D.
             
    
    Output(s):
    x_cen -- Target centroid location.
    c_cum -- The calculated flux sum within the checkbox 
             region.
    
    Example usage:
    
        >> cb_centroid, cb_sum = centroid_1D(psf, cb_center, 
                                 cb_hw)
        
        Find the vector centroid given the checkbox center and 
        halfwidth.
    """
    
    # Collapse input image unto x axis
    vector = np.sum(image, axis=0)
    
    c_sum = 0.0
    xcen = 0.0
        
    for ii in xrange(xpeak - xhw - 1, xpeak + xhw - 1):
        c_sum = c_sum + vector[ii]
        xloc = ii + 1
        xcen = xcen + xloc * vector[ii]
    
    print('(centroid_1D): Sum = ', c_sum)
    
    
    if c_sum == 0:
        print('(centroid_1D): ERROR - divide by zero')
    else:
        xcen = xcen / c_sum
        
    print('(centroid_1D): Centroid = ', xcen-1)
        
    # -1 on both axes, as Python is 0 major    
    return xcen-1, c_sum  
# *************************** centroid_1D ***************************


# *************************** find2D_higher_moments ***************************
def find2D_higher_moments(image, centroid, halfwidths, c_sum):
    """
    Calculate the higher moments of the object in the image.
    
    Find the normalized squared and cubed moments with reference 
    to an origin at the centroid, using the centroid and and sum 
    values calculatated previously.
    
    Keyword arguments:
    image      -- 2 dimensional psf image
    centroid   -- Tuple containing centroid of the object in 
                  image ([x, y] in pixels).
    halfwidths -- The halfwidths in both x and y of the 
                  checkbox from centroid_2D(). (pix)
    c_sum      -- The calculated flux sum within the checkbox 
                 region.
                 
    
    Output(s):
    x_mom, y_mom -- Tuples containing the x and y higher moments 
                    (normalized squared and cubed).
    
    Example usage:
    
        >> x_mom, y_mom = find2D_higher_moments(psf, cb_centroid, 
                          cb_hw, cb_sum)
        
        Find the second and third moments of psf, given the 
        centroid, checkbox halfwidths, and calculated sum.
    """
    
    # Unpack centroid to seperate values
    xcen, ycen = np.floor(centroid)
    xhw, yhw = halfwidths
    
    xmoment2 = 0
    xmoment3 = 0
    ymoment2 = 0
    ymoment3 = 0
        
    # Set up x and y centroid scanning ranges
    x_range = np.array((np.floor(xcen - xhw) - 1, np.ceil(xcen + xhw) - 1))
    y_range = np.array((np.floor(ycen - yhw) - 1, np.ceil(ycen + yhw) - 1))
    
    
    for ii in xrange(np.int(x_range[0]), np.int(x_range[1])):
        for jj in xrange(np.int(y_range[0]), np.int(y_range[1])):
            
            xloc = ii - np.floor(xcen)
            yloc = jj - np.floor(ycen)
            
            xweight = 0
            yweight = 0
            
            xoff = np.abs(ii - xcen)
            yoff = np.abs(jj - ycen)
            
            if xoff <= xhw:
                xweight = 1
            elif xoff > xhw and xoff < (xhw + 1):
                xweight = xhw + 1 - xoff
                
            if yoff <= yhw:
                yweight = 1
            elif yoff > yhw and yoff < (yhw + 1):
                yweight = yhw + 1 - yoff
                
            weight = xweight * yweight
            
            xmoment2 = xmoment2 + xloc**2 * image[jj, ii] * weight
            xmoment3 = xmoment3 + xloc**3 * image[jj, ii] * weight
            ymoment2 = ymoment2 + yloc**2 * image[jj, ii] * weight
            ymoment3 = ymoment3 + yloc**3 * image[jj, ii] * weight
            
    xmoment2 = xmoment2 / c_sum
    xmoment3 = xmoment3 / c_sum
    ymoment2 = ymoment2 / c_sum
    ymoment3 = ymoment3 / c_sum

    # Pack the x and y moments to return to main program
    x_moment = np.array((xmoment2, xmoment3))
    y_moment = np.array((ymoment2, ymoment3))
    
    return x_moment, y_moment   
# *************************** find2D_higher_moments ***************************

# *************************** find1D_higher_moments ***************************
def find1D_higher_moments(image, xcen, xhw, c_sum):
    """
    Calculate the higher moments of the object in the image.
    
    Find the normalized squared and cubed moments with reference 
    to an origin at the centroid, using the centroid and and sum 
    values calculatated previously.
    
    Keyword arguments:
    image -- 2 dimensional psf image
    x_cen -- Target cebtroid location.
    xhw   -- The checkbox halfwidth.
    c_sum -- The calculated flux sum within the checkbox region.
                 
    
    Output(s):
    x_mom -- Vector higher moments.
    
    Example usage:
    
        >> x_mom = find1D_higher_moments(psf, cb_centroid, cb_hw, 
                   cb_sum)
        
        Find the second and third moments of psf, given the 
        vector centroid, checkbox halfwidth, and calculated sum.
    """
    
    # Collapse input image unto x axis
    vector = np.sum(image, axis=0)
    
    xmoment2 = 0.0
    xmoment3 = 0.0
    
    # Set up x and y centroid scanning ranges
    x_range = np.array((np.floor(xcen - xhw) - 1, np.ceil(xcen + xhw) - 1))

    for ii in xrange(np.int(x_range[0]), np.int(x_range[1])): 
        xloc = (ii + 1) - np.floor(xcen)
        
        xweight = 0
        xoff = np.abs(ii - xcen)
        
        if xoff <= xhw:
            xweight = 0
        elif xoff > xhw and xoff < xhw + 1:
            xweight = xhw + 1 - xoff
            
        xmoment2 = xmoment2 + xloc**2 * vector[ii] * xweight
        xmoment3 = xmoment3 + xloc**3 * vector[ii] * xweight
        
    xmoment2 = xmoment2 / c_sum
    xmoment3 = xmoment3 / c_sum
    
    # Pack moments for return to main program
    x_mom = np.array((xmoment2, xmoment3))
    
    return x_mom
# *************************** find1D_higher_moments ***************************



# Print diagnostic load message
print("(jwst_targloc): JWST Target Locate Utilities Version {} loaded!".format(__version__))
