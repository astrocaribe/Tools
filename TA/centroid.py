# Calculate the flux-weighted centroid of a NxN pixel image.
# Input: NxN pixel image grid
# Output: Co-ordinates of centroid
def centroid(grid):
    
    #Diagnostic output
    print "Grid length:", len(grid)
    print "Total:", np.sum(grid)
    print "Average:", np.average(grid)

    # Weights and weighted averages
    weight_i = np.sum(grid, 0)
    weight_j = np.sum(grid, 1)
    weight_avg_i = np.sum(grid, 0) * np.arange(len(grid))
    weight_avg_j = np.sum(grid, 1) * np.arange(len(grid))
        
    # X and Y position of the centroid    
    x = np.sum(weight_avg_i)/np.sum(weight_i)
    y = np.sum(weight_avg_j)/np.sum(weight_j)
    centroid = (x.astype(np.float16), y.astype(np.float16))

    return centroid