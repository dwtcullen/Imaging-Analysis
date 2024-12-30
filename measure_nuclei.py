"""
Nuclei measurement toolbox | David Cullen | Jones Lab 2024
"""

import numpy as np
import image_output
def intensity_sum(regionmask: np.ndarray, 
                  intensity_image: np.ndarray
                 ) -> float:
    '''
    Return the total intensity for all pixels inside the mask.

    Parameters
    ----------
    regionmask : np.ndarray
        A binary mask indicating which pixels belong to the region of interest.
    intensity_image : np.ndarray
        The 2D intensity image corresponding to the regionmask.

    Returns:
    float
        The sum of intensities for all pixels inside the regionmask.
    '''
    return np.sum(intensity_image[regionmask])

def intensity_mean_eroded(regionmask: np.ndarray, 
                          intensity_image: np.ndarray
                         ) -> float:
    '''
    Erode the regionmask, then find the mean intensity of each channel within.
    
    Parameters:
    - regionmask (np.ndarray): The binary mask representing the region of interest.
    - intensity_image (np.ndarray): The multi-channel intensity image.
    
    Returns:
    - float: The mean intensity of the eroded regionmask.
    '''
    eroded = erode(regionmask, erosion = 1)
    eroded_image = intensity_image[eroded]
    return np.mean(eroded_image)

def measure_dist_from_centre(centroid: tuple[int, int],
                             image_centre: float):
    '''
    Basic trig to calculate the distance from the centre of the nucleus to the centre of the image.
    
    Parameters:
    - centroid: Tuple of two elements representing the row and column in which the centre of the nucleus is located in the image.
    -image: ndarray of the current image.
    
    Returns:
    - dist_from_centre: a float representing the distance from the centre of the nucleus to the centre of the image.
    '''
    prop_row, prop_col = centroid
    rows_from_centre = np.abs(prop_row-image_centre)
    cols_from_centre = np.abs(prop_col-image_centre)
    dist_from_centre = np.sqrt(rows_from_centre**2 + cols_from_centre**2)
    return dist_from_centre
    
def get_regionprops(masks, image):
    from skimage.measure import regionprops
    regionprops_image = np.moveaxis(image, 0, -1)
    props = regionprops(masks, intensity_image=regionprops_image, extra_properties=(intensity_sum,))
    return props

def measure_properties(name: tuple[str, str, int], 
                       image: np.ndarray, 
                       regions: list, 
                       results: list,
                       progress_report: str,
                       nuclear_masks,
                       calculate_pearson: bool = False
                      ) -> list[tuple]:
    import matplotlib.pyplot as plt
    
    '''
    Measure properties of nuclei in an image.
    
    Parameters:
    - name: Tuple of three elements representing the well, row, and column of the current image
    - image: numpy ndarray representing the image to be processed
    - regions: list of regionprops objects representing the nuclei in the image
    
    Returns:
    - results: a list of lists where each sublist represents the properties of a single nucleus
    
    The properties included in each sublist are:
    - well: well identifier, a concatenation of the row and column in the form "RowColumn"
    - field: identifier of the field within the well
    - label: identifier of the nucleus within the field
    - area: area of the nucleus, in pixels
    - circularity: measure of circularity, defined as 4 * pi * area / perimeter^2
    - Mean_DAPI: the mean log10 intensity of the nucleus in the DAPI channel
    - Mean_DAPI: the mean log10 intensity of the nucleus in the DAPI channel
    - Mean_Green (optional): the mean log10 intensity of the nucleus in the Green channel, if present
    - Mean_Orange (optional): the mean log10 intensity of the nucleus in the Orange channel, if present
    - Mean_EdU (optional): the mean log10 intensity of the nucleus in the EdU channel, if present
    - Sum_DAPI: the sum log10 intensity of the nucleus in the DAPI channel
    - Sum_Green (optional): the sum log10 intensity of the nucleus in the Green channel, if present
    - Sum_Orange (optional): the sum log10 intensity of the nucleus in the Orange channel, if present
    - Sum_EdU (optional): the sum log10 intensity of the nucleus in the EdU channel, if present
    '''
   
    print(progress_report + f"|Measuring nuclear properties...                                                                                               ", end='\r')
    for props in regions:
        #if np.log10(props.intensity_mean[1])>2.25:
        #    image_output.output_cropped_images(props, image, name)
        circularity = 4 * np.pi * props.area / (props.perimeter**2)
        distance_from_centre = measure_dist_from_centre(props.centroid, image.shape[1]/2)
        row, column, field = name
        result_dict = {
            'Column': column,
            'Row': row
        }
        live_status = 'live'
        
        result = [column,
                  row,
                  row+column,
                  field, 
                  props.label,
                  props.area,
                  circularity,
                  props.centroid,
                  props.centroid[0],
                  props.centroid[1],
                  distance_from_centre, live_status] + [np.log10(i) for i in props.intensity_mean] + [np.log10(i) for i in props.intensity_sum]
        
        #if calculate_pearson:
        #    pearson = pearson(props)
        #    result.append(pearson)
            
        results.append(result)
        
        
def top_5percent_brightest_pixels_mask(image):
    # Get the flattened pixel values
    pixels = image.ravel()
    threshold = np.sort(pixels)[-int(len(pixels)*0.1)]
    result = np.copy(image)
    # Set all pixels below the threshold to zero. Optional: Set all pixels above the threshold to the maximum value, 65535.
    result[result < threshold] = 0
    result[result >= threshold] = 65535
    return result



from itertools import combinations
#def pearson(props):
#import scipy.stats as st
#            dapi_green_pearson, _ = st.pearsonr(top_5percent_brightest_pixels_mask(props.image_intensity[...,0]).ravel(), top_5percent_brightest_pixels_mask(props.image_intensity[...,1]).ravel())
#            dapi_orange_pearson, _ = st.pearsonr(top_5percent_brightest_pixels_mask(props.image_intensity[...,0]).ravel(), top_5percent_brightest_pixels_mask(props.image_intensity[...,2]).ravel())
#            dapi_edu_pearson, _ = st.pearsonr(top_5percent_brightest_pixels_mask(props.image_intensity[...,0]).ravel(), top_5percent_brightest_pixels_mask(props.image_intensity[...,3]).ravel())
#        
#            green_edu_pearson, _ = st.pearsonr(top_5percent_brightest_pixels_mask(props.image_intensity[...,1]).ravel(), top_5percent_brightest_pixels_mask(props.image_intensity[...,3]).ravel())
#            orange_edu_pearson, _ = st.pearsonr(top_5percent_brightest_pixels_mask(props.image_intensity[...,2]).ravel(), top_5percent_brightest_pixels_mask(props.image_intensity[...,3]).ravel())
#            green_orange_pearson, _ = st.pearsonr(top_5percent_brightest_pixels_mask(props.image_intensity[...,1]).ravel(), #top_5percent_brightest_pixels_mask(props.image_intensity[...,2]).ravel())
#        
#            dapi_green_pearson_all, _ = st.pearsonr(props.image_intensity[...,0].ravel(),props.image_intensity[...,1].ravel())
#            dapi_orange_pearson_all, _ = st.pearsonr(props.image_intensity[...,0].ravel(),props.image_intensity[...,2].ravel())
#            dapi_edu_pearson_all, _ = st.pearsonr(props.image_intensity[...,0].ravel(),props.image_intensity[...,3].ravel())
#            
#            green_edu_pearson_all, _ = st.pearsonr(props.image_intensity[...,1].ravel(),props.image_intensity[...,3].ravel())
#            orange_edu_pearson_all, _ = st.pearsonr(props.image_intensity[...,2].ravel(),props.image_intensity[...,3].ravel())
#            green_orange_pearson_all, _ = st.pearsonr(props.image_intensity[...,1].ravel(),props.image_intensity[...,2].ravel())
#
#            result = [column,
#                      column_formatted,
#                      row,
#                      row+column_formatted,
#                      row+column, 
#                      name[2], 
#                      props.label, 
#                      props.area, 
#                      circularity, 
#                      prop_row, 
#                      prop_col, 
#                      dist_from_centre,
#                      dapi_green_pearson,
#                      dapi_orange_pearson,
#                      dapi_edu_pearson,
#                      green_edu_pearson, 
#                      orange_edu_pearson, 
#                      green_orange_pearson, 
#                      dapi_green_pearson_all,
#                      dapi_orange_pearson_all,
#                      dapi_edu_pearson_all,
#                      green_edu_pearson_all, 
#                      orange_edu_pearson_all, 
#                      green_orange_pearson_all] + [np.log10(i) for i in props.intensity_mean] + [np.log10(i) for i in props.intensity_sum]