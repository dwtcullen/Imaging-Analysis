"""
Kinase translocation reporter toolbox | David Cullen | Jones Lab 2024
"""
import numpy as np
def rgb(image: np.ndarray, 
        KTR_channel: str):
    image_height = image.shape[1]
    rgb = np.zeros(shape=(image_height, image_height, 3),dtype=np.uint16)
    
    KTR_channel_dictionary = {'Green': 1, 'Orange': 2}
    
    KTR_channel_number = KTR_channel_dictionary[KTR_channel]
    rgb[:,:,0] = image[0]
    rgb[:,:,1] = image[KTR_channel_number]
    return rgb

from cellpose import models
def segment_cytoplasm(rgb, progress_report, microscope_name, user):
    '''
    Segment cytoplasm in the input RGB image.
   
    Parameters:
    rgb (ndarray): RGB image with shape (rows, cols, channels)
   
    Returns:
    cyto_mask (ndarray): Binary image of cytoplasm segmentation with shape (rows, cols)
    '''

    print(progress_report + "|Segmenting perinuclear rings...                                                                                          ", end='\r')
    # Load the pre-trained model for cytoplasm segmentation
    #model = models.CellposeModel(model_type='ringpose_230402')
    model = models.CellposeModel(model_type=f'/scratch/user/{user}/cellpose/models/ringpose_230402')
    # Evaluate the model on the input RGB image
    if microscope_name == 'Deltavision':
        cyto_masks = model.eval(rgb, diameter=70, channels=[2,1]) 
    else:
        cyto_masks = model.eval(rgb, channels=[2,1])  
    return cyto_masks[0]
    print(progress_report + "|Segmented perinuclear rings!                                                                                          ", end='\r')

from scipy.spatial.distance import cdist
def get_closest_cyto_indices(nuclear_props, 
                             cyto_props):
   
    '''
    This function takes in two lists of region properties of the nuclei and cytoplasm and returns the indices of the closest cytoplasmic regions to each nucleus.

    Parameters:
    - nuclear_props (list): a list of `regionprops` objects, representing the properties of each nucleus in the image
    - cyto_props (list): a list of `regionprops` objects, representing the properties of each cytoplasmic region in the image

    Returns:
    - closest_cyto_indices (numpy.ndarray): an array of shape len(nuclear_props), where each element represents the index of the closest cytoplasmic region to the corresponding nucleus.
    '''
    if len(nuclear_props) == 0 or len(cyto_props) == 0:
        # Return an empty array if there are no cells present in the field of view
        return np.array([])
    # Get the centroid coordinates of each nucleus and cytoplasmic region
    nuclear_coords = [prop.centroid for prop in nuclear_props]
    cyto_coords = [prop.centroid for prop in cyto_props]
    # Calculate the pairwise distances between the nuclei and cytoplasmic regions
    distances = cdist(nuclear_coords, cyto_coords)
    # Get the indices of the closest cytoplasmic region to each nucleus
    closest_cyto_indices = np.argmin(distances, axis=1)

    return closest_cyto_indices

import numpy as np
import cv2
import tifffile as tf
def cytoplasmic_nuclear_ratio(
    nuclear_props: list,
    cyto_props: list,
    nuclear_masks: np.ndarray,
    cyto_masks: np.ndarray,
    intensity_image: np.ndarray,
    name: tuple[str, str, str],
    green_ratio: list,
    orange_ratio: list,
    green_ring_intensity: list,
    orange_ring_intensity: list,
    progress_report: str
    ) -> float:
   
    """
    Given a set of nuclear properties and cytoplasmic properties, this function calculates the ratio of
    the mean intensity of cytoplasmic regions around each nucleus to the mean intensity of the nuclei.

    Parameters
    ----------
    nuclear_props: list
        List of skimage.measure._regionprops._RegionProperties objects representing the properties of nuclei.
    cyto_props: list
        List of skimage.measure._regionprops._RegionProperties objects representing the properties of cytoplasmic regions.
    intensity_image: numpy.ndarray
        Numpy array representing the intensity image.
    cyto_masks: numpy.ndarray
        Binary numpy array representing the cytoplasmic masks.

    Returns
    -------
    dhb: float
        The ratio of cytoplasmic to nuclear mean intensity for each nucleus.
    """
    print(progress_report + "|Measuring cyto/nuc ratios...                                                                                               ", end='\r')
    closest_cyto_indices = get_closest_cyto_indices(nuclear_props, cyto_props)
    dilated_nuclei = cv2.dilate(nuclear_masks, np.ones((5,5), np.uint8))
    ring_masks = np.where(dilated_nuclei==0, cyto_masks, 0)
    for i, nuclear_prop in enumerate(nuclear_props):
        if len(nuclear_props) == 0 or len(cyto_props) == 0:
            green_ratio.append(0)
            orange_ratio.append(0)
            green_ring_intensity.append(0)
            orange_ring_intensity.append(0)
        else:
            closest_cyto_index = closest_cyto_indices[i]
            cyto_prop = cyto_props[closest_cyto_index]
            ring_indices = np.where(ring_masks == cyto_prop.label)
            
            green_image = intensity_image[1]
            orange_image = intensity_image[2]
           
        
            #green_image = intensity_image[:,:,1]
            #orange_image = intensity_image[:,:,2]
       
            green_ring_pixels = green_image[ring_indices]
            orange_ring_pixels = orange_image[ring_indices]

            green_ring_mean_intensity = np.mean(green_ring_pixels)
            orange_ring_mean_intensity = np.mean(orange_ring_pixels)
       
            green_ring_intensity.append(green_ring_mean_intensity)
       
            nuclear_intensity = nuclear_prop.intensity_mean
       
            current_green_ratio = green_ring_mean_intensity / nuclear_intensity[1]
            current_orange_ratio = orange_ring_mean_intensity / nuclear_intensity[2]
            #print(f"Nuclear: {nuclear_intensity}")
            #print(f"Ratio: {ratio}")
            green_ratio.append(current_green_ratio)
            orange_ratio.append(current_orange_ratio)
            orange_ring_intensity.append(orange_ring_mean_intensity)
            #print(f"Nucleus {nuclear_prop.label}, cyt {cyto_prop.label}, at {nuclear_prop.centroid}, {cyto_prop.centroid}")
        print(progress_report + "|Measured Cyto/Nuc ratio!                                                                                      ", end='\r')
            
            
            
def output_ktr_image(name: str, 
                 image: np.ndarray, 
                 masks: np.ndarray,
                 KTR_masks: np.ndarray
                ) -> None:
    '''
    Save example images with nuclear rings.

    Parameters:
        name: str, name of the output file
        meta: Tuple[str, str, str], metadata of the image
        image: np.ndarray, input multi-channel image
        masks: np.ndarray, input 2D masks image

    Returns:
        None, saves the modified image to a file in the Images folder
    '''
    # Compute the inner ring of the masks
    rings = inner_ring(masks, erosion = 1)
    KTR_rings = inner_ring(KTR_masks, erosion = 1)
    # Copy the input image
    output_image = image.copy()
    # Loop through each channel in the output image
    for channel in output_image:
        # Modify the current channel using np.putmask
        np.putmask(channel, rings, rings+10000)
        np.putmask(channel, KTR_rings, rings+10000)
    # Subtract 65536 from the image
    output_image = 65536 - output_image
    # Write the output image to a file
    tf.imwrite(f'Images/{name[0]}{name[1]}_{name[2]}.ome.tif', output_image)