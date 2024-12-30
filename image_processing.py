"""
Image processing toolbox | David Cullen | Jones Lab 2024
"""
import cv2
import tifffile as tf
import numpy as np
def assemble_image(
    group,
    path: str,
    crop: bool = False
    ) -> np.ndarray:
    '''
    Import each file in a group and assemble a multi-channel image.

    Parameters:
        group (GroupType): An object of GroupType, containing a list of strings.
        path (str): A file path to the directory where the image files are located.

    Returns:
        np.ndarray: A multi-channel image as a numpy array.
    '''
    arrays = [tf.imread(str(path) + "/" + each_string) for each_string in group.String]
    if len(arrays) == 3:
        arrays = [arrays[0], arrays[1], np.zeros_like(arrays[0]), arrays[-1]]
    elif len(arrays) == 2:
        arrays = [arrays[0], np.zeros_like(arrays[0]), np.zeros_like(arrays[0]), arrays[-1]]
        
    array = np.array(arrays)   
    if crop:
        crop_size = 400
        array = array[:, crop_size:-crop_size, crop_size:-crop_size]
    return array

def subtract_background(
    image: np.ndarray,
    progress_report: str,
    microscope_name: str,
    CPU_or_GPU: str,
    has_KTR: bool=False):
    import pyclesperanto_prototype as cle
    from skimage import restoration
    '''
    For RPE1 cells, the required radius happens to be 2x the magnification. 20pixels for 10x, 40pixels for 20x.
    
    Parameters:
        image (np.ndarray): A single channel image.
        progress_report (str): A single channel image.
       microscope_name (str): The name of the microscope.
        CPU_or_GPU (str): A string 'CPU' or 'GPU' that represents whether a GPU is being used in this analysis.

    Returns:
        np.ndarray: A binary mask of the segmented nuclei.
    '''
    print(progress_report + "|Subtracting background...                                                                                                             ", end='\r')
    radius_dictionary = {
        'Vader': 20,
        'Maul': 20,
        'Deltavision': 40
    }
    radius = radius_dictionary[microscope_name]
    if CPU_or_GPU == 'CPU':
        if has_KTR:
            image[0] = image[0] - restoration.rolling_ball(image[0], radius=radius)
            image[-1] = image[-1] - restoration.rolling_ball(image[-1], radius=radius)
        else:
            image = np.array([i - restoration.rolling_ball(i, radius=radius) for i in image])
            
    if CPU_or_GPU == 'GPU':
        if has_KTR:
            image[0] = cle.top_hat_sphere(image[0], radius_x=radius, radius_y=radius)
            image[-1] = cle.top_hat_sphere(image[-1], radius_x=radius, radius_y=radius)
        else:
            image = np.array([cle.top_hat_sphere(i, radius_x=radius, radius_y=radius) for i in image])
    print(progress_report + "|Background subtracked.                                                                                                             ", end='\r')
    return image

import cellpose
from cellpose import models
def segment_nuclei(image: np.ndarray,
                   progress_report: str,
                   microscope_name: str,
                   user: str,
                   model_name: str = 'cellpose7'
                  ) -> np.ndarray:
    '''
    Use pre-trained model to segment nuclei in an image. 
    This is the default function, for 10x air.

    Parameters:
        image (np.ndarray): A single channel image.

    Returns:
        np.ndarray: A binary mask of the segmented nuclei.
    '''

    from skimage.segmentation import clear_border
    print(progress_report + "|Segmenting nuclei...                                                                                                             ", end='\r')
    #model = models.CellposeModel(model_type='nuclei')
    #model = models.CellposeModel(model_type='cellpose7')
    model_dict = {
        'cellpose7': models.CellposeModel(model_type=f'/scratch/user/{user}/cellpose/models/cellpose7'),
        'nuclei': models.CellposeModel(model_type='nuclei')
    }
    model = model_dict[model_name]
    if microscope_name == 'Deltavision':
        masks = model.eval(image[0], diameter=60, channels=[0,0])
    else:
        masks = model.eval(image[0], channels=[0,0])
    masks = clear_border(masks[0])      
    return masks

def erode(masks: np.ndarray, 
          erosion: int
         ) -> np.ndarray:
    '''
    Erode the mask.

    Parameters:
        masks (np.ndarray): The masks to erode.
        erosion (int): The size of the erosion to perform.

    Returns:
        np.ndarray: The eroded masks.
    '''
    # Create a square kernel with size (erosion*2+1, erosion*2+1) and data type np.uint8
    kernel = np.ones((erosion*2+1, erosion*2+1), np.uint8)
    # Convert masks to np.ndarray with data type np.uint16
    masks = np.asarray(masks, dtype="uint16")
    # Erode the masks with the kernel using OpenCV's cv2.erode function
    return cv2.erode(masks, kernel)

def inner_ring(masks: np.ndarray, 
               erosion: int
              ) -> np.ndarray:
    '''
    Create a ring on the inside border of the mask.

    Parameters:
        masks (np.ndarray): The masks to create a ring on.
        erosion (int): The size of the erosion to perform.

    Returns:
        np.ndarray: The inside ring of the masks.
    '''
    # Erode the masks
    eroded = erode(masks, erosion)
    # Return the result of original masks minus the eroded masks.
    return masks - eroded

from cv2 import dilate
def outer_ring(masks: np.ndarray, 
               dilation: int, 
               gap: int
              ) -> np.ndarray:
    ''' 
    Create a ring outside the border of the mask. 
    
    Parameters:
        masks: np.ndarray, the binary image masks.
        dilation: int, the size of the dilation.
        gap: int, the size of the gap.

    Returns:
        np.ndarray, the image with the outer ring.
        
    ''' 
    # Define the dilation kernel
    kernel = np.ones((dilation*2+1, dilation*2+1), np.uint8)
    # Define the gap kernel
    gap_kernel = np.ones((gap*2+1, gap*2+1), np.uint8)
    # Convert the masks to a uint8 array
    masks = np.asarray(masks, dtype="uint8")
    masks[masks>1] = 1100
    # Perform a slight dilation on the masks using the dilation kernel
    slight_dilation = cv2.dilate(masks, kernel)
    # Perform a dilation on the masks using the gap kernel
    dilated = cv2.dilate(masks, gap_kernel)
    # Return the difference between the dilated and slight_dilation masks
    return dilated - slight_dilation


    





        
#def decision_tree(name, image, KTR_channel, has_KTR, output_df):
#    
#    decision = {
#        #1 TTTT
#        (True, True, True, True): image_output_TTTT,
#        #2 TTFT
#        (True, True, False, True): image_output_TTFT,
#        #3 TFTT
#        (True, False, True, True): image_output_TFTT,
#        #4 TFFT
#        (True, False, False, True): image_output_TFFT,
#        #5 TTTF
#        (True, True, True, False): image_output_TTTF,
#        #6 FTFT
#        (False, True, False, True): image_output_FTFT,
#        #7 FTFF
#        (False, True, False, False): image_output_FTFF,
#        #8 FFFT
#        (False, False, False, True): image_output_FFFT,
#        #9 FFFF
#        (False, False, False, False): image_output_FFFF
#    }