"""
Propidium Iodide Cell Death Assay Toolbox | David Cullen | Jones Lab 2024
"""

import cellpose
from cellpose import models
def segment_pi(image: np.ndarray,
                   progress_report: str,
                   user: str,
                   model_name: str = '240904_PI'
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
    print(progress_report + "|Segmenting PI...                                                                                                             ", end='\r')
    #model = models.CellposeModel(model_type='nuclei')
    #model = models.CellposeModel(model_type='cellpose7')
    model = models.CellposeModel(model_type=f'/scratch/user/{user}/cellpose/models/{model_name}')
    masks = model.eval(image[1], channels=[0,0])
    masks = clear_border(masks[0])      
    return masks

