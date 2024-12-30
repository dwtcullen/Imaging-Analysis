"""
Image output toolbox | David Cullen | Jones Lab 2024
"""
import numpy as np
import image_processing
import tifffile as tf

def manage_image_output(name: tuple[str, str, int],
                        image: np.ndarray,
                        nuclear_masks: np.ndarray,
                        progress_report: str,
                        string: str,
                        KTR_masks: np.ndarray,
                        has_KTR: bool=False, 
                        output_masks: bool=False, 
                        output_image_stacks: bool=False) -> None:
    '''
    This function manages the decision whether to output images.
    '''
    if output_masks:
        output_mask(name, nuclear_masks, KTR_masks, progress_report, string, has_KTR)
    if output_image_stacks:
        output_image_stack(name, image, nuclear_masks, KTR_masks, progress_report, has_KTR)

def output_mask(name,
                nuclear_masks: np.ndarray, 
                KTR_masks: np.ndarray,
                progress_report: str,
                string: str,
                has_KTR:bool=False
               ) -> None:
    '''
    Save nucear and KTR masks if requested.
    '''
    print(progress_report + "|Saving nuclear masks...                                                                                                             ", end='\r')
    input_string = string.replace('405', 'cpm')
    tf.imwrite(f'Output_Images/Nuclear_Masks/{input_string}', nuclear_masks) 
    #tf.imwrite(f'Output_Images/Nuclear_Masks/
    if has_KTR:
        print(progress_report + "|Saving KTR masks...                                                                                                             ", end='\r')
        tf.imwrite(f'Output_Images/KTR_Masks/{name[0]}{name[1]}_{name[2]}_KTR_masks', KTR_masks)

def output_image_stack(name: str, 
                       image: np.ndarray, 
                       nuclear_masks: np.ndarray,
                       KTR_masks: np.ndarray,
                       progress_report: str, 
                       has_KTR: bool=False
                      ) -> None:
    '''
    Save example images with nuclear and perinuclear rings for image sets containing KTR images.

    Parameters:
        name: str, name of the output file
        meta: Tuple[str, str, str], metadata of the image
        image: np.ndarray, input multi-channel image
        masks: np.ndarray, input 2D masks image

    Returns:
        None, saves the modified image to a file in the Images folder
    '''
    print(progress_report + "|Saving image stacks...                                                                                                             ", end='\r')
    rings = image_processing.inner_ring(nuclear_masks, erosion = 1)
    output_image = image.copy()
    for channel in output_image:
        np.putmask(channel, rings, rings+10000)
        if has_KTR:
            KTR_rings = image_processing.inner_ring(KTR_masks, erosion = 1)
            np.putmask(channel, KTR_rings, rings+10000)
    output_image = 65536 - output_image
    tf.imwrite(f'Output_Images/Stacks/{name[0]}{name[1]}_{name[2]}.ome.tif', output_image)
    
def crop_image(image: np.ndarray, 
               centroid: tuple[int, int], 
               crop_size: int = 220
              ) -> np.ndarray:
    '''
    Crop an image around a given centroid to a specified size. 
    This is used to generate equally-sized cropped images of nuclei.

    Parameters:
    image (numpy.ndarray): The input image to be cropped. 
                           The shape should be (rows, cols, channels).
    centroid (tuple): The centroid (row, col) around which to crop the image.
    crop_size (int, optional): The size (height and width) of the cropped image. 
                               Default to 40 for 10x, 250 for 60x.

    Returns:
    numpy.ndarray: The cropped and padded image with shape (crop_size, crop_size, channels).
    '''
    # Get the shape of the image
    channels, rows, cols = image.shape

    # Calculate the start and end row and column indices for the crop
    start_row = max(0, int(centroid[0] - crop_size / 2))
    end_row = min(rows, int(centroid[0] + crop_size / 2))
    start_col = max(0, int(centroid[1] - crop_size / 2))
    end_col = min(cols, int(centroid[1] + crop_size / 2))
    
    # Crop the image
    cropped_image = image[start_row:end_row, start_col:end_col, :]
    
    # Calculate the padding for each side if the cropped image is not the correct size
    top_padding = bottom_padding = int((crop_size - cropped_image.shape[0]) / 2)
    left_padding = right_padding = int((crop_size - cropped_image.shape[1]) / 2)
    if cropped_image.shape[0] < crop_size:
        bottom_padding = bottom_padding + (crop_size - cropped_image.shape[0]) % 2
    if cropped_image.shape[1] < crop_size:
        right_padding = right_padding + (crop_size - cropped_image.shape[1]) % 2
    
    # Pad the cropped image to the required size
    padded_image = np.pad(cropped_image, ((top_padding, bottom_padding), (left_padding, right_padding), (0, 0)), 'constant')
    
    return padded_image

def output_cropped_images(props, 
                        image: np.ndarray, 
                        name: str, 
                        crop_size: int = 45
                       ):
    """
    This function takes in an image, region properties of cells in the image, the name of the image, a threshold for DAPI intensity, a threshold for EdU intensity, and a crop size, and saves the cropped images of cells that meet the intensity criteria.

    Parameters:
    - image (np.ndarray): a 3D numpy array of shape (height, width, channel), representing the input image
    - name (str): a string representing the name of the image
    - crop_size (int): an integer representing the size of the crop (default is 220)

    Returns:
    - None. Outputs image to folder.
    """
    crop = crop_image(image, props.centroid, crop_size)
    #crop = np.moveaxis(crop, -1, 0)
    #contrasted = np.array([auto_contrast(i) for i in crop])
    #vert = np.vstack((contrasted[0], contrasted[1], contrasted[3]))
        #vert_image = np.vstack((cropped_image[...,0],cropped_image[...,1],cropped_image[...,3]))
        # Save the cropped image
    tf.imwrite(f'Output_Images/{name[0]}{name[1]}_{str(props.label)}_stack.ome.tif', crop)
   
def output_rgb(image: np.ndarray, 
               name: list[str]
              ) -> None:
    '''
    Save the image in RGB format.

    Parameters:
    image: a numpy array with shape (3, 2040, 2040)
        The input image.
    name: a list of 3 strings
        The name of the image.

    Returns:
    None. The outpput image will be saved in rgb folder.
    '''
    # Initialize the RGB image with zeros
    rgb = np.zeros(shape=(2040, 2040, 3), dtype=np.uint16)
    # Set the red channel of the RGB image
    rgb[:, :, 0] = image[0]
    # Set the green channel of the RGB image
    rgb[:, :, 1] = image[1]
    # Save the RGB image to disk
    tf.imsave('Output_Images/' + 'Well_' + name[0]+name[1] + '_Field_' + name[2] + '.tif', rgb)
    
import cv2
import numpy as np
def downsample(image, scale_percent):
    """
    Downsamples an image by a given percentage.

    Parameters:
    - image: Np array.
    - scale_percent (int): percentage to scale the image by.

    Returns:
    - downsampled_image (numpy.ndarray): the downsampled image as a NumPy array.
    """
    # Calculate the new dimensions
    width = int(image.shape[1] * scale_percent / 100)
    height = int(image.shape[0] * scale_percent / 100)
    # Resize the image
    dim = (width, height)
    downsampled_image = cv2.resize(image, dim, interpolation=cv2.INTER_AREA)
    return downsampled_image

import image_processing
def auto_contrast(image, saturation=0.005):
    """
    Adjusts the contrast of an image to maximize the contrast within the image while preserving the relative intensities
    of the features present in the image.

    Parameters:
        image (numpy.ndarray): A 2D or 3D numpy array representing the image.
        saturation (float): The fraction of the lowest and highest pixel intensities to saturate. Default is 0.005.

    Returns:
        numpy.ndarray: A 2D or 3D numpy array with the contrast adjusted.
    """
    # Determine the minimum and maximum pixel intensities
    pmin = np.percentile(image, 100 * saturation)
    pmax = np.percentile(image, 100 * (1 - saturation))
    # Adjust the contrast of the image
    output = np.clip((image - pmin) / (pmax - pmin), 0, 1)
    # Convert the output back to the original data type
    output = output * np.iinfo(image.dtype).max
    output = output.astype(image.dtype)
    return output

def save_example(path: str, 
                 example_well: str,
                 is_4i: bool = False,
                 is_live: bool = False
                ):
    import os
    folder = 'Example'
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    file_list = [entry.name for entry in os.scandir(path) if entry.is_file()]
    metadata = [{'String': file_name, **file_management.extract_file_info(file_name, is_live=is_live, is_4i = is_4i)} for file_name in file_list]
    metadata = pd.DataFrame(metadata)
    metadata = metadata[metadata['Field'] == '5']
    metadata = metadata.groupby(['Row', 'Column'])
    for i, (name, group) in enumerate(metadata):
        well = str(name[0]) + str(name[1])
        if well in example_well:
            array = 65535 - np.array([image_processing.auto_contrast(tf.imread(str(path) + "/" + each_string)) for each_string in group.String])
            altered = []
            for example in array:
                example = np.pad(example, pad_width=3, mode='constant', constant_values=0)
                example = np.pad(example, ((170, 170), (10, 10)), mode='constant', constant_values=65535)
                altered.append(example)
            array1 = np.array(altered)
            h_stack = np.hstack(array1)
            image_processing.add_example_text(h_stack, name)
            cv2.imwrite(f'{folder}/Well_{name[0]}{name[1]}.png', image_processing.downsample(h_stack,50))       

