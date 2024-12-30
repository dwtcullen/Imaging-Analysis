#My python files
import image_processing
import ktr
import timekeeping
import image_output
import file_management
import measure_nuclei

#Necessary packages
import numpy as np

def process_images(
    root_directory: str,
    image_folder: str, 
    output_df_name: str,
    CPU_or_GPU: str,
    microscope_name: str,
    user: str,
    KTR_channel: str,
    model_name: str,
    has_KTR = False,
    output_masks = False,
    output_image_stacks = False,
    output_df = False,
    crop: bool = False,
    calculate_pearson = False):
    
    ''' Call each function to process IF images in folder. '''
    global results
    results = []
    green_ratio = []
    orange_ratio = []
    green_ring_intensity = []
    orange_ring_intensity = []
    KTR_masks = []
    
    file_management.set_directory(root_directory, user)
    metadata = file_management.get_metadata(image_folder, microscope_name, groups = ['Row', 'Column', 'Field'])
    total_groups = len(metadata)
    start_time = timekeeping.current_time()
    file_management.folders(output_masks, has_KTR, output_image_stacks)
    for i, (name, group) in enumerate(metadata):
        string = file_management.get_DAPI_string(group)
        progress_report = timekeeping.progress(i, total_groups, name, start_time)
        image = image_processing.assemble_image(group, image_folder, crop)
        image_output.output_rgb(image, name)
        
