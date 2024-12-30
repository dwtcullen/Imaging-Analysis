#My python files
import image_processing
import ktr
import timekeeping
import image_output
import file_management
import measure_nuclei
#import pi_cell_death

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
        image = image_processing.subtract_background(image, progress_report, microscope_name, CPU_or_GPU, has_KTR)
        nuclear_masks = image_processing.segment_nuclei(image, progress_report, microscope_name, user, model_name)
        if has_KTR:
            KTR_masks = ktr.segment_cytoplasm(ktr.rgb(image, KTR_channel), progress_report, microscope_name, user)
        if output_df:
            nuclear_props = measure_nuclei.get_regionprops(nuclear_masks, image)
            if has_KTR:    
                cyto_props = measure_nuclei.get_regionprops(KTR_masks, image)
                ktr.cytoplasmic_nuclear_ratio(nuclear_props, cyto_props, nuclear_masks, KTR_masks, image, name, green_ratio, orange_ratio, green_ring_intensity, orange_ring_intensity, progress_report)
            measure_nuclei.measure_properties(name, image, nuclear_props, results, progress_report, nuclear_masks, calculate_pearson) 
        image_output.manage_image_output(name, image, nuclear_masks, progress_report, string, KTR_masks, has_KTR, output_masks, output_image_stacks)
    if output_df:
        file_management.write_df(results, image, output_df_name, green_ratio, orange_ratio, green_ring_intensity, orange_ring_intensity, has_KTR)
    else:
        print(f"Image processing completed.                                                                                                               ", end="\r")