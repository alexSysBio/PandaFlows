# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:34:06 2024

@author: alexpapa
"""


import flowio
import numpy as np
import pandas as pd
import os
import skimage
# import LoG_adaptive_image_filter as log
import matplotlib.pyplot as plt


class flow_cytometry_class(object):
    
    def __init__(self, fcs_path, images_folder):
        
        fd = flowio.FlowData(fcs_path)
        
        channels_list = []
        
        for ch in range(1, len(fd.channels)+1):
            
            ch_string = fd.channels[str(ch)]
            channels_list.append(ch_string[list(ch_string.keys())[0]])
        
        print(len(fd.events), 'events')
        print(len(fd.channels), 'channels:')
        print(channels_list)
        print('collecting metadata...')
        self.metadata = fd.text
        
        fd_array = np.reshape(fd.events, (-1, fd.channel_count))
        
        fc_df = pd.DataFrame()
        i = 0 
        for ch in channels_list:
            fc_df[ch] = list(fd_array[:,i])
            i+=1
        
        self.fcs_dataframe = fc_df
        self.channels_list = channels_list
        
        if images_folder != 'none':
            images_list = os.listdir(images_folder)
            print('loading images')
            image_arrays = {}
            for img in images_list:
                img_string = img[:img.find('.tif')]
                # print(img_string, img)
                image_arrays[int(img_string)] = skimage.io.imread(images_folder+'/'+img)
            
            self.image_arrays = image_arrays
            print(len(list(image_arrays.keys())), 'images in total')
        else:
            print(f"no image folder selected: {images_folder}")
    
    def get_image_arrays(self):
        return self.image_arrays
    
    def get_flow_cytometry_dataframe(self):
        return self.fcs_dataframe
    
    def get_metadata(self):
        self.metadata
        

    def segment_cell_images(self, 
                            hard_threshold=775, min_area=35, check_segmentation=False):
        
        # log_params = [1, 1000, 95, 1, 19, -10, 0], 
        # area_threshold = 50,
        # rounds_of_erosion = 2
        # resized_img = skimage.transform.resize(img, (img.shape[0]*resize_factor, img.shape[1]*resize_factor))
        # log_image = log.log_adaptive_filter(img, log_params)[1]
        # log_image = skimage.morphology.remove_small_objects(log_image, area_threshold)
        # log_image = skimage.morphology.convex_hull_object(log_image, connectivity=2)

        cell_masks = {}
        cell_areas = {}
        
        for img_key in self.image_arrays:
            
            img = self.image_arrays[img_key]
            masked_image = img < hard_threshold
            masked_image = skimage.morphology.remove_small_objects(masked_image, min_area)
            masked_image = skimage.morphology.binary_dilation(masked_image)
                
            image_labels = skimage.measure.label(masked_image)
            
            if check_segmentation==True:
                plt.imshow(img)
                plt.colorbar()
                plt.show()
                plt.imshow(masked_image)
                plt.show()
                input('press enter to proceed')
            
            if np.max(image_labels)==1:
                cell_masks[img_key] = masked_image
                cell_areas[img_key] = masked_image[np.nonzero(masked_image)].shape[0]
            elif np.max(image_labels)==0:
                print(img_key, 'no cell segmented')
            elif np.max(image_labels)>1:
                print(img_key, 'more than one objects segmented')
            
        self.cell_masks = cell_masks
        self.cell_areas = cell_areas
        
    
    def get_cell_masks(self):
        return self.cell_masks
    
    def get_cell_areas(self):
        return self.cell_areas
    
    


