# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:34:06 2024

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, HHMI at Stanford University
"""


import flowio
import numpy as np
import pandas as pd
import os
import skimage
import matplotlib.pyplot as plt
import span_histogram_selection as span
import polygon_selection as poly
import pickle_read_save as prs


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
            
        self.channels_list = channels_list
        
        fc_df = fc_df[fc_df>0]
        
        try:
            gate_dict = prs.load_data('stored_gates')
            self.gate_dict = gate_dict
        except FileNotFoundError:
            self.gate_dict = {}

        # fcs_df = ga.apply_gates(dataframe, gate_dict)
        self.fcs_dataframe = fc_df
       
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
        
    
    def reset_gates(self):
        val = input('Warning, this operation will erase all stored gates and remove them from the dataframe: YES/NO: ')
        val = val.upper()
        if val == 'YES':
            fcs_df = self.get_flow_cytometry_dataframe()
            fcs_df = fcs_df.drop(list(self.gate_dict.keys()), axis=1)
            
            self.gate_dict = {}
            prs.save_data(self.gate_dict, 'stored_gates')
            
            self.fcs_dataframe  = fcs_df
            
        elif val == 'NO':
            print(self.get_gates())
            print('Gates have not not been reset.')
            pass

        else:
            self.reset_gates()
    
    def get_gates(self):
        return self.gate_dict
    
    
    def segment_cell_images(self, 
                            hard_threshold=775, min_area=35, check_segmentation=False):

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

    
    
    def histogram_gate(self, variable, log, bin_array):
        
        fcs_df = self.get_flow_cytometry_dataframe()
        gate_dict = self.get_gates()
        
        if log == True:
            # gated_df = gated_df[gated_df[variable]>0]
            x = np.log10(fcs_df[variable])
            log_string = '_log'
        else:
            x = fcs_df[variable]
            log_string = '_lin'
        
        i = 0
        n = 1
        while i==0:
            gate_name = variable+log_string+'_gate_'+str(n)
            if gate_name not in self.gate_dict:
                i+=1
            else:
                n+=1
        
        span.return_selected_ranges(x, bin_array, variable, gate_name)
        
        loaded_gate = prs.load_data(gate_name)
        gate_dict[gate_name] = loaded_gate
        prs.save_data(gate_dict, 'stored_gates')
        self.gate_dict = gate_dict
        os.remove(gate_name)
        print(self.get_gates())
        
        gate_range = gate_dict[gate_name]
        fcs_df[gate_name] = 0
        fcs_df[gate_name] = np.where(x.between(gate_range[0], gate_range[1]), 1, fcs_df[gate_name])
        self.fcs_dataframe = fcs_df
        

def scatter_gate(self, x_variable, y_variable, x_log, y_log):
        
        fcs_df = self.get_flow_cytometry_dataframe()
        gate_dict = self.get_gates()
        
        if x_log == True:
            # gated_df = gated_df[gated_df[variable]>0]
            x = np.log10(fcs_df[x_variable])
            xlog_string = '_log'
        else:
            x = fcs_df[x_variable]
            xlog_string = '_lin'

        if y_log == True:
            # gated_df = gated_df[gated_df[variable]>0]
            y = np.log10(fcs_df[y_variable])
            ylog_string = 'log'
        else:
            y = fcs_df[y_variable]
            ylog_string = 'lin'
        
        log_string = xlog_string+ylog_string
        
        i = 0
        n = 1
        while i==0:
            gate_name = x_variable+'_'+y_variable+log_string+'_gate_'+str(n)
            if gate_name not in self.gate_dict:
                i+=1
            else:
                n+=1
        
        poly.return_markers_in_polygon(x, y, x_variable, y_variable, gate_name, colormap='rainbow')
        
        loaded_gate = prs.load_data(gate_name)
        gate_dict[gate_name] = loaded_gate
        prs.save_data(gate_dict, 'stored_gates')
        self.gate_dict = gate_dict
        os.remove(gate_name)
        print(self.get_gates())
        
        polygon_coordinates = gate_dict[gate_name][0]
        polygon = poly.get_polygon_from_coordinates(polygon_coordinates)
        
        gate_df = pd.DataFrame()
        gate_df['x'] = x
        gate_df['y'] = y
        gate_df['g'] = gate_df.apply(lambda x: poly.get_markers_inside_gate(x.x,x.y, polygon), axis=1)
        
        fcs_df[gate_name] = gate_df.g
        self.fcs_dataframe = fcs_df

    


