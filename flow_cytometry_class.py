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
import access_fcs_fields as af


class flow_cytometry_class(object):
    """_summary_

    Args:
        object (_type_): a class that includes functions for converting .fcs data into Pandas DataFrames,
        segmenting cells from flow cytometry images, and gate histogram or scatter plot data.
    """
    def __init__(self, fcs_path, images_folder, experiment_id):
        
        fd = flowio.FlowData(fcs_path)
        
        self.experiment_id = experiment_id
        
        channels_list = []
        for ch in range(1, len(fd.channels)+1):
            try:
                ch_string = fd.channels[str(ch)]
            except KeyError:
                ch_string = fd.channels[ch]
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
            gate_dict = prs.load_data('stored_gates_'+self.experiment_id)
            self.gate_dict = gate_dict
        except FileNotFoundError:
            gate_dict = {}
            self.gate_dict = gate_dict
        
        print('Applying stored gates to the dataframe...')
        for gt in gate_dict:
            if gate_dict[gt][-1] == 'histogram_gate':
                print(f'Applying histogram gate {gt}')
                gate_range = gate_dict[gt][0]
                x_variable = gate_dict[gt][1]
                log_bool = gate_dict[gt][2]
                x, log_string = af.access_single_column_data(fc_df, x_variable, log_bool)
                fc_df = span.apply_histogram_gate_to_dataframe(x, gate_range, gt, fc_df)
            elif gate_dict[gt][-1] == 'polygon_gate':
                print(f'Applying polygon gate {gt}')
                polygon_coords = gate_dict[gt][0]
                variables = gate_dict[gt][1]
                log_bools = gate_dict[gt][2]
                x, y, log_string = af.access_double_column_data(fc_df, variables, log_bools)
                fc_df = poly.apply_polygon_gate_to_dataframe(x, y, polygon_coords, gt, fc_df)
                
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
    
    
    # ------------------------ all the get functions -----------------------#
    def get_image_arrays(self):
        return self.image_arrays
    
    def get_flow_cytometry_dataframe(self):
        return self.fcs_dataframe
    
    def get_metadata(self):
        self.metadata
        
    def get_cell_masks(self):
        return self.cell_masks

    def get_cell_areas(self):
        return self.cell_areas
    
    def get_gates(self):
        return self.gate_dict

    
    # ------------------------ Image segmentation -----------------------#
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
        
        
    # ------------------------ Gating -----------------------#
    def import_gates(self, gate_path):
        imported_gate_dict = prs.load_data(gate_path)
        val = input('Warning, this operation will overwrite all previously stored gates. Continue?: YES/NO: ')
        val = val.upper()
        i = 0
        while i == 0:
            if val == 'YES':
                self.gate_dict = imported_gate_dict
                prs.save_data(imported_gate_dict, 'stored_gates_'+self.experiment_id)
                i+=1
            elif val == 'NO':
                print('Gate-import aborted by the user')
                pass
                i+=1

    def reset_gates(self):
        val = input('Warning, this operation will erase all stored gates and remove them from the dataframe: YES/NO: ')
        val = val.upper()
        if val == 'YES':
            fcs_df = self.get_flow_cytometry_dataframe()
            fcs_df = fcs_df.drop(list(self.gate_dict.keys()), axis=1)
            
            self.gate_dict = {}
            prs.save_data(self.gate_dict, 'stored_gates_'+self.experiment_id)
            
            self.fcs_dataframe  = fcs_df
            
        elif val == 'NO':
            print(self.get_gates())
            print('Gates have not not been reset.')
            pass

        else:
            self.reset_gates()
            
    
    def histogram_gate(self, variable, log, bin_array):
        
        fcs_df = self.get_flow_cytometry_dataframe()
        gate_dict = self.get_gates()
        
        x, log_string = af.access_single_column_data(fcs_df, variable, log)
        
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
        prs.save_data(gate_dict, 'stored_gates_'+self.experiment_id)
        self.gate_dict = gate_dict
        os.remove(gate_name)
        print(self.get_gates())
        
        gate_range = gate_dict[gate_name][0]
        fcs_df = span.apply_histogram_gate_to_dataframe(x, gate_range, gate_name, fcs_df)
        self.fcs_dataframe = fcs_df
        

    def scatter_gate(self, x_variable, y_variable, x_log, y_log, sample_size=10000):
            
            fcs_df = self.get_flow_cytometry_dataframe()
            gate_dict = self.get_gates()
            
            x, y, log_string = af.access_double_column_data(fcs_df, (x_variable, y_variable), (x_log, y_log))
            
            i = 0
            n = 1
            while i==0:
                gate_name = x_variable+'_'+y_variable+log_string+'_gate_'+str(n)
                if gate_name not in self.gate_dict:
                    i+=1
                else:
                    n+=1
            
            poly.return_markers_in_polygon(x, y, x_variable, y_variable, gate_name, sample_size, colormap='rainbow')
            
            loaded_gate = prs.load_data(gate_name)
            gate_dict[gate_name] = loaded_gate
            prs.save_data(gate_dict, 'stored_gates_'+self.experiment_id)
            self.gate_dict = gate_dict
            os.remove(gate_name)
            print(self.get_gates())
            
            polygon_coordinates = gate_dict[gate_name][0]
            fcs_df = poly.apply_polygon_gate_to_dataframe(x, y, polygon_coordinates, gate_name, fcs_df)
            
            self.fcs_dataframe = fcs_df

        


