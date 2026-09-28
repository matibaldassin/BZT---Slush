'''
MIT License

Copyright (c) 2022 Franco Aquistapace

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

# Import necesary modules
import pandas as pd
import numpy as np
from scipy.special import sph_harm
from sklearn.neighbors import NearestNeighbors

import os
import time
import timeit
import itertools

# Let's streamline the process of opening a file and
# extracting the features
def open_and_extract(fname,features,training=False):
    '''
    Parameters
    ----------
    fname : str
        Sample file path.
    features : list of str
        List with the features of the sample that are going to be used.

    Returns
    -------
    Original and normalized dataframes with the selected features of the 
    sample, ready for training of the SOM or to be classified. Also returns 
    the header lines of the file.

    '''
    
    file = open(fname,'r')
    found_cols = False
    header_lines = []
    line_count = 0
    while not found_cols:
        line = file.readline()
        line_count += 1
        if 'ITEM: ATOMS' in line:
            columns = line.split()
            columns.pop(0) # Pop 'ITEM:'
            columns.pop(0) # Pop 'ATOMS'
            found_cols = True
        header_lines.append(line)
    file.close()
    
    df = pd.read_csv(fname,sep=' ', 
                     skiprows=line_count, names=columns)
    
    
    norm_df = df[features].copy()
    # For the training process we use the min and max values 
    # of every feature from the same sample, and save them to properly
    # rescale the samples to be analized later
    if training:
        for feat in features:
            min_value = norm_df[feat].min()
            max_value = norm_df[feat].max()
            delta = max_value - min_value
            norm_df[feat] = (norm_df[feat] - min_value) / delta
            
            # Save min and max for later
            min_max_features.append([min_value, max_value])
            
    else:
        for i in range(len(features)):
            feat = features[i]
            min_value = min_max_features[i][0]
            max_value = min_max_features[i][1]
            delta = max_value - min_value
            norm_df[feat] = (norm_df[feat] - min_value) / delta
        
    return df, norm_df, header_lines


def export_sample(path, new_df, header_lines, new_cols_str):
    # Save new file with the groups (layers 1 and 2) assigned to each atom
    print('Writing results...')
    new_file = open(path, 'w')

    # Write the header of the file
    for i in range(len(header_lines)):
        line = header_lines[i]
        if i == len(header_lines) - 1:
            line = line.replace('\n', 
                                ' '+new_cols_str+'\n')
        new_file.write(line)

    # Now let's write the new data
    final_string = new_df.to_csv(index=False, sep=' ', 
                             float_format='%s', header=False)
    new_file.write(final_string)


    new_file.close()
    

# Define every auxiliary function here...


def split_sample(og_df,max_atoms):
    '''
    Params:
        og_df : pandas DataFrame
            DataFrame containing at least 'x', 'y' and 'z' columns.
        max_atoms : int
            Maximum atom count per split.
    Output:
        Returns a list of smaller DataFrames with slices of the
        sample.
    '''
    # If total amount of atoms is less than maximum, avoid
    # splitting
    if og_df.shape[0] < max_atoms:
        return [og_df]
    
    # First, we need to get the maximum and minimum x,y,z values
    x_max,x_min = og_df['x'].max(), og_df['x'].min()
    y_max,y_min = og_df['y'].max(), og_df['y'].min()
    z_max,z_min = og_df['z'].max(), og_df['z'].min()
    
    # Order coordinates based on variability
    axes = ['x','y','z']
    axes_var = np.array([x_max-x_min,y_max-y_min,z_max-z_min])
    sorted_axes = [axes[i] for i in np.argsort(axes_var)]
    sorted_var = axes_var[np.argsort(axes_var)]
    
    # Start with two boxes, along first axis
    # If not sufficient, split into more boxes
    done = False
    boxes = 2
    u = sorted_axes[-1] # Prefered axis
    while not done:
        split_dfs = []
        for i in range(boxes):
            start = og_df[u].min() + i*sorted_var[-1]/boxes
            if i == boxes-1:
                split_dfs.append(og_df.loc[og_df[u] >= start])
            else:
                end = og_df[u].min() + (i+1)*sorted_var[-1]/boxes 
                split_dfs.append(og_df.loc[(og_df[u] >= start) & (og_df[u] < end)])
        # Check if condition is satisfied
        done = True
        for j in range(len(split_dfs)):
            if split_dfs[j].shape[0] > max_atoms:
                done = False
        boxes += 1
    
    return split_dfs


def periodic_replicate(df, periodic_dims):
    '''
    Output:
        Generates periodic replicas of the given configurations,
        in all axes.
    '''
    replicas = [df]

    for i, ax in enumerate(['x', 'y', 'z']):
        # Displace all atoms by +-1 box size in
        # the given dimension
        box_size = max(periodic_dims[i]) - min(periodic_dims[i])

        # Init new replicas and displace
        new_df_plus = df.copy()
        new_df_minus = df.copy()
        new_df_plus[ax] = df[ax] + box_size
        new_df_minus[ax] = df[ax] - box_size

        # Save replicas
        replicas.append(new_df_plus)
        replicas.append(new_df_minus)

    return pd.concat(replicas)

    
def slice_padding(og_df, split_df, skin, 
                  periodic_axes, periodic_dims):
    '''
    Params:
        og_df : pandas DataFrame
            DataFrame containing coordinates for the entire
            sample.
        split_df : pandas DataFrame
            DataFrame containing a slice of the sample.
        skin : float
            Distance to use as skin for the padding of the 
            sample split.
        periodic_axes : list of str
            List containing the axes 'x', 'y' or/and 'z'
            that have periodic boundary conditions (PBC).
        periodic_dims : list of tuples
            List containing tuples (u_min,u_max), with the
            dimensions of the simulation box along the 
            axis u with PBC. Must be ordered in correspondence 
            with periodic_axes.
    Output:
        Returns a pandas DataFrame that contains the padding
        of the given slice of the sample.
    '''
    # Generate paddings from original data, for the given
    # slice
    pad_df = og_df.copy()
    
    # Get set of delimiting coordinates
    delim_coords = []
    axes = ['x','y','z']
    for i in range(len(axes)):
        u = axes[i]
        
        # Check if slice is at the limits of the sample,
        # and apply boundary conditions
        u_min, u_max = split_df[u].min(), split_df[u].max()
        if u in periodic_axes:
            box_min, box_max = periodic_dims[periodic_axes.index(u)]
            if np.abs(u_min-box_min) < skin:
                u_min = box_max # Go to the other side of the box
            if np.abs(u_max-box_max) < skin:
                u_max = box_min # Go to the other side of the box
        delim_coords.append([u_min,u_max])
    # Now, let's get the 6 individual paddings
    paddings = []
    positive_faces = []
    negative_faces = []
    
    u_tuples = [['x','y','z'],['y','z','x'],['z','x','y']]
    for n in range(3):
        u1,u2,u3 = u_tuples[n]
        # Faces +-
        face_plus = og_df.loc[
            (og_df[u1] > delim_coords[n][1]) & (og_df[u1] <= delim_coords[n][1]+skin)
        ].copy()
        face_plus = face_plus.loc[
                (face_plus[u2] >= split_df[u2].min()) & (face_plus[u2] <= split_df[u2].max())
            ]
        face_plus = face_plus.loc[
                (face_plus[u3] >= split_df[u3].min()) & (face_plus[u3] <= split_df[u3].max())
            ]
        face_minus = og_df.loc[
            (og_df[u1] < delim_coords[n][0]) & (og_df[u1] >= delim_coords[n][0]-skin)
        ].copy()
        face_minus = face_minus.loc[
                (face_minus[u2] >= split_df[u2].min()) & (face_minus[u2] <= split_df[u2].max())
            ]
        face_minus = face_minus.loc[
                (face_minus[u3] >= split_df[u3].min()) & (face_minus[u3] <= split_df[u3].max())
            ]
            

        positive_faces.append(face_plus.copy())
        negative_faces.append(face_minus.copy())
    
    # Shift faces in order to achieve periodic boundary conditions
    for i in range(3):
        u = axes[i]
        if u in periodic_axes:
            pos_f = positive_faces[i].copy()
            neg_f = negative_faces[i].copy()
            pos_dif = pos_f[u].max() - split_df[u].max()
            neg_dif = neg_f[u].min() - split_df[u].min()
            min_bound, max_bound = periodic_dims[periodic_axes.index(u)]
            if pos_dif < 0:
                pos_f[u] = pos_f[u] - pos_f[u].min() + max_bound
            if neg_dif > 0:
                neg_f[u] = neg_f[u] - neg_f[u].max() + min_bound
            positive_faces[i] = pos_f.copy()
            negative_faces[i] = neg_f.copy()
    
    paddings.extend(positive_faces)
    paddings.extend(negative_faces)
    
    # Now, let's build the remaining padding for the sides
    # of the slice
    for n in range(3):
        u1,u2,u3 = u_tuples[n]
        # In this case, it is easier to define coordinates
        # for each side base on pbc
        start_u1_plus = delim_coords[n][1]
        end_u1_plus = delim_coords[n][1] + skin
        start_u1_minus = delim_coords[n][0]
        end_u1_minus = delim_coords[n][0] - skin
        start_u2_plus = delim_coords[axes.index(u2)][1]
        end_u2_plus = delim_coords[axes.index(u2)][1] + skin
        start_u2_minus = delim_coords[axes.index(u2)][0]
        end_u2_minus = delim_coords[axes.index(u2)][0] - skin
        start_u3, end_u3 = split_df[u3].min(), split_df[u3].max()
        
        if u1 in periodic_axes:
            min_bound_u1, max_bound_u1 = periodic_dims[periodic_axes.index(u1)]
        if u2 in periodic_axes:
            min_bound_u2, max_bound_u2 = periodic_dims[periodic_axes.index(u2)]
        
        # Side 1 plus
        side_1_plus = og_df.loc[
            (og_df[u1] > start_u1_plus) & (og_df[u1] <= end_u1_plus)
        ].copy()
        side_1_plus = side_1_plus.loc[
                (side_1_plus[u2] > start_u2_plus) & (side_1_plus[u2] <= end_u2_plus)
            ]
        side_1_plus = side_1_plus.loc[
                (side_1_plus[u3] >= start_u3) & (side_1_plus[u3] <= end_u3)
            ]
        # Relocate based on pbc
        if side_1_plus[u1].max() < split_df[u1].max() and u1 in periodic_axes:
            side_1_plus[u1] = side_1_plus[u1] - side_1_plus[u1].min() + max_bound_u1
        if side_1_plus[u2].max() < split_df[u2].max() and u2 in periodic_axes:
            side_1_plus[u2] = side_1_plus[u2] - side_1_plus[u2].min() + max_bound_u2
        
        # Side 1 minus
        side_1_minus = og_df.loc[
            (og_df[u1] < start_u1_minus) & (og_df[u1] >= end_u1_minus)
        ].copy()
        side_1_minus = side_1_minus.loc[
                (side_1_minus[u2] > start_u2_plus) & (side_1_minus[u2] <= end_u2_plus)
            ]
        side_1_minus = side_1_minus.loc[
                (side_1_minus[u3] >= start_u3) & (side_1_minus[u3] <= end_u3)
            ]
        # Relocate based on pbc
        if side_1_minus[u1].min() > split_df[u1].min() and u1 in periodic_axes:
            side_1_minus[u1] = side_1_minus[u1] - side_1_minus[u1].max() + min_bound_u1
        if side_1_minus[u2].max() < split_df[u2].max() and u2 in periodic_axes:
            side_1_minus[u2] = side_1_minus[u2] - side_1_minus[u2].min() + max_bound_u2
        
        # Side 2 plus
        side_2_plus = og_df.loc[
            (og_df[u1] > start_u1_plus) & (og_df[u1] <= end_u1_plus)
        ].copy()
        side_2_plus = side_2_plus.loc[
                (side_2_plus[u2] < start_u2_minus) & (side_2_plus[u2] >= end_u2_minus)
            ]
        side_2_plus = side_2_plus.loc[
                (side_2_plus[u3] >= start_u3) & (side_2_plus[u3] <= end_u3)
            ]
        # Relocate based on pbc
        if side_2_plus[u1].max() < split_df[u1].max() and u1 in periodic_axes:
            side_2_plus[u1] = side_2_plus[u1] - side_2_plus[u1].min() + max_bound_u1
        if side_2_plus[u2].min() > split_df[u2].min() and u2 in periodic_axes:
            side_2_plus[u2] = side_2_plus[u2] - side_2_plus[u2].max() + min_bound_u2
        
        # Side 2 minus
        side_2_minus = og_df.loc[
            (og_df[u1] < start_u1_minus) & (og_df[u1] >= end_u1_minus)
        ].copy()
        side_2_minus = side_2_minus.loc[
                (side_2_minus[u2] < start_u2_minus) & (side_2_minus[u2] >= end_u2_minus)
            ]
        side_2_minus = side_2_minus.loc[
                (side_2_minus[u3] >= start_u3) & (side_2_minus[u3] <= end_u3)
            ]
        # Relocate based on pbc
        if side_2_minus[u1].min() > split_df[u1].min() and u1 in periodic_axes:
            side_2_minus[u1] = side_2_minus[u1] - side_2_minus[u1].max() + min_bound_u1
        if side_2_minus[u2].min() > split_df[u2].min() and u2 in periodic_axes:
            side_2_minus[u2] = side_2_minus[u2] - side_2_minus[u2].max() + min_bound_u2
        
        paddings.append(side_1_plus.copy())
        paddings.append(side_1_minus.copy())
        paddings.append(side_2_plus.copy())
        paddings.append(side_2_minus.copy())
        
    
    return pd.concat(paddings)

def get_periodic(header_lines):
    '''
    Params:
        header_lines : list of str
            Last element resulting from the open_and_extract 
            function. Must have info about periodic boundary
            conditions of a given sample.
    Output:
        Returns a tuple (periodic_axes, periodic_dims), 
        where periodic_axes is a list of elements u, with 
        u belonging to ['x','y','z'], indicating which 
        axes are periodic; and periodic_dims is a list of
        tuples with the simulation box bounds along each
        of the given periodic axes.
    '''
    # Go through every line and find the BOX BOUNDS item
    periodic_axes, periodic_dims = [], []
    axes = ['x','y','z']
    for i in range(len(header_lines)):
        line = header_lines[i]
        if 'BOX BOUNDS' in line:
            # Let's extract the periodic boundary conditions
            pbc_string = line.replace('ITEM: BOX BOUNDS','').split()
            # Use only the last 3 tokens (x, y, z boundary conditions)
            pbc_axes = pbc_string[-3:] if len(pbc_string) >= 3 else pbc_string
            for j in range(len(pbc_axes)):
                if 'pp' in pbc_axes[j]:
                    periodic_axes.append(axes[j])
            # Now let's get the periodic dimensions
            for j in range(3):
                if axes[j] in periodic_axes:
                    axis_dims = header_lines[i+j+1].split()
                    periodic_dims.append([float(dim) for dim in axis_dims])
    return (periodic_axes, periodic_dims)

def get_positions(n, data, indexes):
    return list(np.array(data.iloc[indexes[n]][['x','y','z']]).flatten())


def generate_vectors(data,nearest_neighbors):
    '''
    data : pandas DataFrame
        Data containing (at least) the atom positions x, y and z.
    nearest_neighbors : tuple
        Result of a NearestNeighbors neighbor graph search on data, with
        mode='connectivity'.
        
    Returns:
        List of arrays containing neighbor distance vectors for 
        each atom.
    '''
    # Eliminate self connections
    neighbors = nearest_neighbors.astype(bool)
    atom_data = np.array(data[['x','y','z']])
    tot_list = []
    for i in range(data.shape[0]):
        # Get neighbors
        atom_neigh = atom_data[neighbors[i,:]]
        # Calculate distance vectors and turn into np array
        atom_neigh = atom_neigh-atom_data[i,:]
        tot_list.append(atom_neigh)
    
    return tot_list


# Aux function to get cluster sizes for a given data
def get_cluster_sizes(df, cluster_col='cluster_id'):
    '''
    df : pandas DataFrame
        Data containing (at least) the atom cluster_id column.
    cluster_col : str (optional)
        Name for the column representing cluster ids. Default
        value is "cluster_id".
    
    Returns:
        Value counts containing the size of each polar cluster.
    '''
    # Exclude cluster_id=-1
    cluster_data = df.loc[df[cluster_col] > -1, cluster_col]
    # Return cluster sizes
    return cluster_data.value_counts()

