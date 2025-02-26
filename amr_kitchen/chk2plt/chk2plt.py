import os
import re
from tqdm import tqdm
from multiprocessing import Pool
import numpy as np
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import shape_from_header, header_from_indices
from .checkpoint_reader import CheckpointReader


def write_plt_bin_from_chk(args):
    (b_state, b_gradp, b_I_R,
     idxs_state, offsets_gradp, offsets_I_R,
     b_plt, state_field_indices,
     do_gradp, do_species_reactions, floor_massfracs) = args
    """
    Read checkpoint data and write it to a binary file in a plotfile
    ---
    b_state: State vector binary file, the distribution mapping is taken from the state
             vector data so only one checkpoint binary file is read for each function
             call.
    b_gradp: Checkpoint binary files containing the pressure gradient data. A list of files
             for each box in the state vector binary file, following sequential box order
             from the state binary file.
    b_I_R:   Checkpoint binary files containing the species reaction rate data.
             A list of files for each box in the state vector binary file, following
             sequential box order from the state binary file.
    idxs_state: global indices of the boxes in the state vector binary file, in the order
                they are written in the file.
    offsets_gradp: Binary file offsets of the pressure gradient data boxes.
    offsets_I_R: Binary file offsets of the species reaction rate data boxes.
    b_plt: Destination binary file in the plotfile directory where the data is written.
    state_field_indices: Indices of the different fields in the state vector data.
    do_gradp: Flag to include pressure gradient data in the destination plotfile.
    do_species_reactions: Flag to include species net consumption rate in the destination
                          plotfile.
    floor_massfracs: Flag to rescale mass fraction so they sum to 1.0
    """
    offsets_plt = []
    mins_plt = []
    maxs_plt = []

    with open(b_plt, 'wb') as bp:
        with open(b_state, 'rb') as bs:
            bid = 0
            while True:
                try:
                    state_shape = shape_from_header(bs.readline().decode('ascii'))
                except:
                    break
                data = np.fromfile(bs, 'float64', np.prod(state_shape))
                # Remove ghosts cells
                box_indices = idxs_state[bid]
                box_shape = box_indices[1] - box_indices[0] + 1
                n_ghosts = (state_shape[:-1] - box_shape)//2
                data = data.reshape(state_shape, order='F')
                data = data[n_ghosts[0]:-n_ghosts[0],
                            n_ghosts[1]:-n_ghosts[1],
                            n_ghosts[2]:-n_ghosts[2], :]
                # Normalizing species mass fractions
                if floor_massfracs:
                    Y_start = state_field_indices['Y_start']
                    Y_end = state_field_indices['Y_end']
                    Y_sum = np.sum(data[..., Y_start:Y_end], axis=-1)
                    data[..., Y_start:Y_end] /= Y_sum[..., np.newaxis]
                # Read pressure gradient data
                if do_gradp:
                    with open(b_gradp[bid], 'rb') as bg:
                        bg.seek(offsets_gradp[bid])
                        shape_gradp = shape_from_header(bg.readline().decode('ascii'))
                        data_gradp = np.fromfile(bg, 'float64', np.prod(shape_gradp))
                    data_gradp = data_gradp.reshape(shape_gradp, order='F')
                    data = np.concatenate([data, data_gradp], axis=-1)
                # Read species reaction rate data
                if do_species_reactions:
                    with open(b_I_R[bid], 'rb') as br:
                        br.seek(offsets_I_R[bid])
                        shape_I_R = shape_from_header(br.readline().decode('ascii'))
                        data_I_R = np.fromfile(br, 'float64', np.prod(shape_I_R))
                    data = np.concatenate([data, data_I_R], axis=-1)

                hw = header_from_indices(box_indices[0],
                                         box_indices[1],
                                         data.shape[-1])
                offsets_plt.append(bp.tell())
                mins_plt.append(np.min(data, axis=(0, 1, 2)))
                maxs_plt.append(np.max(data, axis=(0, 1, 2)))
                bp.write(hw)
                bp.write(data.flatten(order='F').tobytes())
                bid += 1
    return offsets_plt, mins_plt, maxs_plt


class chk2plt(CheckpointReader):
    """
    Class inheriting from the checkpoint reader to
    convert checkpoint data to a plotfile
    """

    def __init__(self, chkdir, target_plotfile=None, species=None,
                 gradp=True, species_reactions=False, floor_massfracs=True, pltdir=None):
        """
        Constructor which does the conversion from checkpoint to plotfile
        chkdir: Directory of the checkpoint
        target_plotfile: Plotfile to use to obtain the species names
        species: list of the expected species in the checkpoint
                 (if a target plotfile is not specified)
        gradp: flag to include pressure gradient data in the plotfile
        species_reactions: flag to include species reactions in the plotfile
        floor_massfracs: flag to normalize species mass fractions to 1.0
        pltdir: User defined plotfile output directory default is the checkpoint
                step with the plotfile prefix
        """
        super().__init__(chkdir, maxmins=False)
        self.do_gradp = gradp
        self.do_species_reactions = species_reactions
        self.floor_massfracs = floor_massfracs
        if pltdir is None:
            root, chkbase = os.path.split(chkdir)
            self.pltdir = os.path.join(root, chkbase.replace('chk', 'plt'))
        else:
            self.pltdir = pltdir
        if target_plotfile is not None:
            self.pck_ref = PlotfileCooker(target_plotfile, header_only=True)
            species_data = [f for f in self.pck_ref.fields if re.search(r'Y\(.*\)', f) is not None]
            # Fall back to reaction rates
            if len(species_data) == 0:
                species_data = [f for f in self.pck_ref.fields if re.search(r'I_R\(.*\)', f) is not None]
                self.species = [y.replace('I_R(', '')[:-1] for y in species_data]
            else:
                self.species = [y.replace('Y(', '')[:-1] for y in species_data]
        elif len(species) > 0:
            self.species = species
        else:
            raise ValueError(("Checkpoints do not contain species names "
                              "so they must be user defined"))
        # Number of fields in the output
        self.nfields_out = self.nfields['state']
        if self.do_gradp:
            self.nfields_out += self.nfields['gradp']
        if self.do_species_reactions:
            self.nfields_out += self.nfields['I_R']
        # Fields in the output
        fields = ['x_velocity',
                  'y_velocity',
                  'z_velocity',
                  'density']
        for sp in self.species:
            fields.append(f'Y({sp})')
        for f in ['rhoh', 'temp', 'RhoRT']:
            fields.append(f)
        if self.do_gradp:
            for f in ['gradpx', 'gradpy', 'gradpz']:
                fields.append(f)
        if self.do_species_reactions:
            for sp in self.species:
                fields.append(f'I_R({sp})')
        self.fields_out = fields

        if len(self.fields_out) != self.nfields_out:
            raise ValueError(('The number of fields found in the checkpoint does'
                              ' not match the target plotfile the species in each'
                              ' probably do not match'))
        self.convert()


    def convert(self):
        """
        Main function converting the checkpoint to plotfile
        """
        # Make the plotfile structure
        os.makedirs(os.path.join(os.getcwd(), self.pltdir), exist_ok=True)
        for level in range(self.max_level + 1):
            level_dir = f"Level_{level}"
            os.makedirs(os.path.join(os.getcwd(), self.pltdir, level_dir), exist_ok=True)

        # Use the state vector binary files as the reference box distribution
        for level in range(self.max_level + 1):
            # Level root directories
            lv_chk_root = os.path.join(self.chkdir, f"Level_{level}")
            lv_plt_root = os.path.join(self.pltdir, f"Level_{level}")
            # Data subsets binary files
            all_bin_gradp = self.boxes[level]['gradp_paths']
            all_bin_I_R = self.boxes[level]['I_R_paths']
            # indices of all boxes at the current level
            all_box_indices = self.boxes[level]['indices']
            # Empty arrays to store the level header infos
            all_offsets_plt = np.empty(self.nboxes[level])
            all_mins_plt = np.empty((self.nboxes[level], self.nfields_out))
            all_maxs_plt = np.empty((self.nboxes[level], self.nfields_out))
            all_binfiles_plt = np.empty(self.nboxes[level], dtype=object)
            # Box indices at the current level
            lv_bids = np.arange(self.nboxes[level])

            # divide by unique state vector data binary file
            state_bin_box_ids = [] #ids of the boxes in each state binary file
            mp_args = []
            # For each unique state vector binary file
            state_bin_unique = np.unique(self.boxes[level]['state_paths'])
            for state_bin in state_bin_unique:
                # Box indices in the current state vector binary file
                bids = lv_bids[self.boxes[level]['state_paths'] == state_bin]
                # Offsets of the boxes in the binary file
                offsets_state = self.boxes[level]['state_offsets'][bids]
                # Get the offsets order so we can read the file sequentially
                offsets_state_order = np.argsort(offsets_state)
                # Reorder the current box indices
                bids = bids[offsets_state_order]
                state_bin_box_ids.append(bids)
                # Target binary file in the plotfile matching the state vector binary
                bin_path_plt = os.path.join(lv_plt_root, state_bin.replace('state', 'Cell'))
                # Store the target plotfile binary files in the empty array
                all_binfiles_plt[bids] = state_bin.replace('state', 'Cell')
                # Single multiprocessing call
                mp_call = [os.path.join(lv_chk_root, state_bin), # state binary file path
                           # paths to each box in the gradp binary files
                           [os.path.join(lv_chk_root, b_gradp) for b_gradp in all_bin_gradp[bids]],
                           # paths to each box in the I_R binary files
                           [os.path.join(lv_chk_root, b_I_R) for b_I_R in all_bin_I_R[bids]],
                           all_box_indices[bids], # indices for the current boxes
                           self.boxes[level]['gradp_offsets'][bids], # offsets in the gradp binary files
                           self.boxes[level]['I_R_offsets'][bids], # offsets in the I_R binary files
                           bin_path_plt, # Path to the target plotfile binary file
                           self.state_field_indices, # Field indices in the state vector
                           self.do_gradp, # gradp flag
                           self.do_species_reactions, # species mass fraction flag
                           self.floor_massfracs] # Floor mass fracs flag
                mp_args.append(mp_call)
            # Call the multiprocessing
            out = tqdm(Pool().imap(write_plt_bin_from_chk, mp_args), total = len(mp_args))

            for (offsets, mins, maxs), bid in zip(out, state_bin_box_ids):
                all_offsets_plt[bid] = offsets
                all_mins_plt[bid, :] = mins
                all_maxs_plt[bid, :] = maxs

            self.write_level_header(level, all_binfiles_plt, all_offsets_plt,
                                    all_mins_plt, all_maxs_plt)
        self.write_global_header()

    def write_level_header(self, level, all_binfiles_plt, all_offsets_plt,
                           all_mins_plt, all_maxs_plt):
        # Write the Level header
        with open(os.path.join(self.pltdir, f'Level_{level}', 'Cell_H'), 'w') as lhead:
            lhead.write('1\n')
            lhead.write('1\n')
            lhead.write(f'{self.nfields_out}\n')
            lhead.write('0\n')
            lhead.write(f'({self.nboxes[level]} 0\n')
            for bidx in self.boxes[level]['indices']:
                lhead.write((f'(({bidx[0][0]},{bidx[0][1]},{bidx[0][2]})'
                             f' ({bidx[1][0]},{bidx[1][1]},{bidx[1][2]})'
                             f' (0,0,0))\n'))
            lhead.write(')\n')
            lhead.write(f'{self.nboxes[level]}\n')
            for bfile, offset in zip(all_binfiles_plt, all_offsets_plt):
                lhead.write(f'FabOnDisk: {bfile} {int(offset)}\n')
            lhead.write('\n')
            lhead.write(f'{self.nboxes[level]},{self.nfields_out}\n')
            for minvals in all_mins_plt:
                lhead.write(','.join([f"{m:.16e}" for m in minvals]) + ',\n')
            lhead.write('\n')
            lhead.write(f'{self.nboxes[level]},{self.nfields_out}\n')
            for maxvals in all_maxs_plt:
                lhead.write(','.join([f"{m:.16e}" for m in maxvals]) + ',\n')

    def write_global_header(self):
        # Write the global header
        with open(os.path.join(self.pltdir, 'Header'), 'w') as ghead:
            ghead.write('HyperCLaw-V1.1\n')
            ghead.write(f'{self.nfields_out}\n')
            for f in self.fields_out:
                ghead.write(f'{f}\n')
            ghead.write(f'{len(self.geo_hi)}\n')
            ghead.write(f'{self.time}\n')
            ghead.write(f'{self.max_level}\n')
            ghead.write(' '.join([f'{p}' for p in self.geo_lo]) + '\n')
            ghead.write(' '.join([f'{p}' for p in self.geo_hi]) + '\n')
            # Assume refinement factors are always 2
            ghead.write(' '.join(['2' for _ in range(self.max_level)]) + '\n')
            # Grid sizes
            gridsize_str = []
            for level in range(self.max_level + 1):
                gs_str = [f'{gs -1}' for gs in self.grid_sizes[level]]
                gridsize_str.append('((0,0,0) (' + ','.join(gs_str) + ') (0,0,0))')
            ghead.write(' '.join(gridsize_str) + '\n')
            # Step numbers
            ghead.write(' '.join([f'{self.step_number}' for _ in range(self.max_level + 1)]) + '\n')
            # grid resolutions
            for dx_vals in self.dx:
                ghead.write(' '.join([f'{dx}' for dx in dx_vals]) + '\n')
            ghead.write('0\n')
            ghead.write('0\n')
            for level in range(self.max_level + 1):
                ghead.write(f"{level} {self.nboxes[level]} {self.time}\n")
                ghead.write(f"{self.step_number}\n")
                for box in self.compute_boxes_bounds(level):
                    for c in range(3):
                        ghead.write(f"{box[0][c]} {box[1][c]}\n")
                ghead.write(f'Level_{level}/Cell\n')

