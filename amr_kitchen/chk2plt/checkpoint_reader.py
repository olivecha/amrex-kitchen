import os
import numpy as np
from amr_kitchen.utils import shape_from_header


class CheckpointReader(object):
    """
    A class to read PeleLMeX checkpoint metadata
    This takes for granted the variables PeleLMeX puts
    in its checkpoints so it would probably not work
    for other solvers/codes without adjustments
    """
    # Indices of the variables in the state vector
    state_field_indices = {'x_velocity':0,
                           'y_velocity':1,
                           'z_velocity':2,
                           'density':3,
                           'Y_start':4,
                           'Y_end':-3,
                           'rhoh':-3,
                           'temp':-2,
                           'RhoRT':-1}
    # With data subset in the checkpoint has ghost
    # cells in the binary files
    data_has_ghost = {'I_R': False,
                      'divU': True,
                      'gradp': False,
                      'p': True,
                      'state': True}

    def __init__(self, chkdir, maxmins=False):
        """
        Read the checkpoint metadata and store it in a class
        chkdir: path to the checkpoint directory
        maxmins: flag to read max/min values in the cell headers
        """
        self.chkdir = chkdir
        # Read the checkpoint header
        with open(os.path.join(self.chkdir, 'Header')) as hfile:
            self.version = hfile.readline().replace('\n', '')
            self.max_level = int(hfile.readline().replace('\n', ''))
            self.step_number = int(hfile.readline().replace('\n', ''))
            # Sometimes there is an int there
            value = float(hfile.readline())
            # Check if its an integer
            if (value % 1 == 0):
                # Read the time on the next line
                self.time = float(hfile.readline().replace('\n', ''))
            else:
                # Other time its the step time
                self.time = value
            # TODO: find out what those are in PeleLMeX
            self.small_float1 = float(hfile.readline().replace('\n', ''))
            self.small_float2 = float(hfile.readline().replace('\n', ''))
            self.geo_lo = np.array([float(n) for n in hfile.readline().split()])
            self.geo_hi = np.array([float(n) for n in hfile.readline().split()])
            self.boxes = []
            # Read the mesh boxing data
            for lv in range(self.max_level + 1):
                level_boxes = {'indices':[]}
                nboxes, nfaces = [int(n) for n in hfile.readline().replace('(', '').split()]
                for _ in range(nboxes):
                    b_start, b_end, b_faces = hfile.readline().split()
                    b_start = b_start.replace('(', '').replace(')', '')
                    b_start = [int(n) for n in b_start.split(',')]
                    b_end = b_end.replace('(', '').replace(')', '')
                    b_end = [int(n) for n in b_end.split(',')]
                    level_boxes['indices'].append([b_start, b_end])
                hfile.readline()
                level_boxes['indices'] = np.array(level_boxes['indices'])
                self.boxes.append(level_boxes)
            self.pressure = float(hfile.readline().replace('\n', ''))
            # Normally the coordinate system is here
            value = float(hfile.readline())
            if (value % 1 == 0):
                self.sys_coord = int(value)
                # Then there is a zero after
                assert 0 == int(hfile.readline())
                # To store the typ vals
                field_counter = 0
                self.typvals = []
            else:
                # Sometimes the coord system is not here
                field_counter = 1
                self.typvals = [value]
            # Read the typical values
            for l in hfile:
                self.typvals.append(float(l))
                field_counter += 1
            self.nfields = {'typval': field_counter}
        # Compute grid sizes

        self.grid_sizes = [np.max(self.boxes[0]['indices'][:, 1, :], axis=0) + 1]
        for lv in range(1, self.max_level + 1):
            self.grid_sizes.append(self.grid_sizes[-1] * 2)
        self.nboxes = np.array([self.boxes[lv]['indices'].shape[0] for lv in range(self.max_level + 1)])
        # Compute grid resolution
        self.domain = self.geo_hi - self.geo_lo
        self.dx = np.array([self.domain / self.grid_sizes[lv] for lv in range(self.max_level + 1)])
        # Compute uniform grids
        self.box_centers = []
        for lv in range(self.max_level + 1):
            grids = [np.linspace(self.geo_lo[coord] + self.dx[lv][coord]/2, 
                                 self.geo_hi[coord] - self.dx[lv][coord]/2,
                                 self.grid_sizes[lv][coord]) for coord in range(3)]
            self.box_centers.append(grids)

        # Read level fields headers
        for lv in range(self.max_level + 1):
            for sub_data in ['I_R', 'divU', 'gradp', 'p', 'state']:
                (self.nfields[f'{sub_data}'],
                 self.boxes[lv][f'{sub_data}_paths'],
                 self.boxes[lv][f'{sub_data}_offsets'],
                 self.boxes[lv][f'{sub_data}_mins'],
                 self.boxes[lv][f'{sub_data}_maxs']) = self.read_level_header(lv, f'{sub_data}_H', maxmins=maxmins)

    def compute_boxes_bounds(self, level):
        """
        Compute the boxes bounds at an amr level to write in the
        plotfile header
        """
        level_indices = self.boxes[level]['indices']
        boxes = np.empty_like(level_indices, dtype=float)
        boxes[:, :, 0] = self.box_centers[level][0][level_indices[:, :, 0]]
        boxes[:, :, 1] = self.box_centers[level][1][level_indices[:, :, 1]]
        boxes[:, :, 2] = self.box_centers[level][2][level_indices[:, :, 2]]
        boxes[:, 0, :] -= self.dx[level][0]/2
        boxes[:, 1, :] += self.dx[level][0]/2
        return boxes

    def read_level_header(self, level, hfilename, maxmins=False):
        """
        Read a level header and return paths, offsets and minmaxes
        """
        paths = []
        offsets = []
        min_vals = []
        max_vals = []

        lvhead_path = os.path.join(self.chkdir,
                                   f'Level_{level}',
                                   hfilename)
        with open(lvhead_path) as hfile:
            hfile.readline()
            hfile.readline()
            n_fields = int(hfile.readline())
            hfile.readline()
            # Skip the box indices
            hfile.readline()
            for _ in range(self.nboxes[level]):
                hfile.readline()
            hfile.readline()
            hfile.readline()
            for _ in range(self.nboxes[level]):
                _, box_path, box_offset = hfile.readline().split()
                paths.append(box_path)
                offsets.append(int(box_offset))
            if maxmins:
                hfile.readline()
                hfile.readline()
                for _ in range(self.nboxes[level]):
                    box_mins = np.array(hfile.readline().split(',')[:-1], dtype=float)
                    min_vals.append(box_mins)
                hfile.readline()
                hfile.readline()
                for _ in range(self.nboxes[level]):
                    box_maxs = np.array(hfile.readline().split(',')[:-1], dtype=float)
                    max_vals.append(box_maxs)
        return n_fields, np.array(paths), np.array(offsets), np.array(min_vals), np.array(max_vals)

    def read_box(self, bid, prefix, lv):
        """
        Read a single FAB from the checkpoint data
        """
        box_file = self.boxes[lv][f'{prefix}_paths'][bid]
        box_file_path = os.path.join(self.chkdir, f'Level_{lv}', box_file)
        with open(box_file_path, 'rb') as bfile:
            bfile.seek(self.boxes[lv][f'{prefix}_offsets'][bid])
            shape = shape_from_header(bfile.readline().decode('ascii'))
            data = np.fromfile(bfile, 'float64', np.prod(shape))
        # Remove ghosts cells
        if self.data_has_ghost[prefix]:
            box_indices = self.boxes[lv]['indices'][bid]
            box_shape = box_indices[1] - box_indices[0] + 1
            n_ghosts = (shape[:-1] - box_shape)//2
            data = data.reshape(shape, order='F')
            data = data[n_ghosts[0]:-n_ghosts[0],
                        n_ghosts[1]:-n_ghosts[1],
                        n_ghosts[2]:-n_ghosts[2], :]
            return data
        else:
            return data.reshape(shape, order='F')
