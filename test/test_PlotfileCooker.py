import os
import unittest
import warnings
import numpy as np

from amr_kitchen import PlotfileCooker, DTYPE
from amr_kitchen import plotfile_cooker
from amr_kitchen.plotfile_cooker import (region_overlap, expand_array_nd,
                                         zoom_array_nd)

class TestSliceData(unittest.TestCase):
    pfile2d = "test_assets/example_plt_2d"
    pfile3d = "test_assets/example_plt_3d"

    def test_load2d(self):
        for Lv in [0, 1]:
            hdr = PlotfileCooker(self.pfile2d, limit_level=Lv)
            self.assertIsInstance(hdr, PlotfileCooker)
            self.assertEqual(hdr.ndims, 2)
            self.assertEqual(hdr.limit_level, Lv)

    def test_load3d(self):
        for Lv in [0, 1, 2]:
            hdr = PlotfileCooker(self.pfile3d, limit_level=Lv)
            self.assertIsInstance(hdr, PlotfileCooker)
            self.assertEqual(hdr.ndims, 3)
            self.assertEqual(hdr.limit_level, Lv)

    def test_minmaxs2d(self):
        hdr = PlotfileCooker(self.pfile2d, maxmins=True)
        self.assertIsInstance(hdr, PlotfileCooker)
        self.assertTrue('mins' in hdr.cells[0])
        self.assertTrue('maxs' in hdr.cells[0])

    def test_minsmaxs3d(self):
        hdr = PlotfileCooker(self.pfile3d, maxmins=True)
        self.assertIsInstance(hdr, PlotfileCooker)
        self.assertTrue('mins' in hdr.cells[0])
        self.assertTrue('maxs' in hdr.cells[0])

    def test_field_slice_offset(self):
        # Slicing the fields must select the same fields as slicing the
        # array of every field, including when the slice does not start
        # at the first field and when it has a step
        for pfile in [self.pfile2d, self.pfile3d]:
            hdr = PlotfileCooker(pfile, serial=True)
            nfields = len(hdr.fields)
            full = hdr[:][0][0]
            self.assertEqual(full.shape[-1], nfields)
            for key in [slice(None, 3), slice(3, 6), slice(0, 6, 2),
                        slice(3, nfields, 3), slice(nfields - 2, None),
                        slice(2, None)]:
                data = hdr[key][0][0]
                self.assertTrue(np.array_equal(data, full[..., key]),
                                (pfile, key))
            # The fields are read in the order they are stored so a
            # backward slice is not supported
            with self.assertRaises(ValueError):
                hdr[6:2:-1][0]

    def test_bybinfile_iterator2d(self):
        hdr = PlotfileCooker(self.pfile2d)
        for lv in range(hdr.limit_level + 1):
            for bfname, offsets, indexes in hdr.bybinfile(lv):
                with open(bfname) as bf:
                    for idx, ofs in zip(indexes, offsets):
                        shape = [idx[1][i] - idx[0][i] + 1 for i in range(hdr.ndims)]
                        shape.append(len(hdr.fields))
                        bf.seek(ofs)
                        arr = np.fromfile(bf, DTYPE, np.prod(shape))

    def test_bybinfile_iterator3d(self):
        hdr = PlotfileCooker(self.pfile3d)
        for lv in range(hdr.limit_level + 1):
            for bfname, offsets, indexes in hdr.bybinfile(lv):
                with open(bfname) as bf:
                    for idx, ofs in zip(indexes, offsets):
                        shape = [idx[1][i] - idx[0][i] + 1 for i in range(hdr.ndims)]
                        shape.append(len(hdr.fields))
                        bf.seek(ofs)
                        arr = np.fromfile(bf, DTYPE, np.prod(shape))


class TestCoordinateSlice(unittest.TestCase):
    """
    Slicing the grid coordinates ('x', 'y', 'z') through
    PlotfileCooker.__getitem__. The box_points method is used as the
    ground truth for the expected coordinate arrays.
    """
    pfile2d = "test_assets/example_plt_2d"
    pfile3d = "test_assets/example_plt_3d"

    def test_single_coord_box(self):
        # A single coordinate field is returned as the bare tiled 3D array
        hdr = PlotfileCooker(self.pfile3d)
        for lv in range(hdr.limit_level + 1):
            for bid in range(len(hdr.cells[lv]['indexes'])):
                points = hdr.box_points(lv, bid)
                for axis, name in enumerate(['x', 'y', 'z']):
                    data = hdr[name][lv][bid]
                    self.assertEqual(data.ndim, hdr.ndims)
                    self.assertTrue(np.array_equal(data, points[axis]))

    def test_integer_index_convention(self):
        # Coordinates are the three last fields when indexing by integer
        hdr = PlotfileCooker(self.pfile3d)
        nf = len(hdr.fields)
        lv, bid = 1, 3
        points = hdr.box_points(lv, bid)
        for axis in range(3):
            data = hdr[nf + axis][lv][bid]
            self.assertTrue(np.array_equal(data, points[axis]))
        # Index beyond the coordinates raises
        with self.assertRaises(IndexError):
            hdr[nf + 3][lv][bid]

    def test_mixed_field_and_coord(self):
        # Coordinates can be mixed with data fields, order is preserved
        hdr = PlotfileCooker(self.pfile3d)
        lv, bid = 1, 3
        points = hdr.box_points(lv, bid)
        data = hdr[['z', 'temp', 'x']][lv][bid]
        self.assertEqual(data.shape[-1], 3)
        self.assertTrue(np.array_equal(data[..., 0], points[2]))
        self.assertTrue(np.array_equal(data[..., 1], hdr['temp'][lv][bid]))
        self.assertTrue(np.array_equal(data[..., 2], points[0]))

    def test_multi_field_order_with_coord(self):
        # Two data fields in descending index order plus a coordinate: the
        # requested column order must be preserved (not the sorted order).
        # In this plotfile temp=5 and density=3, so the requested field
        # indices [5, 3] differ from their sorted order [3, 5].
        hdr = PlotfileCooker(self.pfile3d)
        lv, bid = 1, 3
        data = hdr[['temp', 'density', 'z']][lv][bid]
        self.assertEqual(data.shape[-1], 3)
        self.assertTrue(np.array_equal(data[..., 0], hdr['temp'][lv][bid]))
        self.assertTrue(np.array_equal(data[..., 1], hdr['density'][lv][bid]))
        self.assertTrue(np.array_equal(data[..., 2], hdr.box_points(lv, bid)[2]))

    def test_negative_box_index(self):
        # Negative box indices are supported on the coordinate path
        hdr = PlotfileCooker(self.pfile3d)
        lv = 1
        last = len(hdr.cells[lv]['indexes']) - 1
        self.assertTrue(np.array_equal(hdr['x'][lv][-1],
                                       hdr.box_points(lv, last)[0]))
        mixed = hdr[['temp', 'x']][lv][-1]
        self.assertTrue(np.array_equal(mixed[..., 0], hdr['temp'][lv][last]))
        self.assertTrue(np.array_equal(mixed[..., 1],
                                       hdr.box_points(lv, last)[0]))

    def test_iter_method_coord_not_implemented(self):
        # The on-the-fly .iter() data iterator does not support coordinates
        hdr = PlotfileCooker(self.pfile3d)
        with self.assertRaises(NotImplementedError):
            hdr[['temp', 'x']][0].iter([0])

    def test_coord_only_multi(self):
        # Multiple coordinates with no data field
        hdr = PlotfileCooker(self.pfile3d)
        lv, bid = 2, 10
        points = hdr.box_points(lv, bid)
        data = hdr[['x', 'z']][lv][bid]
        self.assertEqual(data.shape[-1], 2)
        self.assertTrue(np.array_equal(data[..., 0], points[0]))
        self.assertTrue(np.array_equal(data[..., 1], points[2]))

    def test_coord_writable(self):
        # The tiled coordinate array must be writable (for masking, etc.)
        hdr = PlotfileCooker(self.pfile3d)
        data = hdr['x'][1][0]
        data[0, 0, 0] = -999.0
        self.assertEqual(data[0, 0, 0], -999.0)

    def test_iter_read_order_alignment(self):
        # Iterating yields coordinates in the binary file read order, and
        # the data field column stays aligned with the coordinate column
        hdr = PlotfileCooker(self.pfile3d)
        for lv in range(hdr.limit_level + 1):
            read_order = [int(b) for b in hdr.binfile_read_order(lv)]
            # Coordinate only iteration
            for n, cbox in enumerate(hdr['y'][lv]):
                expected = hdr.box_points(lv, read_order[n])[1]
                self.assertTrue(np.array_equal(cbox, expected))
            # Mixed data + coordinate iteration
            for n, mbox in enumerate(hdr[['density', 'z']][lv]):
                bid = read_order[n]
                self.assertTrue(np.array_equal(mbox[..., 0],
                                               hdr['density'][lv][bid]))
                self.assertTrue(np.array_equal(mbox[..., 1],
                                               hdr.box_points(lv, bid)[2]))

    def test_box_slice_and_mask(self):
        # Slices and boolean masks of boxes return lists of combined arrays
        hdr = PlotfileCooker(self.pfile3d)
        lv = 1
        # Contiguous slice
        sliced = hdr[['density', 'y']][lv][1:4]
        for off, bid in enumerate(range(1, 4)):
            self.assertTrue(np.array_equal(sliced[off][..., 0],
                                           hdr['density'][lv][bid]))
            self.assertTrue(np.array_equal(sliced[off][..., 1],
                                           hdr.box_points(lv, bid)[1]))
        # Boolean mask
        mask = np.zeros(len(hdr.cells[lv]['indexes']), dtype=bool)
        mask[2] = True
        mask[5] = True
        masked = hdr['x'][lv][mask]
        for cbox, bid in zip(masked, np.flatnonzero(mask)):
            self.assertTrue(np.array_equal(cbox,
                                           hdr.box_points(lv, int(bid))[0]))

    def test_amr_masked_iter_with_coord(self):
        # amr_masked_iter must yield the masked coordinate data in the
        # binary file read order, with the masks staying aligned. We
        # reproduce the expected sequence and compare it to the iterator.
        hdr = PlotfileCooker(self.pfile3d)
        masks = hdr.get_amr_masks()
        expected = []
        for lv in range(hdr.limit_level + 1):
            read_order = [int(b) for b in hdr.binfile_read_order(lv)]
            for n, bid in enumerate(read_order):
                xbox = hdr.box_points(lv, bid)[0]
                if lv < hdr.max_level:
                    mask = masks[lv][n]
                    if mask.any():
                        expected.append(xbox[mask])
                else:
                    expected.append(xbox.ravel(order='F'))
        produced = list(hdr.amr_masked_iter('x'))
        self.assertEqual(len(produced), len(expected))
        for got, exp in zip(produced, expected):
            self.assertTrue(np.array_equal(got, exp))

    def test_unchanged_fullslice(self):
        # PlotfileCooker[:] must keep returning only the data fields
        hdr = PlotfileCooker(self.pfile3d)
        data = hdr[:][0][0]
        self.assertEqual(data.shape[-1], len(hdr.fields))

    def test_2d_coords(self):
        # x and y are available in 2D, z raises a clear error
        hdr = PlotfileCooker(self.pfile2d)
        data = hdr[['x', 'y']][0][0]
        self.assertEqual(data.shape[-1], 2)
        with self.assertRaises(IndexError):
            hdr['z']


class TestRegionHelpers(unittest.TestCase):
    """
    Unit tests of the functions used to place the AMR box data in a
    region slice. These do not read any plotfile data: the box data is
    replaced by the global cell indices so that the expected array is
    known exactly.
    """

    def test_region_overlap_against_numpy(self):
        # Tile the axis with boxes of a given shape at a coarser level,
        # then rebuild any slice of the fine axis from the boxes and
        # compare it with slicing the fine axis directly
        size = 48
        full = np.arange(size, dtype='float64')
        for factor, box_shape in [(1, 6), (1, 7), (2, 4), (2, 5), (4, 3)]:
            # Boxes tiling the coarse axis (the last one is truncated)
            n_coarse = size // factor
            bounds = [(lo, min(lo + box_shape, n_coarse) - 1)
                      for lo in range(0, n_coarse, box_shape)]
            for start in [0, 1, 5, 23, 47]:
                for stop in [0, 6, 24, 47, 48]:
                    for step in [1, 2, 3, 5, 7]:
                        expected = full[start:stop:step]
                        out = np.full(len(expected), np.nan)
                        for lo, hi in bounds:
                            # Box bounds on the fine axis (half open)
                            blo, bhi = lo * factor, (hi + 1) * factor
                            overlap = region_overlap(blo, bhi, start, stop,
                                                     step)
                            if overlap is None:
                                continue
                            out[overlap[0]] = full[blo:bhi][overlap[1]]
                        self.assertTrue(np.array_equal(out, expected),
                                        (factor, box_shape, start, stop, step))

    def test_region_overlap_no_selected_cell(self):
        # A box can overlap the region without holding a selected cell
        # when the step is larger than its width
        self.assertIsNone(region_overlap(1, 3, 0, 10, 4))
        self.assertIsNotNone(region_overlap(1, 3, 0, 10, 2))
        # No overlap at all
        self.assertIsNone(region_overlap(10, 20, 0, 10, 1))
        self.assertIsNone(region_overlap(0, 10, 10, 20, 1))

    def test_expand_array_nd(self):
        # The 3D expansion must match amr_kitchen.utils.expand_array3d
        arr3d = np.arange(2 * 3 * 4, dtype='float64').reshape(2, 3, 4)
        expanded = expand_array_nd(arr3d, 3)
        self.assertEqual(expanded.shape, (6, 9, 12))
        for index in np.ndindex(expanded.shape):
            coarse = tuple(i // 3 for i in index)
            self.assertEqual(expanded[index], arr3d[coarse])
        # And work in 2D as well
        arr2d = np.arange(3 * 5, dtype='float64').reshape(3, 5)
        expanded = expand_array_nd(arr2d, 2)
        self.assertEqual(expanded.shape, (6, 10))
        for index in np.ndindex(expanded.shape):
            self.assertEqual(expanded[index],
                             arr2d[tuple(i // 2 for i in index)])
        # A factor of one is a no-op
        self.assertTrue(np.array_equal(expand_array_nd(arr3d, 1), arr3d))

    def test_zoom_array_nd_shape_and_bounds(self):
        # The zoom must subdivide each cell (grid_mode) instead of
        # stretching the array between its end points
        arr = np.arange(3 * 3 * 3, dtype='float64').reshape(3, 3, 3)
        zoomed = zoom_array_nd(arr, 4)
        self.assertEqual(zoomed.shape, (12, 12, 12))
        # Linear interpolation of a monotonic array stays in range
        self.assertGreaterEqual(zoomed.min(), arr.min())
        self.assertLessEqual(zoomed.max(), arr.max())
        # A constant box stays constant everywhere: with the scipy
        # default mode zeros would be pulled in at the boundaries
        constant = zoom_array_nd(np.full((3, 4, 5), 7.0), 3)
        self.assertEqual(constant.shape, (9, 12, 15))
        self.assertTrue(np.allclose(constant, 7.0))

    def test_zoom_array_nd_nans(self):
        # NaNs must not propagate through the whole array, but every
        # cell influenced by a NaN input cell must be NaN
        arr = np.zeros((4, 4, 4))
        arr[1, 1, 1] = np.nan
        zoomed = zoom_array_nd(arr, 2)
        self.assertEqual(zoomed.shape, (8, 8, 8))
        self.assertTrue(np.isnan(zoomed).any())
        self.assertFalse(np.isnan(zoomed).all())
        # The cells of the NaN input cell are NaN
        self.assertTrue(np.isnan(zoomed[2:4, 2:4, 2:4]).all())
        # And no value was invented where the data was zero
        self.assertTrue(np.all(zoomed[~np.isnan(zoomed)] == 0.0))
        # An all NaN box stays all NaN with the expanded shape
        all_nan = zoom_array_nd(np.full((2, 3, 4), np.nan), 2)
        self.assertEqual(all_nan.shape, (4, 6, 8))
        self.assertTrue(np.isnan(all_nan).all())


class TestRegionSlice(unittest.TestCase):
    """
    Slicing a region of the data on the uniform grid at the finest level
    with PlotfileCooker[fields][x, y, z]
    """
    pfile2d = "test_assets/example_plt_2d"
    pfile3d = "test_assets/example_plt_3d"
    # 3D plotfile where the finest level only covers half the domain,
    # used to test that coarse level data fills the rest
    pfile3deb = "test_assets/plt_eb_3d"

    def uniform_array(self, pck, field):
        """
        The whole uniform grid at the finest level, assembled from the
        AMR boxes of every level from the coarsest to the finest
        """
        level = pck.limit_level
        shape = tuple(int(size) for size in pck.grid_sizes[level])
        out = np.full(shape, np.nan)
        for lv in range(level + 1):
            factor = 2 ** (level - lv)
            for bid in range(len(pck.cells[lv]['indexes'])):
                lo, hi = pck.cells[lv]['indexes'][bid]
                data = pck[field][lv][bid]
                if factor > 1:
                    data = expand_array_nd(data, factor)
                out[tuple(slice(int(lo[d]) * factor, (int(hi[d]) + 1) * factor)
                          for d in range(pck.ndims))] = data
        return out

    def test_region_matches_box_data(self):
        # The whole domain must match the data assembled box by box, and
        # both the Ellipsis and the full slices select it
        pck = PlotfileCooker(self.pfile3d, serial=True)
        expected = self.uniform_array(pck, 'temp')
        self.assertEqual(expected.shape,
                         tuple(pck.grid_sizes[pck.limit_level]))
        for key in [(slice(None), slice(None), slice(None)), Ellipsis]:
            data = pck['temp'][key]
            self.assertEqual(data.shape, expected.shape)
            self.assertTrue(np.array_equal(data, expected))

    def test_region_slices_like_numpy(self):
        # Every supported index must select the same cells as numpy
        # indexing of the whole uniform array
        pck = PlotfileCooker(self.pfile3d, serial=True)
        expected = self.uniform_array(pck, 'temp')
        keys = [(slice(1, 5), slice(1, 5), slice(1, 5)),
                (slice(None, None, 2), slice(3, None), slice(None, 7)),
                (slice(3, None, 3), slice(None, None, 5), slice(2, 30, 4)),
                (slice(-8, None), slice(None, -8), slice(4, 6)),
                (slice(5, 2), slice(None), slice(None)),
                (slice(None), Ellipsis),
                (Ellipsis, slice(0, 4))]
        for key in keys:
            data = pck['temp'][key]
            self.assertEqual(data.shape, expected[key].shape, key)
            self.assertTrue(np.array_equal(data, expected[key]), key)

    def test_region_integer_index_drops_axis(self):
        # An axis indexed with an integer is dropped, like numpy does
        pck = PlotfileCooker(self.pfile3d, serial=True)
        expected = self.uniform_array(pck, 'temp')
        for key in [(slice(None), slice(None), 7),
                    (3, slice(None), slice(None)),
                    (-1, slice(-4, None), slice(None, None, 2)),
                    (2, 3, slice(None)),
                    (2, 3, 4)]:
            data = pck['temp'][key]
            self.assertEqual(data.shape, expected[key].shape, key)
            self.assertTrue(np.array_equal(data, expected[key]), key)

    def test_region_negative_step(self):
        # Negative steps reverse the axis
        pck = PlotfileCooker(self.pfile3d, serial=True)
        expected = self.uniform_array(pck, 'temp')
        for key in [(slice(None, None, -1), slice(None), slice(None)),
                    (slice(None, None, -1),) * 3,
                    (slice(20, 4, -3), slice(None, None, -2), slice(None))]:
            data = pck['temp'][key]
            self.assertEqual(data.shape, expected[key].shape, key)
            self.assertTrue(np.array_equal(data, expected[key]), key)

    def test_region_field_axis(self):
        # A single field drops the trailing field axis, multiple fields
        # keep it in the requested order
        pck = PlotfileCooker(self.pfile3d, serial=True)
        key = (slice(0, 8), slice(4, 12), slice(None, None, 4))
        temp = pck['temp'][key]
        self.assertEqual(temp.ndim, 3)
        density = pck['density'][key]
        # The requested order is preserved (temp is field 5, density is 3)
        both = pck[['temp', 'density']][key]
        self.assertEqual(both.shape, temp.shape + (2,))
        self.assertTrue(np.array_equal(both[..., 0], temp))
        self.assertTrue(np.array_equal(both[..., 1], density))
        # A single field in a list keeps the field axis
        single = pck[['temp']][key]
        self.assertEqual(single.shape, temp.shape + (1,))
        self.assertTrue(np.array_equal(single[..., 0], temp))
        # A field slice selects the fields of the slice
        sliced = pck[3:6][key]
        self.assertEqual(sliced.shape, temp.shape + (3,))
        self.assertTrue(np.array_equal(sliced[..., 2], temp))
        # A single field index behaves like the field name
        self.assertTrue(np.array_equal(pck[5][key], temp))

    def test_region_coords(self):
        # The 'x', 'y' and 'z' virtual fields hold the cell center
        # coordinates of the finest level grid
        pck = PlotfileCooker(self.pfile3d, serial=True)
        level = pck.limit_level
        key = (slice(2, 20, 3), slice(None, None, -2), 5)
        for axis, name in enumerate(['x', 'y', 'z']):
            data = pck[name][key]
            grid = pck.grids[level][axis]
            # Global indices of the cells selected along the axis
            indices = np.atleast_1d(
                np.arange(pck.grid_sizes[level][axis])[key[axis]])
            if axis == 2:
                # The z axis is indexed with an integer so it was dropped,
                # the coordinate is constant over the two other axes
                self.assertEqual(data.ndim, 2)
                self.assertTrue(np.allclose(data, grid[indices[0]]))
            else:
                # The coordinate only varies along its own axis
                self.assertEqual(data.shape[axis], len(indices))
                for n, index in enumerate(indices):
                    self.assertTrue(np.allclose(np.moveaxis(data, axis, 0)[n],
                                                grid[index]))
        # Coordinates can be mixed with data fields, in order
        mixed = pck[['temp', 'z']][key]
        self.assertEqual(mixed.shape[-1], 2)
        self.assertTrue(np.array_equal(mixed[..., 0], pck['temp'][key]))
        self.assertTrue(np.allclose(mixed[..., 1], pck['z'][key]))

    def test_region_coarse_level_fill_3d(self):
        # Cells which the finest level does not cover must hold the value
        # of the coarse cell containing them. The expected value is read
        # from the level 0 boxes without using the expansion function.
        pck = PlotfileCooker(self.pfile3deb, serial=True)
        self.assertEqual(pck.limit_level, 1)
        region = pck['density'][...]
        shape = region.shape
        # Cells covered by the finest level
        covered = np.zeros(shape, dtype=bool)
        for lo, hi in pck.cells[1]['indexes']:
            covered[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1] = True
        # The test is only meaningful if the coverage is partial
        self.assertTrue(0 < covered.mean() < 1)
        # Level 0 data on its own grid
        coarse = np.full(tuple(int(s) for s in pck.grid_sizes[0]), np.nan)
        for bid, (lo, hi) in enumerate(pck.cells[0]['indexes']):
            coarse[lo[0]:hi[0] + 1,
                   lo[1]:hi[1] + 1,
                   lo[2]:hi[2] + 1] = pck['density'][0][bid]
        # Coarse cell holding each cell of the finest level grid
        grids = np.meshgrid(*[np.arange(size) // 2 for size in shape],
                            indexing='ij')
        expected = coarse[grids[0], grids[1], grids[2]]
        self.assertFalse(np.isnan(expected).any())
        self.assertTrue(np.array_equal(region[~covered], expected[~covered]))
        # And the covered cells hold the finest level data
        self.assertFalse(np.isnan(region).any())

    def test_region_2d(self):
        # Region slices of a 2D plotfile take two indices
        pck = PlotfileCooker(self.pfile2d, serial=True)
        expected = self.uniform_array(pck, 'temp')
        for key in [(slice(None), slice(None)),
                    Ellipsis,
                    (slice(100, 140, 3), slice(None, None, -4)),
                    (64, slice(None))]:
            data = pck['temp'][key]
            self.assertEqual(data.shape, expected[key].shape, key)
            self.assertTrue(np.array_equal(data, expected[key]), key)
        # Three indices are not valid for a 2D plotfile
        with self.assertRaises(IndexError):
            pck['temp'][0:4, 0:4, 0:4]

    def test_region_limit_level(self):
        # The region is sliced on the grid of the finest level which was
        # read from the plotfile, not of the finest level of the plotfile
        pck = PlotfileCooker(self.pfile3d, limit_level=1, serial=True)
        data = pck['temp'][...]
        self.assertEqual(data.shape, tuple(pck.grid_sizes[1]))
        self.assertTrue(np.array_equal(data, self.uniform_array(pck, 'temp')))
        # Indices are bound to that grid
        with self.assertRaises(IndexError):
            pck['temp'][20, 0, 0]

    def test_region_dtype(self):
        # The data type of the plotfile is preserved, except when a
        # coordinate is requested (the grids are stored in float64)
        pck = PlotfileCooker(self.pfile3d, serial=True)
        box_dtype = pck['temp'][pck.limit_level][0].dtype
        self.assertEqual(pck['temp'][0:4, 0:4, 0:4].dtype, box_dtype)
        self.assertEqual(pck[['temp', 'x']][0:4, 0:4, 0:4].dtype,
                         np.dtype('float64'))
        self.assertEqual(pck['y'][0:4, 0:4, 0:4].dtype, np.dtype('float64'))

    def test_region_empty_selection(self):
        # An empty selection returns an empty array of the right shape
        pck = PlotfileCooker(self.pfile3d, serial=True)
        self.assertEqual(pck['temp'][4:4, :, :].shape, (0, 32, 32))
        self.assertEqual(pck[['temp', 'density']][8:2, :, 0:4].shape,
                         (0, 32, 4, 2))

    def test_region_step_skips_unread_boxes(self):
        # A box overlapping the region can hold no selected cell when the
        # step is larger than the box: it must not be read at all
        pck = PlotfileCooker(self.pfile3d, serial=True)
        original = plotfile_cooker.mp_read_box_single_field
        read_boxes = []

        def counting_read(args):
            read_boxes.append(args[1])
            return original(args)

        try:
            plotfile_cooker.mp_read_box_single_field = counting_read
            pck['temp'][:, :, :]
            n_whole_domain = len(read_boxes)
            read_boxes.clear()
            pck['temp'][0:32:16, :, :]
            n_stepped = len(read_boxes)
        finally:
            plotfile_cooker.mp_read_box_single_field = original
        # The whole domain reads every box of every level
        self.assertEqual(n_whole_domain,
                         sum(len(pck.cells[lv]['indexes'])
                             for lv in range(pck.limit_level + 1)))
        # Only the boxes holding one of the two selected x planes
        expected = 0
        for lv in range(pck.limit_level + 1):
            factor = 2 ** (pck.limit_level - lv)
            for lo, hi in pck.cells[lv]['indexes']:
                covered = range(int(lo[0]) * factor, (int(hi[0]) + 1) * factor)
                if any(plane in covered for plane in range(0, 32, 16)):
                    expected += 1
        self.assertLess(expected, n_whole_domain)
        self.assertEqual(n_stepped, expected)

    def test_region_serial_and_parallel(self):
        # Both the serial and the multiprocessing box readers must give
        # the same result
        key = (slice(4, 20), slice(None), 3)
        serial = PlotfileCooker(self.pfile3d, serial=True)['temp'][key]
        parallel = PlotfileCooker(self.pfile3d)['temp'][key]
        self.assertTrue(np.array_equal(serial, parallel))

    def test_region_level_selection_unchanged(self):
        # A single integer still selects an AMR level
        pck = PlotfileCooker(self.pfile3d, serial=True)
        self.assertTrue(np.array_equal(pck['temp'][1][0],
                                       pck['temp'][1][0]))
        self.assertEqual(pck['temp'][-1][0].ndim, 3)
        with self.assertRaises(ValueError):
            pck['temp'][10]
        # A bare slice is not a valid AMR level
        with self.assertRaises(TypeError):
            pck['temp'][0:2]

    def test_region_index_errors(self):
        pck = PlotfileCooker(self.pfile3d, serial=True)
        # Wrong number of indices
        with self.assertRaises(IndexError):
            pck['temp'][0:4, 0:4]
        with self.assertRaises(IndexError):
            pck['temp'][0:4, 0:4, 0:4, 0:4]
        # More than one Ellipsis
        with self.assertRaises(IndexError):
            pck['temp'][..., 0:4, ...]
        # Out of bounds integer index
        with self.assertRaises(IndexError):
            pck['temp'][32, 0, 0]
        with self.assertRaises(IndexError):
            pck['temp'][0, -33, 0]
        # Unsupported index type
        with self.assertRaises(TypeError):
            pck['temp'][[0, 1], 0:4, 0:4]

    def test_region_size_warning(self):
        # A warning is emitted when the output array is larger than
        # REGION_SIZE_WARNING
        pck = PlotfileCooker(self.pfile3d, serial=True)
        threshold = plotfile_cooker.REGION_SIZE_WARNING
        try:
            plotfile_cooker.REGION_SIZE_WARNING = 0
            with self.assertWarns(UserWarning):
                pck['temp'][0:2, 0:2, 0:2]
        finally:
            plotfile_cooker.REGION_SIZE_WARNING = threshold
        # And none below the threshold
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            pck['temp'][0:2, 0:2, 0:2]


class TestInterpMode(unittest.TestCase):
    """
    PlotfileCooker.set_interp_mode, defining how coarse AMR level data is
    expanded onto the finest level grid in region slices
    """
    pfile3deb = "test_assets/plt_eb_3d"

    def test_default_mode_is_repeat(self):
        pck = PlotfileCooker(self.pfile3deb, serial=True)
        self.assertIs(pck.interp_mode, expand_array_nd)
        repeated = pck['density'][...]
        pck.set_interp_mode('repeat')
        self.assertIs(pck.interp_mode, expand_array_nd)
        self.assertTrue(np.array_equal(pck['density'][...], repeated))

    def test_zoom_mode(self):
        # The interpolated data must stay within the range of the
        # plotfile data and differ from the repeated data
        pck = PlotfileCooker(self.pfile3deb, serial=True)
        repeated = pck['density'][...]
        pck.set_interp_mode('zoom')
        self.assertIs(pck.interp_mode, zoom_array_nd)
        zoomed = pck['density'][...]
        self.assertEqual(zoomed.shape, repeated.shape)
        self.assertFalse(np.isnan(zoomed).any())
        self.assertGreaterEqual(zoomed.min(), repeated.min())
        self.assertLessEqual(zoomed.max(), repeated.max())
        # The cells covered by the finest level are read directly so they
        # are identical, the coarse cells are interpolated
        covered = np.zeros(zoomed.shape, dtype=bool)
        for lo, hi in pck.cells[1]['indexes']:
            covered[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1] = True
        self.assertTrue(np.array_equal(zoomed[covered], repeated[covered]))
        self.assertFalse(np.array_equal(zoomed[~covered], repeated[~covered]))

    def test_custom_mode(self):
        # A user supplied function is called once per field and per box
        # below the finest level
        pck = PlotfileCooker(self.pfile3deb, serial=True)
        calls = []

        def expand(arr, factor):
            calls.append((arr.shape, factor))
            return expand_array_nd(arr, factor)

        pck.set_interp_mode(expand)
        self.assertIs(pck.interp_mode, expand)
        custom = pck['density'][...]
        self.assertTrue(len(calls) > 0)
        self.assertTrue(all(factor == 2 for _, factor in calls))
        pck.set_interp_mode('repeat')
        self.assertTrue(np.array_equal(custom, pck['density'][...]))

    def test_custom_mode_bad_shape(self):
        # An expansion function returning the wrong shape is reported
        pck = PlotfileCooker(self.pfile3deb, serial=True)
        pck.set_interp_mode(lambda arr, factor: arr)
        with self.assertRaises(ValueError):
            pck['density'][...]

    def test_unknown_mode(self):
        pck = PlotfileCooker(self.pfile3deb, serial=True)
        with self.assertRaises(ValueError):
            pck.set_interp_mode('linear')

    def test_mode_not_used_at_finest_level(self):
        # The finest level data is read as is, so a region of the finest
        # level only does not call the expansion function
        pck = PlotfileCooker(self.pfile3deb, limit_level=0, serial=True)

        def fail(arr, factor):
            raise AssertionError("the expansion function was called")

        pck.set_interp_mode(fail)
        data = pck['density'][...]
        self.assertEqual(data.shape, tuple(pck.grid_sizes[0]))
