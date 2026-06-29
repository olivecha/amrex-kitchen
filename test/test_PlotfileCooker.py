import os
import unittest
import numpy as np

from amr_kitchen import PlotfileCooker, DTYPE

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

