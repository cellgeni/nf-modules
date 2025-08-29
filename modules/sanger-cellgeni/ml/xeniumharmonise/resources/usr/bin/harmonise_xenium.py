#!/usr/bin/env python3

import pyarrow.lib as _lib
import uuid
import tiledb
import numpy as np
import os
# from dask.distributed import Client
import dask
import dask.array as da
import dask.dataframe as dd
import tifffile
import zarr
import json
import logging
import datetime
import pandas as pd
import palom
from itertools import repeat
from tqdm.contrib.concurrent import thread_map
import sqlite3
from pprint import pprint
import scipy.sparse as sp
import fire

# default run configuration

DASK_SCHEDULER_ADDRESS='tcp://farm22-head2:40883'
CONCURRENCY=4
DATA_ROOT='/lustre/scratch126/cellgen/team361/projects/histology_to_gene_expression/workspace/data'

#def easyloader(sample_bundle_uri: str, harmonised_dataset_uri: str):
    #"""
    ###Automatically identify the xenium_bundle_uri and hematoxylin_eosin_image_uri in the sample_bundle_uri
    ##and call the main function with these parameters.
    #"""
    # find the xenium bundle and h&e image
    #lfiles = os.listdir(sample_bundle_uri)
    #xenium_bundle_uri = os.path.join(sample_bundle_uri, [i for i in lfiles if 'output-X' in i][0])
    #hematoxylin_eosin_image_uri = os.path.join(sample_bundle_uri, [i for i in lfiles if '.ndpi' in i][0])

    # call the main function
    #main(xenium_bundle_uri, hematoxylin_eosin_image_uri, harmonised_dataset_uri)


def main(xenium_bundle_uri: str, hematoxylin_eosin_image_uri: str, harmonised_dataset_uri: str):
    """Harmonise a xenium dataset with a h&e image

    Generate a harmonised store for a xenium dataset with a h&e image. The h&e image is assumed to be
    generated from the same slide as the xenium dataset, and must have been converted to `.ndpi`.While we end up
    re-implementing logic in the `spatialdata_io.xenium` library, this gives us much more flexibility. In particular,
    spatialdata only supports zarr as a storage engine, and therefore does not have good support for sparse arrays.
    """

    # create a tiledb group

    tiledb.group_create(harmonised_dataset_uri)

    # define function for chunked read and write

    def empty_dask_array_from_tiff(tiff_uri, series=None, chunks=None):
        """This is a utility function to create an empty dask array with the correct
        dimensions and shape, such that `read_block_from_tiff` can be run with `map_blocks`
        to lazily load chunks from the tiff file."""

        # suppress warning from tifffile
        # https://github.com/scverse/spatialdata-io/blob/7451288357b9d74886989f86e8d1166d386630e4/src/spatialdata_io/readers/xenium.py#L313
        class IgnoreSpecificMessage(logging.Filter):
            def filter(self, record: logging.LogRecord) -> bool:
                if "OME series cannot read multi-file pyramids" in record.getMessage():
                    return False
                return True

        logger = tifffile.logger()
        logger.addFilter(IgnoreSpecificMessage())

        tifffile_store = tifffile.imread(tiff_uri, aszarr=True)
        zarr_store = zarr.open(tifffile_store, 'r')

        if series is None:
            zarr_array = zarr_store
        else:
            zarr_array = zarr_store[series]

        if chunks is None:
            chunks = zarr_array.chunks

        dask_array = da.empty(
            shape=zarr_array.shape,
            chunks=chunks,
            dtype=zarr_array.dtype
        )
        return dask_array

    def read_block_from_tiff(block, tiff_uri, series=None, block_info=None):

        # suppress warning from tifffile
        # https://github.com/scverse/spatialdata-io/blob/7451288357b9d74886989f86e8d1166d386630e4/src/spatialdata_io/readers/xenium.py#L313
        class IgnoreSpecificMessage(logging.Filter):
            def filter(self, record: logging.LogRecord) -> bool:
                if "OME series cannot read multi-file pyramids" in record.getMessage():
                    return False
                return True

        logger = tifffile.logger()
        logger.addFilter(IgnoreSpecificMessage())

        if block_info is None:
            return block
        slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
        tifffile_store = tifffile.imread(tiff_uri, aszarr=True)
        zarr_store = zarr.open(tifffile_store, 'r')
        if series is not None:
            zarr_array = zarr_store[series]
        else:
            zarr_array = zarr_store
        return zarr_array[slices]

    def write_block_to_tiledb(block, tiledb_uri, attribute="px_value", block_info=None): # required due to a bug in dask.array.to_tiledb()
        if block_info is None:
            return dask.array.array([[[0]]])
        slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
        with tiledb.DenseArray(tiledb_uri, mode='w') as tiledb_array:
            tiledb_array[slices] = {attribute: block}
        return dask.array.array([[[0]]])

    def write_block_multiple_attributes_to_tiledb(block, tiledb_uri, attribute_dimension, attributes, block_info=None):
        """Unpack array to attributes on a given axis; expected to act over all attributes"""
        if block_info is None:
            return dask.array.array([[[0]]])
        slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
        value_slices = tuple([value_slice for dim, value_slice in enumerate(slices) if dim != attribute_dimension])
        if len(attributes) != block.shape[attribute_dimension]:
            raise ValueError('incorrect number of attributes')
        unstacked_block = np.moveaxis(block, attribute_dimension, 0)
        block_dict = {attributes[i]: unstacked_block[i] for i in range(len(attributes))}
        with tiledb.DenseArray(tiledb_uri, mode='w') as tiledb_array:
            tiledb_array[value_slices] = block_dict
        return dask.array.array([[[0]]])

    # tiledb.libtiledb.TileDBError: [TileDB::C++API] Error: Failed to create context
    # occurs with `to_tiledb()`

    # load the immunofluorescence image shape and chunk size

    print('processing immunoflourescence')

    immunofluorescence_raw_uri = os.path.join(xenium_bundle_uri, 'morphology.ome.tif')
    immunofluorescence = empty_dask_array_from_tiff(immunofluorescence_raw_uri, series='0', chunks=(1, 16384, 16384))

    # create an array for the immunofluorescence z stack

    immunofluorescence_uri = os.path.join(harmonised_dataset_uri, 'immunofluorescence')

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="z_position", domain=(0, immunofluorescence.shape[0] - 1), tile=1, dtype=np.int32),
            tiledb.Dim(name="px_rows", domain=(0, immunofluorescence.shape[1] - 1), tile=4096, dtype=np.int32),
            tiledb.Dim(name="px_cols", domain=(0, immunofluorescence.shape[2] - 1), tile=4096, dtype=np.int32)
    ),
        attrs=[
            tiledb.Attr(name="px_value", dtype=immunofluorescence.dtype, filters=[tiledb.ZstdFilter()])
        ]
    )
    tiledb.Array.create(immunofluorescence_uri, schema)

    # add data to the array

    (immunofluorescence
     .map_blocks(read_block_from_tiff, immunofluorescence_raw_uri, series=0)
     .map_blocks(write_block_to_tiledb, immunofluorescence_uri)
     .compute()) #FIXME this now seems to work with tiledb to_tiledb, so rewrite and remove the custom function

    # load the immunofluorescence autofocus projection images

    print('processing immunoflourescence focus')

    immunofluorescence_focus_filenames = sorted(os.listdir(os.path.join(xenium_bundle_uri, 'morphology_focus')))
    immunofluorescence_focus_raw_uri = os.path.join(xenium_bundle_uri, 'morphology_focus', immunofluorescence_focus_filenames[0])

    immunofluorescence_focus = empty_dask_array_from_tiff(immunofluorescence_focus_raw_uri, chunks=(1, 16384, 16384))

    # create an array for the immunofluorescence autofocus projection images

    immunofluorescence_focus_uri = os.path.join(harmonised_dataset_uri, 'immunofluorescence_focus')

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="stain", domain=(0, immunofluorescence_focus.shape[0] - 1), tile=1, dtype=np.int32),
            tiledb.Dim(name="px_rows", domain=(0, immunofluorescence_focus.shape[1] - 1), tile=4096, dtype=np.int32),
            tiledb.Dim(name="px_cols", domain=(0, immunofluorescence_focus.shape[2] - 1), tile=4096, dtype=np.int32)
        ),
        attrs=[
            tiledb.Attr(name="px_value", dtype=immunofluorescence_focus.dtype, filters=[tiledb.ZstdFilter()])
        ]
    )
    tiledb.Array.create(immunofluorescence_focus_uri, schema)

    # add data to the array

    (immunofluorescence_focus
     .map_blocks(read_block_from_tiff, immunofluorescence_focus_raw_uri)
     .map_blocks(write_block_to_tiledb, immunofluorescence_focus_uri)
     .compute())

    # load the cell and nucleus segmentation masks

    print('processing segmentation masks')

    cells_zarr_group = zarr.open_group(store=zarr.storage.ZipStore(os.path.join(xenium_bundle_uri, 'cells.zarr.zip'), mode='r'))

    cell_segmentation_masks = da.from_zarr(cells_zarr_group['masks'][1])
    nucleus_segmentation_masks = da.from_zarr(cells_zarr_group['masks'][0])

    # create arrays for the cell and nucleus segmentation masks

    cell_segmentation_mask_uri = os.path.join(harmonised_dataset_uri, 'cell_segmentation_mask')
    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="y", domain=(0, cell_segmentation_masks.shape[0] - 1), tile=4096, dtype=np.int32),
            tiledb.Dim(name="x", domain=(0, cell_segmentation_masks.shape[1] - 1), tile=4096, dtype=np.int32)
        ),
        attrs=[
            tiledb.Attr(name="instance_label_id", dtype=cell_segmentation_masks.dtype, filters=[tiledb.ZstdFilter()])
        ]
    )
    tiledb.Array.create(cell_segmentation_mask_uri, schema)

    nucleus_segmentation_mask_uri = os.path.join(harmonised_dataset_uri, 'nuclei_segmentation_mask')
    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="y", domain=(0, nucleus_segmentation_masks.shape[0] - 1), tile=4096, dtype=np.int32),
            tiledb.Dim(name="x", domain=(0, nucleus_segmentation_masks.shape[1] - 1), tile=4096, dtype=np.int32)
        ),
        attrs=[
            tiledb.Attr(name="instance_label_id", dtype=nucleus_segmentation_masks.dtype, filters=[tiledb.ZstdFilter()])
        ]
    )
    tiledb.Array.create(nucleus_segmentation_mask_uri, schema)

    # add data to the array

    cell_segmentation_masks.rechunk(chunks=(16384, 16384)).map_blocks(write_block_to_tiledb, cell_segmentation_mask_uri, "instance_label_id").compute()
    nucleus_segmentation_masks.rechunk(chunks=(16384, 16384)).map_blocks(write_block_to_tiledb, nucleus_segmentation_mask_uri, "instance_label_id").compute()

    # load the cell and nucleus segmentation information

    print('processing boundaries')

    cell_boundaries = dd.read_parquet(os.path.join(xenium_bundle_uri, 'cell_boundaries.parquet'))[['cell_id', 'label_id']].drop_duplicates()
    nucleus_boundaries = dd.read_parquet(os.path.join(xenium_bundle_uri, 'nucleus_boundaries.parquet'))[['cell_id', 'label_id']].drop_duplicates()

    # create arrays for the cell and nucleus segmentation information

    cell_segmentation_uri = os.path.join(harmonised_dataset_uri, 'cell_segmentation')
    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="cell_id", domain=(None, None), tile=64, dtype='ascii'),
        ),
        attrs=[
            tiledb.Attr(name="mask_instance_label_id", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
        ],
        sparse=True,
        allows_duplicates=False # each cell can only have a single segmentation
    )
    tiledb.Array.create(cell_segmentation_uri, schema)

    nucleus_segmentation_uri = os.path.join(harmonised_dataset_uri, 'nucleus_segmentation')
    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="cell_id", domain=(None, None), tile=64, dtype='ascii'),
        ),
        attrs=[
            tiledb.Attr(name="mask_instance_label_id", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
        ],
        sparse=True,
        allows_duplicates=True # each cell can have multiple nucleus segmentations
    )
    tiledb.Array.create(nucleus_segmentation_uri, schema)

    # write the data

    def write_partition(partition, tiledb_uri, partition_info=None):
        if partition_info is None:
            return 1
        with tiledb.open(tiledb_uri, 'w') as tiledb_array:
            tiledb_array[
                partition['cell_id'],
            ] = {
                'mask_instance_label_id': partition['label_id'],
            }
        return 1

    cell_boundaries.map_partitions(write_partition, cell_segmentation_uri).compute()
    nucleus_boundaries.map_partitions(write_partition, nucleus_segmentation_uri).compute()

    # load the cell summary information

    print('processing cell summary')

    cell_summary = dd.read_parquet(os.path.join(xenium_bundle_uri, 'cells.parquet')).reset_index()

    # get the domain size from the immunofluorescence

    immunofluorescence_resolution_um_per_px = 0.2125
    y_domain = (0, immunofluorescence_focus.shape[1] * immunofluorescence_resolution_um_per_px)
    x_domain = (0, immunofluorescence_focus.shape[2] * immunofluorescence_resolution_um_per_px)

    # create an array for the cell summary information

    cell_summary_uri = os.path.join(harmonised_dataset_uri, 'cell_summary')
    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="cell_index", domain=(0, cell_summary['index'].max().compute()), tile=64, dtype=np.int64),
            tiledb.Dim(name="cell_id", domain=(None, None), tile=64, dtype='ascii'),
            tiledb.Dim(name="y_centroid", domain=y_domain, tile=256, dtype=np.float64),
            tiledb.Dim(name="x_centroid", domain=x_domain, tile=256, dtype=np.float64),
        ),
        attrs=[
            tiledb.Attr(name="transcript_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="control_probe_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="genomic_control_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="control_codeword_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="unassigned_codeword_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="deprecated_codeword_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="total_counts", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="cell_area", dtype=np.float64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="nucleus_area", dtype=np.float64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="nucleus_count", dtype=np.int64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="segmentation_method", dtype=np.dtype('U'), filters=[tiledb.ZstdFilter()]),
        ],
        sparse=True,
        allows_duplicates=False
    )
    tiledb.Array.create(cell_summary_uri, schema)

    # add data to the array

    def write_partition(partition, tiledb_uri, partition_info=None):
        if partition_info is None:
            return 1
        with tiledb.open(tiledb_uri, 'w') as tiledb_array:
            tiledb_array[
                partition['index'],
                partition['cell_id'],
                partition['y_centroid'],
                partition['x_centroid']
            ] = {
                'transcript_counts': partition['transcript_counts'],
                'control_probe_counts': partition['control_probe_counts'],
                'genomic_control_counts': partition['genomic_control_counts'],
                'control_codeword_counts': partition['control_codeword_counts'],
                'unassigned_codeword_counts': partition['unassigned_codeword_counts'],
                'deprecated_codeword_counts': partition['deprecated_codeword_counts'],
                'total_counts': partition['total_counts'],
                'cell_area': partition['cell_area'],
                'nucleus_area': partition['nucleus_area'],
                'nucleus_count': partition['nucleus_count'],
                'segmentation_method': partition['segmentation_method'],
            }
        return 1

    cell_summary.map_partitions(write_partition, cell_summary_uri).compute()

    # import cell feature matrix

    print('processing cell features')

    cell_features_group = zarr.open_group(store=zarr.storage.ZipStore(os.path.join(xenium_bundle_uri, 'cell_feature_matrix.zarr.zip'), mode='r'))

    feature_ids = cell_features_group['cell_features'].attrs['feature_ids'] #TODO record this information
    feature_keys = cell_features_group['cell_features'].attrs['feature_keys'] #TODO check genes are the same as in the transcripts table
    feature_types = cell_features_group['cell_features'].attrs['feature_types']
    gene_indices = [x == 'gene' for x in feature_types]

    cell_features_cell_id = da.from_zarr(cell_features_group['cell_features']['cell_id']) # identical to `cells_cell_id`

    cell_features_data = da.from_zarr(cell_features_group['cell_features']['data'])
    cell_features_indices = da.from_zarr(cell_features_group['cell_features']['indices'])
    cell_features_indptr = da.from_zarr(cell_features_group['cell_features']['indptr'])

    # create a sparse array to store the data

    cell_features_uri = os.path.join(harmonised_dataset_uri, 'cell_features')

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="gene_index", domain=(0, 5000), tile=64, dtype=np.int32),
            tiledb.Dim(name="cell_index", domain=(0, cell_features_cell_id.shape[0] - 1), tile=64, dtype=np.int32),
        ),
        attrs=[
            tiledb.Attr(name="count", dtype=np.uint32, filters=[tiledb.ZstdFilter()]),
        ],
        sparse=True,
        allows_duplicates=False
    )
    tiledb.Array.create(cell_features_uri, schema)

    # add data to the array

    @dask.delayed()
    def construct_csr_matrix(data, indices, indptr, shape):
        return sp.csr_matrix((data, indices, indptr), shape=shape)

    def write_csr_matrix_to_tiledb(block, tiledb_uri, block_info=None):
        if block_info is None:
            return dask.array.array([[0]])

        block_origin = [start for start, stop in block_info[None]['array-location']]

        # block is a scipy.csr_matrix

        block_coo = block.tocoo()

        # write with a block offset

        with tiledb.open(tiledb_uri, 'w') as tiledb_array:
            tiledb_array[block_coo.row + block_origin[0], block_coo.col + block_origin[1]] = block_coo.data

        return dask.array.array([[0]])

    cell_features_shape = (len(feature_ids), cell_features_cell_id.shape[0])
    delayed_csr_matrix = construct_csr_matrix(cell_features_data, cell_features_indices, cell_features_indptr, cell_features_shape)
    dask_sparse_array = da.from_delayed(delayed_csr_matrix, shape=cell_features_shape, dtype=cell_features_data.dtype)
    dask_sparse_array_genes_only = dask_sparse_array[gene_indices, :]
    dask_sparse_array_genes_only.rechunk((1024, -1)).map_blocks(write_csr_matrix_to_tiledb, cell_features_uri).compute()

    # load the h&e image

    print('processing he image')

    he_image = empty_dask_array_from_tiff(hematoxylin_eosin_image_uri, series='0', chunks=(16384, 16384, 3))

    # create an array for the h&e image

    he_image_uri = os.path.join(harmonised_dataset_uri, 'hematoxylin_eosin_image')

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="channel", domain=(0, he_image.shape[2] - 1), tile=1, dtype=np.int32),
            tiledb.Dim(name="px_rows", domain=(0, he_image.shape[0] - 1), tile=4096, dtype=np.int32),
            tiledb.Dim(name="px_cols", domain=(0, he_image.shape[1] - 1), tile=4096, dtype=np.int32)
        ),
        attrs=[
            tiledb.Attr(name="px_value", dtype=he_image.dtype, filters=[tiledb.ZstdFilter()])
        ]
    )
    tiledb.Array.create(he_image_uri, schema)

    # add data to the array

    (he_image
     .map_blocks(read_block_from_tiff, hematoxylin_eosin_image_uri, series='0')
     .transpose([2, 0, 1])
     .map_blocks(write_block_to_tiledb, he_image_uri)
     .compute())

    # align the h&e image

    print('processing aligned he image')

    immunofluorescence_image_dapi = da.from_tiledb(immunofluorescence_focus_uri, attribute="px_value", chunks=(1, 8192, 8192))[0, :, :]
    he_image_green_channel = da.from_tiledb(he_image_uri, attribute="px_value", chunks=(1, 8192, 8192))[1, :, :]
    he_image = da.from_tiledb(he_image_uri, attribute="px_value", chunks=(1, 8192, 8192))[:, :, :]

    THUMBNAIL = 5

    immunofluorescence_image_dapi_thumbnail = immunofluorescence_image_dapi[::2 ** THUMBNAIL, ::2 ** THUMBNAIL]
    he_image_green_channel_thumbnail = he_image_green_channel[::2 ** THUMBNAIL, ::2 ** THUMBNAIL]

    palom_aligner = palom.align.Aligner(
        ref_img=immunofluorescence_image_dapi,
        moving_img=he_image_green_channel,
        ref_thumbnail=immunofluorescence_image_dapi_thumbnail,
        moving_thumbnail=he_image_green_channel_thumbnail,
        ref_thumbnail_down_factor=2 ** THUMBNAIL,
        moving_thumbnail_down_factor=2 ** THUMBNAIL
    )

    palom_aligner.coarse_register_affine(n_keypoints=4000)
    palom_aligner.compute_shifts() # this step has a high memory requirement in a 'finalize' dask job on large input images
    palom_aligner.constrain_shifts()

    aligned_he_image = palom.align.block_affine_transformed_moving_img(
        ref_img=immunofluorescence_image_dapi,
        moving_img=he_image,
        mxs=palom_aligner.block_affine_matrices_da
    )

    # create an array for the aligned h&e image

    aligned_he_image_uri = os.path.join(harmonised_dataset_uri, 'hematoxylin_eosin_image_aligned')

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="channel", domain=(0, aligned_he_image.shape[0] - 1), tile=1, dtype=np.int32),
            tiledb.Dim(name="px_rows", domain=(0, aligned_he_image.shape[1] - 1), tile=4096, dtype=np.int32),
            tiledb.Dim(name="px_cols", domain=(0, aligned_he_image.shape[2] - 1), tile=4096, dtype=np.int32)
        ),
        attrs=[
            tiledb.Attr(name="px_value", dtype=aligned_he_image.dtype, filters=[tiledb.ZstdFilter()])
        ]
    )
    tiledb.Array.create(aligned_he_image_uri, schema)

    # add data to the array

    (aligned_he_image
     .rechunk(chunks=(3, 16384, 16384))
     .map_blocks(write_block_to_tiledb, aligned_he_image_uri)
     .compute())

    # load the transcripts

    print('processing transcripts')

    # get the domain size for transcripts from the immunofluorescence
    immunofluorescence_resolution_um_per_px = 0.2125
    y_domain = (0, immunofluorescence_focus.shape[1] * immunofluorescence_resolution_um_per_px)
    x_domain = (0, immunofluorescence_focus.shape[2] * immunofluorescence_resolution_um_per_px)

    transcripts = (dd.read_parquet(os.path.join(xenium_bundle_uri, 'transcripts.parquet'))
                   .query(f'y_location > {y_domain[0]}')
                   .query(f'y_location < {y_domain[1]}')
                   .query(f'x_location > {x_domain[0]}')
                   .query(f'x_location < {x_domain[1]}')) # there's a strange issue where transcript positions can be negative or occur just outside of the imaging domain

    # load the gene panel information

    #transcripts_zarr_store = zarr.open(os.path.join(xenium_bundle_uri, 'transcripts.zarr.zip'), 'r')
    #transcripts_attrs = dict(transcripts_zarr_store.attrs)

    with open(os.path.join(xenium_bundle_uri, 'gene_panel.json')) as gene_panel_json:
        gene_panel = json.load(gene_panel_json)

    gene_info = [x for x in gene_panel['payload']['targets'] if x['type']['descriptor'] == 'gene']
    gene_ids = {i: x['type']['data']['id'] for i, x in enumerate(gene_info)}
    gene_names = {i: x['type']['data']['name'] for i, x in enumerate(gene_info)}
    gene_codewords = {i: x['codewords'] for i, x in enumerate(gene_info)}
    codeword_gene_index = {v: k for k, l in gene_codewords.items() for v in l}

    #TODO save this information in an array

    # select only codewords with gene targets

    permitted_codewords = sorted([x for v in gene_codewords.values() for x in v])
    gene_transcripts = transcripts[transcripts['codeword_index'].isin(permitted_codewords)]

    # add the gene index

    gene_transcripts_with_gene_index = gene_transcripts.assign(
        gene_index=gene_transcripts['codeword_index'].map(codeword_gene_index, meta=pd.Series(dtype="int32"))
    )

    # create a sparse array to store the data

    transcripts_uri=os.path.join(harmonised_dataset_uri, 'transcripts')

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name="z", domain=(-10000, 10000), tile=256, dtype=np.float32),
            tiledb.Dim(name="y", domain=y_domain, tile=256, dtype=np.float32),
            tiledb.Dim(name="x", domain=x_domain, tile=256, dtype=np.float32),
            tiledb.Dim(name="gene_index", domain=(0, len(gene_codewords)-1), tile=1, dtype=np.int32),
        ),
        attrs=[
            tiledb.Attr(name="transcript_id", dtype=np.uint64, filters=[tiledb.ZstdFilter()]),
            tiledb.Attr(name="qv", dtype=np.float32, filters=[tiledb.ZstdFilter()]),
        ],
        sparse=True,
        allows_duplicates=True
    )
    tiledb.Array.create(transcripts_uri, schema)

    # add data to the array

    def write_partition(partition, tiledb_uri, partition_info=None):
        if partition_info is None:
            return 1
        with tiledb.open(tiledb_uri, 'w') as tiledb_array:
            tiledb_array[
                partition['z_location'],
                partition['y_location'],
                partition['x_location'],
                partition['gene_index']
            ] = {
                'transcript_id': partition['transcript_id'],
                'qv': partition['qv']
            }
        return 1

    (gene_transcripts_with_gene_index
     .repartition(partition_size="1024MB")
     .map_partitions(write_partition, transcripts_uri).compute())

    # add the arrays to the group

    with tiledb.Group(harmonised_dataset_uri, "w") as group_root:
        group_root.add(immunofluorescence_uri, "immunofluorescence")
        group_root.add(immunofluorescence_focus_uri, "immunofluorescence_focus")
        group_root.add(cell_segmentation_mask_uri, "cell_segmentation_masks")
        group_root.add(cell_segmentation_uri, "cell_segmentation")
        group_root.add(nucleus_segmentation_mask_uri, "nucleus_segmentation_masks")
        group_root.add(nucleus_segmentation_uri, "nucleus_segmentation")
        group_root.add(cell_summary_uri, "cell_summary")
        group_root.add(cell_features_uri, "cell_features")
        group_root.add(he_image_uri, "hematoxylin_eosin_image")
        group_root.add(aligned_he_image_uri, "hematoxylin_eosin_image_aligned")
        group_root.add(transcripts_uri, "transcripts")


def write_record(section_label, uri):

    with sqlite3.connect(os.path.join(DATA_ROOT, 'metadata.db')) as conn:
        cursor = conn.cursor()

        cursor.execute('''
        INSERT INTO harmonised_datasets (section, uri) VALUES ((SELECT id FROM sections WHERE label=?), ?)
        ''', [section_label, uri])

        conn.commit()


def run_and_write_record(kwargs):
    run_kwargs = kwargs['run_kwargs']
    record_kwargs = kwargs['record_kwargs']

    print("starting run")
    pprint(run_kwargs)
    try:
        main(**run_kwargs)
    except Exception as e:
        raise RuntimeError(f'Error processing sample with configuration {run_kwargs}') from e

    print("recording run")
    pprint(record_kwargs)
    write_record(**record_kwargs)


# if __name__ == "__main__":

#     # configure tiledb
#     cfg = tiledb.Ctx().config()
#     cfg.update({'py.init_buffer_bytes': 1024 ** 2 * 50})
#     tiledb.default_ctx(cfg)

#     # connect to a dask cluster for distributed compute
#     client = Client(DASK_SCHEDULER_ADDRESS)

#     # query the available sections and datasets
#     with sqlite3.connect(os.path.join(DATA_ROOT, 'metadata.db')) as conn:
#         cursor = conn.cursor()
#         cursor.execute('''
#         SELECT sections.label, technologies.name, datasets.uri
#         FROM sections
#             JOIN datasets on sections.id = datasets.section
#             JOIN technologies ON datasets.technology = technologies.id
#             LEFT JOIN harmonised_datasets ON sections.id = harmonised_datasets.section
#         WHERE
#             harmonised_datasets.section IS NULL
#         ''')
#         response = cursor.fetchall()
#         datasets = pd.DataFrame(response, columns=['section_label', 'technology', 'uri']).pivot(index=['section_label'], columns='technology', values='uri')
#         datasets_dict = datasets.transpose().to_dict()
#         print("will harmonise datasets:")
#         pprint(datasets_dict)

#     # create output directory and error if the output already exists
#     os.makedirs(os.path.join(DATA_ROOT, 'harmonised'), exist_ok=True)

#     # generate output array groups
#     uuids = [str(uuid.uuid4()) for x in datasets_dict]
#     arguments = [{
#         "run_kwargs": {
#             "xenium_bundle_uri": v['Xenium In Situ'],
#             "hematoxylin_eosin_image_uri": v['NanoZoomer'],
#             "harmonised_dataset_uri": os.path.join(DATA_ROOT, 'harmonised', u),
#         },
#         "record_kwargs": {
#             "section_label": k,
#             "uri": os.path.join('harmonised', u),
#         }
#     } for (k, v), u in zip(datasets_dict.items(), uuids)]

#     print("will run with arguments")
#     pprint(arguments)

#     # start processing in parallel
#     print(f"starting execution pool with {CONCURRENCY} workers")

#     result = thread_map(
#         run_and_write_record,
#         arguments,
#         max_workers=CONCURRENCY
#     )

#     print("processing complete")

if __name__ == '__main__':
  fire.Fire()