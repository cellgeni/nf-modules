#!/usr/bin/env python3

import dask.array as da
import tiledb
import numpy as np
from numpy import dtype
import dask
import itertools
# from dask.distributed import Client
import datetime
import os
from itertools import repeat
# from tqdm.contrib.concurrent import thread_map
# import sqlite3
from pprint import pprint
import json
import uuid
import zarr
import fire

# default run configuration

DASK_SCHEDULER_ADDRESS='tcp://farm22-head2:40883'
CONCURRENCY=4
DATA_ROOT='/lustre/scratch126/cellgen/team361/projects/histology_to_gene_expression/workspace/data'

dask.config.set(scheduler = 'single-threaded')

# feature definition

FEATURE_TYPE_NAME='Spatially binned transcript count tiles'
VERSION = '1.1.0'


def generate_data_array(harmonised_dataset_absolute_path, data_array_absolute_path, bin_size_um_per_bin=0.5):
    """Generate spatially binned transcript counts across the entire capture area

    Generate a sparse array of shape (gene_index, h, w) binned counts of transcripts, at the specified bin size. When
    the bin size is equal to the resolution of the rescaled morphology image, positions in these arrays should
    correspond with one another.
        """

    with tiledb.Group(harmonised_dataset_absolute_path, 'r') as group_root:
        transcripts_uri = group_root['transcripts'].uri

    # compute the scale factor and new domain

    scale_factor = 1 / bin_size_um_per_bin

    with tiledb.open(transcripts_uri, 'r') as transcripts_array:
        domains = {dim.name: dim.domain for dim in transcripts_array.schema.domain}


    y_domain_binned = (np.floor(domains['y'][0] * scale_factor).astype(np.uint32), np.floor(domains['y'][1] * scale_factor).astype(np.uint32) - 1)
    x_domain_binned = (np.floor(domains['x'][0] * scale_factor).astype(np.uint32), np.floor(domains['x'][1] * scale_factor).astype(np.uint32) - 1)
    output_array_shape = (
        domains['gene_index'][1] + 1,
        np.floor(domains['y'][1] * scale_factor).astype(np.uint32),
        np.floor(domains['x'][1] * scale_factor).astype(np.uint32)
    )


    # chunks = (256, 4096, 4096)
    # shards = (256, 4096 * 4, 4096 * 4)
    chunks = (64, 512, 512)
    shards = (2048, 2048, 2048)
    shape = tuple([int(i)for i in output_array_shape])
    dtype = np.uint8

    zarr_array = zarr.create_array(
        store=data_array_absolute_path,
        # mode="w",
        shape=shape,
        chunks=chunks,
        shards=shards,
        dtype=dtype
    )

    processing_chunks = tuple([max(512, chunk) for chunk in chunks])



    # # create an array for the binned transcripts

    # schema = tiledb.ArraySchema(
    #     domain=tiledb.Domain(
    #         tiledb.Dim(name="gene_index", domain=(0, 5000), tile=64, dtype=np.int32),
    #         tiledb.Dim(name="y", domain=y_domain_binned, tile=256, dtype=np.int32),
    #         tiledb.Dim(name="x", domain=x_domain_binned, tile=256, dtype=np.int32),
    #     ),
    #     attrs=[
    #         tiledb.Attr(name="count", dtype=np.uint8, filters=[tiledb.ZstdFilter()]),
    #     ],
    # )
    # tiledb.Array.create(data_array_absolute_path, schema)

    # add data to the array

    def load_tiledb_as_dense_array(block, scale_factor, array_uri, dtype=np.uint8, block_info=None):
        if block_info is None:
            return block  # map_blocks expects to get an array back
        dense_array_slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
        block_shape = tuple(stop - start for start, stop in block_info[None]['array-location'])
        block_origin = tuple(start for start, stop in block_info[None]['array-location'])

        # convert coordinates to the transcripts sparse array
        sparse_array_slices = tuple([
            dense_array_slices[0],
            slice(dense_array_slices[1].start / scale_factor, dense_array_slices[1].stop / scale_factor),
            slice(dense_array_slices[2].start / scale_factor, dense_array_slices[2].stop / scale_factor)
        ])

        # tiledb query on the transcripts sparse array
        with tiledb.open(array_uri, 'r') as transcripts_array:
            result = transcripts_array[:, sparse_array_slices[1], sparse_array_slices[2], sparse_array_slices[0]] # be careful here, the coordinates are switched and the z coordinate is aggregated

        # convert query results to rescaled coordinates and reposition relative to the block origin
        gene_indexes = result['gene_index'] - block_origin[0]
        ys = (result['y'] * scale_factor) - block_origin[1]
        xs = (result['x'] * scale_factor) - block_origin[2]

        # for tiledb queries "with a floating-point array domain, index bounds are inclusive"
        # we therefore need to adjust the result to be non-inclusive on the upper bound
        # we do this on the query results rather than the transformed results, as this is less
        # likely to incorrectly exclude transcripts outside of the desired bounds

        within_bound_idx = [(y < sparse_array_slices[1].stop) & (x < sparse_array_slices[2].stop) for y, x in zip(result['y'], result['x'])]

        gene_indexes = gene_indexes[within_bound_idx]
        ys = ys[within_bound_idx]
        xs = xs[within_bound_idx]

        # bin the transcripts
        y_bins = np.floor(ys).astype(np.int32)
        x_bins = np.floor(xs).astype(np.int32)
        bins = np.array([gene_indexes, y_bins, x_bins])
        unique_bins, bin_counts = np.unique(bins, axis=1, return_counts=True)

        # convert to a dense array
        dense_array = np.zeros(shape=block_shape, dtype=dtype)
        dense_array[
            unique_bins[0],
            unique_bins[1],
            unique_bins[2]
        ] = bin_counts

        # return the dense array
        return dense_array



    def write_block(block, zarr_array, block_info=None):
        if block_info is None:
            return da.array([[[0]]])
        slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
        zarr_array[slices] = block
        return da.array([[[0]]])
    # def write_block(block, tiledb_uri, block_info=None):
    #     if block_info is None:
    #         return da.array([[[0]]]) # map_blocks expects to get an array back
    #     slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
    #     with tiledb.DenseArray(tiledb_uri, mode='w') as tiledb_array:
    #         tiledb_array[slices] = {"count": block}
    #     return da.array([[[0]]])

    # result = rescaled_image.map_blocks(write_block, zarr_array)
    # result.compute()

    result = (da.zeros(shape=output_array_shape, chunks=(256, 4096, 4096), dtype=np.uint8) # chunks used for processing
     .map_blocks(load_tiledb_as_dense_array, scale_factor, transcripts_uri)
     .map_blocks(write_block, zarr_array))

    result.compute()

    meta = {
        'name': 'spatially_binned_expression',
        'display_name': FEATURE_TYPE_NAME,
        'version': VERSION,
        'datas_source': harmonised_dataset_absolute_path,
        'parameters': {
            'bin_size_um_per_bin': bin_size_um_per_bin,
            'shape': shape,
            'dtype': str(dtype),
            'chunks': chunks,
            'shards': shards
        }
    }

    with open(f'{data_array_absolute_path}.json', 'w') as file:
        json.dump(meta, file)



# def run_and_write_record(kwargs):
#     generate_data_array_kwargs = kwargs['generate_data_array_kwargs']
#     generate_index_array_kwargs = kwargs['generate_index_array_kwargs']
#     record_kwargs = kwargs['record_kwargs']

#     print("starting data array generation")
#     pprint(generate_data_array_kwargs)
#     try:
#         generate_data_array(**generate_data_array_kwargs)
#     except Exception as e:
#         raise RuntimeError(f'Error generating data array for sample with configuration {generate_data_array_kwargs}') from e

#     print("starting index array generation")
#     pprint(generate_index_array_kwargs)
#     try:
#         generate_idx_array(**generate_index_array_kwargs)
#     except Exception as e:
#         raise RuntimeError(f'Error generating index array for sample with configuration {generate_index_array_kwargs}') from e

#     print("recording run")
#     pprint(record_kwargs)
#     write_record(**record_kwargs)


# if __name__ == "__main__":

#     # configure tiledb
#     cfg = tiledb.Ctx().config()
#     cfg.update({'py.init_buffer_bytes': 1024 ** 2 * 50})
#     tiledb.default_ctx(cfg)

#     # connect to a dask cluster for distributed compute
#     client = Client(DASK_SCHEDULER_ADDRESS)

#     # query the available run configurations
#     with sqlite3.connect(os.path.join(DATA_ROOT, 'metadata.db')) as conn:
#         cursor = conn.cursor()
#         cursor.execute('''
#             SELECT region_sets.id, feature_flavors.id, feature_flavors.parameters, harmonised_datasets.uri 
#             FROM feature_types 
#                 JOIN feature_flavors ON feature_types.id = feature_flavors.type
#                 JOIN region_flavors ON feature_flavors.generated_for = region_flavors.id
#                 JOIN region_sets ON region_flavors.id = region_sets.flavor
#                 JOIN sections ON region_sets.section = sections.id
#                 JOIN harmonised_datasets ON sections.id = harmonised_datasets.section
#                 LEFT JOIN feature_sets 
#                     ON feature_flavors.id = feature_sets.flavor 
#                     AND region_sets.id = feature_sets.generated_for
#             WHERE feature_types.name=?
#                 AND feature_sets.id IS NULL
#             ''', [FEATURE_TYPE_NAME])
#         response = cursor.fetchall()
#         run_configurations = [
#             {"region_set_id": a, "feature_flavor_id": b, "parameters": json.loads(c), "harmonised_dataset_uri": d} for
#             a, b, c, d in response]
#         print("will run with configuration:")
#         pprint(run_configurations)

#     # create output directory
#     os.makedirs(os.path.join(DATA_ROOT, 'data_arrays'), exist_ok=True)
#     os.makedirs(os.path.join(DATA_ROOT, 'idx_arrays'), exist_ok=True)

#     # generate output arrays
#     data_array_uuids = [str(uuid.uuid4()) for x in run_configurations]
#     index_array_uuids = [str(uuid.uuid4()) for x in run_configurations]
#     arguments = [{
#         "generate_data_array_kwargs": {
#             "harmonised_dataset_absolute_path": os.path.join(DATA_ROOT, r['harmonised_dataset_uri']),
#             "data_array_absolute_path": os.path.join(DATA_ROOT, 'data_arrays', data_array_uuid),
#             "bin_size_um_per_bin": r['parameters']['bin_size_um_per_bin'],
#         },
#         "generate_index_array_kwargs": {
#             "data_array_absolute_path": os.path.join(DATA_ROOT, 'data_arrays', data_array_uuid),
#             "idx_array_absolute_path": os.path.join(DATA_ROOT, 'idx_arrays', index_array_uuid),
#             "tile_size_px": r['parameters']['tile_size_px'],
#             "stride_size_px": r['parameters']['stride_size_px'],
#             "exclude_out_of_tissue": r['parameters']['exclude_out_of_tissue']
#         },
#         "record_kwargs": {
#             "feature_flavor_id": r['feature_flavor_id'],
#             "region_set_id": r['region_set_id'],
#             "data_array_uri": os.path.join('data_arrays', data_array_uuid),
#             "index_array_uri": os.path.join('idx_arrays', index_array_uuid)
#         }
#     } for r, data_array_uuid, index_array_uuid in zip(run_configurations, data_array_uuids, index_array_uuids)]

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