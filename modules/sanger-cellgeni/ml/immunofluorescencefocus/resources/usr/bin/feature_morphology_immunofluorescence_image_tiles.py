#!/usr/bin/env python3

import dask.array as da
import numpy as np
from dask_image.ndinterp import affine_transform
import tiledb
# from dask.distributed import Client
import datetime
import os
from itertools import repeat
# from tqdm.contrib.concurrent import thread_map
import dask
import json
# import sqlite3
from pprint import pprint
import uuid
import zarr
import fire

# default run configuration

DASK_SCHEDULER_ADDRESS='tcp://farm22-head2:40883'
CONCURRENCY=4
DATA_ROOT='/lustre/scratch126/cellgen/team361/projects/histology_to_gene_expression/workspace/data'

dask.config.set(scheduler = 'single-threaded')

# feature definition

FEATURE_TYPE_NAME='Morphology immunofluorescence image tiles'
VERSION = '1.0.1'



def generate_data_array(harmonised_dataset_absolute_path, data_array_absolute_path, tile_size_px=256, resolution_um_per_px=0.5, clipping_percentile=98):
    """Generate normalized intensity values for an entire immunofluorescence image

    Generate a large array of shape (c, h, w) with normalized and clipped intensity values, representing the entire
    slide. This feature generation includes image rescaling and re-chunking. We assume the immunofluorescence image is of
    resolution 0.2125 um / px, as described in the 10x documentation.
    """


    immunofluorescence_focus_uri = os.path.join(harmonised_dataset_absolute_path, 'immunofluorescence_focus')
    image_array = da.from_tiledb(immunofluorescence_focus_uri, attribute='px_value', chunks=(1, 8192, 8192))
   


    # rescale the image

    scale_factor = 0.2125 / resolution_um_per_px
    

    transform_matrix = np.array([
        [1, 0, 0, 0],
        [0, 1 / scale_factor, 0, 0],
        [0, 0, 1 / scale_factor, 0],
        [0, 0, 0, 1]
    ])

    rescaled_image = affine_transform(
        image_array,
        transform_matrix,
        output_shape=(4, int(image_array.shape[1] * scale_factor), int(image_array.shape[2] * scale_factor)),
        output_chunks=(1, 8192, 8192)
    )

    # normalize and clip the intensity values

    percentiles = da.percentile(rescaled_image, clipping_percentile, axis=(1, 2))
    

    clipped = da.stack([
        da.clip(rescaled_image[i], 0, percentiles[i])
        for i in range(len(rescaled_image))
    ])

    max_values = clipped.max(axis=(1, 2))
    

    normalized = da.stack([
        clipped[i] / max_values[i]
        for i in range(len(clipped))
    ])

    normalized_uint8 = (normalized * 255).astype(np.uint8)


    # chunks = (1, tile_size_px, tile_size_px)
    # shards = (1, tile_size_px * 4, tile_size_px * 4)
    chunks = (4, tile_size_px, tile_size_px)
    shards = (4, tile_size_px * 4, tile_size_px * 4)
    shape = normalized_uint8.shape
    dtype = normalized_uint8.dtype

    zarr_array = zarr.create_array(
        store=data_array_absolute_path,
        # mode="w",
        shape=shape,
        chunks=chunks,
        shards=shards,
        dtype=dtype
    )

    processing_chunks = tuple([max(512, chunk) for chunk in chunks])

    def write_block(block, zarr_array, block_info=None):
        if block_info is None:
            return da.array([[[0]]])
        slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
        zarr_array[slices] = block
        return da.array([[[0]]])

    result = normalized_uint8.map_blocks(write_block, zarr_array)
    result.compute()
    meta = {
        'name': 'immunofluorescence',
        'display_name': FEATURE_TYPE_NAME,
        'version': VERSION,
        'datas_source': harmonised_dataset_absolute_path,
        'parameters': {
            'resolution_um_per_px': resolution_um_per_px,
            "clipping_percentile": clipping_percentile,
            'tile_size_px': tile_size_px,
            'shape': shape,
            'dtype': str(dtype),
            'chunks': chunks,
            'shards': shards,
            'pixel_size_um': 0.2125
        }
    }

    with open(f'{data_array_absolute_path}.json', 'w') as file:
        json.dump(meta, file)
    # # create the output array

    # schema = tiledb.ArraySchema(
    #    domain=tiledb.Domain(
    #        tiledb.Dim(name="stain", domain=(0, normalized_uint8.shape[0] - 1), tile=1, dtype=np.int32),
    #        tiledb.Dim(name="px_rows", domain=(0, normalized_uint8.shape[1] - 1), tile=tile_size_px, dtype=np.int32),
    #        tiledb.Dim(name="px_cols", domain=(0, normalized_uint8.shape[2] - 1), tile=tile_size_px, dtype=np.int32)
    #    ),
    #    attrs=[
    #        tiledb.Attr(name="px_value", dtype=normalized_uint8.dtype, filters=[tiledb.ZstdFilter()])
    #    ]
    # )
    # tiledb.Array.create(data_array_absolute_path, schema)

    # # add data to the array

    # def write_block(block, tiledb_uri, block_info=None):
    #     if block_info is None:
    #         return da.array([[[1]]]) # map_blocks expects to get an array back
    #     slices = tuple(slice(start, stop) for start, stop in block_info[None]['array-location'])
    #     with tiledb.DenseArray(tiledb_uri, mode='w') as tiledb_array:
    #         tiledb_array[slices] = {"px_value": block}
    #     return da.array([[[1]]])

    # normalized_uint8.map_blocks(write_block, data_array_absolute_path).compute()



# def write_record(feature_flavor_id, region_set_id, data_array_uri, index_array_uri):

#     with sqlite3.connect(os.path.join(DATA_ROOT, 'metadata.db')) as conn:
#         cursor = conn.cursor()

#         cursor.execute('''
#         INSERT INTO feature_sets (flavor, generated_for, data_array_uri, idx_array_uri) VALUES (?, ?, ?, ?)
#         ''', [feature_flavor_id, region_set_id, data_array_uri, index_array_uri])

#         conn.commit()


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
#             "tile_size_px": r['parameters']['tile_size_px'],
#             "resolution_um_per_px": r['parameters']['resolution_um_per_px'],
#             "clipping_percentile": r['parameters']['clipping_percentile']
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