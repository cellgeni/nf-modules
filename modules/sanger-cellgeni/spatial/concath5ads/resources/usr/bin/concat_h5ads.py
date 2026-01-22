#!/usr/bin/env python3
import click
import scanpy as sc

@click.command()
@click.version_option(version='0.1.0')  # Add version option here
@click.argument('input_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('-o', '--output', required=True, type=click.Path(), help='Output h5ad file path')
@click.option('--axis', type=int, default=0, help='Axis to concatenate (0=observations, 1=variables)')
def concat_h5ads(input_files, output, axis):
    """Concatenate multiple h5ad files into a single h5ad file."""
    
    if not input_files:
        raise click.UsageError("At least one input file is required")
    
    click.echo(f"Loading {len(input_files)} file(s)...")
    adata_list = [sc.read_h5ad(f) for f in input_files]
    
    click.echo("Concatenating...")
    adata_concat = sc.concat(adata_list, axis=axis, join='outer')
    
    click.echo(f"Writing output to {output}")
    adata_concat.write_h5ad(output)
    
    click.echo("Done!")


if __name__ == '__main__':
    concat_h5ads()