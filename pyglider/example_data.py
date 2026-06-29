import tarfile

import pooch


def get_example_data(outdir='./'):
    """
    Get example data sets and configuration files.

    Downloads ``example-data.tar.gz`` from the pyglider documentation site and
    extracts it into *outdir*, creating ``outdir/example-data/``.

    Parameters
    ----------
    outdir : str or Path, optional
        Location to extract the example files into.  They will be put at
        ``outdir/example-data/``.  Default is to unpack in the
        current directory.
    """
    tarball = pooch.retrieve(
        'https://pyglider.readthedocs.io/en/stable/example-data.tar.gz',
        known_hash=None,
    )

    with tarfile.open(tarball, 'r:gz') as tf:
        tf.extractall(outdir)


__all__ = ['get_example_data']
