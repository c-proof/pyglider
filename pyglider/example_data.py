from zipfile import ZipFile
import pooch

def get_example_data(outdir='./'):
    """
    Get example data sets and configuration files

    Parameters
    ----------
    outdir : str or Path, optional
        Location to extract the example files into.  They will be put at
        ``outdir/pyglider-example-data/``.  Default is to unpack in the
        current directory.
    """
    zipfile = pooch.retrieve("https://github.com/c-proof/pyglider-example-data/archive/refs/heads/main.zip",
                             known_hash=None)

    with ZipFile(zipfile, 'r') as zipObj:
        # Extract all the contents of zip file in outdir
        zipObj.extractall(outdir)

__all__ = ['get_example_data']