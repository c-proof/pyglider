from zipfile import ZipFile
import pooch

def get_example_data(outdir='./'):
    """
    Get example data sets and configuration files

    Parameters
    ----------
    outdir : str or Path, optional
        Location to extract the example files into.  They will be put at
        ``outdir/example-data/``.  Default is to unpack in the
        current directory.
    """
    zipfile = pooch.retrieve("https://cproof.uvic.ca/pyglider-example-data/pyglider-example-data.zip",
                             known_hash='5643a5301530e8dd60060a357cd9ed88eb1e84d761710c2a4013bc3c1817a859')

    with ZipFile(zipfile, 'r') as zipObj:
        # Extract all the contents of zip file in outdir
        zipObj.extractall(outdir)

__all__ = ['get_example_data']