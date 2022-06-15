# from pathlib import Path
#
# library_dir = Path(__file__).parent.parent.absolute()

def pytest_sessionstart(session):
    pass
    # os.mkdir(library_dir / '_example_data')
    # pyexamp.get_example_data(library_dir / '_example_data')


def pytest_sessionfinish(session, exitstatus):
    pass