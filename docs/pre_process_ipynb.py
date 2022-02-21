
import sys
import glob
import json
import pathlib
import shutil


def hide_cells(notebook, tags=None):
    """
    Finds the tags in each cell and removes it.
    Returns dict without 'hide' tagged cells.
    """
    if tags is None:
        tags = ['hide']

    clean = []
    for cell in notebook['cells']:
        if not set(tags).intersection(cell['metadata'].get('tags', list())):
            clean.append(cell)

    notebook['cells'] = clean

    return notebook


def change_kernel(notebook):
    """
    Vanillafy the kernelspec.
    """
    new_kernelspec = {
        "display_name": "Python 3 (ipykernel)",
        "language": "python",
        "name": "python3",
    }
    notebook['metadata']['kernelspec'].update(new_kernelspec)
    return notebook


def main(path):
    """
    Process the IPYNB files in path, save in place (side-effect).
    """
    fnames = glob.glob(path.strip('/') + '/[!_]*.ipynb')  # Not files with underscore.
    outpath = pathlib.Path('userguide')
    if outpath.exists():
        shutil.rmtree(outpath)
    outpath.mkdir(exist_ok=True)

    for fname in fnames:
        with open(fname, encoding='utf-8') as f:
            notebook = json.loads(f.read())

        notebook = change_kernel(notebook)
        notebook = hide_cells(notebook)
        filepart = pathlib.Path(fname).name

        with open(outpath / filepart, 'w') as f:
            _ = f.write(json.dumps(notebook))

    return


if __name__ == '__main__':
    print(sys.argv[1])
    _ = main(sys.argv[1])