import os
import shutil
from tqdm import tqdm


def normalize_path(path: str) -> str:
    """Remove mlcube-specific components from the given path

    Args:
        path (str): mlcube path

    Returns:
        str: normalized path
    """
    # for this specific problem, we know that all paths start with `/mlcube_io*`
    # and that this pattern won't change, shrink or grow. We can therefore write a
    # simple, specific solution
    if path.startswith("/mlcube_io"):
        return path[12:]

    # In case the path has already been normalized
    return path

def update_row_with_dict(df, d, idx):
    for key in d.keys():
        df.loc[idx, key] = d.get(key)


def get_id_tp(index: str):
    return index.split("|")


def set_files_read_only(path):
    for root, dirs, files in os.walk(path):
        for file_name in files:
            file_path = os.path.join(root, file_name)
            os.chmod(file_path, 0o444)  # Set read-only permission for files

        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            set_files_read_only(
                dir_path
            )  # Recursively call the function for subdirectories


def cleanup_storage(remove_folders):
    for folder in remove_folders:
        shutil.rmtree(folder, ignore_errors=True)


class MockTqdm(tqdm):
    def __getattr__(self, attr):
        return lambda *args, **kwargs: None
