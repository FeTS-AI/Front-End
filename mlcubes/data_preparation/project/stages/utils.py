import os
from tqdm import tqdm


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


class MockTqdm(tqdm):
    def __getattr__(self, attr):
        return lambda *args, **kwargs: None
