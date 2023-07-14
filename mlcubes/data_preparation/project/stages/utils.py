import os
from tqdm import tqdm


def update_row_with_dict(df, d, idx):
    for key in d.keys():
        df.loc[idx, key] = d.get(key)


def get_id_tp(index: str):
    return index.split("|")

class MockTqdm(tqdm):
    def __getattr__(self, attr):
        return lambda *args, **kwargs: None