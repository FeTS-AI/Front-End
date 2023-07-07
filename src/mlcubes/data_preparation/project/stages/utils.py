import os


def update_row_with_dict(df, d, idx):
    for key in d.keys():
        df.loc[idx, key] = d.get(key)


def get_id_tp(index: str):
    return index.split("|")
