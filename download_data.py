#!/usr/bin/env python
from pathlib import Path

import requests
import tqdm

def download_data(url):
    this_file_path = Path(__file__).resolve()
    data_dir = this_file_path.parent / "data"
    data_dir.mkdir(exist_ok=True)

    file_list = ["debug_graph.dot", "debug_graph.pt"]
    for file in file_list:
        response = requests.get(url + file, stream=True)
        total_size_in_bytes = int(response.headers.get('content-length', 0))
        block_size = 1024 * 1024  # 1 MB
        out_file = data_dir / file
        with out_file.open('wb') as f:
            for data in tqdm.tqdm(response.iter_content(chunk_size=block_size),
                                total=total_size_in_bytes//block_size, 
                                unit='MB', unit_scale=True,
                                desc=out_file.name):
                f.write(data)

    print("Downloaded data to", data_dir)

if __name__ == "__main__":
    url = "https://portal.nersc.gov/project/m3443/xju/walk_through_debug/"
    download_data(url)
