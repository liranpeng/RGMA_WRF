import requests

token = "/global/homes/h/heroplr/.edl_token"

headers = {
    "Authorization": f"Bearer {token}"
}

with open("/pscratch/sd/h/heroplr/wrf_scratch/MERRA2_Data/missing_data.txt") as f:
    urls = f.read().splitlines()

for url in urls:
    filename = url.split("/")[-1]
    print(f"Downloading {filename}...")

    with requests.get(url, headers=headers, stream=True) as r:
        r.raise_for_status()
        with open(filename, "wb") as f_out:
            for chunk in r.iter_content(chunk_size=8192):
                f_out.write(chunk)
