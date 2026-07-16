# =============================================================================
# Download surface GRIB2 files per timestep from ECMWF open-data (IFS 0.25 deg)
# Native output frequency: 0-144 h every 3 h, 150-360 h every 6 h
# Full documentation: https://github.com/ecmwf/ecmwf-opendata
# Developed by: @ICPAC
# =============================================================================
import os
from time import sleep
from ecmwf.opendata import Client
from datetime import datetime


# forecast initialization
forecast_date = datetime.now().strftime("%Y%m%d")
time = "00"

output_dir = f"./ecmwf_forecasts_{forecast_date}/"
os.makedirs(output_dir, exist_ok=True)

# base variables available at every step
base_params = ["2t",         # 2 metre temperature (instantaneous)
               "2d",         # 2 metre dewpoint temperature
               "10u",        # 10 metre U wind component
               "10v",        # 10 metre V wind component
               "sp",         # Surface pressure
               "ssrd",       # Surface short-wave radiation downwards (accumulated)
               ]


def params_for_step(step):

    if step == 0:
        return base_params
    elif step <= 144:
        return base_params + ["mx2t3", "mn2t3"]
    else:
        return base_params + ["mx2t6", "mn2t6"]


# initialize ECMWF client
client = Client(
    source='ecmwf',  # open data source
    model="ifs",
    resol="0p25",
    preserve_request_order=False,
    infer_stream_keyword=True,
)

steps_3h = list(range(0, 145, 3))    # 0, 3, 6, ..., 144
steps_6h = list(range(150, 361, 6))  # 150, 156, ..., 360
forecast_steps = steps_3h + steps_6h

max_retries = 3

for step in forecast_steps:
    step_str = f"{step:03d}"

    surf_file = os.path.join(
        output_dir, f"ECMWF_sfc_{forecast_date}{time}_{step_str}.grib2")

    # skip if file already exists
    if os.path.exists(surf_file):
        print(f"File already exists: {step_str}h, skipping to the next download...")
        continue

    print(f"Downloading forecast data at {step_str}h...")

    # download surface data (with retries for transient network failures)
    for attempt in range(1, max_retries + 1):
        try:
            client.retrieve(
                type="fc",
                stream="oper",
                levtype="sfc",
                date=forecast_date,
                time=time,
                step=step,
                param=params_for_step(step),
                target=surf_file
            )
            break
        except Exception as e:
            print(f"  Attempt {attempt}/{max_retries} failed for step {step_str}h: {e}")
            # remove partial file so the next attempt / rerun starts clean
            if os.path.exists(surf_file):
                os.remove(surf_file)
            if attempt < max_retries:
                sleep(30)
            else:
                print(f"  Giving up on step {step_str}h - rerun the script to retry.")

print(f"\n All GRIB2 files are ready in {output_dir} for preprocessing.")
