import os
from pathlib import Path
if 'ESMFMKFILE' not in os.environ:
    os.environ['ESMFMKFILE'] = str(Path(os.__file__).parent.parent / 'esmf.mk')

import xscen as xs
from xscen import CONFIG
from dask import config as dskconf
from dask.distributed import Client
import clisops.core.subset as subset
import xarray as xr
from PET_calculations_canada import mix_pet
import glob
import logging
import warnings

# based on xscen workflow templates: https://xscen.readthedocs.io/en/latest/templates.html
xs.load_config(
    "paths.yml", "config.yml", verbose=(__name__ == "__main__"), reset=True
)

if "logging" in CONFIG:
    logger = logging.getLogger("pet_poly")

if __name__ == "__main__":

    warnings.filterwarnings("ignore")

    # set dask  configuration
    daskkws = CONFIG["dask"].get("client", {})
    dskconf.set(**{k: v for k, v in CONFIG["dask"].items() if k != "client"})

    # copy config to the top of the log file
    if "logging" in CONFIG and "file" in CONFIG["logging"]["handlers"]:
        f1 = open(CONFIG["logging"]["handlers"]["file"]["filename"], "a+")
        f2 = open("config.yml", 'r')
        f1.write(f2.read())
        f1.close()
        f2.close()

    # initialize Project Catalog (only do this once, if the file doesn't already exist)
    if not os.path.exists(CONFIG["pcat"]):
        pcat = xs.ProjectCatalog(
            CONFIG["pcat"], **CONFIG['create']
        )
    else:
        # load project catalog
        pcat = xs.ProjectCatalog(CONFIG["pcat"])

    # --- EXTRACT & REGRID ---
    if "extract" in CONFIG["tasks"]:

        # filter catalog for data that we want
        cat = xs.search_data_catalogs(**CONFIG["extract"]["simulation"]["search_data_catalogs"])

        cur = {"id": CONFIG['regrid']['target'],
                "xrfreq": "D",
                "processing_level": "extracted",
                "variable": "tasmax"
                }
        if not pcat.exists_in_cat(**cur):
            # load target grid
            with (
                Client(**CONFIG["extract"]["reconstruction"]["terraclimate"]["dask"], **daskkws),
                xs.measure_time(name="extracted ds_target", logger=logger)
            ):
                ds_target = xs.extract_dataset(catalog=cat[CONFIG['regrid']['target']],
                                               variables_and_freqs={'tasmax': "D"},
                                               **CONFIG['extract']["simulation"]["extract_dataset"])
                xs.save_and_update(ds_target['D'], pcat, **CONFIG['save'])
        else:
            ds_target = pcat.search(**cur).to_dask().sel(time=slice('1976', '1978')).load()

        # iterate on types to extract (reconstruction, simulation)
        for source_type, type_dict in CONFIG["extract"].items():
            if source_type == "reconstruction":
                if "reconstruction" in CONFIG["tasks"]:
                    with (xs.measure_time(name="reconstruction extraction / regridding", logger=logger)):
                        for rec_type, rec_dict in type_dict.items():
                            if rec_type=="cru_ts":
                                ds = xr.open_dataset(rec_dict["open_ds"])
                            else:
                                paths = glob.glob(rec_dict["open_ds"])
                                ds = xr.open_mfdataset(paths)
                            ds_can = subset.subset_shape(ds, **rec_dict["subset_shape"])
                            ds_can = ds_can.drop('crs')

                            ds_regrid = xs.regrid_dataset(
                                ds=ds_can,
                                ds_grid=ds_target,
                                **CONFIG['regrid']['regrid_dataset']
                            )
                            ds_regrid = xs.clean_up(ds_regrid, add_attrs=rec_dict["add_att"])
                            xs.save_and_update(ds_regrid, pcat, **CONFIG['save'])
            else:

                # iterate over ids from the search
                for ds_id, dc in cat.items():
                    # attrs of current iteration that are relevant now
                    cur = {
                        "id": ds_id,
                        "xrfreq": "D",
                        "processing_level": "extracted",
                        "variable": ["pet_hg", "pet_da"],
                    }
                    # check if steps was already done
                    if not pcat.exists_in_cat(**cur):
                        with (
                            Client(**type_dict["dask"], **daskkws),
                            xs.measure_time(name=f" {cur}", logger=logger)
                        ):
                            # create dataset from sub-catalog with right domain and periods
                            ds_dict = xs.extract_dataset(
                                catalog=dc,
                                **type_dict["extract_dataset"],
                            )

                            # iterate over the different datasets/frequencies
                            for key_freq, ds in ds_dict.items():
                                if dc.unique('bias_adjust_project')[0] != 'ESPO-G6-E5L':
                                    ds = xs.regrid_dataset(
                                        ds=ds,
                                        ds_grid=ds_target,
                                        **CONFIG['regrid']['regrid_dataset']
                                    )
                                ind = mix_pet(ds) # calculate PET
                                # save to zarr
                                xs.save_and_update(ind, pcat, **CONFIG['save'])

# --- Diagnostics ---
