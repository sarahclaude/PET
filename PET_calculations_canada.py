from xclim.core.units import convert_units_to
from xclim.indices.helpers import extraterrestrial_solar_radiation, _gather_lat
import xclim as xc
import xarray as xr

# test PET here and add them to xclim later
def doggers_allen(tasmax, tasmin, pr, tas=None, lat=None):
    if lat is None:
        lat = _gather_lat(tasmin if tas is None else tas)

    tasmin = convert_units_to(tasmin, "degC")
    tasmax = convert_units_to(tasmax, "degC")
    pr = convert_units_to(pr, "mm/day")
    if tas is None:
        tas = (tasmin + tasmax) / 2
    else:
        tas = convert_units_to(tas, "degC")

    ra = extraterrestrial_solar_radiation(tasmin.time, lat)
    ra = convert_units_to(ra, "MJ m-2 d-1")

    #dooger and allen formula
    out = 0.0013 * 0.408 * ra * (tas + 17)  * ((tasmax - tasmin) - (0.0123 * pr)) ** 0.76
    out = out.clip(0)
    return out

# calculate HG and DA PET and put them together
def mix_pet(ds):
    hg_pet = xc.indicators.atmos.potential_evapotranspiration(tasmin=ds['tasmin'], tasmax=ds['tasmax'], method='HG85')
    hg_pet.name = 'pet_hg'
    hg_pet = hg_pet * 86400
    hg_pet.attrs["units"] = "mm/day"
    da_pet = doggers_allen(tasmin=ds['tasmin'], tasmax=ds['tasmax'], pr=ds['pr'])
    da_pet.name = 'pet_da'
    da_pet.attrs = hg_pet.attrs

    for k, v in da_pet.attrs.items():
        if "hg85" in v:
            da_pet.attrs[k] = v.replace("hg85", "doggersallen")
    da_pet.attrs["cell_methods"] = "tasmin: time: minimum within days tasmax: time: maximum within days, pr: time: sum within days"
    pet = xr.merge([hg_pet, da_pet])
    pet.attrs = ds.attrs
    pet.attrs["cat:variable"] = ["pet_hg", "pet_da"]
    return pet


