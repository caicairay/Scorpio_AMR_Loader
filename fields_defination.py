import numpy as np
from yt.units import dimensions

def velocity_field(comp):
    def _velocity(field, data):
        return data["stream", "moment_%d" % comp]/data["stream","density"]
    return _velocity

def _setup_velocity_fields(ds):
    us = ds.unit_system
    for comp in range(1,4):
        ds.add_field(("gas","velocity_%s" % comp),
                     sampling_type="cell",
                     function=velocity_field(comp),
                     units = us["velocity"])

def setup_fluid_fields(ds):
    us = ds.unit_system
    _setup_velocity_fields(ds)
    def _kinetic_energy_density(field, data):
        vx=data["gas","velocity_1"]
        vy=data["gas","velocity_2"]
        v2 = vx**2+vy**2
        return 0.5*data["stream","density"]*v2
    ds.add_field(("gas", "kinetic_energy_density"), function=_kinetic_energy_density,
                 units=us["density"]*us["velocity"]**2,
                 dimensions=dimensions.density*dimensions.velocity**2,
                 sampling_type="cell")
    def _magnetic_energy_density(field, data):
        emag = 0.5 * data['stream', 'magnetic_1']**2
        emag += 0.5 * data['stream', 'magnetic_2']**2
#        emag /= 4 * np.pi
        return emag
    ds.add_field(('gas', 'magnetic_energy_density'), function=_magnetic_energy_density,
                 units=us["density"] * us["velocity"] ** 2,
                 dimensions=dimensions.density * dimensions.velocity ** 2,
                 sampling_type='cell')
    def _full_thermal_pressure_MHD(field, data):
        pthermal = (1.666667 - 1)*(data['stream','energy_density']-data['gas','kinetic_energy_density']
                    -data['gas','magnetic_energy_density'])
        return pthermal
    ds.add_field(('gas', 'thermal_pressure'), function=_full_thermal_pressure_MHD,
                 units=us['density']*us['velocity']**2,
                 dimensions=dimensions.density*dimensions.velocity**2,
                 sampling_type='cell')
