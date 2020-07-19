import numpy as np
from yt.units import dimensions
"""
TODO:
    - missing dimensions
    - alias
    - parameters, like gamma etc
    - thermal pressure
"""
code_density = "code_mass / code_length**3"
code_moment = "code_mass / code_length**2 / code_time"
code_pressure = "code_mass / code_length / code_time**2"

direction_alias= {
    "cartesian": ("x","y","z")
    }
def velocity_field(comp):
    def _velocity(field, data):
        return data["stream", "mom%d" % comp]/data["stream","density"]
    return _velocity
def _setup_velocity_fields(ds):
    us = ds.unit_system
    for comp in range(3):
        ds.add_field(("gas","velocity_%s" % direction_alias["cartesian"][comp]),
                     sampling_type="cell",
                     function=velocity_field(comp+1),
                     units = us["velocity"])

def magnetic_field(comp):
    def _magnetic(field, data):
        return data["stream", "bl%d" % comp]+data["stream","br%d" % comp]
    return _magnetic
def _setup_magnetic_fields(ds):
    for comp in range(3):
        ds.add_field(("gas","magnetic_field_%s" % direction_alias["cartesian"][comp]),
                     sampling_type="cell",
                     function=magnetic_field(comp+1),
                     units = "code_magnetic")

def energy_density(field,data):
    return data["stream","ene"]
def _setup_energy_density(ds):
    ds.add_field(
        ("gas","energy_density"),
        sampling_type="cell",
        function=energy_density,
        units=code_pressure
        )

def setup_fluid_fields(ds):
    us = ds.unit_system
    _setup_velocity_fields(ds)
    _setup_magnetic_fields(ds)
    _setup_energy_density(ds)
    def _kinetic_energy_density(field, data):
        v2 = 0
        for comp in range(3):
            v2 += data["gas","velocity_%s" % direction_alias["cartesian"][comp]]**2
        return 0.5*data["gas","density"]*v2
    ds.add_field(("gas", "kinetic_energy_density"), function=_kinetic_energy_density,
                 units=us["density"]*us["velocity"]**2,
                 dimensions=dimensions.density*dimensions.velocity**2,
                 sampling_type="cell")
    def _magnetic_energy_density(field, data):
        emag =0
        for comp in range(3):
            emag += 0.5 * data['gas', 'magnetic_field_%s'
                    % direction_alias["cartesian"][comp]]**2
#        emag /= 4 * np.pi
        return emag
    ds.add_field(('gas', 'magnetic_energy_density'), function=_magnetic_energy_density,
                 units=us["density"] * us["velocity"] ** 2,
                 dimensions=dimensions.density * dimensions.velocity ** 2,
                 sampling_type='cell')
    # def _full_thermal_pressure_MHD(field, data):
    #     pthermal = (1.666667 - 1)*(data['stream','energy_density']-data['gas','kinetic_energy_density']
    #                 -data['gas','magnetic_energy_density'])
    #     return pthermal
    # ds.add_field(('gas', 'thermal_pressure'), function=_full_thermal_pressure_MHD,
    #              units=us['density']*us['velocity']**2,
    #              dimensions=dimensions.density*dimensions.velocity**2,
    #              sampling_type='cell')
