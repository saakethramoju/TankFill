import math
import fullplot as fplt

from thermoprop import *
from fullflow import *


# ---- Units ---- #

PSIA_TO_PA = 6894.76
IN_TO_M = 0.0254
IN2_TO_M2 = 0.00064516
FT_TO_M = 0.3048
L_TO_M3 = 1/1000
MINUTE_TO_SECONDS = 60


# ---- Configuration ---- #

DEWAR_PRESSURE = 50.0 * PSIA_TO_PA
AMBIENT_PRESSURE = 14.67 * PSIA_TO_PA
AMBIENT_TEMPERATURE = 298.15

PASSTHROUGH_LENGTH = 1.0
PASSTHROUGH_DIAMETER = 0.0254

TANK_DIAMETER = 9.5 * IN_TO_M
TANK_VOLUME = 90 * L_TO_M3

FILL_SYSTEM_CDA = 0.00429 * IN2_TO_M2
VENT_VALVE_CDA = 3.3116e-5

GRAVITY = 9.80665


# ---- Solver ---- #

DT = 0.5 * MINUTE_TO_SECONDS
T_FINAL = 45 * MINUTE_TO_SECONDS

TARGET_TANK_FILL_PERCENT = 95.0

VERBOSE = True
STATISTICS = True

FILENAME = "lox_tank_fill"









# ---- Geometry ---- #

PASSTHROUGH_AREA = math.pi * PASSTHROUGH_DIAMETER**2 / 4
PASSTHROUGH_VOLUME = PASSTHROUGH_AREA * PASSTHROUGH_LENGTH

TANK_AREA = math.pi * TANK_DIAMETER**2 / 4
TANK_HEIGHT = TANK_VOLUME / TANK_AREA

TOTAL_VOLUME = PASSTHROUGH_VOLUME + TANK_VOLUME

TARGET_TANK_FILL_FRACTION = TARGET_TANK_FILL_PERCENT / 100.0

if not 0.0 < TARGET_TANK_FILL_PERCENT < 100.0:
    raise ValueError("TARGET_TANK_FILL_PERCENT must be between 0 and 100.")


# ---- Helper Functions ---- #

def tank_geometry(liquid_volume):
    if liquid_volume <= PASSTHROUGH_VOLUME:
        liquid_height = liquid_volume / PASSTHROUGH_AREA
        tank_liquid_volume = 0.0
    else:
        tank_liquid_volume = liquid_volume - PASSTHROUGH_VOLUME
        liquid_height = PASSTHROUGH_LENGTH + tank_liquid_volume / TANK_AREA

    tank_fill_fraction = tank_liquid_volume / TANK_VOLUME

    return {
        "liquid_height": liquid_height,
        "tank_liquid_volume": tank_liquid_volume,
        "tank_fill_fraction": tank_fill_fraction,
    }


def flash_split(raw_quality):
    quality = max(0.0, min(1.0, raw_quality))
    return {"quality": quality}


def positive_flow(mass_flow):
    return {"mass_flow": max(0.0, mass_flow)}


# ---- Custom Components ---- #

class LiquidInventory(Component):
    def __init__(
        self,
        name,
        network,
        mass,
        mass_flow_in,
    ):
        self.mass_dot = 0.0
        self.setup()

    def evaluate_states(self):
        self.mass_dot = self.mass_flow_in.value

    @property
    def dynamics(self):
        return [(self.mass, self.mass_dot)]


# ---- Network ---- #

LOXFill = Network("LOX Tank Fill")


# ---- States ---- #

liquid_mass = State(0.0, bounds=(0.0, None))

ullage_pressure = State(AMBIENT_PRESSURE, bounds=(1.0, None))
ullage_temperature = State(AMBIENT_TEMPERATURE)

ullage_air_mass_fraction = State(1.0, bounds=(0.0, 1.0))
ullage_gox_mass_fraction = State(0.0, bounds=(0.0, 1.0))

ullage_composition = {"air": ullage_air_mass_fraction, "O2": ullage_gox_mass_fraction}
gox_composition = {"air": 0.0, "O2": 1.0}


# ---- Dewar LOX ---- #

DewarLOX = Lookup(
    "Dewar LOX",
    LOXFill,
    Fluid,
    "LOX",
    pressure=DEWAR_PRESSURE,
    quality=0.0,
)


# ---- Tank Liquid Reference ---- #

LOXTankLiquid = Lookup(
    "LOX Tank Liquid",
    LOXFill,
    Fluid,
    "LOX",
    pressure=ullage_pressure,
    quality=0.0,
)

liquid_volume = liquid_mass / LOXTankLiquid.density
gas_volume = TOTAL_VOLUME - liquid_volume


# ---- Tank Geometry ---- #

TankGeometry = Lookup(
    "Tank Geometry",
    LOXFill,
    tank_geometry,
    liquid_volume,
    outputs=("liquid_height", "tank_liquid_volume", "tank_fill_fraction"),
)

liquid_height = TankGeometry.liquid_height
tank_liquid_volume = TankGeometry.tank_liquid_volume
tank_fill_fraction = TankGeometry.tank_fill_fraction

liquid_pressure = ullage_pressure + LOXTankLiquid.density * GRAVITY * liquid_height


# ---- Flash Properties ---- #

FlashLOX = Lookup(
    "Flash LOX",
    LOXFill,
    Fluid,
    "LOX",
    pressure=liquid_pressure,
    quality=0.0,
)

FlashGOX = Lookup(
    "Flash GOX",
    LOXFill,
    Fluid,
    "LOX",
    pressure=liquid_pressure,
    quality=1.0,
)

raw_flash_quality = (DewarLOX.enthalpy - FlashLOX.enthalpy) / (FlashGOX.enthalpy - FlashLOX.enthalpy)

FlashSplit = Lookup(
    "Flash Split",
    LOXFill,
    flash_split,
    raw_flash_quality,
    outputs=("quality",),
)

flash_quality = FlashSplit.quality


# ---- Combined Fill System ---- #

FillSystem = DischargeCoefficient(
    "Combined Fill System",
    LOXFill,
    upstream_pressure=DewarLOX.pressure,
    downstream_pressure=liquid_pressure,
    density=DewarLOX.density,
    discharge_coefficient=1.0,
    cross_sectional_area=FILL_SYSTEM_CDA,
)

gox_mass_flow_in = FillSystem.mass_flow * flash_quality
lox_mass_flow_in = FillSystem.mass_flow - gox_mass_flow_in


# ---- Liquid Inventory ---- #

LOXInventory = LiquidInventory(
    "LOX Inventory",
    LOXFill,
    mass=liquid_mass,
    mass_flow_in=lox_mass_flow_in,
)


# ---- Ullage Gas Properties ---- #

UllageAir = Lookup(
    "Ullage Air",
    LOXFill,
    IdealGas,
    "air",
    pressure=ullage_pressure,
    temperature=ullage_temperature,
)

UllageGOX = Lookup(
    "Ullage GOX",
    LOXFill,
    IdealGas,
    "O2",
    pressure=ullage_pressure,
    temperature=ullage_temperature,
)

ullage_gas_constant = ullage_air_mass_fraction * UllageAir.gas_constant + ullage_gox_mass_fraction * UllageGOX.gas_constant
ullage_cp = ullage_air_mass_fraction * UllageAir.specific_heat_cp + ullage_gox_mass_fraction * UllageGOX.specific_heat_cp
ullage_cv = ullage_air_mass_fraction * UllageAir.specific_heat_cv + ullage_gox_mass_fraction * UllageGOX.specific_heat_cv
ullage_gamma = ullage_cp / ullage_cv

ullage_density = ullage_pressure / (ullage_gas_constant * ullage_temperature)


# ---- Vent Valve ---- #

VentValve = CompressibleOrifice(
    "LOX Tank Vent Valve",
    LOXFill,
    upstream_total_pressure=ullage_pressure,
    upstream_total_temperature=ullage_temperature,
    downstream_pressure=AMBIENT_PRESSURE,
    discharge_coefficient=1.0,
    cross_sectional_area=VENT_VALVE_CDA,
    gas_constant=ullage_gas_constant,
    specific_heat_ratio=ullage_gamma,
)

VentOutflow = Lookup(
    "Vent Outflow",
    LOXFill,
    positive_flow,
    VentValve.mass_flow,
    outputs=("mass_flow",),
)

vent_mass_flow = VentOutflow.mass_flow


# ---- Ullage Volume ---- #

UllageVolume = Volume(
    "Ullage Volume",
    LOXFill,
    pressure=ullage_pressure,
    volume=gas_volume,
    density=ullage_density,
    mass_flow_in=gox_mass_flow_in,
    mass_flow_out=vent_mass_flow,
)


# ---- Ullage Composition ---- #

UllageComposition = Composition(
    "Ullage Composition",
    LOXFill,
    inlets=[(gox_mass_flow_in, gox_composition)],
    outlets=[(vent_mass_flow, None)],
    solve=ullage_composition,
    mass=UllageVolume.mass,
)


# ---- Fill Redline ---- #

TankFillLimit = fplt.Trace(
    x=[0.0, T_FINAL],
    y=[TARGET_TANK_FILL_FRACTION, TARGET_TANK_FILL_FRACTION],
    name="Tank Fill Limit",
    role="redline",
)

TankFillSensor = Sensor(
    "Tank Fill Sensor",
    LOXFill,
    reading=tank_fill_fraction,
    conditions=TankFillLimit,
)

FillSequence = Sequence("Fill Sequence", LOXFill)
FillSequence.abort(condition=(TankFillSensor, "Tank Fill Limit"), message="Target tank fill percentage reached.")


# ---- Tracking ---- #

LOXFill.track("Dewar Pressure [Pa]", DewarLOX.pressure)
LOXFill.track("Dewar Temperature [K]", DewarLOX.temperature)
LOXFill.track("Dewar Density [kg/m^3]", DewarLOX.density)

LOXFill.track("Fill System Mass Flow [kg/s]", FillSystem.mass_flow)
LOXFill.track("Incoming LOX Mass Flow [kg/s]", lox_mass_flow_in)
LOXFill.track("Incoming GOX Mass Flow [kg/s]", gox_mass_flow_in)

LOXFill.track("Flash Quality [-]", flash_quality)
LOXFill.track("Flash LOX Temperature [K]", FlashLOX.temperature)
LOXFill.track("Flash GOX Temperature [K]", FlashGOX.temperature)

LOXFill.track("LOX Mass [kg]", liquid_mass)
LOXFill.track("Total Liquid Volume [m^3]", liquid_volume)
LOXFill.track("Tank Liquid Volume [m^3]", tank_liquid_volume)
LOXFill.track("Tank Liquid Volume [L]", tank_liquid_volume / L_TO_M3)
LOXFill.track("Tank Fill Fraction [-]", tank_fill_fraction)
LOXFill.track("Tank Fill Percent [%]", tank_fill_fraction * 100.0)

LOXFill.track("Liquid Height [m]", liquid_height)
LOXFill.track("Ullage Pressure [Pa]", ullage_pressure)
LOXFill.track("Liquid Bottom Pressure [Pa]", liquid_pressure)

LOXFill.track("Tank LOX Temperature [K]", LOXTankLiquid.temperature)
LOXFill.track("Tank LOX Density [kg/m^3]", LOXTankLiquid.density)

LOXFill.track("Ullage Volume [m^3]", gas_volume)
LOXFill.track("Ullage Mass [kg]", UllageVolume.mass)
LOXFill.track("Ullage Temperature [K]", ullage_temperature)
LOXFill.track("Ullage Density [kg/m^3]", ullage_density)

LOXFill.track("Ullage Air Mass Fraction [-]", ullage_air_mass_fraction)
LOXFill.track("Ullage GOX Mass Fraction [-]", ullage_gox_mass_fraction)
LOXFill.track("Ullage Composition Sum [-]", ullage_air_mass_fraction + ullage_gox_mass_fraction)

LOXFill.track("Ullage Gas Constant [J/kg/K]", ullage_gas_constant)
LOXFill.track("Ullage Cp [J/kg/K]", ullage_cp)
LOXFill.track("Ullage Cv [J/kg/K]", ullage_cv)
LOXFill.track("Ullage Gamma [-]", ullage_gamma)

LOXFill.track("Vent Mass Flow [kg/s]", vent_mass_flow)
LOXFill.track("Vent Choked [-]", VentValve.choked)


# ---- Solve ---- #

Transient(LOXFill).solve(
    dt=DT,
    t_final=T_FINAL,
    filename=FILENAME,
    verbose=VERBOSE,
    statistics=STATISTICS,
)