from objects.compartment import *

# Subclasses (inheritances) of the class compartment add extra attributes to the compatment that define the type of compartment

UTOPIA_surfaceSea_water_compartments = ["Ocean_Surface_Water", "Coast_Surface_Water"]

UTOPIA_water_compartments = [
    "Ocean_Mixed_Water",
    "Ocean_Column_Water",
    "Coast_Column_Water",
    "Surface_Freshwater",
    "Bulk_Freshwater",
]

UTOPIA_deep_soil_compartments = [
    "Beaches_Deep_Soil",
    "Background_Soil",
    "Impacted_Soil",
]

UTOPIA_soil_surface_compartments = [
    "Beaches_Soil_Surface",
    "Background_Soil_Surface",
    "Impacted_Soil_Surface",
]

UTOPIA_sediment_compartment = [
    "Sediment_Freshwater",
    "Sediment_Ocean",
    "Sediment_Coast",
]

UTOPIA_air_compartments = ["Air"]


# listOfProcessess=['k_discorporation','k_fragmentation','k_heteroaggregation',
# 'k_heteroaggregate_breackup','k_biofouling', 'k_defouling', 'k_advective_transport',
#  'k_settling', 'k_rising', "k_sea_spray_aerosol", "k_wind_trasport",
#  "k_wet_depossition", "k_dry_deposition", "k_sediment_resuspension",
#  "k_burial", "k_percolation", "k_runoff_transport", "k_tillage"]


class compartment_water(Compartment):
    # added new attributes relative to the water compartments such as SPM concentration, flow velocity etc.
    def __init__(
        self,
        Cname,
        SPM_mgL,
        waterFlow_m3_s,
        T_K,
        G,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
        flowVelocity_m_s=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.SPM_mgL = SPM_mgL
        self.flowVelocity_m_s = flowVelocity_m_s
        self.waterFlow_m3_s = waterFlow_m3_s
        self.T_K = T_K
        self.G = G  # Shear rate (G, in s−1)
        self.processess = [
            "discorporation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
            "mixing",
        ]
        # if waterFlow_m3_s == "nan":
        #     waterFlow_m3_s = self.flowVelocity_m * self.Cdepth_m * self.Cwidth_m
        # else:
        #     pass


class compartment_surfaceSea_water(Compartment):
    # added new processess to the list of processess. new attributes that migth be needed to this processess should be added here
    def __init__(
        self,
        Cname,
        SPM_mgL,
        waterFlow_m3_s,
        T_K,
        G,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
        flowVelocity_m_s=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.SPM_mgL = SPM_mgL
        self.flowVelocity_m_s = flowVelocity_m_s
        self.waterFlow_m3_s = waterFlow_m3_s
        self.T_K = T_K
        self.G = G  # Shear rate (G, in s−1)
        self.processess = [
            "discorporation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
            "mixing",
            "sea_spray_aerosol",
            "beaching",
        ]


class compartment_sediment(Compartment):
    def __init__(
        self,
        Cname,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.processess = [
            "discorporation",
            "fragmentation",
            "sediment_resuspension",
            "burial",
        ]


class compartment_soil_surface(Compartment):
    def __init__(
        self,
        Cname,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )

        self.processess = [
            "discorporation",
            "fragmentation",
            "runoff_transport",
            "percolation",
            "soil_air_resuspension",
            "soil_convection",
        ]
        # self.earthworm_density_in_m3 = earthworm_density_in_m3
        # self.Qrunoff_m3 = Qrunoff_m3


class compartment_deep_soil(Compartment):
    def __init__(
        self,
        Cname,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.processess = [
            "discorporation",
            "fragmentation",
            "sequestration_deep_soils",
            "soil_convection",
        ]


# retention_in_soil (straining?) of the particles in soil following heteroaggregation with geocolloids?
# shall we also include heteroaggregation/heteroaggegrate break-up processess in the soil compartment? In SimpleBox for Nano they do account for aggregation and attachment

# Difference between retention in soil and sequestration deep soil: sequestrations deep soil is like burial in deep sediments (elemination process-->out of the system)


class compartment_air(Compartment):
    def __init__(
        self,
        Cname,
        T_K=None,
        wind_speed_m_s=None,
        I_rainfall_mm=None,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
        flowVelocity_m_s=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.T_K = T_K
        self.wind_speed_m_s = wind_speed_m_s
        self.I_rainfall_mm = I_rainfall_mm
        self.flowVelocity_m_s = flowVelocity_m_s
        self.processess = [
            "discorporation",
            "fragmentation",
            "wind_trasport",
            "dry_deposition",
            "wet_deposition",
        ]
        # shall we also include heteroaggregation/heteroaggegrate break-up processess in the air compartment?
