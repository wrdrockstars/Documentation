from netCDF4 import Dataset

# Open the netCDF file
file_path = r"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\TEST\AWT_MODIS_FSC_20010101.NC"
var_name = 'fSCA'
with Dataset(file_path, 'r') as dataset:
    # List all variables
    for var_name in dataset.variables:
        var = dataset.variables[var_name]
        print(f"Variable: {var_name}")
        
        # Check for compression
        if hasattr(var, 'filters'):
            filters = var.filters()
            if filters['zlib']:
                print(f"  Compression: Enabled")
                print(f"  Compression Level: {filters['complevel']}")
            else:
                print(f"  Compression: Not Enabled")
        else:
            print(f"  Compression: Not Enabled")
