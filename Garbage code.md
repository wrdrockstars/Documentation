The following code, which was popularly known as Udta_teer_3, was modified and added a few lines. These lines became imperative because while running the whole code the volatile memory (RAM) of the LINUX server was getting occupied (>90%) which was causing errors in the process and made the task very heavy for the system to execute. The sole reason, as we anticipated, could be the formation of unwanted cache or temporary files in the loop when it opens an NC file in data variable using xarray. So to overcome this general/popular problem of handling and processing large data files, the following few lines were added (shown separately with spaces).

        
        data.close()
        del data, clipped_data, fsca, fsca_binary
        # Explicit garbage collection
        gc.collect()
