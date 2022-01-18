# Visualise
tar_visnetwork()

# Run all models
#tar_make()

# In parallel - use background process
tar_make_future(workers = 3)#, callr_function = callr::r_bg)
tar_make_future(workers = 3, callr_function = callr::r_bg)

#tar_meta(fields = warnings)
