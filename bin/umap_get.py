import pandas as pd
import umap
import os

def umap_process(covid_dist):

    u = umap.UMAP(metric = "precomputed", random_state = 2020)
    covid_umap = u.fit_transform(covid_dist)
    covid_umap_df = pd.DataFrame(data = covid_umap)
    return(covid_umap_df)


