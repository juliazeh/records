#!/usr/bin/env/python

"a package for pulling occurrenct data from GBIF"

import pandas as pd
import requests
import numpy as np


class Records:
    """
Create a Records class instance with GBIF occurrence records stored 
    in a pandas DataFrame for a queried taxon between a range of years. 
    Parameters:
    -----------
    q: str
        Query taxonomic name. 
    interval: tuple
        Range of years to return results for. Should be (min, max) tuple.
    Attributes:
    -----------
    params: The parameter dictionary to filter GBIF search.
    df: Pandas DataFrame with returned records.
    sdf: A view of the 'df' DataFrame selecting only three relevant columns.
    """
    def __init__(self, q=None, interval=None, **kwargs):

        self.params = {
            'q': q,
            'year': ",".join([str(i) for i in interval]),
            'basisOfRecord': "PRESERVED_SPECIMEN",
            'hasCoordinate': "true",
            'hasGeospatialIssue': "false",
            "country": "US",
            "offset": "0",
            "limit": "300",
        }

        self.params.update(kwargs)
        self.df = self._get_all_records()

    
    @property
    def sdf(self):
        """
        Return a copy of the current .df dataframe selecting only the three
        most generally relevant columns: species, year, and stateProvince. 
        This is only meant for viewing and will raise a warning if you try to
        modify it since it is a copy, and thus you would be setting values on 
        a selection of a selection. See pandas docs in the warning for detalis.
        """
        return self.df[["species", "year", "country", "stateProvince"]]

    def _get_all_records(self):
        """
        Obtains records for given taxon name query and interval using params
        dictionary and stores the results in a dataframe.
        """
        baseurl = "http://api.gbif.org/v1/occurrence/search?"
        data = []

        while 1:
            # make request and store results
            res = requests.get(
                url=baseurl,
                params=self.params,
            )
            
            # check for errors
            res.raise_for_status()

            # increment counter
            self.params["offset"] = str(int(self.params["offset"]) + 300)

            # concatenate data
            idata = res.json()
            data += idata["results"]

            # stop when end of record is reached
            if idata["endOfRecords"]:
                break

        return pd.DataFrame(data)

class Epochs:
    """
    Returns an Epochs class instance that includes GBIF occurrence records 
    and stores an extra label with each query that includes the interval 
    (epoch) during which those records were collected, and returns all 
    records in a sorted pandas dataframe. 
    Parameters:
    -----------
    q: str
        Query taxonomic name. 
    start: int
        Earliest year from which to search for records. 
    end: int
        Latest year from which to search for records. 
    epochsize: int
        of years to return results for. Should be (min, max) tuple.
    Attributes:
    -----------
    df: Pandas DataFrame with returned records.
    sdf: A view of the 'df' DataFrame selecting only four relevant columns.
    """ 
    def __init__(self, q, start, end, epochsize, **kwargs):
        
        # make range of epochs
        epochs = range(start, end, epochsize)

        # get Record objects across the epoch range
        rdicts = {
            i: Records(q, (i, i + epochsize), **kwargs) for i in epochs
        }

        # add epoch to each dataframe
        for epoch in rdicts:
            rdicts[epoch].df["epoch"] = epoch

        # if rdicts, then build dataframe, otherwise skip it. 
        if rdicts:

            # concatenate all dataframes into one
            self.df = pd.concat([i.df for i in rdicts.values()])

            # sort values by year, and reset index without keeping old index
            self.df = (
                self.df
                .sort_values(by="year")
                .reset_index(drop=True)
                )
        else:
            self.df = pd.DataFrame([])


    @property
    def sdf(self):
        """
        Return a copy of the current .df dataframe selecting only the three
        most generally relevant columns: species, year, epoch, country, 
        and stateProvince. 
        This is only meant for viewing and will raise a warning if you try to
        modify it since it is a copy, and thus you would be setting values on 
        a selection of a selection. See pandas docs in the warning for detalis.
        """
        return self.df[
            ["species", "year", "epoch", "country", "stateProvince"]]



    def simpsons_diversity(self, by):
        """
        Calculates simpon's diversity index: the probability that any two 
        sampled individuals are the same species. Enter a key for groupby
        as a list of single or multiple keys.
        Parameters
        -----------
        by: str
            A column name used by .groupby to group samples prior to 
            calculating simpson's diversity. For example, enter 
        """
        # group on 'by' keyword, and exclude records missing data for 'by'.
        # and then count species in each group and calculate simp's div.
        data = (
            self.df[self.df[by].notna() & self.df.species.notna()]
            .groupby(by)
            .species
            .apply(calculate_simpsons_diversity)
            )

        # set zero values of simpson's diversity to nan
        data[data == 0] = np.nan
        return data




# utility functions
def load_epochs_from_csv(filepath):
    # init an empty epoch instance
    ep = Epochs("", 0, 0, 1)

    # load existing dataframe to instance's .df attribute
    ep.df = pd.read_csv(filepath, index_col=0)
    return ep


def calculate_simpsons_diversity(arr):
    "internal function to calculate simpson's diversity"
    simps = 0.
    for taxon in arr.unique():
        # proportion of individuals that are this species
        simps += (np.sum(arr == taxon) / arr.shape[0]) ** 2
    return 1. - simps
