#!/usr/bin/env/python

import pandas as pd
import requests


class Records:
    """
    Returns a dataframe with all records from GBIF for given taxon and interval.
    Format: Records(q=taxon, interval=(years))
    """
    def __init__(self, q=None, interval=None):

        self.q = q
        self.interval = interval
        self.params = {
            "q": self.q,
            "year": self.interval,
            "offset": "0"
        }
        self.df = self._get_all_records()

    def _get_all_records(self):
        """
        Obtains records for given taxon name query and interval using params
        dictionary and stores the results in a dataframe.
        """
        baseurl = "http://api.gbif.org/v1/occurrence/search?"
        start = 0
        data = []

        while 1:
            # make request and store results
            res = requests.get(
                url=baseurl,
                params=self.params,
            )
            # increment counter
            self.params["offset"] = str(int(self.params["offset"]) + 300)

            # concatenate data
            idata = res.json()
            data += idata["results"]

            # stop when end of record is reached
            if idata["endOfRecords"]:
                break

        return pd.DataFrame(data)
