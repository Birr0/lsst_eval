from db import query
from astropy.time import Time
from tqdm import tqdm
import pandas as pd
import os

def survey_years():
    one_year_timespan = query(
        '''
        select min(midPointMjdTai), max(midPointMjdTai)
        from dp03_catalogs_1yr.DiaSource as DiaSource
        '''
    )
    ten_year_timespan = query(
        '''
        select min(midPointMjdTai), max(midPointMjdTai)
        from dp03_catalogs_10yr.DiaSource as DiaSource
    ''')

    catalog_delta = one_year_timespan["max"].data - one_year_timespan["min"].data # Time epoch for first year catalog

    year_dates = []
    year_dates.append(Time([one_year_timespan["min"].data, \
                    one_year_timespan["max"].data], format="mjd"))

    for i in range(1, 9): # Iterate through all DR - split into years
        if i == 8:
            year_dates.append(Time([one_year_timespan["min"].data  + (i)*catalog_delta, \
                    ten_year_timespan["max"].data], format="mjd"))
        else:
            year_dates.append(Time([one_year_timespan["min"].data  + (i)*catalog_delta, \
                        one_year_timespan["max"].data + (i)*catalog_delta], format="mjd"))
    
    return year_dates

def get_yearly_photometry():
    year_dates = survey_years()
    for i, date in tqdm(enumerate(year_dates)):
        if not os.path.isdir("./DR"):
            os.mkdir(f"./DR/")
            os.mkdir(f"./DR/DR{i + 1}")
        
        for j in tqdm(range(1, 13)):
            delta = ((date[1][0].value - date[0][0].value) / 100) # Split into months
            start_date = date[0][0].value + (j - 1) * delta
            end_date = date[0][0].value + j * delta 

            df = query(
                stmt = f'''
                select mag, magErr, band, phaseAngle, DiaSource.ssObjectId, midPointMjdTai
                from dp03_catalogs_10yr.DiaSource as DiaSource
                join dp03_catalogs_10yr.SSSource as sss on sss.ssObjectId = DiaSource.ssObjectId
                where midPointMjdTai >= {start_date} and midPointMjdTai <= {end_date}
                LIMIT 10000
            '''
            ).to_table().to_pandas()

            df.to_csv(f"./DR/DR{i + 1}/{j}_photometry_{start_date}-{end_date}")
            del df, start_date, end_date, delta # Clear memory
        break

#get_yearly_photometry()