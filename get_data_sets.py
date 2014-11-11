import urllib
from random import randint
from time import sleep

#testfile = urllib.URLopener()
#testfile.retrieve("" \
#                  ,"Data/")


#GLOBAL TEMPERATURE SERIES
testfile = urllib.URLopener()
testfile.retrieve("http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.3.0.0.annual_ns_avg.txt" \
                  ,"Data/HadCRUT.4.3.0.0.annual_ns_avg.txt")

#testfile = urllib.URLopener()
#testfile.retrieve("http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt" \
#                  ,"Data/GLB.Ts+dSST.txt")


testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.ncdc.noaa.gov/pub/data/mlost/operational/products/aravg.mon.land_ocean.90S.90N.v3.5.4.201409.asc" \
                  ,"Data/aravg.mon.land_ocean.90S.90N.v3.5.4.201409.asc")

testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.ncdc.noaa.gov/pub/data/mlost/operational/products/aravg.ann.land_ocean.90S.90N.v3.5.4.201409.asc" \
                  ,"Data/aravg.ann.land_ocean.90S.90N.v3.5.4.201409.asc")



#CPC Arctic Oscillation
testfile = urllib.URLopener()
testfile.retrieve("http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii" \
                  ,"Data/monthly.ao.index.b50.current.ascii")


#CPC North Atlantic Oscillation Index
testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/nao_index.tim" \
                  ,"Data/nao_index.tim")


#CPC Southern Oscillation Index
testfile = urllib.URLopener()
testfile.retrieve("http://www.cpc.ncep.noaa.gov/data/indices/soi" \
                  ,"Data/soi")

#CPC Antarctic Oscillation Index
testfile = urllib.URLopener()
testfile.retrieve("http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii" \
                  ,"Data/monthly.aao.index.b79.current.ascii")

#CPC Nino 1
testfile = urllib.URLopener()
testfile.retrieve("http://www.esrl.noaa.gov/psd/data/correlation/nina1.data" \
                  ,"Data/nina1.data")

#CPC nino 4
testfile = urllib.URLopener()
testfile.retrieve("http://www.esrl.noaa.gov/psd/data/correlation/nina4.data" \
                  ,"Data/nina4.data")

#CPC nino 3.4
testfile = urllib.URLopener()
testfile.retrieve("http://www.esrl.noaa.gov/psd/data/correlation/nina34.data" \
                  ,"Data/nina34.data")

#CPC AMO
testfile = urllib.URLopener()
testfile.retrieve("http://www.esrl.noaa.gov/psd/data/correlation/amon.us.long.data" \
                  ,"Data/amon.us.long.data")
