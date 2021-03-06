import urllib
from random import randint
from time import sleep

#testfile = urllib.URLopener()
#testfile.retrieve("" \
#                  ,"Data/")


ncdc_version = "v3.5.4.201410"
hadsst_version = "3.1.1.0"
crutem_version = "4.3.0.0"
hadcrut_version = "4.3.0.0"
cowtan_and_way_version = "v2_0_0"
uah_version = "5.6"

###################
## CPC ONI / Nino indices
###################
testfile = urllib.URLopener()
testfile.retrieve("http://www.cpc.ncep.noaa.gov/data/indices/3mth.nino34.81-10.ascii.txt" \
                  ,"Data/3mth.nino34.81-10.ascii.txt")

testfile = urllib.URLopener()
testfile.retrieve("http://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt" \
                  ,"Data/oni.ascii.txt")

testfile = urllib.URLopener()
testfile.retrieve("http://www.cpc.ncep.noaa.gov/data/indices/ersst3b.nino.mth.81-10.ascii" \
                  ,"Data/ersst3b.nino.mth.81-10.ascii")



assert False

#################
## Cowtan and Way
#################
testfile = urllib.URLopener()
testfile.retrieve("http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_annual_"+cowtan_and_way_version+".txt" \
                  ,"Data/had4_krig_annual_"+cowtan_and_way_version+".txt")


testfile = urllib.URLopener()
testfile.retrieve("http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_"+cowtan_and_way_version+".txt" \
                  ,"Data/had4_krig_"+cowtan_and_way_version+".txt")

#Hybrid version from 1979
testfile = urllib.URLopener()
testfile.retrieve("http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_short_uah_"+cowtan_and_way_version+".txt" \
                  ,"Data/had4_short_uah_"+cowtan_and_way_version+".txt")

#Hybrid annual 
testfile = urllib.URLopener()
testfile.retrieve("http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_short_uah_annual_"+cowtan_and_way_version+".txt" \
                  ,"Data/had4_short_uah_annual_"+cowtan_and_way_version+".txt")


assert False

################
## NCDC
################
#NCDC LSAT
testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.ncdc.noaa.gov/pub/data/mlost/operational/products/aravg.ann.land.90S.90N."+ncdc_version+".asc" \
                  ,"Data/aravg.ann.land.90S.90N."+ncdc_version+".asc")

#NCDC SST
testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.ncdc.noaa.gov/pub/data/mlost/operational/products/aravg.ann.ocean.90S.90N."+ncdc_version+".asc" \
                  ,"Data/aravg.ann.ocean.90S.90N."+ncdc_version+".asc")

#NCDC Monthly GMT
testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.ncdc.noaa.gov/pub/data/mlost/operational/products/aravg.mon.land_ocean.90S.90N."+ncdc_version+".asc" \
                  ,"Data/aravg.mon.land_ocean.90S.90N."+ncdc_version+".asc")

#NCDC Annual GMT
testfile = urllib.URLopener()
testfile.retrieve("ftp://ftp.ncdc.noaa.gov/pub/data/mlost/operational/products/aravg.ann.land_ocean.90S.90N."+ncdc_version+".asc" \
                  ,"Data/aravg.ann.land_ocean.90S.90N."+ncdc_version+".asc")

assert False

##############
## JISAO PDO
##############

testfile = urllib.URLopener()
testfile.retrieve("http://jisao.washington.edu/pdo/PDO.latest" \
                  ,"Data/PDO.latest")

##############
## UAH MSU
##############

testfile = urllib.URLopener()
testfile.retrieve("http://www.nsstc.uah.edu/data/msu/t2lt/uahncdc_lt_"+uah_version+".txt" \
                  ,"Data/uahncdc_lt_"+uah_version+".txt")


########################
## Hadley Centre / UEA
########################
#CRUTEM4
testfile = urllib.URLopener()
testfile.retrieve("http://www.metoffice.gov.uk/hadobs/crutem4/data/diagnostics/global/nh+sh/CRUTEM."+crutem_version+".global_n+s" \
                  ,"Data/CRUTEM."+crutem_version+".global_n+s")

#HadSST3
testfile = urllib.URLopener()
testfile.retrieve("http://www.metoffice.gov.uk/hadobs/hadsst3/data/HadSST."+hadsst_version+"/diagnostics/HadSST."+hadsst_version+"_annual_globe_ts.txt" \
                  ,"Data/HadSST."+hadsst_version+"_annual_globe_ts.txt")

#HadCRUT 4 GLOBAL TEMPERATURE SERIES
testfile = urllib.URLopener()
testfile.retrieve("http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT."+hadcrut_version+".annual_ns_avg.txt" \
                  ,"Data/HadCRUT."+hadcrut_version+".annual_ns_avg.txt")

testfile = urllib.URLopener()
testfile.retrieve("http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT."+hadcrut_version+".monthly_ns_avg.txt" \
                  ,"Data/HadCRUT."+hadcrut_version+".monthly_ns_avg.txt")


                  
assert False


#GISTEMP
#Do not use: returns 403 Forbidden error
#testfile = urllib.URLopener()
#testfile.retrieve("http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt" \
#                  ,"Data/GLB.Ts+dSST.txt")







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
