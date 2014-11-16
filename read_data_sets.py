import numpy as np
from time_series import *

def read_jisao_pdo():
    f=open("Data/PDO.latest",'r')

    for i in range(1,32):
        f.readline()

    pdo_year = []
    pdo_month = []
    pdo_data = []

    for line in f:
        line = line.strip()
        if line != "":
            columns = line.split()
            year = columns[0].split('*')
            for m in range(1,len(columns)):
                pdo_year.append(float(year[0]))
                pdo_month.append(float(m))
                pdo_data.append(float(columns[m]))
        else:
            break
            
    f.close()
    
    pdo_ts = monthly_time_series(pdo_year, pdo_month, pdo_data)

    return pdo_ts

def read_uah(version="5.6"):
    f=open("Data/uahncdc_lt_"+version+".txt",'r')

    uah_year = []
    uah_month = []
    uah_data = []

    line = f.readline()
    readon = 1

    for line in f:
        line = line.strip()
        columns = line.split()
        if readon == 1:
            if columns[0] == "Year":
                readon = 0

            if readon:
                uah_year.append(int(columns[0]))
                uah_month.append(int(columns[1]))
                uah_data.append(float(columns[2]))

    f.close()
    
    uah_ts = monthly_time_series(uah_year, uah_month, uah_data)

    return uah_ts

def read_hadley_monthly(filename):
    f = open(filename, 'r')
    hadcrut_year = []
    hadcrut_anom = []
    hadcrut_month = []

    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        ym = columns[0].split('/')
        hadcrut_year.append(float(ym[0]))
        hadcrut_month.append(float(ym[1]))
        hadcrut_anom.append(float(columns[1]))
    
    f.close()
    had_ts = monthly_time_series(hadcrut_year,
                                 hadcrut_month,
                                 hadcrut_anom)

    return had_ts
    
def read_hadcrut4_monthly(version):
    return read_hadley_monthly('Data/HadCRUT.'+version+'.monthly_ns_avg.txt')

def read_hadley(filename):
#read hadley format annual data sets and make annual time series out of them
    f = open(filename, 'r')

    hadcrut_year = []
    hadcrut_anom = []
    hadcrut_upper_unc = []
    hadcrut_lower_unc = []

    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        hadcrut_year.append(float(columns[0]))
        hadcrut_anom.append(float(columns[1]))
        hadcrut_upper_unc.append(float(columns[11]))
        hadcrut_lower_unc.append(float(columns[10]))
    
    f.close()
    had_ts = time_series(hadcrut_year,
                         hadcrut_anom,
                         hadcrut_lower_unc,
                         hadcrut_upper_unc)

    return had_ts


def read_hadsst3(version):
    return read_hadley('Data/HadSST.'+version+'_annual_globe_ts.txt')

def read_crutem4(version):
    return read_hadley('Data/CRUTEM.'+version+'.global_n+s')

def read_hadcrut4(version):
    return read_hadley('Data/HadCRUT.'+version+'.annual_ns_avg.txt')


def read_ncdc_format_monthly(filename):
    f = open(filename, 'r')
    ncdc_year = []
    ncdc_month = []
    ncdc_anom = []

    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        ncdc_year.append(float(columns[0]))
        ncdc_month.append(float(columns[1]) )
        ncdc_anom.append(float(columns[2]))

    f.close()

    ncdc_ts = monthly_time_series(ncdc_year,
                                  ncdc_month,
                                  ncdc_anom)

    return ncdc_ts
    
def read_ncdc_monthly(version):
    return read_ncdc_format_monthly('Data/aravg.mon.land_ocean.90S.90N.'+version+'.asc')

def read_ncdc_format(filename):
    f = open(filename, 'r')
    ncdc_year = []
    ncdc_anom = []
    ncdc_lounc = []
    ncdc_hiunc = []

    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        ncdc_year.append(float(columns[0]))
        ncdc_anom.append(float(columns[1]))
        ncdc_lounc.append(float(columns[1]) - 2 * np.sqrt(float(columns[2])) )
        ncdc_hiunc.append(float(columns[1]) + 2 * np.sqrt(float(columns[2])) )

    f.close()

    ncdc_ts = time_series(ncdc_year,
                          ncdc_anom,
                          ncdc_lounc,
                          ncdc_hiunc)

    return ncdc_ts


def read_ncdc_sst(version):
    return read_ncdc_format('Data/aravg.ann.ocean.90S.90N.'+version+'.asc')

def read_ncdc_lsat(version):
    return read_ncdc_format('Data/aravg.ann.land.90S.90N.'+version+'.asc')

def read_ncdc(version):
    return read_ncdc_format('Data/aravg.ann.land_ocean.90S.90N.'+version+'.asc')


def read_giss_block_monthly(f, block_length, giss_year, giss_month, giss_anom):
#burn first two pointless lines of each block
    line = f.readline()
    line = f.readline()
    for j in range(1,block_length+1):
        line = f.readline()
        columns = line.split()
        for i in range(1,13):
            giss_year.append(float(columns[0]))
            giss_month.append(float(i))
            giss_anom.append(float(columns[i])/100.)
            
    return (giss_year, giss_month, giss_anom)

def read_giss_monthly():
#    f = open('Data/GLB.Ts+dSST.txt','r')
    f = open('Data/GLB.TsERSST.GHCN.CL.PA.txt','r')

    giss_year = []
    giss_month = []
    giss_anom = []
  
#read header information and discard
    for i in range(1,7):
        f.readline()

#read first block of 21 years
    giss_year, giss_month, giss_anom = read_giss_block_monthly(f, 21, giss_year, giss_month, giss_anom)
#read five blocks of 20 years   
    for i in range(1,6):
        giss_year, giss_month, giss_anom = read_giss_block_monthly(f, 20, giss_year, giss_month, giss_anom)
#final block has less than 20
    giss_year, giss_month, giss_anom = read_giss_block_monthly(f, 13, giss_year, giss_month, giss_anom)

#final year is incomplete so calculate from monthlies
    g = f.readline()
    columns = g.split()

    for i in range(1,13):
        if columns[i] != '****':
            giss_year.append(float(columns[0]))
            giss_month.append(float(i))
            giss_anom.append(float(columns[i])/100.)
    
    giss_ts = monthly_time_series(giss_year,
                                  giss_month,
                                  giss_anom)

    return giss_ts

def read_giss_block(f, block_length, giss_year, giss_anom):
    line = f.readline()
    line = f.readline()
    for j in range(1,block_length+1):
        line = f.readline()
        columns = line.split()
        giss_year.append(float(columns[0]))
        giss_anom.append(float(columns[13])/100.)
    return (giss_year,giss_anom)


def read_giss():
#    f = open('Data/GLB.Ts+dSST.txt','r')
    f = open('Data/GLB.TsERSST.GHCN.CL.PA.txt','r')
    
    giss_year = []
    giss_anom = []

#read header information and discard
    for i in range(1,7):
        f.readline()

#read first block of 21 years
    giss_year, giss_anom = read_giss_block(f, 21, giss_year, giss_anom)
#read five blocks of 20 years   
    for i in range(1,6):
        giss_year, giss_anom = read_giss_block(f, 20, giss_year, giss_anom)
#final block has less than 20
    giss_year, giss_anom = read_giss_block(f, 13, giss_year, giss_anom)

#final year is incomplete so calculate from monthlies
    g = f.readline()
    columns = g.split()
    giss_year.append(float(columns[0]))
    final_year = []
    for i in range(1,13):
        if columns[i] != '****':
            final_year.append(float(columns[i])/100.)
    giss_anom.append(np.mean(final_year))

    giss_lounc = []
    giss_hiunc = []

    for i in range(0,len(giss_year)):
        if giss_year[i] >= 1880 and giss_year[i] <= 1900:
            giss_lounc.append(giss_anom[i]-0.08)
            giss_hiunc.append(giss_anom[i]+0.08)
        elif giss_year[i] >= 1901 and giss_year[i] <= 1950:
            giss_lounc.append(giss_anom[i]-0.05)
            giss_hiunc.append(giss_anom[i]+0.05)
        elif giss_year[i] >= 1951 and giss_year[i] <= 3000:
            giss_lounc.append(giss_anom[i]-0.05)
            giss_hiunc.append(giss_anom[i]+0.05)
    
    giss_ts = time_series(giss_year,
                          giss_anom,
                          giss_lounc,
                          giss_hiunc)

    return giss_ts


def read_cowtan_and_way(version):
    f = open('Data/had4_krig_annual_'+version+'.txt','r')
    
    cw_year = []
    cw_anom = []
    cw_lounc = []
    cw_hiunc = []

    for line in f:
        columns = line.split()
        cw_year.append(int(columns[0]))
        cw_anom.append(float(columns[1]))
        cw_lounc.append(float(columns[1]) - 2* float(columns[2]))
        cw_hiunc.append(float(columns[1]) + 2* float(columns[2]))
        
    cw_ts = time_series(cw_year,
                          cw_anom,
                          cw_lounc,
                          cw_hiunc)

    return cw_ts

def read_cowtan_and_way_hybrid(version):
    f = open('Data/had4_short_uah_annual_'+version+'.txt','r')
    
    cw_year = []
    cw_anom = []
    cw_lounc = []
    cw_hiunc = []

    for line in f:
        columns = line.split()
        cw_year.append(int(columns[0]))
        cw_anom.append(float(columns[1]))
        cw_lounc.append(float(columns[1]) - 2* float(columns[2]))
        cw_hiunc.append(float(columns[1]) + 2* float(columns[2]))
        
    cw_ts = time_series(cw_year,
                          cw_anom,
                          cw_lounc,
                          cw_hiunc)

    return cw_ts
