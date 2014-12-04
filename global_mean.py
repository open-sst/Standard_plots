import numpy as np
#import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from time_series import *
from read_data_sets import *

def combine_series(series1, series2, series3):
    #average together three annual time series and return an annual time series
    #hardcoded to do this from 1880 to 2014. It takes the uncertainty estimates
    #from the first of the three series.
    comb_year = []
    comb_anom = []
    comb_upper_unc = []
    comb_lower_unc = []

    for year in range(1880,2015):
            
        ind1 = series1.times.index(year)
        ind2 = series2.times.index(year)
        ind3 = series3.times.index(year)

        av = np.mean([series1.data[ind1], series2.data[ind2], series3.data[ind3] ])

        comb_year.append(year)
        comb_anom.append(av)
        comb_upper_unc.append(av+series1.hiunc[ind1]-series1.data[ind1])
        comb_lower_unc.append(av+series1.lounc[ind1]-series1.data[ind1])

        combined = time_series(comb_year,
                               comb_anom,
                               comb_lower_unc,
                               comb_upper_unc)


    return combined




###########################
###########################
###########################
###########################
##
##      MAIN STRIP
##
###########################
###########################
###########################
###########################


ncdc_version = "v3.5.4.201410"
hadsst_version = "3.1.1.0"
crutem_version = "4.3.0.0"
hadcrut_version = "4.3.0.0"
cowtan_and_way_version = "v2_0_0"


print("GLOBAL AVERAGE TEMPERATURES")

era_monthly = read_era_interim()
era_ts = era_monthly.annualise(10)
era_ts.add_name("ERA")
era_ts.add_offset(0.3)

cwh_ts = read_cowtan_and_way_hybrid_monthly(cowtan_and_way_version)
cwh_ts = cwh_ts.annualise()
cwh_ts.add_name("Cowtan and Way Hybrid")
cwh_ts.add_offset(0.3)

cw_ts = read_cowtan_and_way_monthly(cowtan_and_way_version)
cw_ts = cw_ts.annualise()
cw_ts.add_name("Cowtan and Way")

had_ts = read_hadcrut4(hadcrut_version)
had_ts.add_name("HadCRUT."+hadcrut_version)

ncdc_ts = read_ncdc(ncdc_version)
ncdc_ts.rebaseline(1961,1990)
ncdc_ts.add_name("MLOST")

giss_ts = read_giss()
giss_ts.rebaseline(1961,1990)
giss_ts.add_name("GISTEMP")

#combine HadCRUT, NCDC and GISS to get the WMO combined series
#which uses HadCRUT uncertainty estimates
combined_ts = combine_series(had_ts, ncdc_ts, giss_ts)
combined_ts.add_name("Combined")



#plot the global temperatures with the non-linear trend
#estimate provided by Prof David Stephenson of the University
#of Exeter, as used in the IPCC report
combined_trend_ts = read_stephenson_trends()

combined_ts.plot_ts('Black')
combined_trend_ts.plot_ts_with_unc('Red','Pink')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.axis((1875,2016,-0.6,0.7))
plt.xticks(range(1880,2020,20), fontsize = 16)
plt.yticks(np.arange(-0.6,0.8,0.2), fontsize = 16)
#plt.text(1900,0.0,"DRAFT",fontsize=100,color="AliceBlue",alpha=0.5)
plt.text(1940,-0.58,"Analysis provided by Prof. David Stephenson, University of Exeter",fontsize=10)

plt.show()

#print out the whole of the combined time series
combined_ts.print_ts()

#plot variety of different views
combined_ts.modded_skyscraper_diagram()
combined_ts.plot_ranking_diagram_mono(7,1)
ninos = read_oni_nino_categories()
combined_ts.plot_skyscraper_diagram(ninos)
combined_ts.plot_ts_highlight_nino(ninos)

#make diagram showing all three data sets
had_ts.plot_ts_with_unc('Black','AliceBlue')
giss_ts.plot_ts('Red')
ncdc_ts.plot_ts('DeepSkyBlue')
cw_ts.plot_ts('Green')
cwh_ts.plot_ts('Pink')
era_ts.plot_ts('Gold')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.legend(loc='best')
plt.xticks(range(1860,2020,20), fontsize = 16)
plt.yticks(np.arange(-0.8,0.8,0.2), fontsize = 16)
#plt.savefig('gmt.png', bbox_inches='tight')
plt.show()


#print the running mean of the combined time series
print("Running mean")
combined_ts.print_running_mean(10)

#print the 5 warmest years in the combined and other data sets
combined_ts.print_ordered_ts(5)

era_ts.print_ordered_ts(5)
had_ts.print_ordered_ts(5)
ncdc_ts.print_ordered_ts(5)
giss_ts.print_ordered_ts(5)
cw_ts.print_ordered_ts(5)
cwh_ts.print_ordered_ts(5)


print("")
print("Global average SST")

ncdc_sst_ts = read_ncdc_sst(ncdc_version)
ncdc_sst_ts.rebaseline(1961,1990)
ncdc_sst_ts.add_name("ERSSTv3")
ncdc_sst_ts.print_ordered_ts(5)


hadsst_ts = read_hadsst3(hadsst_version)
latest_year = np.mean([0.342, 0.314, 0.347, 0.478, 0.477, 0.563, 0.551, 0.644, 0.574, 0.529])
hadsst_ts.add_year(2014,latest_year,latest_year-0.1,latest_year+0.1)
hadsst_ts.add_name("HadSST."+hadsst_version)
hadsst_ts.print_ordered_ts(5)


hadsst_ts.plot_ts_with_unc('Black','AliceBlue')
ncdc_sst_ts.plot_ts_with_unc('DeepSkyBlue','Yellow')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.legend(loc='best')
plt.show()


print("")
print("Global average LSAT")

ncdc_lsat_ts = read_ncdc_lsat(ncdc_version)
ncdc_lsat_ts.rebaseline(1961,1990)
ncdc_lsat_ts.add_name("NCDC LSAT")
ncdc_lsat_ts.print_ordered_ts(10)

print("")

crutem4_ts = read_crutem4(crutem_version)
crutem4_ts.add_name("CRUTEM."+crutem_version)
crutem4_ts.print_ordered_ts(10)


crutem4_ts.plot_ts_with_unc('Black','AliceBlue')
ncdc_lsat_ts.plot_ts_with_unc('DeepSkyBlue','Yellow')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.legend(loc='best')
plt.show()


print("")
print("Year to date comparisons")

era_monthly = read_era_interim()
era_ts = era_monthly.annualise(10)
era_ts.add_name("ERA")

had_monthly = read_hadcrut4_monthly(hadcrut_version)
had_ts = had_monthly.annualise(10)
had_ts.add_name("HadCRUT")

ncdc_monthly = read_ncdc_monthly(ncdc_version)
ncdc_monthly.rebaseline(1961,1990)
ncdc_ts = ncdc_monthly.annualise(10)
ncdc_ts.add_name("NCDC")

giss_monthly = read_giss_monthly()
giss_monthly.rebaseline(1961,1990)
giss_ts = giss_monthly.annualise(10)
giss_ts.add_name("GISTEMP")

combined_ts = combine_series(had_ts, ncdc_ts, giss_ts)
combined_ts.add_name("Combined")

era_ts.print_ordered_ts(5)
had_ts.print_ordered_ts(5)
ncdc_ts.print_ordered_ts(5)
giss_ts.print_ordered_ts(5)
combined_ts.print_ordered_ts(5)
