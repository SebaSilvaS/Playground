# Call libraries
suppressMessages({
  library(xcms); library(tidyverse); library(data.table); library(MSnbase)
  library(CAMERA); library(glycanPredict); library(viridis)
})

ppm_to_mz = function(mz, noise){
  ppm = mz / 1000000 * noise
  return(ppm)
}

Dulce_fetch = function(data, cwp=NULL, pdp=NULL, 
                 return_everything=T){
  
  # Defining a default "Ungrouped" category if no grouping is defined.
  if (is.null(sampleGroups)){
    sampleGroups = rep("Ungrouped", nrow(data))
  }
  
  # Fetch 
  processed_data = findChromPeaks(data, param=cwp) %>% groupChromPeaks(param=pdp)
  
  message("Peaks picked and grouped!")
  
  if (return_everything){
    peaks_data = as.data.frame(processed_data@msFeatureData[["chromPeaks"]])
    features_data = as.data.frame(processed_data@msFeatureData[["featureDefinitions"]])
    
    return(list(data=processed_data, 
                peaks=peaks_data, 
                features=features_data))
  }
  return(processed_data)
}


foo = function(x, list){
  x = x+1
  
  
}


Dulce_toCAMERA = function(data, names=NULL, classes="Unclassified"){
  
  if (class(data)!="XCMSnExp"){stop("Dulce error: 'data' object is not from 'XCMSnExp' class.")}
  data_converted = as(data, "xcmsSet")
  
  if (is.null(names)){
    names = sprintf("sample_%03d", 1:nrow(data))
  }
  sampnames(data_converted) = names
  sampclass(data_converted) = classes
  
  message("Dulce note: Do not worry, 'sampclass' and 'sampnames' have already been set.")
  return(data_converted)
}

Dulce_annotate = function(data, isotopes=T, adducts=T, 
                          perfwhm=0.5, mzabs=0.01, cor_eic_th=0.75,
                          polarity=NULL){
  
  data = xsAnnotate(data) %>% groupFWHM(perfwhm = perfwhm)
  
  if (isotopes){
    data = data %>% findIsotopes(mzabs=mzabs) %>% 
      groupCorr(cor_eic_th=cor_eic_th) 
  }
  
  if (adducts){
    data = data %>% findAdducts(polarity=polarity)
  }
  
  return(data)
}


Dulce_trimIsotopes = function(data, rtmin=0, rtmax=Inf){
  
  library(data.table)
  if (class(data)!="xsAnnotate"){stop("Dulce error: 'data' object is not from 'xsAnnotate' class.")}
  
  data = getPeaklist(data)
  
  data = data %>% filter(between(rt, rtmin, rtmax))
  
  data_isotopes = data %>% filter(isotopes!="") %>% 
    mutate(isotope_group=sub(x=isotopes, pat="\\[M.*", rep="")) %>% 
    group_by(isotope_group) %>% 
    mutate(isotopes = paste(isotopes, collapse=",")) %>% 
    distinct(isotope_group, .keep_all=T)
  
  data = data %>% filter(isotopes=="") %>% bind_rows(data_isotopes) %>% 
    dplyr::select(-isotope_group)
}

Dulce_predict = function(data, dp1=1, dp2=6, ESI_mode="neg", 
                         scan_range1=0, scan_range2=2000, 
                         pent_option=1, label="none",
                         modifications=c("deoxy", "unsaturated", "sulphate"),
                         noise=0.001, noise_unit="mz"){
  
  if (!all(c("mzmin","mzmax") %in% colnames(data))){
    stop("Dulce error: check if 'mzmin' and 'mzmax' are columns in the 'data' object. They are needed to overlap the predictions.")
  }
  
  predicted = predictGlycans(dp1=dp1, dp2=dp2, ESI_mode=ESI_mode, 
                             scan_range1=scan_range1, scan_range2=scan_range2,
                             pent_option=pent_option, modifications=modifications, 
                             label = label) %>% 
    pivot_longer(-c(name, dp, mass, formula), names_to = "ion", values_to = "mz") %>% 
    drop_na()
  
  if (noise_unit=="mz"){
    predicted = predicted %>% mutate(mzmin=mz-noise,
                                     mzmax=mz+noise)
  } else if (noise_unit=="ppm"){
    predicted = predicted %>% mutate(mzmin=mz-ppm_to_mz(mz, noise),
                                     mzmax=mz+ppm_to_mz(mz, noise))
  } else {
    stop("Dulce error: noise_unit value is not recognized. Please choose between 'mz' and 'ppm'")
  }
  
  setDT(predicted)
  setDT(data)
  setkey(predicted, mzmin, mzmax)
  predicted = foverlaps(data,predicted) %>% drop_na(name)
  
  return(predicted)
}

Dulce_doAll = function(data,
                       ppm=10, peakwidth=c(5,50), snthresh=10, prefilter=c(3,1000), 
                       sampleGroups=NULL, fetch_noise=5000, binSize=0.01, bw=5, minFraction=0.5,
                       names=NULL, classes="Unclassified",
                       isotopes=T, adducts=T, perfwhm=0.5, mzabs=0.01, 
                       cor_eic_th=0.75, polarity=NULL,
                       rtmin=0, rtmax=Inf,
                       dp1=1, dp2=6, ESI_mode="neg", scan_range1=0, scan_range2=2000, 
                       pent_option=1, label="none",
                       modifications=c("deoxy", "unsaturated", "sulphate"),
                       predict_noise=0.001, noise_unit="mz",
                       return_everything=F){
  
  data_fetched = Dulce_fetch(data, ppm=ppm, peakwidth=peakwidth, snthresh=snthresh, prefilter=prefilter, 
                     sampleGroups=sampleGroups, noise=fetch_noise, binSize=binSize, bw=bw, 
                     minFraction=minFraction, return_everything=return_everything)
  
  if (return_everything){
    data_peaks = data_fetched$peaks
    data_features = data_fetched$features
    data_fetched = data_fetched$data
  }
  
  data_converted = Dulce_toCAMERA(data_fetched, names=names, classes=classes)
  data_annotated = Dulce_annotate(data_converted, isotopes=isotopes, adducts=adducts, 
                                  perfwhm=perfwhm, mzabs=mzabs, cor_eic_th=cor_eic_th, polarity=polarity)
  data_trimmed = Dulce_trimIsotopes(data_annotated, rtmin=rtmin, rtmax=rtmax)
  data_predicted = Dulce_predict(data_trimmed)
  
  if (return_everything){
    return(list(data_fetched=data_fetched, 
                data_peaks=data_peaks, 
                data_features=data_features,
                data_converted=data_converted,
                data_annotated=data_annotated,
                data_trimmed=data_trimmed,
                data_predicted=data_predicted))
  }
  
  return(data_predicted)
}


# Sample code -------------------------------------------------------------
file.paths = dir(path="mzML_example2", all.files=F, full.names=T)
pheno.data = data.frame(name=basename(file.paths) %>% 
                          gsub(pat="MS31_20220203_|_\\d{2}.mzML|_100xdilute", 
                               rep=""),
                        sampletype=basename(file.paths) %>% 
                          gsub(pat=".*fulldigest.*",
                               rep="digest") %>% 
                          gsub(pat=".*solventblank.*",
                               rep="solvent_blank"))


pheno.data = pheno.data[c(1,2),]
file.paths = file.paths[c(1,2)]

data = readMSData(files=file.paths, 
                  pdata=new("NAnnotatedDataFrame", pheno.data),
                  mode="onDisk")

# Filter by only MS data (not MS/MS)
data = data[data@featureData@data$msLevel==1]

register(SerialParam())

lista = Dulce_fetch(data, sampleGroups=data$sampletype, minFraction=0.4)

data_converted = Dulce_toCAMERA(list$data)

data_annotated = Dulce_annotate(data_converted, polarity="negative")

data_formargot = getPeaklist(data_annotated)

data_trimmed = Dulce_trimIsotopes(data_annotated, rtmin=300, rtmax=2500)


data_predicted = Dulce_predict(data_trimmed)


data_predicted = Dulce_doAll(data)


a = Dulce_fetch(data)







