# Call libraries
suppressMessages({
  library(xcms); library(tidyverse); library(data.table); library(MSnbase)
  library(CAMERA); library(glycanPredict); library(viridis)
})

ppm_to_mz = function(mz, noise){
  ppm = mz / 1000000 * noise
  return(ppm)
}

Dulce_fetch = function(data, cwp=NULL, pdp=NULL, return_everything=F){
  
  # If no param objects defined, default is used.
  if (is.null(cwp)){cwp = CentWaveParam()}
  else if (!is.null(cwp) & class(cwp)!="CentWaveParam"){
    stop("Dulce error: 'cwp' is not from 'CentWaveParam'. Check for ?CentWaveParam or set it as NULL to use default arguments.")}
  
  if (is.null(pdp)){pdp = PeakDensityParam(sampleGroups=rep("Ungrouped", nrow(data)))}
  else if (!is.null(pdp) & class(pdp)!="PeakDensityParam"){
    stop("Dulce error: 'pdp' is not from 'PeakDensityParam'. Check for ?PeakDensityParam or set it as NULL to use default arguments.")}
  
  # Defining a default "Ungrouped" category if no grouping is defined.
  if (length(pdp@sampleGroups)==0){
    pdp@sampleGroups = rep("Ungrouped", nrow(data))
    message("Dulce warning: No sample groups were defined. One group under the name of 'Ungrouped' will be created.")}
  
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

Dulce_isotAddu = function(data, isotopes=T, adducts=T, 
                          perfwhm=0.5, mzabs=0.01, cor_eic_th=0.75,
                          polarity=NULL){
  
  if (is.null(polarity)){
    stop("Dulce error: NULL polarity? Check it twice.")
  }
  
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
  
  if (class(data)!="xsAnnotate"){stop("Dulce error: 'data' object is not from 'xsAnnotate' class.")}
  
  data = getPeaklist(data) %>% filter(between(rt, rtmin, rtmax))
  
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
                         ppm=NULL, mzabs=NULL){
  
  if (!all(c("mzmin","mzmax") %in% colnames(data))){
    stop("Dulce error: check if 'mzmin' and 'mzmax' are columns in the 'data' object. They are needed to overlap the predictions.")
  }
  
  if (is.null(ppm) & is.null(mzabs)){stop("Dulce error: 'ppm' and 'mzabs' have NULL values. Please specify one.")}
  if ((!is.null(ppm) & !is.null(mzabs))){
    message("Dulce warning: 'ppm' and 'mzabs' were specified (it is one or the other, not both). Using only 'ppm'")
    mzabs = NULL
  }
  
  predicted = predictGlycans(dp1=dp1, dp2=dp2, ESI_mode=ESI_mode, 
                             scan_range1=scan_range1, scan_range2=scan_range2,
                             pent_option=pent_option, modifications=modifications, 
                             label = label) %>% 
    pivot_longer(-c(name, dp, mass, formula), names_to = "ion", values_to = "mz") %>% 
    drop_na()
  
  if (!is.null(mzabs)){
    predicted = predicted %>% mutate(mzmin=mz-mzabs,
                                     mzmax=mz+mzabs)
  } else if (!is.null(ppm)){
    predicted = predicted %>% mutate(mzmin=mz-ppm_to_mz(mz, ppm),
                                     mzmax=mz+ppm_to_mz(mz, ppm))
  } else {
    stop("Dulce error: noise_unit value is not recognized. Please choose between 'mz' and 'ppm'")
  }
  
  setDT(predicted)
  setDT(data)
  setkey(predicted, mzmin, mzmax)
  predicted = foverlaps(data,predicted) %>% drop_na(name)
  
  return(predicted)
}

Dulce_doAll = function(data, cwp=NULL, pdp=NULL,
                       names=NULL, classes="Unclassified",
                       isotopes=T, adducts=T, perfwhm=0.5, mzabs_isotAdd=0.01, 
                       cor_eic_th=0.75, polarity=NULL,
                       rtmin=0, rtmax=Inf,
                       dp1=1, dp2=6, ESI_mode="neg", scan_range1=0, scan_range2=2000, 
                       pent_option=1, label="none",
                       modifications=c("deoxy", "unsaturated", "sulphate"),
                       ppm=NULL, mzabs=NULL,
                       return_everything=F){
  
  data_fetched = Dulce_fetch(data, cwp=cwp, pdp=pdp, return_everything=return_everything)
  
  if (return_everything){
    data_peaks = data_fetched$peaks
    data_features = data_fetched$features
    data_fetched = data_fetched$data
  }
  
  data_converted = Dulce_toCAMERA(data_fetched, names=names, classes=classes)
  data_annotated = Dulce_isotAddu(data_converted, isotopes=isotopes, adducts=adducts, 
                                  perfwhm=perfwhm, mzabs=mzabs_isotAdd, cor_eic_th=cor_eic_th, polarity=polarity)
  data_trimmed = Dulce_trimIsotopes(data_annotated, rtmin=rtmin, rtmax=rtmax)
  data_predicted = Dulce_predict(data_trimmed, dp1=dp1, dp2=dp2, ESI_mode=ESI_mode, 
                                 scan_range1=scan_range1, scan_range2=scan_range2, 
                                 pent_option=pent_option, label=label,
                                 modifications=modifications,
                                 ppm=ppm, mzabs=mzabs)
  
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
file.paths = dir(path="mzML_example3", all.files=F, full.names=T)
pheno.data = data.frame(name=basename(file.paths) %>% 
                          gsub(pat="MS31_20220203_|_\\d{2}.mzML|_100xdilute", 
                               rep=""),
                        sampletype1=basename(file.paths) %>% 
                          gsub(pat=".*blank.*",
                               rep="blank") %>% 
                          gsub(pat=".*lam.*",
                               rep="laminarin") %>% 
                          gsub(pat=".*yeastmannan.*",
                               rep="yeastmannan"),
                        sampletype2=basename(file.paths) %>% 
                          gsub(pat=".*blank.*",
                               rep="blank") %>% 
                          gsub(pat=".*lam_omix.*|.*gh76.*",
                               rep="positive control") %>% 
                          gsub(pat=".*fitdog.*",
                               rep="sample"),
                        rep=basename(file.paths) %>% 
                          gsub(pat=".*rep2.*|.*-2-.*",
                               rep="B") %>% 
                          gsub(pat=".*rep1.*|.*-1-.*",
                               rep="A"))
pheno.data = pheno.data[c(1,2),]
file.paths = file.paths[c(1,2)]
data = readMSData(files=file.paths, 
                  pdata=new("NAnnotatedDataFrame", pheno.data),
                  mode="onDisk")

# Filter by only MS data (not MS/MS)
data = data[data@featureData@data$msLevel==1]

register(SerialParam())

data_processed = Dulce_doAll(data, polarity="negative", ppm=10)







