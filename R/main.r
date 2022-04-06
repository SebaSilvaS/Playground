# Call libraries
suppressMessages({
  library(xcms); library(tidyverse); library(data.table); library(MSnbase)
  library(CAMERA); library(glycanPredict); library(viridis); library(gridExtra)
})

devtools::install_github("margotbligh/glycanPredict")

ppm_to_mz = function(mz, noise){
  ppm = mz / 1000000 * noise
  return(ppm)
}

predictGlycansParams = setClass("predictGlycansParams",
                                slots=list(dp1="numeric", 
                                           dp2="numeric", 
                                           ESI_mode="character", 
                                           scan_range1="numeric", 
                                           scan_range2="numeric", 
                                           pent_option="numeric", 
                                           label="character",
                                           modifications="character"),
                                prototype=prototype(dp1=1, 
                                                    dp2=6, 
                                                    ESI_mode="neg", 
                                                    scan_range1=0, 
                                                    scan_range2=2000, 
                                                    pent_option=1, 
                                                    label="none",
                                                    modifications=c("deoxy", "unsaturated", "sulphate")))

a = CentWaveParam()
ppm(a)

Dulce_fetch = function(data, cwp=NULL, pdp=NULL, return_everything=F){
  
  # If no param objects defined, default is used.
  if (is.null(cwp)){cwp = CentWaveParam()}
  else if (!is.null(cwp) & class(cwp)!="CentWaveParam"){
    stop("Dulce error: 'cwp' is not an 'CentWaveParam' object. Check for ?CentWaveParam or set it as NULL to use default arguments.")}
  
  if (is.null(pdp)){pdp = PeakDensityParam(sampleGroups=rep("Ungrouped", nrow(data)))}
  else if (!is.null(pdp) & class(pdp)!="PeakDensityParam"){
    stop("Dulce error: 'pdp' is not an 'PeakDensityParam' object. Check for ?PeakDensityParam or set it as NULL to use default arguments.")}
  
  # Defining a default "Ungrouped" category if no grouping is defined.
  if (length(pdp@sampleGroups)==0){
    pdp@sampleGroups = rep("Ungrouped", nrow(data))
    message("Dulce warning: No sample groups were defined. .,name of 'Ungrouped' will be created.")}
  
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
  suppressMessages({data_converted = as(data, "xcmsSet")})
  
  if (is.null(names)){
    names = sprintf("sample_%03d", 1:nrow(data))
  }
  sampnames(data_converted) = names
  sampclass(data_converted) = classes

  return(data_converted)
}

Dulce_find = function(data, isotopes=T, adducts=T, 
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

Dulce_annotate = function(data, pgp=NULL,
                           ppm=NULL, mzabs=NULL){
    
  if (!all(c("mzmin","mzmax") %in% colnames(data))){
    stop("Dulce error: check if 'mzmin' and 'mzmax' are columns in the 'data' object. They are needed to overlap the predictions.")
  }
  
  if (is.null(ppm) & is.null(mzabs)){stop("Dulce error: 'ppm' and 'mzabs' have NULL values. Please specify one.")}
  if ((!is.null(ppm) & !is.null(mzabs))){
    message("Dulce warning: 'ppm' and 'mzabs' were specified (it is one or the other, not both). Using only 'ppm'")
    mzabs = NULL
  }
  
  if (is.null(pgp)){pgp = predictGlycansParams()}
  else if (!is.null(pgp) & class(pgp)!="predictGlycansParams"){
    stop("Dulce error: 'cwp' is not an 'predictGlycansParams' object. Check for ?predictGlycansParams or set it as NULL to use default arguments.")}
  
  
  predicted = predictGlycans(dp1=pgp@dp1, dp2=pgp@dp2, ESI_mode=pgp@ESI_mode, 
                             scan_range1=pgp@scan_range1, scan_range2=pgp@scan_range2,
                             pent_option=pgp@pent_option, modifications=pgp@modifications, 
                             label=pgp@label) %>% 
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
                       pgp=NULL,
                       ppm=NULL, mzabs=NULL,
                       return_everything=F){
  
  data_fetched = Dulce_fetch(data, cwp=cwp, pdp=pdp, return_everything=return_everything)
  
  if (return_everything){
    data_peaks = data_fetched$peaks
    data_features = data_fetched$features
    data_fetched = data_fetched$data
  }
  
  data_converted = Dulce_toCAMERA(data_fetched, names=names, classes=classes)
  data_annotated = Dulce_find(data_converted, isotopes=isotopes, adducts=adducts, 
                                  perfwhm=perfwhm, mzabs=mzabs_isotAdd, cor_eic_th=cor_eic_th, polarity=polarity)
  data_trimmed = Dulce_trimIsotopes(data_annotated, rtmin=rtmin, rtmax=rtmax)
  data_predicted = Dulce_annotate(data_trimmed, pgp=pgp,
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

Dulce_display = function(data, peak=NULL, rtime_window=10, ppm_window=25, rt_scale=1, mz_scale=10){
  
  # Checks
  
  # Peak limits
  mz_min = peak$mzmin
  mz_max = peak$mzmax
  rt_min = peak$rtmin
  rt_max = peak$rtmax
  
  # Plot limits
  mz_window = ppm_to_mz(peak$mz, ppm_window) * mz_scale
  mz_bottom = mz_min - mz_window 
  mz_top = mz_max + mz_window   
  rt_bottom = rt_min - rtime_window * rt_scale
  rt_top = rt_max + rtime_window * rt_scale
  
  # Gathering data
  data_header = header(data)
  spectra_window = data_header %>% filter(between(retentionTime, rt_bottom, rt_top))
  
  data = data[data@featureData@data$spIdx %in% spectra_window$spectrum]
  
  neighbourhood = data.frame(mz = unlist(mz(data), use.names=F),
                             intensity = unlist(intensity(data), use.names=F),
                             rt = rep(rtime(data), lapply(mz(data), length))) %>% 
                  filter(between(mz, mz_bottom, mz_top))
  
  rt_mz = ggplot(neighbourhood, aes(x=rt, y=mz, col=intensity)) +
    geom_point(size=2, alpha=0.8) +
    scale_color_viridis(option="magma") +
    ggplot2::annotate(geom="rect",
                      ymin=mz_min,
                      ymax=mz_max,
                      xmin=rt_min, 
                      xmax=rt_max,
                      alpha=0.1, fill="red", col="red") + 
    theme_bw() + theme(legend.position = "none") + 
    ylim(c(mz_bottom, mz_top)) + xlim(c(rt_bottom, rt_top)) +  
    labs(x="Retention time (seconds)",
         y="M/z")
  
  rt_intensity = ggplot(neighbourhood %>% filter(between(mz, mz_min, mz_max)), 
                        aes(x=rt, y=intensity)) + 
    geom_line() + 
    ggplot2::annotate(geom="rect", 
                      xmin=rt_min, 
                      xmax=rt_max,
                      ymin=-Inf, 
                      ymax=Inf, 
                      alpha=0.1, fill="red", col="red")+
    xlim(c(rt_bottom, rt_top)) + 
    theme_bw() + labs(x=NULL) + theme(legend.position="none")
  
  mz_intensity = ggplot(neighbourhood, aes(x=round(mz,3), y=intensity)) + 
    geom_col() + 
    ggplot2::annotate(geom="rect", 
                      xmin=mz_min,
                      xmax=mz_max,
                      ymin=-Inf, 
                      ymax=Inf, 
                      alpha=0.1, fill="red", col="red")+
    theme_bw() + coord_flip() + labs(x=NULL) + 
    xlim(c(mz_bottom, mz_top))
  
  empty = ggplot() + theme_void() + ggplot2::annotate(geom="text", x=0, y=0, label=paste0("peak ",peak$sn))
  
  plot = arrangeGrob(rt_intensity, empty, rt_mz, mz_intensity, 
                     ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  return(plot)
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

data = data %>% Dulce_fetch(return_everything = T)
peak = data$peaks[1,]

ggsave(filename=paste0("C:/Users/ssilva/ownCloud/Lab rotations/R with Margot/Playground/R/", peak$sn, ".png"),
       Dulce_display(data$data, peak),
       width=4, height=3, unit)

mz(data$data)
register(SerialParam())

data_processed = Dulce_doAll(data, polarity="negative", ppm=10)







