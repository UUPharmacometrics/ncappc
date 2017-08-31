pop_hist <- function(obsdata,
                     smeanData, 
                     svarData,
                     title=NULL,
                     quant=c(0.025,0.975)) 
{
  ## create data
  sim_mean_dat <- smeanData
  sim_mean_dat$type <- "Mean"
  
  sim_var_dat <- svarData
  sim_var_dat$type <- "Variance"
  
  sim_dat <- rbind(sim_mean_dat,sim_var_dat)
  
  long_dat <- sim_dat %>% tidyr::gather(-type,key = "var",value = "value")
  long_dat <- long_dat %>% dplyr::mutate(var2=paste0(var," :: ",type))
  
  vars <- unique(long_dat$var2)
  
  dfs <- list()
  for (v in vars){
    dfs[[length(dfs)+1]] <- long_dat %>% dplyr::filter(var2==v)
  }
  
  #bw <- function(x) diff(unname(quantile(as.numeric(x),c(0.005,0.985))))/(2*IQR(as.numeric(x)))/length(as.numeric(x))^(1/3)
  
  bw <- function(x) 2 * IQR(x) / (length(x)^(1/3))
  
  stuff<- ggplot(data=long_dat, mapping=aes(x=value) )
  stuff <- stuff + facet_wrap( ~ var2 , scale="free")
  
  for(i in 1:length(dfs)){
    stuff<- stuff + geom_histogram(data=dfs[[i]], binwidth=bw(dfs[[i]]$value),colour="black", fill="white")
  }
  
  long_sum <- obsdata %>% dplyr::summarise_all(mean,na.rm=TRUE) %>% 
    tidyr::gather(key = "var",value = "value") %>% 
    dplyr::mutate(type="Mean") %>% 
    dplyr::mutate(var2=paste0(var," :: ",type))
  
  long_sum_2 <- obsdata %>% dplyr::summarise_all(var,na.rm=TRUE) %>% 
    tidyr::gather(key = "var",value = "value") %>% 
    dplyr::mutate(type="Variance") %>% 
    dplyr::mutate(var2=paste0(var," :: ",type))
  
  long_sum <- rbind(long_sum,long_sum_2)
  
  stuff <- stuff + 
    geom_vline(data=long_sum, aes(xintercept=value,colour="Observed"),
               linetype="solid", size=1 ) 
  
  pct_dat_sum <- long_dat %>% dplyr::group_by(var2) %>% 
    dplyr::summarise(low=quantile(value,probs=min(quant)),
              high=quantile(value,probs=max(quant))) %>% 
    tidyr::gather(low,high,key = "loc",value = "value")
  
  med_dat_sum <- long_dat %>% dplyr::group_by(var2) %>% 
    dplyr::summarise(median=median(value,na.rm=TRUE)) 
  
  stuff <- stuff + geom_vline(data=pct_dat_sum,aes(xintercept=value,
                                                   color="Simulated (2.5%, 50%, 97.5%)"),   
                              linetype="dotted", size=1) +
    geom_vline(data=med_dat_sum,aes(xintercept=median,
                                    color="Simulated (2.5%, 50%, 97.5%)"),   
               linetype="solid", size=1)
  
  
  
  stuff <- stuff + theme(legend.position="bottom") +
    scale_colour_manual(name='', values=c('Observed'='red', 
                                          'Simulated (2.5%, 50%, 97.5%)'='grey40')) +
    theme(legend.position="bottom") 
  #guides(colour = guide_legend(override.aes = list(linetype=c(1,0)
  #                                                 , shape=c(NA, 16))))
  
  if(!is.null(title)){
    Label <-title
    stuff <- stuff + ggtitle(Label)
  }
  
  return(stuff)
  
}
