hist_mean_var <- function(obs_data, sim_data,
                          param, 
                          strat_vars=NULL) {
  
  NSIM <- ID <- figlbl <- NULL
  
  # get stratification variable
  strats <- NULL
  if(!is.null(strat_vars)) strats <- rlang::syms(strat_vars)
  
  # summarize the simulated data along stratifications if present
  sum_dat <- sim_data
  if(!is.null(strat_vars)){
    for(i in 1:length(strats)){
      sum_dat <- sum_dat %>% dplyr::rename_all(~sub(paste0('STRAT',i),strats[[i]],.x))
    }
  }
  sum_dat <- sum_dat %>%  dplyr::select(param,NSIM,!!!strats) %>% dplyr::group_by(NSIM) 
  if(!is.null(strats)) sum_dat <- sum_dat %>% dplyr::group_by(!!!strats,add=TRUE)
  sum_dat <- sum_dat %>% dplyr::summarise_all(dplyr::funs("mean","median","var"),na.rm=T) %>% 
    dplyr::ungroup()
  
  # get data to look as it should
  #obs_data <- outData 
  if(!is.null(strat_vars)){
    for(i in 1:length(strats)){
      obs_data <- obs_data %>% dplyr::rename_all(~sub(paste0('STRAT',i),strats[[i]],.x))
    }
  }
  obs_data <- obs_data %>%  dplyr::filter(ID!="") %>% dplyr::select(param,!!!strats) 
  if(!is.null(strats)) obs_data <- obs_data %>% dplyr::group_by(!!!strats)
  obs_data <- obs_data %>% tidyr::nest()
  
  mean_data <- sum_dat %>%  dplyr::select(dplyr::matches(".*\\_mean$"),!!!strats) %>% 
    setNames(., sub("_mean$", "", names(.)))
  if(!is.null(strats)) mean_data <- mean_data %>% dplyr::group_by(!!!strats)
  mean_data <- mean_data %>% tidyr::nest()
  
  var_data <- sum_dat %>%  dplyr::select(dplyr::matches(".*\\var$"),!!!strats) %>% 
    setNames(., sub("_var$", "", names(.)))
  if(!is.null(strats)) var_data <- var_data %>% dplyr::group_by(!!!strats)
  var_data <- var_data %>% tidyr::nest()
  
  if(!is.null(strats)){
    join_vars <- strats %>% paste()
    plot_data <- mean_data %>%  dplyr::rename(mean_data=data) %>% 
      dplyr::full_join(var_data,by=join_vars) %>% dplyr::rename(var_data=data) %>% 
      dplyr::full_join(obs_data,by=join_vars) %>% dplyr::rename(obs_data=data)
    
    str_lab <- c()
    for(i in 1:length(strats)) {
      plot_data <- plot_data %>% dplyr::mutate(!!paste0("lab",i):=paste(paste0(strats[[i]]),paste0(!!strats[[i]]), sep = '='))
      str_lab <- c(str_lab,paste0("lab",i))
    }
    str_lab_var <- rlang::syms(str_lab)
    plot_data <- plot_data %>% dplyr::mutate("figlbl"=paste(!!!str_lab_var,sep=","))
  } else {
    plot_data <- mean_data %>% dplyr::rename(mean_data=data)
    plot_data$var_data <- var_data$data
    plot_data$obs_data <- obs_data$data
    plot_data$figlbl      <- list(NULL)
  }
  
  plot_data <- plot_data %>% dplyr::mutate(plot=purrr::pmap(list(obs_data,mean_data,var_data,figlbl),hist_mean_var_plot))
  
  pop_hist_list <- plot_data$plot
  return(pop_hist_list)
  #   
  #   
  #   pop_hist_list <- list()
  #   
  #   # Population histogram (case=1)
  #   if (case == 1){
  #     smeanData   <- data.frame()
  #     smedianData <- data.frame()
  #     svarData <- data.frame()
  #     #for (i in 1:length(lasdf)){
  #     for (i in 1:nsim){
  #       #tmdf        <- subset(data.frame(lasdf[[i]]), select=param)
  #       tmdf <- dasdf %>% dplyr::filter(NSIM==i) %>% dplyr::select(param)
  #       meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
  #       smeanData   <- rbind(smeanData, meanPrm)
  #       medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
  #       smedianData <- rbind(smedianData, medianPrm)
  #       varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
  #       svarData <- rbind(svarData, varPrm)
  #     }
  #     obsdata     <- subset(outData, select=param, ID!="")
  #     figlbl      <- NULL
  #     histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
  #     gdr         <- histpopgrob$gdr
  #     mylegend    <- histpopgrob$legend
  #     lheight     <- histpopgrob$lheight
  #     pop_hist_tmp <- pop_hist(obsdata,
  #                              smeanData, 
  #                              svarData,
  #                              title=figlbl)
  #     if (printOut){
  #       fl <- paste0(usrdir,"/PopMean")
  #       if (figFormat=="tiff"){
  #         eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
  #       }else{
  #         eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
  #       }
  #       suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #       dev.off()
  #     }
  #     suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #     ggr <- grid.grab()
  #     popplot[[length(popplot)+1]] <- ggr
  #     pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
  #   }
  #   
  #   # Population histogram (case=2)
  #   if (case == 2){
  #     for (s1 in 1:npopStr1){
  #       if (nrow(dasdf[dasdf$STRAT1==popStr1[s1],]) == 0) next
  #       smeanData   <- data.frame()
  #       smedianData <- data.frame()
  #       svarData <- data.frame()
  #       for (i in 1:nsim){
  #       #for (i in 1:length(lasdf)){
  #         #tmdf        <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==popStr1[s1])
  #         tmdf <- dasdf %>% dplyr::filter(NSIM==i,STRAT1==popStr1[s1]) %>% dplyr::select(param)
  #         meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
  #         smeanData   <- rbind(smeanData, meanPrm)
  #         medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
  #         smedianData <- rbind(smedianData, medianPrm)
  #         varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
  #         svarData <- rbind(svarData, varPrm)
  #       }
  #       obsdata     <- subset(outData, select=param, ID!="" & STRAT1==popStr1[s1])
  #       figlbl      <- paste0(popStrNm1,"-",popStr1[s1])
  #       histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
  #       gdr         <- histpopgrob$gdr
  #       mylegend    <- histpopgrob$legend
  #       lheight     <- histpopgrob$lheight
  #       pop_hist_tmp <- pop_hist(obsdata,
  #                                smeanData, 
  #                                svarData,
  #                                title=figlbl)
  #       if (printOut){
  #         fl <- paste0(usrdir,"/PopMean_",figlbl)
  #         if (figFormat=="tiff"){
  #           eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
  #         }else{
  #           eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
  #         }
  #         suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #         dev.off()
  #       }
  #       suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #       ggr <- grid.grab()
  #       popplot[[length(popplot)+1]] <- ggr
  #       pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
  #     }
  #   }
  #   
  #   # Population histogram (case=3)
  #   if (case == 3){
  #     for (s1 in 1:npopStr1){
  #       for (s2 in 1:npopStr2){
  #         if (nrow(dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2],]) == 0) next
  #         smeanData   <- data.frame()
  #         smedianData <- data.frame()
  #         svarData <- data.frame()
  #         for (i in 1:nsim){
  #         #for (i in 1:length(lasdf)){
  #           #tmdf        <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
  #           tmdf <- dasdf %>% dplyr::filter(NSIM==i,STRAT1==popStr1[s1],STRAT2==popStr2[s2]) %>% dplyr::select(param)
  #           meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
  #           smeanData   <- rbind(smeanData, meanPrm)
  #           medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
  #           smedianData <- rbind(smedianData, medianPrm)
  #           varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
  #           svarData <- rbind(svarData, varPrm)
  #         }
  #         obsdata     <- subset(outData, select=param, ID!="" & STRAT1==popStr1[s1] & STRAT2==popStr2[s2])
  #         figlbl      <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2])
  #         histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
  #         gdr         <- histpopgrob$gdr
  #         mylegend    <- histpopgrob$legend
  #         lheight     <- histpopgrob$lheight
  #         pop_hist_tmp <- pop_hist(obsdata,
  #                                  smeanData, 
  #                                  svarData,
  #                                  title=figlbl)
  #         if (printOut){
  #           fl <- paste0(usrdir,"/PopMean_",figlbl)
  #           if (figFormat=="tiff"){
  #             eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
  #           }else{
  #             eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
  #           }
  #           suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #           dev.off()
  #         }
  #         suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #         ggr <- grid.grab()
  #         popplot[[length(popplot)+1]] <- ggr
  #         pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
  #       }
  #     }
  #   }
  #   
  #   # Population histogram (case=4)
  #   if (case == 4){
  #     for (s1 in 1:npopStr1){
  #       for (s2 in 1:npopStr2){
  #         for (s3 in 1:npopStr3){
  #           if (nrow(dasdf[dasdf$STRAT1==popStr1[s1] & dasdf$STRAT2==popStr2[s2] & dasdf$STRAT3==popStr3[s3],]) == 0) next
  #           smeanData   <- data.frame()
  #           smedianData <- data.frame()
  #           svarData <- data.frame()
  #           for (i in 1:nsim){
  #           # for (i in 1:length(lasdf)){
  #             #tmdf        <- subset(data.frame(lasdf[[i]]), select=param, STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
  #             tmdf <- dasdf %>% dplyr::filter(NSIM==i,STRAT1==popStr1[s1],STRAT2==popStr2[s2],STRAT3==popStr3[s3]) %>% dplyr::select(param)
  #             meanPrm     <- as.data.frame(lapply(tmdf, FUN=function(x) mean(as.numeric(x[!is.na(x)]))))
  #             smeanData   <- rbind(smeanData, meanPrm)
  #             medianPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) median(as.numeric(x[!is.na(x)]))))
  #             smedianData <- rbind(smedianData, medianPrm)
  #             varPrm   <- as.data.frame(lapply(tmdf, FUN=function(x) var(as.numeric(x[!is.na(x)]))))
  #             svarData <- rbind(svarData, varPrm)
  #           }
  #           obsdata     <- subset(outData, select=param, ID!="" & STRAT1==popStr1[s1] & STRAT2==popStr2[s2] & STRAT3==popStr3[s3])
  #           figlbl      <- paste0(popStrNm1,"-",popStr1[s1],"_",popStrNm2,"-",popStr2[s2],"_",popStrNm3,"-",popStr3[s3])
  #           histpopgrob <- histpop.plot(obsdata=obsdata,simdata=smedianData,figlbl=figlbl,param=param,cunit=cunit,tunit=tunit,spread=spread)
  #           gdr         <- histpopgrob$gdr
  #           mylegend    <- histpopgrob$legend
  #           lheight     <- histpopgrob$lheight
  #           pop_hist_tmp <- pop_hist(obsdata,
  #                                    smeanData, 
  #                                    svarData,
  #                                    title=figlbl)
  #           if (printOut){
  #             fl <- paste0(usrdir,"/PopMean_",figlbl)
  #             if (figFormat=="tiff"){
  #               eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",compression=\"lzw\",res=120)")))
  #             }else{
  #               eval(parse(text=paste0(figFormat,"(file=\"",fl,".",figFormat,"\",height=",hth,",width=",wth,",units=\"cm\",res=120)")))
  #             }
  #             suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #             dev.off()
  #           }
  #           suppressMessages(suppressWarnings(grid.arrange(gdr,mylegend, heights=unit.c(unit(1, "npc") - lheight, lheight))))
  #           ggr <- grid.grab()
  #           popplot[[length(popplot)+1]] <- ggr
  #           pop_hist_list[[length(pop_hist_list)+1]] <- pop_hist_tmp
  #         }
  #       }
  #     }
  #   }
}