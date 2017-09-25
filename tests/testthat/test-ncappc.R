context("Check that NCA calculations are working")

test_that("Check if ncappc is available",{
  expect_true(exists("ncappc"), TRUE)
  ncappc.db <- data.frame(TIME=c(0:10),CONC=c(0,1,2,3.5,3.6,3.2,2.4,1.6,1.3,1.1,0.8))
  
  TIME <- ncappc.db$TIME
  CONC <- ncappc.db$CONC
  TIME.positive <- ncappc.db[ncappc.db$TIME >= 0,"TIME"]
  CONC.positive <- ncappc.db[ncappc.db$CONC >= 0,"CONC"]
  
  expect_equal(sum(TIME),sum(TIME.positive))
  expect_equal(sum(CONC),sum(CONC.positive))

})


test_that("Extravascular single dose data is computed corectly",{
  exvas_sd_dat <- readr::read_table("data_exvascsd.txt",skip = 1) # Can move this to data
  
  exvas_sd_true_result <- readr::read_csv("results_exvas_sd_true.csv")
  
  # transpose the true results
  exvas_sd_true_result <- exvas_sd_true_result %>%
    tidyr::gather(var, val, 2:ncol(.)) %>%
    tidyr::spread_(names(.)[1], "val")
  
  #### test for linear method 
  out <- ncappc(obsFile=exvas_sd_dat,
                psnOut=TRUE,
                method="linear", 
                evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- exvas_sd_true_result %>% 
    dplyr::filter(var=="Linear_trapezoidal_WinNonlin") %>% 
    dplyr::rename("ID"=var) %>% dplyr::mutate(ID=1)
  
  # are all the names for the true result in the computed result?
  expect_length(names(true_vals)[!(names(true_vals) %in% names(calc_vals))],0)
  
  # select the names we can compare and order the same for both
  calc_vals_sub <- calc_vals %>% dplyr::select(names(calc_vals)[(names(calc_vals) %in% names(true_vals))]) %>% 
    dplyr::select(names(true_vals))
  
  # values the same?
  expect_equal(data.frame(calc_vals_sub),data.frame(true_vals),tolerance=0.005)
  
  #### test for log-linear method 
  out <- ncappc(obsFile=exvas_sd_dat,
                psnOut=TRUE,
                method="linearup-logdown", 
                evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- exvas_sd_true_result %>% 
    dplyr::filter(var=="Linear_log_trapezoidal_WinNonlin") %>% 
    dplyr::rename("ID"=var) %>% dplyr::mutate(ID=1)
  
  # are all the names for the true result in the computed result?
  expect_length(names(true_vals)[!(names(true_vals) %in% names(calc_vals))],0)
  
  # select the names we can compare and order the same for both
  calc_vals_sub <- calc_vals %>% dplyr::select(names(calc_vals)[(names(calc_vals) %in% names(true_vals))]) %>% 
    dplyr::select(names(true_vals))
  
  # values the same?
  expect_equal(data.frame(calc_vals_sub),data.frame(true_vals),tolerance=0.005)

})

test_that("Extravascular multiple dose data is computed corectly",{
  exvas_md_dat <- readr::read_table("data_exvascmd.txt",skip = 1) # Can move this to data
  
  exvas_md_dat[9,"TIME"] <- 96
  exvas_md_true_result <- readr::read_csv("results_exvas_md_true.csv")
  
  exvas_md_true_result <- exvas_md_true_result %>%
    tidyr::gather(var, val, 2:ncol(.)) %>%
    tidyr::spread_(names(.)[1], "val")
  
  #### test for linear method 
  out <- ncappc(obsFile=exvas_md_dat,
                psnOut=TRUE,
                method="linear", 
                doseType = "ss",
                #timeNmObs = "TAD",
                doseTime = 96,
                Tau = 12,
                evid = TRUE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- exvas_md_true_result %>% 
    dplyr::filter(var=="Linear_trapezoidal_WinNonlin") %>% 
    dplyr::rename("ID"=var) %>% dplyr::mutate(ID=1)
  
  # are all the names for the true result in the computed result?
  expect_length(names(true_vals)[!(names(true_vals) %in% names(calc_vals))],0)
  
  # select the names we can compare and order the same for both
  calc_vals_sub <- calc_vals %>% dplyr::select(names(calc_vals)[(names(calc_vals) %in% names(true_vals))]) %>% 
    dplyr::select(names(true_vals))
  
  # values the same?
  expect_equal(data.frame(calc_vals_sub),data.frame(true_vals),tolerance=0.005)
  
  #### test for log-linear method 
  out <- ncappc(obsFile=exvas_md_dat,
                psnOut=TRUE,
                method="linearup-logdown", 
                doseType = "ss",
                #timeNmObs = "TAD",
                doseTime = 96,
                Tau = 12,
                evid = TRUE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- exvas_md_true_result %>% 
    dplyr::filter(var=="Linear_log_trapezoidal_WinNonlin") %>% 
    dplyr::rename("ID"=var) %>% dplyr::mutate(ID=1)
  
  # are all the names for the true result in the computed result?
  expect_length(names(true_vals)[!(names(true_vals) %in% names(calc_vals))],0)
  
  # select the names we can compare and order the same for both
  calc_vals_sub <- calc_vals %>% dplyr::select(names(calc_vals)[(names(calc_vals) %in% names(true_vals))]) %>% 
    dplyr::select(names(true_vals))
  
  # values the same?
  expect_equal(data.frame(calc_vals_sub),data.frame(true_vals),tolerance=0.005)
  
})
