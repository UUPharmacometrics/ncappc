context("Check that NCA calculations are working")

test_that("Extravascular single dose data is computed corectly",{
  
  exvas_sd_dat <- readr::read_table("data_exvas_sd.txt",skip = 1) # Can move this to data
  
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
  
  exvas_md_dat <- readr::read_table("data_exvas_md.txt",skip = 1) # Can move this to data
  
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

test_that("IV single dose data is computed corectly",{
  
  iv_sd_dat <- readr::read_table("data_iv_sd.txt",skip = 1) # Can move this to data
  
  iv_sd_true_result <- readr::read_csv("results_iv_sd_true.csv")
  
  iv_sd_true_result <- iv_sd_true_result %>%
    tidyr::gather(var, val, 2:ncol(.)) %>%
    tidyr::spread_(names(.)[1], "val")
  
  #### test for linear method 
  out <- ncappc(obsFile=iv_sd_dat,
                psnOut=TRUE,
                method="linear",
                adminType = "iv-bolus",
                #doseType = "ss",
                #timeNmObs = "TAD",
                #doseTime = 96,
                #Tau = 12,
                evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_sd_true_result %>% 
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
  out <- ncappc(obsFile=iv_sd_dat,
                psnOut=TRUE,
                method="linearup-logdown", 
                adminType = "iv-bolus",
                #doseType = "ss",
                #timeNmObs = "TAD",
                #doseTime = 96,
                #Tau = 12,
                evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
 
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_sd_true_result %>% 
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

test_that("IV multiple dose data is computed corectly",{
  
  iv_md_dat <- readr::read_table("data_iv_md.txt",skip = 1) # Can move this to data
  iv_md_dat[9,"TIME"] <- 96
  
  iv_md_true_result <- readr::read_csv("results_iv_md_true.csv")
  
  iv_md_true_result <- iv_md_true_result %>%
    tidyr::gather(var, val, 2:ncol(.)) %>%
    tidyr::spread_(names(.)[1], "val")
  
  #### test for linear method 
  out <- ncappc(obsFile=iv_md_dat,
                psnOut=TRUE,
                method="linear",
                adminType = "iv-bolus",
                doseType = "ss",
                #timeNmObs = "TAD",
                doseTime = 96,
                backExtrp = T,
                Tau = 12,
                #evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_md_true_result %>% 
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
  out <- ncappc(obsFile=iv_md_dat,
                psnOut=TRUE,
                method="linearup-logdown", 
                adminType = "iv-bolus",
                doseType = "ss",
                #timeNmObs = "TAD",
                doseTime = 96,
                backExtrp = T,
                Tau = 12,
                #evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_md_true_result %>% 
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


test_that("IV infusion single dose data is computed corectly",{
  
  iv_inf_sd_dat <- readr::read_table("data_iv_inf_sd.txt",skip = 1) # Can move this to data
  
  iv_inf_sd_true_result <- readr::read_csv("results_iv_inf_sd_true.csv")
  
  iv_inf_sd_true_result <- iv_inf_sd_true_result %>%
    tidyr::gather(var, val, 2:ncol(.)) %>%
    tidyr::spread_(names(.)[1], "val")
  
  #### test for linear method 
  out <- ncappc(obsFile=iv_inf_sd_dat,
                psnOut=TRUE,
                method="linear",
                adminType = "iv-infusion",
                #doseType = "ss",
                #timeNmObs = "TAD",
                #doseTime = 96,
                #Tau = 12,
                evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_inf_sd_true_result %>% 
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
  out <- ncappc(obsFile=iv_inf_sd_dat,
                psnOut=TRUE,
                method="linearup-logdown", 
                adminType = "iv-infusion",
                #doseType = "ss",
                #timeNmObs = "TAD",
                #doseTime = 96,
                #Tau = 12,
                #TI=2,
                evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_inf_sd_true_result %>% 
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

test_that("IV infusion multiple dose data is computed corectly",{
  
  iv_inf_md_dat <- readr::read_table("data_iv_inf_md.txt",skip = 1) # Can move this to data
  iv_inf_md_dat[9,"TIME"] <- 96
  
  iv_inf_md_true_result <- readr::read_csv("results_iv_inf_md_true.csv")
  
  iv_inf_md_true_result <- iv_inf_md_true_result %>%
    tidyr::gather(var, val, 2:ncol(.)) %>%
    tidyr::spread_(names(.)[1], "val")
  
  #### test for linear method 
  out <- ncappc(obsFile=iv_inf_md_dat,
                psnOut=TRUE,
                method="linear",
                adminType = "iv-infusion",
                doseType = "ss",
                #timeNmObs = "TAD",
                doseTime = 96,
                backExtrp = T,
                Tau = 12,
                TI=2,
                #evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_inf_md_true_result %>% 
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
  out <- ncappc(obsFile=iv_inf_md_dat,
                psnOut=TRUE,
                method="linearup-logdown", 
                adminType = "iv-infusion",
                doseType = "ss",
                #timeNmObs = "TAD",
                doseTime = 96,
                backExtrp = T,
                Tau = 12,
                TI=2,
                #evid = FALSE, 
                onlyNCA = T,
                extrapolate = T,
                printOut = F,
                noPlot = T)
  
  
  calc_vals <- out$ncaOutput %>% tibble::as_tibble() %>% 
    dplyr::mutate_all(dplyr::funs(as.numeric(as.character(.))))
  true_vals <- iv_inf_md_true_result %>% 
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

