library(car)
library(MASS)
library(tidyverse)
library(openxlsx)

dat <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/updatedata/dat.xlsx")
var_to_test <- colnames(dat %>%
                          dplyr::select(-annee, -region, -dplyr::contains("AMR_")))

j <- "log_AMC_ville_penicillines"
Y <- "AMR_EHPAD_FQ_R"
d <- dat[, c(Y, "annee", "region", j)]
d <- d[complete.cases(d), ]
x <- d[[j]]

# fonction boxcox
transf_boxcox <- function(x){
  if (is.null(x) || length(x) == 0) return(NULL)
  if (any(x <= 0, na.rm = TRUE)) return(rep(NA, length(x)))
  bc <- MASS::boxcox(lm(x ~ 1), plotit = FALSE)
  lambda <- bc$x[which.max(bc$y)]
  if (abs(lambda) < 1e-6){
    x_boxcox <- log(x)
  } else {
    x_boxcox <- (x^lambda - 1) / lambda
  }
  return(x_boxcox)
}

hist(transf_boxcox(x))

# transformation yeo-johnson
transf_yeo_johnson <- function(x){
yj <- car::powerTransform(x, family="yjPower")
x_yj <- car::yjPower(x, coef(yj))
return(x_yj)}

# transformation INT rank-based
transf_INT <- function(x){
  n <- length(x)
  x_int <- qnorm((rank(x, na.last="keep") - 0.5) / n)
  return(x_int)}

# RMSE
rmse <- function(mod, y){
  n <- length(y)
  py <- predict(mod)
  if (length(py) != length(y)){
    stop("Y et predict(mod) doivent être de la même longueur")
  }
  return(sqrt((1/n)*sum((y-py)^2)))
}

# RRMSE
rrmse <- function(mod, y){
  n <- length(y)
  py <- predict(mod)
  if (length(py) != length(y)){
    stop("Y et predict(mod) doivent être de la même longueur")
  }
  rmse <- sqrt((1/n)*sum((y-py)^2))
  rrmse <- rmse/sd(y)
  return(rrmse)
}

results_model2 <- function(var_to_test, Y, 
                           save_excel = FALSE,
                           criteria = "RMSE",
                           digits = 2){

Tab <- data.frame(
    variable = character(),
    model_choosen = character(),
    forme = character(),
    transformation = character(),
    betaX = numeric(),
    p.value_betaX = numeric(),
    betaX2 = numeric(),
    p.value_betaX2 = numeric(),
    phi = numeric(),
    sigma_ind = numeric(),
    sigma_res = numeric(),
    stringsAsFactors = FALSE
  )

Tab2 <- data.frame(
  variable = character(),
  model = character(),
  transformation = character(),
  p.value_DW_res = numeric(),
  p.value_JB_res = numeric(),
  p.value_BP_res = numeric(),
  p.value_JB_random_res = numeric(),
  respect_assumption = character(),
  AIC = numeric(),
  RMSE = numeric()
)

#i <- 1 
for (j in var_to_test){
  #print(j)
  d <- dat[, c(Y, "annee", "region", j)]
  d <- d[complete.cases(d), ]
  
  if (nrow(d) < 3) next
  if (sd(d[[j]]) == 0) next
  if (sd(d[[Y]]) == 0) next
  
  # transformation et hypothèses
  if (startsWith(j, "log_")){d[[j]] <- exp(d[[j]])}
  transfos <- list(
    none = function(x) x,
    log = function(x) log(x),
    sqrt = function(x) sqrt(x),
    sinhlog = function(x) log(sinh(x)),
    boxcox = function(x) transf_boxcox(x),
    yeojohnson = function(x) transf_yeo_johnson(x),
    INT = function(x) transf_INT(x)
  )
  
  # les formulas
  form1 <- as.formula(paste0(Y, " ~ X"))
  form2 <- as.formula(paste0(Y, " ~ X + I(X^2)"))
  
  found_model <- FALSE
  
  # stocker les modèles
  valid_models <- list()
  model_stats <- data.frame()
  model_param <- data.frame()
  
  ### LMM
  for (t in names(transfos)) {
      Xtmp <- tryCatch(
        transfos[[t]](d[[j]]),
        error = function(e) NULL
      )
      
      if (is.null(Xtmp)) next
      if (length(Xtmp) == 0) next
      if (all(is.na(Xtmp))) next
      if (any(!is.finite(Xtmp))) next   
      if (sd(Xtmp, na.rm = TRUE) == 0) next
      
      d$X <- scale(Xtmp, center = TRUE, scale = FALSE)
      
      if (var(d$X, na.rm = TRUE) < 1e-10) next 
      mm1 <- model.matrix(form1, data = d)
      if (kappa(t(mm1) %*% mm1, exact = TRUE) > 1e12) next  
      mm2 <- model.matrix(form2, data = d)
      if (kappa(t(mm2) %*% mm2, exact = TRUE) > 1e12) next
    
    out.lmm1 <- tryCatch(nlme::lme(
      form1,
      data = d,
      random = ~ 1|region,
      correlation = corAR1(form = ~ annee | region),
      method = "REML"),
      error = function(e) NULL
      )
    
    if(is.null(out.lmm1)) next
    
    out.lmm2 <- tryCatch(nlme::lme(
      form2,
      data = d, 
      random = ~ 1|region,
      correlation = corAR1(form = ~ annee | region),
      method = "REML"),
      error = function(e) NULL
      )
    
    if(is.null(out.lmm2)) next
    
    # test du modèle
    aov <- anova(out.lmm2, out.lmm1)
    
    #choix de modèle
    mod <- out.lmm1
    if (aov$`p-value`[2] < 0.05) {
      mod <- out.lmm2
    }
    
    ## Hypothèses
    # test de Jarque-Bera pour la normalité résiduelle
    jb <- tseries::jarque.bera.test(resid(mod))
    jb.p <- jb$p.value
    
    # test de breush-pagan pour l'homoscédasticité
    res <- resid(mod, type = "pearson")
    fit <- fitted(mod)
    bp <- lmtest::bptest(res ~ fit)
    bp.p <- bp$p.value
    
    # durbin-watson pour l'indep des résidus
    res <- resid(mod)
    dw <- lmtest::dwtest(res ~ 1)
    dw.p <- dw$p.value
    
    # Test de normalité des residus des random effects
    resid.intercept <- ranef(mod)[["(Intercept)"]]
    jb.rand <- tseries::jarque.bera.test(resid.intercept)
    jb.rand.p <- jb.rand$p.value
    
    if (bp.p > 0.05 & jb.p > 0.05 & dw.p > 0.05 & jb.rand.p > 0.05){
      #message(j, "Transformation retenue (LMM) : ", t)
      best_model.1 <- mod
      model_type <- "LMM"
      found_model <- TRUE
      
      coef_Tab <- summary(best_model.1)$tTable[, c("Value", "p-value")]
      varce <- VarCorr(best_model.1)
      df1.LMM <- data.frame(
        variable = j,
        model_choosen = "LMM",
        forme = ifelse(identical(best_model.1, out.lmm2), "quadratique", "lineaire"),
        transformation = t,
        betaX = coef_Tab[2,1],
        p.value_betaX = coef_Tab[2,2],
        betaX2 = ifelse(identical(best_model.1, out.lmm2), coef_Tab[3,1], NA),
        p.value_betaX2 = ifelse(identical(best_model.1, out.lmm2), coef_Tab[3,2], NA),
        phi = coef(best_model.1$modelStruct$corStruct, unconstrained = FALSE),
        sigma_ind = as.numeric(varce[1,2]),
        sigma_res = as.numeric(varce[2,2]),
        stringsAsFactors = FALSE
      )
      
      df2.LMM <- data.frame(
        variable = j,
        model = "LMM",
        transformation = t,
        p.value_DW_res = dw.p,
        p.value_JB_res = jb.p,
        p.value_BP_res = bp.p,
        p.value_JB_random_res = jb.rand.p,
        respect_assumption = ifelse(bp.p > 0.05 & jb.p > 0.05 & dw.p > 0.05 & jb.rand.p > 0.05, "Yes", "No"),
        AIC = stats::AIC(mod),
        RMSE = rmse(mod, d[[j]])
      )
      
      model_stats <- rbind(model_stats, df2.LMM)
      model_param <- rbind(model_param, df1.LMM)
    }
  
  }
  
  ### GLS
  for (t in names(transfos)) {
    Xtmp <- tryCatch(
      transfos[[t]](d[[j]]),
      error = function(e) NULL
    )
    
    if (is.null(Xtmp)) next
    if (length(Xtmp) == 0) next
    if (all(is.na(Xtmp))) next
    if (any(!is.finite(Xtmp))) next  
    if (sd(Xtmp, na.rm = TRUE) == 0) next
    
    d$X <- scale(Xtmp, center = TRUE, scale = FALSE)
    
    if (var(d$X, na.rm = TRUE) < 1e-10) next         
    mm1 <- model.matrix(form1, data = d)
    if (kappa(t(mm1) %*% mm1, exact = TRUE) > 1e12) next  
    mm2 <- model.matrix(form2, data = d)
    if (kappa(t(mm2) %*% mm2, exact = TRUE) > 1e12) next
      
      out.gls1 <- tryCatch(nlme::gls(
        form1,
        data = d,
        correlation = corAR1(form = ~ annee | region)),
        error = function(e) NULL)
      
      if(is.null(out.gls1)) next
      
      out.gls2 <- tryCatch(nlme::gls(
        form2,
        data = d, 
        correlation = corAR1(form = ~ annee | region)),
        error = function(e) NULL)
      
      if (is.null(out.gls2)) next
      
      # test du modèle
      aov <- anova(out.gls2, out.gls1)
      
      #choix de modèle
      mod <- out.gls1
      if (aov$`p-value`[2] < 0.05) {
        mod <- out.gls2
      }
      
      ## Hypothèses
      # test de Jarque-Bera pour la normalité résiduelle
      jb <- tseries::jarque.bera.test(resid(mod))
      jb.p <- jb$p.value
      
      # test de breush-pagan pour l'homoscédasticité
      res <- resid(mod, type = "pearson")
      fit <- fitted(mod)
      bp <- lmtest::bptest(res ~ fit)
      bp.p <- bp$p.value
      
      # durbin-watson pour l'indep des résidus
      res <- resid(mod)
      dw <- lmtest::dwtest(res ~ 1)
      dw.p <- dw$p.value
      
      if (bp.p > 0.05 & jb.p > 0.05 & dw.p > 0.05){
        best_model.2 <- mod
        model_type <- "GLS"
        found_model <- TRUE
        
        coef_Tab <- summary(best_model.2)$tTable[, c("Value", "p-value")]
        df1.GLS <- data.frame(
          variable = j,
          model_choosen = "GLS",
          forme = ifelse(identical(best_model.2, out.gls2), "quadratique", "lineaire"),
          transformation = t,
          betaX = coef_Tab[2,1],
          p.value_betaX = coef_Tab[2,2],
          betaX2 = ifelse(identical(best_model.2, out.gls2), coef_Tab[3,1], NA),
          p.value_betaX2 = ifelse(identical(best_model.2, out.gls2), coef_Tab[3,2], NA),
          phi = coef(best_model.2$modelStruct$corStruct, unconstrained = FALSE),
          sigma_ind = NA,
          sigma_res = as.numeric(summary(mod)$sigma),
          stringsAsFactors = FALSE
        )
        
        df2.GLS <- data.frame(
          variable = j,
          model = "GLS",
          transformation = t,
          p.value_DW_res = dw.p,
          p.value_JB_res = jb.p,
          p.value_BP_res = bp.p,
          p.value_JB_random_res = NA,
          respect_assumption = ifelse(bp.p > 0.05 & jb.p > 0.05 & dw.p > 0.05, "Yes", "No"),
          AIC = stats::AIC(mod),
          RMSE = rmse(mod, d[[j]])
        )
        
        model_stats <- rbind(model_stats, df2.GLS)
        model_param <- rbind(model_param, df1.GLS)
      }
    
   }

  
  ### GEE
  d$X <- scale(d[[j]], center = TRUE, scale = FALSE)
  if (var(d$X, na.rm = TRUE) < 1e-10) next         # Fix 3
  mm1 <- model.matrix(form1, data = d)
  if (kappa(t(mm1) %*% mm1, exact = TRUE) > 1e12) next  # Fix 1
  mm2 <- model.matrix(form2, data = d)
  if (kappa(t(mm2) %*% mm2, exact = TRUE) > 1e12) next
      
  out.gee1 <- tryCatch(geepack::geeglm(
        form1,
        id = region, 
        #waves = annee, 
        data = d, 
        corstr = "ar1", 
        std.err = "fij"),
        error = function(e) NULL
        )
  if(is.null(out.gee1)) next
      
    out.gee2 <- tryCatch(geepack::geeglm(
        form2,
        id = region, 
        #waves = annee, 
        data = d, 
        corstr = "ar1", 
        std.err = "fij"),
        error = function(e) NULL)
  
  if (is.null(out.gee2)) next
      
      # test du modèle
      aov <- anova(out.gee1, out.gee2, test = "Wald")
      
      #choix de modèle
      mod <- out.gee1
      if (aov$`P(>|Chi|)` < 0.05) {
        mod <- out.gee2
      }
      
      ## Hypothèses
      # durbin-watson pour l'indep des résidus
      res <- resid(mod)
      dw <- lmtest::dwtest(res ~ 1)
      dw.p <- dw$p.value
      
      if (dw.p < 0.05){
        #message(j, "Transformation retenue (LMM) : ", t)
        best_model.3 <- mod
        model_type <- "GEE"
        
        coef_Tab <- as.data.frame(summary(mod)$coefficients[, c("Estimate", "Pr(>|W|)")])
        #print(coef_Tab$`Pr(>|W|)`)
        df1.GEE <- data.frame(
          variable = j,
          model_choosen = "GEE",
          forme = ifelse(identical(best_model.3, out.gee2), "quadratique", "lineaire"),
          transformation = "none",
          betaX = coef_Tab[2,1],
          p.value_betaX = coef_Tab[2,2],
          betaX2 = ifelse(identical(best_model.3, out.gee2), coef_Tab[3,1], NA),
          p.value_betaX2 = ifelse(identical(best_model.3, out.gee2), coef_Tab[3,2], NA),
          phi = mod$geese$alpha,
          sigma_ind = NA,
          sigma_res = NA,
          stringsAsFactors = FALSE
        )
        
        df2.GEE <- data.frame(
          variable = j,
          model = "GEE",
          transformation = "none",
          p.value_DW_res = dw.p,
          p.value_JB_res = NA,
          p.value_BP_res = NA,
          p.value_JB_random_res = NA,
          respect_assumption = ifelse(dw.p < 0.05, "Yes", "No"),
          AIC = geepack::QIC(best_model.3)[1],
          RMSE = rmse(best_model.3, d[[j]])
        )
        
        model_stats <- rbind(model_stats, df2.GEE)
        model_param <- rbind(model_param, df1.GEE)
      }
      
      if (criteria == "RMSE"){
        best_index <- which.min(model_stats$RMSE)
        df2 <- model_stats[best_index,]
        df1 <- model_param[best_index,]
        Tab2 <- rbind.data.frame(Tab2, df2)
        Tab <- rbind.data.frame(Tab, df1)
      }
      else {
        best_index <- which.min(model_stats$AIC)
        df2 <- model_stats[best_index,]
        df1 <- model_param[best_index,]
        Tab2 <- rbind.data.frame(Tab2, df2)
        Tab <- rbind.data.frame(Tab, df1)
      }
    }
  
  return(list(Tab, Tab2))
}

     
# test
res2 <- results_model2(var_to_test = var_to_test, Y = "AMR_ville_FQ_R")
a <- res2[[1]]
View(a)
resb <- results_model2(var_to_test = var_to_test, Y = "AMR_EHPAD_FQ_R")
b <- resb[[1]]
View(b)

setdiff(var_to_test, res2[[1]]$variable)

quality <- res2[[2]]
View(quality)
boxplot(quality$RMSE)

quality %>% 
  ggplot(aes(x = 1, y = RMSE)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
  geom_point(aes(group = variable),
             width = 0.1,
             show.legend = FALSE, 
             position = position_dodge(0.9)) +
  geom_text(data = quality %>% dplyr::filter(RMSE > 20), 
            aes(group = variable, label = variable),
            hjust = -0.1,
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  theme_bw()
  
summary(quality$RMSE)  
  
  
