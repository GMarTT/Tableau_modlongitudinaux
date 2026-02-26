Y <- "AMR_EHPAD_FQ_R"

results_model <- function(var_to_test, Y, model, save_excel = FALSE, print_common_variables = FALSE){
  
  #--------------------
  # GEE
  #--------------------
  
  Tab <- data.frame(
    variable = character(),
    forme = character(),
    betaX = numeric(),
    p.value_betaX = numeric(),
    betaX2 = numeric(),
    p.value_betaX2 = numeric(),
    phi = numeric(),
    p.value_DW_res = numeric(),
    stringsAsFactors = FALSE
  )
  i <- 1
  for (j in var_to_test) {
    Tab[i, "variable"] <- j
    d <- dat2[, c(Y, "annee", "region", j)]
    d <- d[complete.cases(d), ]
    
    if (nrow(d) < 3) next
    if (sd(d[[j]]) == 0) next
    if (sd(d[[Y]]) == 0) next
    
    # GEE
    # centrage de la variable
    d$Xc <- scale(d[[j]], center = TRUE, scale = FALSE)
    
    # les formulas
    form1 <- as.formula(paste(Y, "~ Xc"))
    form2 <- as.formula(paste(Y, "~ Xc + I(Xc^2)"))
    
    out.gee1 <- geepack::geeglm(
      form1,
      id = region, 
      #waves = annee, 
      data = d, 
      corstr = "ar1", 
      std.err = "fij")
    
    out.gee2 <- geepack::geeglm(
      form2,
      id = region, 
      #waves = annee, 
      data = d, 
      corstr = "ar1", 
      std.err = "fij")
    
    # test du modèle
    mod <- out.gee1
    Tab[i, "forme"] <- "lineaire"
    if (!j %in% c("ETP_prod_agri", "nb_chef_exploitation_agri")) {
      aov <- anova(out.gee1, out.gee2, test = "Wald")
      if (aov$`P(>|Chi|)` < 0.05) {
        mod <- out.gee2
        Tab[i, "forme"] <- "quadratique"
      }}
    
    # extraire les betas et leur p-value 
    coef_tab <- summary(mod)$coefficients[, c("Estimate", "Pr(>|W|)")]  
    Tab[i, "betaX"] <- round(coef_tab[2,1],2)
    Tab[i, "p.value_betaX"] <- round(coef_tab[2,2],2)
    if (identical(mod, out.gee2)){
      Tab[i, "betaX2"] <- round(coef_tab[3,1],2)
      Tab[i, "p.value_betaX2"] <- round(coef_tab[3,2],2)
    }
    Tab[i, "phi"] <- round(mod$geese$alpha, 2)
    # test de durbin-watson pour autocorrélation des résidus
    dw <- lmtest::dwtest(mod)
    dwp <- ifelse(dw$p.value < 0.001, "<0.001", round(dw$p.value,2))
    Tab[i, "p.value_DW_res"] <- dwp
    i <- i+1  
  }
  
  Tab_show <- Tab %>% 
    dplyr::filter(p.value_betaX < 0.05)
  
  #----------------
  # LMM
  #----------------
  
  Tab2 <- data.frame(
    variable = character(),
    random_effect = character(),
    forme = character(),
    betaX = numeric(),
    p.value_betaX = numeric(),
    betaX2 = numeric(),
    p.value_betaX2 = numeric(),
    sigma1 = numeric(),
    sigma2 = numeric(),
    phi = numeric(),
    stringsAsFactors = FALSE
  )
  i <- 1
  for (j in var_to_test) {
    Tab2[i, "variable"] <- j
    d <- dat2[, c(Y, "annee", "region", j)]
    d <- d[complete.cases(d), ]
    
    if (nrow(d) < 3) next
    if (sd(d[[j]]) == 0) next
    if (sd(d[[Y]]) == 0) next
    
    # centrage de la variable
    d$Xc <- scale(d[[j]], center = TRUE, scale = FALSE)
    
    # test de l'effet aléatoire
    form.lmm <- as.formula(paste(Y, "~ Xc + (1|region)"))
    
    out.lm <- lm(
      form1,
      data = d)
    
    out.lmm <- lme4::lmer(form.lmm,
                          data = d, REML = FALSE)
    
    aov <- anova(out.lmm, out.lm)
    Tab2[i, "random_effect"] <- "Yes"
    if (aov$`Pr(>Chisq)`[2] >= 0.05) {
      mod <- out.lm
      Tab2[i, "random_effect"] <- "No"
    }
    
    # Choix entre linéaire et polynomial 
    
    out.lmm1 <- nlme::lme(
      form1,
      data = d,
      random = ~ 1|region,
      correlation = corAR1(form = ~ annee | region),
      method = "REML")
    
    out.lmm2 <- nlme::lme(
      form2,
      data = d, 
      random = ~ 1|region,
      correlation = corAR1(form = ~ annee | region),
      method = "REML")
    
    # test du modèle
    mod <- out.lmm1
    Tab2[i, "forme"] <- "lineaire"
    #if (!j %in% c("ETP_prod_agri",   "nb_chef_exploitation_agri")) {
    aov <- anova(out.lmm2, out.lmm1)
    if (aov$`p-value`[2] < 0.05) {
      mod <- out.lmm2
      Tab2[i, "forme"] <- "quadratique"
    }
    #}
    
    # extraire les betas et leur p-value 
    if (inherits(mod, "lme")) {
      coef_Tab2 <- summary(mod)$tTable[, c("Value", "p-value")]
    } else if (inherits(mod, "lm")) {
      coef_Tab2 <- summary(mod)$coefficients[, c("Estimate", "Pr(>|t|)")]
    } 
    Tab2[i, "betaX"] <- round(coef_Tab2[2,1],2)
    Tab2[i, "p.value_betaX"] <- round(coef_Tab2[2,2],2)
    if (identical(mod, out.lmm2)){
      Tab2[i, "betaX2"] <- round(coef_Tab2[3,1],2)
      Tab2[i, "p.value_betaX2"] <- round(coef_Tab2[3,2],2)
    }
    Tab2[i, "phi"] <- round(coef(mod$modelStruct$corStruct, unconstrained = FALSE), 2)
    varce <- VarCorr(mod)
    Tab2[i, "sigma1"] <- round(as.numeric(varce[1,2]), 2)
    Tab2[i, "sigma2"] <- round(as.numeric(varce[2,2]), 2)
    i <- i+1  
  }
  
  Tab2_show <- Tab2 %>% 
    dplyr::filter(p.value_betaX < 0.05)
  
  #----------------------
  # GLS
  #----------------------
  
  Tab3 <- data.frame(
    variable = character(),
    random_effect = character(),
    forme = character(),
    betaX = numeric(),
    p.value_betaX = numeric(),
    betaX2 = numeric(),
    p.value_betaX2 = numeric(),
    sigma = numeric(),
    phi = numeric(),
    stringsAsFactors = FALSE
  )
  i <- 1
  for (j in var_to_test) {
    Tab3[i, "variable"] <- j
    d <- dat2[, c(Y, "annee", "region", j)]
    d <- d[complete.cases(d), ]
    
    if (nrow(d) < 3) next
    if (sd(d[[j]]) == 0) next
    if (sd(d[[Y]]) == 0) next
    
    # GEE
    # centrage de la variable
    d$Xc <- scale(d[[j]], center = TRUE, scale = FALSE)
    
    # test de l'effet aléatoire
    out.lm <- lm(
      form1,
      data = d)
    
    out.lmm <- lme4::lmer(form.lmm,
                          data = d, REML = FALSE)
    
    aov <- anova(out.lmm, out.lm)
    Tab3[i, "random_effect"] <- "Yes"
    if (aov$`Pr(>Chisq)`[2] >= 0.05) {
      mod <- out.lm
      Tab3[i, "random_effect"] <- "No"
    }
    
    # Choix entre linéaire et polynomial 
    
    out.gls1 <- nlme::gls(
      form1,
      data = d,
      correlation = corAR1(form = ~ annee | region))
    
    out.gls2 <- nlme::gls(
      form2,
      data = d, 
      correlation = corAR1(form = ~ annee | region))
    
    # test du modèle
    mod <- out.gls1
    Tab3[i, "forme"] <- "lineaire"
    #if (!j %in% c("ETP_prod_agri",   "nb_chef_exploitation_agri")) {
    aov <- anova(out.gls1, out.gls2)
    if (aov$`p-value`[2] < 0.05) {
      mod <- out.gls2
      Tab3[i, "forme"] <- "quadratique"
    }
    #}
    
    # extraire les betas et leur p-value 
    if (inherits(mod, "gls")) {
      coef_Tab3 <- summary(mod)$tTable[, c("Value", "p-value")]
    } else if (inherits(mod, "lm")) {
      coef_Tab3 <- summary(mod)$coefficients[, c("Estimate", "Pr(>|t|)")]
    } 
    Tab3[i, "betaX"] <- round(coef_Tab3[2,1],2)
    Tab3[i, "p.value_betaX"] <- round(coef_Tab3[2,2],2)
    if (identical(mod, out.gls2)){
      Tab3[i, "betaX2"] <- round(coef_Tab3[3,1],2)
      Tab3[i, "p.value_betaX2"] <- round(coef_Tab3[3,2],2)
    }
    Tab3[i, "phi"] <- round(coef(mod$modelStruct$corStruct, unconstrained = FALSE), 2)
    Tab3[i, "sigma"] <- round(as.numeric(summary(mod)$sigma), 2)
    i <- i+1  
  }
  
  Tab3_show <- Tab3 %>% 
    dplyr::filter(p.value_betaX < 0.05)
  
  # afficher les variables significatives en commun dans les trois approches
  if (print_common_variables == TRUE){
    print(Reduce(intersect, list(Tab_show$variable,Tab2_show$variable,Tab3_show$variable)))}
  
  # sauvegarder les datasets dans excel
  if (save_excel == TRUE){
    FonctionsUtiles::save_datasets_to_excel(datasets = list("GEE" = Tab, "LMM" = Tab2, "GLS" = Tab3),
                                            file_name = paste0("sel_var_Y_", Y,".xlsx"))}
  
  # affichage selon le modèle
  if (model == "GEE") return(DT::datatable(Tab_show))
  if (model == "LMM") return(DT::datatable(Tab2_show))
  if (model == "GLS") return(DT::datatable(Tab3_show))
  
}

# test
results_model(var_to_test = var_to_test, Y = "AMR_EHPAD_FQ_R", model = "GEE")

