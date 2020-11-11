#3.0 Fully generic function for 3 cancer types, PSA included, results for preprint paper

library(tidyr)
library(dplyr)
library(formattable)
library(flexsurv)  # log-logistic and gompertz distributions
library(ggplot2)   # plotting survival curves
library(gtools)

##  ADAPTED FROM THE MODEL OF DEGELING ET AL, DETAILS AS FOLLOWS:
##  This code was written to perform the analyses presented in the manuscript:
##
##  An inverse stage-shift model to estimate the excess mortality and health 
##  economic impact of delayed access to cancer services due to the COVID-19 pandemic
##
##  by Koen Degeling, Nancy N Baxter, Jon Emery, Fanny Franchini, Peter Gibbs,
##  G Bruce Mann, Grant McArthur, Benjamin J Solomon, Maarten J IJzerman
##
##  For more information or when using/adapting the code, please refer to/cite the
##  pre-print version of the manuscript available from MedRxiv:
##  https://www.medrxiv.org/content/10.1101/2020.05.30.20117630v1
##
##  This code was written by Dr Koen Degeling, Cancer Health Services Research,
##  Centre for Cancer Research & Centre for Health Policy, Faculty of Medicine,
##  Dentistry and Health Sciences, The University of Melbourne, Melbourne, Australia
##
##  The code was reviewed by Dr Fanny Franchini, Cancer Health Services Research,
##  Centre for Cancer Research & Centre for Health Policy, Faculty of Medicine,
##  Dentistry and Health Sciences, The University of Melbourne, Melbourne, Australia
##
##  For questions or other enquiries, please contact Koen Degeling by email:
##  koen.degeling@unimelb.edu.au
##

### 1. INITIALIZATION ----

# Clear environment
# rm(list = ls()); gc();

# Loading packages
#library(dplyr);     # handling data.frames
#library(flexsurv);  # log-logistic and gompertz distributions
#library(ggplot2);   # plotting survival curves

# Global parameters
days_in_year <- 365.25;

### 2. FUNCTIONS TO PERFORM THE ANALYSIS ----

# In this section custom functions are defined to perform the analysis. 
# The 'fun_fit_distrs' function fits parameteric distributions based on quantiles and 
# corresponding percentiles.
# The 'fun_plot_distrs' function plots the survival data and fitted distributions for 
# graphical inspection. 
# The 'fun_integrate_distrs' function determines the expected survival in life years 
# of a certain time horizon for a distribution fitted by the 'fun_fit_distrs' function 
# by integrating the corresponding survival curve.


fun_fit_distrs <- function(q, p) {
  
  # This function fits a series of parametric distributions to quantiles and corresponding 
  # percentiles. It does so through nonlinear least squares regression on the parameters on 
  # log-scale. It fits the exponential, Gamma, Gompertz, log-logistic, log-Normal and
  # Weibull distributions. It return the fitted distribution parameters and sum of the
  # squared residuals in a list.
  #
  # INPUTS:
  # - p           vector of percentiles for which q provided quantiles
  # - q           vector of quantiles corresponding to the percentiles provided to p
  #
  # USES:
  # - pgompertz   function for the Gompertz distribution from the flexsurv package
  # - pllogis     function for the log-logistic distribution from the flexsurv package
  #
  # OUTPUTS:
  # - list        a list with two named objects:
  #               - data      
  #               - distrs    list named according to the distribution names: exponential, 
  #                           gamma, gompertz, loglogistic, lognormal, and weibull. Each 
  #                           list item contains a list of three ojects:
  #                           - pars          the estimated parameters on log scale
  #                           - convergence   logical indicating convergence
  #                           - SSM           the sum of the squared residuals
  
  # Merge the two vectors 'p' and 'q' into a data.frame for use in the call of the 
  # nonlinear least squares regression. Also sort
  df_fit <- data.frame(p = p, q = q);
  
  # The start value for the exponential distribution is estimated based on the latest
  # survival data point. To easily extract that data point, we sort the data first.
  df_fit <- arrange(df_fit, q);
  
  # Estimating the mortality rate according to the standard rate-probability formula
  start_q <- tail(df_fit, 1)$q;
  start_p <- tail(df_fit, 1)$p;
  start_rate <- -(1/start_q)*log(start_p);
  
  # The estimated rate is used as start value for the regression. Note that parameters
  # are estimated on log scale.
  fit_exponential <- nls(formula = p ~ pexp(q = q, rate = exp(log_rate), lower.tail = FALSE),
                         start = list(log_rate = log(start_rate)),
                         data = df_fit); 
  
  # The Weibull distribution is fitted subsequently, because start values for its
  # parameters can be defined based on the rate of the exponential distribution and its
  # parameters themselves are useful to define start values for the other distributions.
  # The start values of the Weibull distribution used basically resemble the fitted
  # exponential distribution.
  fit_weibull <- nls(formula = p ~ pweibull(q = q, shape = exp(log_shape), scale = exp(log_scale), lower.tail = FALSE),
                     start = list(log_shape = log(1), log_scale = -coef(fit_exponential)["log_rate"]),
                     data = df_fit);
  
  # The rate of exponential distribution is used to define the start value for the rate
  # parameter of the Gamma distribution, and the start value for the shape parameters is
  # based on that of the Weibull distribution.
  fit_gamma <- nls(formula = p ~ pgamma(q = q, shape = exp(log_shape), rate = exp(log_rate), lower.tail = FALSE),
                   start = list(log_shape = coef(fit_weibull)["log_shape"], log_rate = coef(fit_exponential)["log_rate"]),
                   data = df_fit);
  
  # The Gompertz distribution is the most difficult one to fit. Because it is most similar
  # to the Gamma distribution, the fitted parameters of the Gamma distribution are used as
  # start values. If fitting the Gompertz is unsuccessful using built-in function, it is
  # performed in a custom way using optim directly.
  fit_gompertz <- tryCatch({nls(formula = p ~ pgompertz(q = q, shape = exp(log_shape), rate = exp(log_rate), lower.tail = FALSE),
                                start = list(log_shape = coef(fit_gamma)["log_shape"], log_rate = coef(fit_gamma)["log_rate"]),
                                data = df_fit)}, error = function(e) NULL);
  if(is.null(fit_gompertz)) {
    fit_gompertz <- optim(
      par = c(coef(fit_gamma)["log_shape"], coef(fit_gamma)["log_rate"]),
      fn = function(x, p, q) sum((p - pgompertz(q = q, shape = exp(x["log_shape"]), rate = exp(x["log_rate"]), lower.tail = FALSE))^2),
      p = df_fit$p,
      q = df_fit$q);
  }
  
  # For the log-logistic distribution, start values for the parameters are based on the 
  # parameters of the Weibull distribution.
  fit_loglogistic <- nls(formula = p ~ pllogis(q = q, shape = exp(log_shape), scale = exp(log_scale), lower.tail = FALSE),
                         start = list(log_shape = coef(fit_weibull)["log_shape"], log_scale = coef(fit_weibull)["log_scale"]),
                         data = df_fit);
  
  # For the log-Normal distribution, start values for both parameters are based on the scale
  # parameter of the Weibull distribution. Note that this is the only distribution that does
  # not require the parameters to be transformed back to the original scale.
  fit_lognormal <- tryCatch({nls(formula = p ~ plnorm(q = q, meanlog = log_mean, sdlog = log_sd, lower.tail = FALSE),
                       start = list(log_mean = coef(fit_weibull)["log_scale"], log_sd = coef(fit_weibull)["log_scale"]),
                       data = df_fit)}, error = function(e) NULL);
  if(is.null(fit_lognormal)) {
    fit_lognormal <- nls(formula = p ~ pexp(q = q, rate = exp(log_rate), lower.tail = FALSE),
        start = list(log_rate = log(start_rate)),
        data = df_fit); 
  }
  
  # Define a list 'list_out' that contains the data used for fitting as argument 'data' 
  # and within argument 'distrs' for each distribution its parameters, convergence and 
  # sum of the squared residuals
  list_out <- list(
    
    # The data.frame used for fitting
    data = df_fit,
    
    # List including all distributions
    distrs = list(
      
      # Exponential
      exponential = list(
        pars = coefficients(fit_exponential),
        convergence = fit_exponential$convInfo$isConv,
        SSR = sum(residuals(fit_exponential)^2)
      ),
      
      # Gamma
      gamma = list(
        pars = coefficients(fit_gamma),
        convergence = fit_gamma$convInfo$isConv,
        SSR = sum(residuals(fit_gamma)^2)
      ),
      
      # Gompertz
      gompertz = if(class(fit_gompertz) == "nls") {
        list(
          pars = coefficients(fit_gompertz),
          convergence = fit_gompertz$convInfo$isConv,
          SSR = sum(residuals(fit_gompertz)^2)
        )
      } else {
        gompertz = list(
          pars = fit_gompertz$par,    
          convergence = (fit_gompertz$convergence == 0),
          SSR = fit_gompertz$value
        )
      },
      
      # Log-normal
      lognormal = list(
        pars = coefficients(fit_lognormal),
        convergence = fit_lognormal$convInfo$isConv,
        SSR = sum(residuals(fit_lognormal)^2)
      ),
      
      # Log-logistic
      loglogistic = list(
        pars = coefficients(fit_loglogistic),
        convergence = fit_loglogistic$convInfo$isConv,
        SSR = sum(residuals(fit_loglogistic)^2)
      ),
      
      # Weibull
      weibull = list(
        pars = coefficients(fit_weibull),
        convergence = fit_weibull$convInfo$isConv,
        SSR = sum(residuals(fit_weibull)^2)
      )
      
    )
  );
  
  # Return the list
  return(list_out);
  
};


fun_plot_distrs <- function(list_out, time_horizon) {
  
  # This function plots the survival data and fitted parametric distributions based on the
  # list_out object returned by the fun_fit_distrs function.
  #
  # INPUTS:
  # - list_out        list returned by the fun_fit_distrs function, containing the data and
  #                   fitted distributions
  # - time_horizon    the time horizon for which the parametric distributions are to be
  #                   plotted
  #
  # USES:
  # - ggplot          plotting functions from the ggplot2 package
  #
  # OUTPUTS:
  # - plot            the survival plot
  
  # Extracting the survival data and distributions from the list_out object.
  df_fit <- list_out$data;
  distrs <- list_out$distrs;
  
  # Setting some general parameters for the plot. For the type and width of the lines, there
  # will be a difference between the actual survival data (obs = observed) and the simulated
  # data according to the parametric distributions (sim = simulated).
  x_values <- seq(from = 0, to = time_horizon, by = 1/12);
  line_width <- c(obs = 1.5, sim = 1);
  line_type <- c(obs = 1, sim = 2);
  colour_blind_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7");
  
  # Plot the distributions in order from generally worst to best fitting distributions, so 
  # that the latter are plotted over the former
  ggplot() +
    
    # General plotting parameters
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = colour_blind_palette) + 
    labs(title = "Observed vs. simulated survival",
         x = "Time in Years",
         y = "Survival Probability") +
    theme(legend.title = element_blank()) + 
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14,face = "bold"),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14, face = "bold")) + 
    
    # Observed data
    geom_line(aes(x = df_fit$q, 
                  y = df_fit$p, 
                  colour = "Data"),
              size = line_width["obs"],
              linetype = line_type["obs"]) + 
    
    # Exponential
    geom_line(aes(x = x_values, 
                  y = pexp(q = x_values, 
                           rate = exp(distrs$exponential$pars["log_rate"]),
                           lower.tail = FALSE),
                  colour = "Exponential"),
              size = line_width["sim"],
              linetype = line_type["sim"]) + 
    
    # Log-logistic
    geom_line(aes(x = x_values, 
                  y = pllogis(q = x_values, 
                              shape = exp(distrs$loglogistic$pars["log_shape"]),
                              scale = exp(distrs$loglogistic$pars["log_scale"]),
                              lower.tail = FALSE),
                  colour = "Log-logistic"),
              size = line_width["sim"],
              linetype = line_type["sim"]) + 
    
    # Log-normal
    geom_line(aes(x = x_values, 
                  y = plnorm(q = x_values, 
                             meanlog = distrs$lognormal$pars["log_mean"],
                             sdlog = distrs$lognormal$pars["log_sd"],
                             lower.tail = FALSE),
                  colour = "Log-normal"),
              size = line_width["sim"],
              linetype = line_type["sim"]) + 
    
    # Gamma
    geom_line(aes(x = x_values, 
                  y = pgamma(q = x_values, 
                             shape = exp(distrs$gamma$pars["log_shape"]),
                             rate = exp(distrs$gamma$pars["log_rate"]),
                             lower.tail = FALSE),
                  colour = "Gamma"),
              size = line_width["sim"],
              linetype = line_type["sim"]) +
    
    # Weibull
    geom_line(aes(x = x_values, 
                  y = pweibull(q = x_values, 
                               shape = exp(distrs$weibull$pars["log_shape"]),
                               scale = exp(distrs$weibull$pars["log_scale"]),
                               lower.tail = FALSE),
                  colour = "Weibull"),
              size = line_width["sim"],
              linetype = line_type["sim"]) +  
    
    # Gompertz
    geom_line(aes(x = x_values, 
                  y = pgompertz(q = x_values, 
                                shape = exp(distrs$gompertz$pars["log_shape"]),
                                rate = exp(distrs$gompertz$pars["log_rate"]),
                                lower.tail = FALSE),
                  colour = "Gompertz"),
              size = line_width["sim"],
              linetype = line_type["sim"]);
  
};


fun_integrate_dist <- function(dist, log_pars, time_horizon) { ## Koen, I have edited the function goal here
  
  # This function returns the expected survival in life years over a chosen time horizon 
  # for a distribution fitted using the 'fun_fit_distrs' function, by integrating the 
  # corresponding survival curve.
  #
  # //delete// This function plots the survival data and fitted parametric distributions based on the
  # list_out objected returned by the fun_fit_distrs function. //delete//
  #
  # INPUTS:
  # - dist            character identifying the distribution type
  # - log_pars        vector with the parameters of the distribution on log scale
  # - time_horizon    the time horizon for which the expected survival is to be calculated
  #
  # USES:
  # - pgompertz       function for the Gompertz distribution from the flexsurv package
  # - pllogis         function for the log-logistic distribution from the flexsurv package
  #
  # OUTPUTS:
  # - out             expected survival in life years over the time horizon 
  
  # Initialize the out object to NULL so it can be easily be checked whether the integration
  # was successful.
  out <- NULL;
  
  # Exponential
  if(dist == "exponential") {
    out <- integrate(
      f = pexp,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      rate = exp(log_pars["log_rate"])
    )$value;
  }
  
  # Gamma
  if(dist == "gamma") {
    out <- integrate(
      f = pgamma,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      rate = exp(log_pars["log_rate"])
    )$value;
  }
  
  # Gompertz
  if(dist == "gompertz") {
    out <- integrate(
      f = pgompertz,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      rate = exp(log_pars["log_rate"])
    )$value;
  }
  
  # Log-logistic
  if(dist == "loglogistic") {
    out <- integrate(
      f = pllogis,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      scale = exp(log_pars["log_scale"])
    )$value;
  }
  
  # Log-normal
  if(dist == "lognormal") {
    out <- integrate(
      f = plnorm,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      meanlog = log_pars["log_mean"],
      sdlog = log_pars["log_sd"]
    )$value;
  }
  
  # Weibull
  if(dist == "weibull") {
    out <- integrate(
      f = pweibull,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      scale = exp(log_pars["log_scale"])
    )$value;
  }
  
  # Check whether integration successful
  if(is.null(out)) stop("No integration has been performed. Perhaps the distribution is not supported?");
  
  # The time_horizon is included in the name of the object for convenience
  names(out) <- paste0("time_horizon_", time_horizon);
  
  # Return the output
  return(out);
  
};


####STAGE SHIFT MODELLING DEFINED AS A FUNCTION####

#############stage_shift function##########
##INPUTS: TTI hazard ratio (scalar), stage-specific survival data (data frame 5 years (columns) by 5 stages (rows)), total annual incidence (scalar), incidence for each stage (vector), whether to use within-stage TTI HR (1=use)

##OUTPUTS: List: [1] Table with excess mortality over 5 years, LY lost over 5 years, LY lost over 10 years, for both 3 and 6 month delays, [2] Predicted stage distribution, [3] predicted survival no delay, [4] predicted survival for 3 month delay, [5] predited survival for 3m and 6m delay w TTI HR applied

#Model function:

stage_shift <- function(HR_TTI,df_survival_stage,tot_incidence,stage_numbers, within_stage) {
  # #testing
  # HR_TTI <- Base_HR_can_delay_I
  # df_survival_stage <- Base_df_survival_can
  # tot_incidence <- Base_n_incidence_can
  # stage_numbers <- Base_v_stage_numbers_can
  
  ### 3. STAGE SHIFT MODELLING ----
  
  # This section contains the actual stage shift modelling for each of the cancer types
  # considered. Stage shifts considered are: from stage I to II, III or IV, and from II to III or IV. Transitions to III/IV are assumed
  # to result in maintaining the observed III/IV ratio. No transitions between stage III and IV are included. Stage III and IV are 
  # assumed to present as normal.
  
  ## 3.1 SPECIFIC CANCER INPUT DATA----
  
  ## DATA USED FOR MODELLING
  
  # The stage-specific survival data:
  df_survival_can <- df_survival_stage;
  
  # The predicted incidence for 2020
  n_incidence_can <- tot_incidence;
  
  # The distribution of the disease stage at treatment initiation:
  v_stage_numbers_can <- stage_numbers
  v_stage_proportions_can <- v_stage_numbers_can / sum(v_stage_numbers_can);
  
  # The costs of treating a specific stage of disease have been approximated from the
  # publication by Goldsbury et al. (2018):
  v_costs_can <- c(I = 50699, II = 67069, III = 98632, IV = 105934, unknown = 100860);
  
  
  ## ESTIMATING TIME TO STAGE PROGRESSION (TTSP)
  
  # For can cancer, TTSP is estimated by approximating how long a delay needs to be for the
  # 5-year survival of stage I patients to match that of stage II patients based on the hazard
  # ratio of a surgical delay on overall survival.
  
  # The required hazard ratio for a stage I patient to match the 5-year survival of a stage II 
  # patients is calculated as the ratio of the mortality rates. The 5-year survival of 100% 
  # for stage I patients is problematic, because it suggests no events, which does not allow 
  # a rate to be calculated. Therefore, we assume the 5-year survival probability to be 99%.
  #No longer 100% survival for stage I in data so correction no longer needed.
  #v2.1 converted to transitions between all stages, including transitioning up multiple stages, note that
  # transitions to stage 3/4 are handled so that stage 3/4 ratio remains the same. Assumes presentations at 3
  # are less afffected by disruption as more likely to be urgent/emergency.
  
  mortality_rate_can_I <- -(1/5) * log(filter(df_survival_can, disease_stage == "I")$survival_5y);
  mortality_rate_can_II <- -(1/5) * log(filter(df_survival_can, disease_stage == "II")$survival_5y);
  mortality_rate_can_III <- -(1/5) * log(filter(df_survival_can, disease_stage == "III")$survival_5y);
  mortality_rate_can_IV <- -(1/5) * log(filter(df_survival_can, disease_stage == "IV")$survival_5y);
  
  HR_can_required_I_II <- mortality_rate_can_II / mortality_rate_can_I;
  HR_can_required_I_III <- mortality_rate_can_III / mortality_rate_can_I;
  HR_can_required_I_IV <- mortality_rate_can_IV / mortality_rate_can_I;
  HR_can_required_II_III <- mortality_rate_can_III / mortality_rate_can_II;
  HR_can_required_II_IV <- mortality_rate_can_IV / mortality_rate_can_II;
  HR_can_required_III_IV <- mortality_rate_can_IV / mortality_rate_can_III;
  
  # Hazard ratio of surgical delays on overall survival for stage I patients:
  # - HR on overall survival is 1.018 for a 1-week delay (Khorana, 2019).
  # Can add different HR for each stage progression.
  HR_can_delay_I <- HR_TTI;
  HR_can_delay_II <- HR_TTI;
  HR_can_delay_III <- HR_TTI;
  
  t_can_delay_I <- 7 / days_in_year * 12;     # in months
  
  # The HRs are converted to regression coefficients by transferring to log scale to model
  # the coefficient as a linear function of time
  coef_can_delay_I <- log(HR_can_delay_I) / t_can_delay_I;
  coef_can_delay_II <- log(HR_can_delay_II) / t_can_delay_I;
  coef_can_delay_III <- log(HR_can_delay_III) / t_can_delay_I;
  
  coef_can_required_I_II <- log(HR_can_required_I_II);
  coef_can_required_I_III <- log(HR_can_required_I_III);
  coef_can_required_I_IV <- log(HR_can_required_I_IV);
  coef_can_required_II_III <- log(HR_can_required_II_III);
  coef_can_required_II_IV <- log(HR_can_required_II_IV);
  coef_can_required_III_IV <- log(HR_can_required_III_IV);
  
  # TTSP (in months) can now be calculated
  TTSP_can_I_II <- coef_can_required_I_II/ coef_can_delay_I;
  TTSP_can_I_III <- coef_can_required_I_III / coef_can_delay_I;
  TTSP_can_I_IV <- coef_can_required_I_IV / coef_can_delay_I;
  TTSP_can_II_III <- coef_can_required_II_III / coef_can_delay_II;
  TTSP_can_II_IV <- coef_can_required_II_IV / coef_can_delay_II;
  TTSP_can_III_IV <- coef_can_required_III_IV / coef_can_delay_III;
  
  # proportion of stage 3 in 3/4 at baseline
  stage_III_IV_prop <- v_stage_numbers_can[3]/(v_stage_numbers_can[3]+v_stage_numbers_can[4])
  
  ## CALCULATING THE STAGE SHIFT
  
  # Based on the TTSP, the proportion of patients who are expected to start treatment at a
  # more advanced stage can be calculated. Since only the expected TTSP is known, and not its
  # distribution, the only distribution that can be used is the exponential distribution for
  # which the rate parameter is defined as: 1 / expected value.
  
  #Calculate transition probabilities sequentially from smallest to largest
  # i.e. To progress to stage III need to also progress to stage II
  
  # Proportion of patients that will progress following a 3-month delay, adjusted for the fact
  # that the service disruption will only be for 3 months, i.e. 25% of 2020 patients:
  prop_progress_can_I_IV_3mo <- pexp(q = 3, rate = 1 / TTSP_can_I_IV)* (3/12);
  prop_progress_can_I_III_3mo <- pexp(q = 3, rate = 1 / TTSP_can_I_III)* (3/12) - prop_progress_can_I_IV_3mo;
  prop_I_3to4 <- (prop_progress_can_I_IV_3mo + prop_progress_can_I_III_3mo);
  prop_progress_can_I_III_3mo <- prop_I_3to4*stage_III_IV_prop
  prop_progress_can_I_IV_3mo <- prop_I_3to4*(1 - stage_III_IV_prop)
  prop_progress_can_I_II_3mo <-  pexp(q = 3, rate = 1 / TTSP_can_I_II) * (3/12) - (prop_progress_can_I_IV_3mo + prop_progress_can_I_III_3mo);
  prop_progress_can_II_IV_3mo <- pexp(q = 3, rate = 1 / TTSP_can_II_IV)* (3/12);
  prop_progress_can_II_III_3mo <- pexp(q = 3, rate = 1 / TTSP_can_II_III)*(3/12) - prop_progress_can_II_IV_3mo;
  prop_II_3to4 <- (prop_progress_can_II_IV_3mo + prop_progress_can_II_III_3mo);
  prop_progress_can_II_III_3mo <- prop_II_3to4*stage_III_IV_prop
  prop_progress_can_II_IV_3mo <- prop_II_3to4*(1 - stage_III_IV_prop)
  prop_progress_can_III_IV_3mo <- 0 #pexp(q = 3, rate = 1 / TTSP_can_III_IV)* (3/12); assume stage 3 still present
  
  # Proportion of patients that will progress following a 6-month delay, adjusted for the fact
  # that the service disruption will only be for 6 months, i.e. 50% of 2020 patients:
  prop_progress_can_I_IV_6mo <- pexp(q = 6, rate = 1 / TTSP_can_I_IV)* (6/12);
  prop_progress_can_I_III_6mo <- pexp(q = 6, rate = 1 / TTSP_can_I_III)*(6/12) - prop_progress_can_I_IV_6mo;
  prop_I_3to4 <- (prop_progress_can_I_IV_6mo + prop_progress_can_I_III_6mo);
  prop_progress_can_I_III_6mo <- prop_I_3to4*stage_III_IV_prop
  prop_progress_can_I_IV_6mo <- prop_I_3to4*(1 - stage_III_IV_prop)
  prop_progress_can_I_II_6mo <-  pexp(q = 6, rate = 1 / TTSP_can_I_II) * (6/12) - (prop_progress_can_I_IV_6mo + prop_progress_can_I_III_6mo);
  prop_progress_can_II_IV_6mo <- pexp(q = 6, rate = 1 / TTSP_can_II_IV)* (6/12);
  prop_progress_can_II_III_6mo <- pexp(q = 6, rate = 1 / TTSP_can_II_III)*(6/12) - prop_progress_can_II_IV_6mo;
  prop_II_3to4 <- (prop_progress_can_II_IV_6mo + prop_progress_can_II_III_6mo);
  prop_progress_can_II_III_6mo <- prop_II_3to4*stage_III_IV_prop
  prop_progress_can_II_IV_6mo <- prop_II_3to4*(1 - stage_III_IV_prop)
  prop_progress_can_III_IV_6mo <- 0 #pexp(q = 6, rate = 1 / TTSP_can_III_IV)* (6/12); assume stage 3 still present
  
  
  # Based on these proportions, the distribution of disease stage at treatment initiation
  # following a 3 month or 6 month service disruption and delay can be calculated:
  v_stage_proportions_can_3mo <- c(
    I = v_stage_proportions_can["I"] * (1 - (prop_progress_can_I_II_3mo + prop_progress_can_I_III_3mo + prop_progress_can_I_IV_3mo)),
    II = v_stage_proportions_can["II"] + (v_stage_proportions_can["I"] * prop_progress_can_I_II_3mo) - (v_stage_proportions_can["II"] * prop_progress_can_II_III_3mo) - (v_stage_proportions_can["II"] * prop_progress_can_II_IV_3mo),
    III = v_stage_proportions_can["III"] + (v_stage_proportions_can["I"] * prop_progress_can_I_III_3mo) + (v_stage_proportions_can["II"] * prop_progress_can_II_III_3mo) - (v_stage_proportions_can["III"] * prop_progress_can_III_IV_3mo),
    IV = v_stage_proportions_can["IV"] + (v_stage_proportions_can["I"] * prop_progress_can_I_IV_3mo) + (v_stage_proportions_can["II"] * prop_progress_can_II_IV_3mo) + (v_stage_proportions_can["III"] * prop_progress_can_III_IV_3mo),
    unknown = v_stage_proportions_can["unknown"]
  );
  
  v_stage_proportions_can_6mo <- c(
    I = v_stage_proportions_can["I"] * (1 - (prop_progress_can_I_II_6mo + prop_progress_can_I_III_6mo + prop_progress_can_I_IV_6mo)),
    II = v_stage_proportions_can["II"] + (v_stage_proportions_can["I"] * prop_progress_can_I_II_6mo) - (v_stage_proportions_can["II"] * prop_progress_can_II_III_6mo) - (v_stage_proportions_can["II"] * prop_progress_can_II_IV_6mo),
    III = v_stage_proportions_can["III"] + (v_stage_proportions_can["I"] * prop_progress_can_I_III_6mo) + (v_stage_proportions_can["II"] * prop_progress_can_II_III_6mo) - (v_stage_proportions_can["III"] * prop_progress_can_III_IV_6mo),
    IV = v_stage_proportions_can["IV"] + (v_stage_proportions_can["I"] * prop_progress_can_I_IV_6mo) + (v_stage_proportions_can["II"] * prop_progress_can_II_IV_6mo) + (v_stage_proportions_can["III"] * prop_progress_can_III_IV_6mo),
    unknown = v_stage_proportions_can["unknown"]
  );
  
  #Add additional delay hazard - assumes proportional hazards
  # - determine if this is to be used based on function input
  if(within_stage != 1){HR_can_delay_I <- 1}
  #HR for 3 & 6 month delay periods (12 and 24 weeks)
  HR_can_delay_3m <- exp(log(HR_can_delay_I)*(3*4.345))
  HR_can_delay_6m <- exp(log(HR_can_delay_I)*(6*4.345))
  #Calculate the survival probabilities with this delay added
  #3 month
  df_survival_can_delay_3m <- df_survival_can %>%
    select(contains("survival")) # select columns
  df_survival_can_delay_3m <- df_survival_can_delay_3m ^ HR_can_delay_3m #apply HR
  df_survival_can_delay_3m <- bind_cols(select(df_survival_can,contains("stage")),df_survival_can_delay_3m) # remake matrix
  #6 month
  df_survival_can_delay_6m <- df_survival_can %>%
    select(contains("survival")) # select columns
  df_survival_can_delay_6m <- df_survival_can_delay_6m ^ HR_can_delay_6m #apply HR
  df_survival_can_delay_6m <- bind_cols(select(df_survival_can,contains("stage")),df_survival_can_delay_6m) # remake matrix
  
  # Checking whether the distribution are still appropriate, i.e. sum up to 1
  #sum(v_stage_proportions_can_3mo);
  #sum(v_stage_proportions_can_6mo);
  
  
  ## ESTIMATING THE IMPACT OF THE STAGE SHIFTS
  
  # According to each distribution of stages at treatment initiation, the combined survival
  # is calculated by weighting the year-specific probabilities
  
  # Extracting the survival data from the data.frame with stage-specific data
  m_survival_can <- select(df_survival_can, contains("survival"));
  m_survival_can_3m <- select(df_survival_can_delay_3m, contains("survival"));
  m_survival_can_6m <- select(df_survival_can_delay_6m, contains("survival"));
  
  # Calculating the weighted survival
  v_survival_can_base <- apply(m_survival_can, 2, function(col) sum(v_stage_proportions_can * col));
  v_survival_can_3mo <- apply(m_survival_can_3m, 2, function(col) sum(v_stage_proportions_can_3mo * col));
  v_survival_can_6mo <- apply(m_survival_can_6m, 2, function(col) sum(v_stage_proportions_can_6mo * col));
  
  # To allow for extrapolation and integration to calculate expected survival in terms of 
  # life years, a series of distributions are fitted to the data: exponential, Gamma, Gompertz
  # log-Normal, log-logistic, and Weibull distributions. These are fitted based on the quantiles
  # and percentiles of the weighted survival using the custom fun_fit_distrs function. The most
  # appropriate distribution is selected by graphical inspection using the 'fun_plot_distrs'
  # function and based on the sum of the squared residuals of the non-linear least squares 
  # regression models used to estimate the parameters on log scale.
  
  # Define a factor with the quantiles corresponding to the weighted percentiles
  v_q <- 0:5;
  
  # Selecting the distribution: baseline
  # - Based on that beyond 5 years background mortality is more likely than cancer-specific
  #   mortality, and Gompertz distributions are known to represent background mortality well,
  #   the Gompertz is selected, providing the most conservative estimate of survival.
  l_distr_can_base <- fun_fit_distrs(q = v_q, p = v_survival_can_base); 
  fun_plot_distrs(list_out = l_distr_can_base, time_horizon = 10);
  sapply(l_distr_can_base$distrs, function(l) l$SSR);
  dist_can_base <- "exponential";
  
  # Selecting the distribution: 3 month disruption and delay
  # - Given that the weighted survival for the 3 month disruption and delay scenario is almost
  #   identical to that of the baseline scenario, and no differences between survival due to
  #   different parametric distributions should be introduced, the Gompertz is selected.
  l_distr_can_3mo <- fun_fit_distrs(q = v_q, p = v_survival_can_3mo); 
  fun_plot_distrs(list_out = l_distr_can_3mo, time_horizon = 10);
  sapply(l_distr_can_3mo$distrs, function(l) l$SSR);
  dist_can_3mo <- "exponential";
  
  # Selecting the distribution: 6 month disruption and delay
  # - Same as for the 3 month disruption and delay
  l_distr_can_6mo <- fun_fit_distrs(q = v_q, p = v_survival_can_6mo); 
  fun_plot_distrs(list_out = l_distr_can_6mo, time_horizon = 10);
  sapply(l_distr_can_6mo$distrs, function(l) l$SSR);
  dist_can_6mo <- "exponential";
  
  # Obtaining the expected survival in life years for two time horizons: 5 and 10 years
  # - Because the built-in integrate function used within the fun_integrate_dist function
  #   is not vectorized by default, a loop approach is taken using sapply. Alternatively,
  #   the fun_integrate_dist function can be vectorized
  v_time_horizons <- c(5, 10);
  
  v_ly_can_base <- sapply(v_time_horizons, function(time_horizon) {
    fun_integrate_dist(dist = dist_can_base,
                       log_pars = l_distr_can_base$distrs[[dist_can_base]]$pars,
                       time_horizon = time_horizon)
  });
  
  v_ly_can_3mo <- sapply(v_time_horizons, function(time_horizon) {
    fun_integrate_dist(dist = dist_can_3mo,
                       log_pars = l_distr_can_3mo$distrs[[dist_can_3mo]]$pars,
                       time_horizon = time_horizon)
  });
  
  v_ly_can_6mo <- sapply(v_time_horizons, function(time_horizon) {
    fun_integrate_dist(dist = dist_can_6mo,
                       log_pars = l_distr_can_6mo$distrs[[dist_can_6mo]]$pars,
                       time_horizon = time_horizon)
  });
  
  # The expected costs of treating the 2020 population are estimated by weighting the
  # stage-specific healthcare costs by the distribution of disease stage at treatment
  # initatiation according to the different scenarios
  # - No differences in costs because we only consider stage I -> stage II shifts and
  #   the treatment costs are expected to be similar
  costs_can_base <- sum(v_stage_proportions_can * v_costs_can);
  costs_can_3mo <- sum(v_stage_proportions_can_3mo * v_costs_can);
  costs_can_6mo <- sum(v_stage_proportions_can_6mo * v_costs_can);
  
    #(1 - v_survival_can_base[6])*n_incidence_can # Total deaths from annual cohort
  
  ### MAKE OUTPUT TABLE###
  #Excess mortality, excess mortality as % of expected mortality, LY lost over 5 year horizon & 10 year horizon, 3 month and 6 month delays from disruption
  Output_table <- data.frame(
    Duration = c("3 months","6 months"),
    excess_mortality = c((v_survival_can_base[6] - v_survival_can_3mo[6]) * n_incidence_can,(v_survival_can_base[6] - v_survival_can_6mo[6]) * n_incidence_can),
    excess_mortality_p = c((v_survival_can_base[6] - v_survival_can_3mo[6])/(1-v_survival_can_base[6]),(v_survival_can_base[6] - v_survival_can_6mo[6])/(1-v_survival_can_base[6]))*100,
    life_years_lost_5 = c(((v_ly_can_base[1] - v_ly_can_3mo[1]) * n_incidence_can),((v_ly_can_base[1] - v_ly_can_6mo[1]) * n_incidence_can)),
    life_years_lost_10 = c(((v_ly_can_base[2] - v_ly_can_3mo[2]) * n_incidence_can),((v_ly_can_base[2] - v_ly_can_6mo[2]) * n_incidence_can))
  )
  
  #Table of stage distributions
  stage_dist_out <- as.data.frame(rbind(v_stage_proportions_can, v_stage_proportions_can_3mo,v_stage_proportions_can_6mo))
  
  #Matrix of weighted survival estimates in first 5 years - to allow calculation of excess mortality for other years
  v_survival_tab <- rbind(v_survival_can_base, v_survival_can_3mo, v_survival_can_6mo)
  
  #List containing all outputs
  Output_ls <- list(Output_table,stage_dist_out,df_survival_can,df_survival_can_delay_3m,df_survival_can_delay_6m, v_survival_tab)
  
  return(Output_ls)
}

############################################################################################################################

####LUNG BASE CASE######
#Base case parameters - lung
Base_HR_lung_delay_I <- 1.032

# The stage-specific survival data for colorectal cancer is taken from: - Now UK Data from England (ONS) "Cancer Survival in England: adults diagnosed between 2013 and 2017 and followed up to 2018"
Base_df_survival_lung <- data.frame(
  disease_stage = c("I", "II", "III", "IV", "unknown"),
  survival_0y = c(1.000, 1.000, 1.000, 1.000, 1.000),
  survival_1y = c(0.852, 0.687, 0.460, 0.177, 0.519),
  survival_2y = c(0.752, 0.561, 0.325, 0.110, 0.453),
  survival_3y = c(0.663, 0.459, 0.230, 0.068, 0.395),
  survival_4y = c(0.585, 0.375, 0.163, 0.042, 0.344),
  survival_5y = c(0.516, 0.306, 0.115, 0.026, 0.300)
);


# The predicted incidence for 2020 is also taken from the average of 2017/2018 in Scottish cancer registry:
Base_n_incidence_lung <-9243/2;
# The distribution of the disease stage at treatment initiation is based on Scottish registry 2017-2018 sum:
Base_v_stage_numbers_lung <- c(I = 1658, II = 688, III = 1995, IV = 4235, unknown = 667);

#Use stage-shift function

Output_lung_base <- stage_shift(Base_HR_lung_delay_I,Base_df_survival_lung,Base_n_incidence_lung,Base_v_stage_numbers_lung,0)

stage_dist_table <- Output_lung_base[[2]]*100
row.names(stage_dist_table) <- c("No delay","3 months","6 months")
formattable(stage_dist_table, align=c("c","c"),digits = 3)

Mortality_tab <- Output_lung_base[[1]]
#Rounding
Mortality_tab$excess_mortality <- round(Mortality_tab$excess_mortality,2)
Mortality_tab$excess_mortality_p <- round(Mortality_tab$excess_mortality_p,2)
Mortality_tab$life_years_lost_5 <- round(Mortality_tab$life_years_lost_5,2)
Mortality_tab$life_years_lost_10 <- round(Mortality_tab$life_years_lost_10,2)
#Use formattable function to produce the output
formattable(Mortality_tab,col.names = c("Duration of disruption","Excess mortality","Excess mortality %","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"))

###Colorectal base case###
#Base case parameters - colorectal

Base_HR_colorectal_delay_I <- 1.005
# Stage specific suvival from SCAN (Lothian only?) data 2009-2014 followed to 2019.
Base_df_survival_colorectal <- data.frame(
  disease_stage = c("I", "II", "III", "IV", "unknown"),
  survival_0y = c(1.000, 1.000, 1.000, 1.000, 1.000),
  survival_1y = c(0.979, 0.937, 0.902, 0.451, 0.722),
  survival_2y = c(0.95, 0.892, 0.807, 0.24, 0.58),
  survival_3y = c(0.924, 0.845, 0.723, 0.149, 0.49),
  survival_4y = c(0.876, 0.792, 0.664, 0.109, 0.427),
  survival_5y = c(0.844, 0.739, 0.619, 0.093, 0.391)
);


# The predicted incidence for 2020 is also taken from the average of 2017/2018 in Scottish cancer registry:
Base_n_incidence_colorectal <-6732/2;
# The distribution of the disease stage at treatment initiation is based on Scottish registry 2017-2018 sum:
Base_v_stage_numbers_colorectal <- c(I = 1105, II = 1612, III = 1681, IV = 1604, unknown = 730);

#Apply stage-shift function
Output_colorectal_base <- stage_shift(Base_HR_colorectal_delay_I,Base_df_survival_colorectal,Base_n_incidence_colorectal,Base_v_stage_numbers_colorectal,0)

stage_dist_table <- Output_colorectal_base[[2]]*100
row.names(stage_dist_table) <- c("No delay","3 months","6 months")
formattable(stage_dist_table, align=c("c","c"),digits = 3)

Mortality_tab <- Output_colorectal_base[[1]]
#Rounding
Mortality_tab$excess_mortality <- round(Mortality_tab$excess_mortality,2)
Mortality_tab$excess_mortality_p <- round(Mortality_tab$excess_mortality_p,2)
Mortality_tab$life_years_lost_5 <- round(Mortality_tab$life_years_lost_5,2)
Mortality_tab$life_years_lost_10 <- round(Mortality_tab$life_years_lost_10,2)
#Use formattable function to produce the output
formattable(Mortality_tab,col.names = c("Duration of disruption","Excess mortality","Excess mortality %","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"))

####Breast Cancer#####
#Base case parameters - breast
Base_HR_breast_delay_I <- 1.018
# The stage-specific survival data for breast cancer is taken from Scottish registry data using data from diagnoses 2010-2018:
Base_df_survival_breast <- data.frame(
  disease_stage = c("I", "II", "III", "IV", "unknown"),
  survival_0y = c(1.000, 1.000, 1.000, 1.000, 1.000),
  survival_1y = c(0.995, 0.989, 0.962, 0.929, 0.994),
  survival_2y = c(0.986, 0.964, 0.880, 0.789, 0.973),
  survival_3y = c(0.973, 0.935, 0.816, 0.739, 0.938),
  survival_4y = c(0.959, 0.906, 0.769, 0.676, 0.899),
  survival_5y = c(0.943, 0.879, 0.721, 0.65, 0.845)
);
# The predicted incidence for 2020 is also taken from the average of 2017/2018 in Scottish cancer registry:
Base_n_incidence_breast <- 8814/2;
# The distribution of the disease stage at treatment initiation is based on Scottish registry 2017-2018 sum:
Base_v_stage_numbers_breast <- c(I = 3561, II = 4042, III = 664, IV = 450, unknown = 94);

#Apply stage-shift function
Output_breast_base <- stage_shift(Base_HR_breast_delay_I,Base_df_survival_breast,Base_n_incidence_breast,Base_v_stage_numbers_breast,0)

stage_dist_table <- Output_breast_base[[2]]*100
row.names(stage_dist_table) <- c("No delay","3 months","6 months")
formattable(stage_dist_table, align=c("c","c"),digits = 3)

Mortality_tab <- Output_breast_base[[1]]
#Rounding
Mortality_tab$excess_mortality <- round(Mortality_tab$excess_mortality,2)
Mortality_tab$excess_mortality_p <- round(Mortality_tab$excess_mortality_p,2)
Mortality_tab$life_years_lost_5 <- round(Mortality_tab$life_years_lost_5,2)
Mortality_tab$life_years_lost_10 <- round(Mortality_tab$life_years_lost_10,2)
#Use formattable function to produce the output
formattable(Mortality_tab,col.names = c("Duration of disruption","Excess mortality","Excess mortality %","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"))

#################Exess Mortality at different time points###################
#Lung
Excess_tab_lung <- as.data.frame(Output_lung_base[[6]])
Excess_tab_lung[4,] <- (Excess_tab_lung[1,] - Excess_tab_lung[2,]) * Base_n_incidence_lung
Excess_tab_lung[5,] <- (Excess_tab_lung[1,] - Excess_tab_lung[3,]) * Base_n_incidence_lung
Excess_tab_lung <- cbind(Mortality_tab$Duration, Excess_tab_lung[4:5,])
formattable(Excess_tab_lung,col.names = c("Duration of disruption","Year 0","Year 1","Year 2", "Year 3", "Year 4", "Year 5"), align=c("l","c"), row.names = FALSE, digits = 3)

#Colorectal
Excess_tab_colorectal <- as.data.frame(Output_colorectal_base[[6]])
Excess_tab_colorectal[4,] <- (Excess_tab_colorectal[1,] - Excess_tab_colorectal[2,]) * Base_n_incidence_colorectal
Excess_tab_colorectal[5,] <- (Excess_tab_colorectal[1,] - Excess_tab_colorectal[3,]) * Base_n_incidence_colorectal
Excess_tab_colorectal <- cbind(Mortality_tab$Duration, Excess_tab_colorectal[4:5,])
formattable(Excess_tab_colorectal,col.names = c("Duration of disruption","Year 0","Year 1","Year 2", "Year 3", "Year 4", "Year 5"), align=c("l","c"), row.names = FALSE, digits = 3)

#Breast
Excess_tab_breast <- as.data.frame(Output_breast_base[[6]])
Excess_tab_breast[4,] <- (Excess_tab_breast[1,] - Excess_tab_breast[2,]) * Base_n_incidence_breast
Excess_tab_breast[5,] <- (Excess_tab_breast[1,] - Excess_tab_breast[3,]) * Base_n_incidence_breast
Excess_tab_breast <- cbind(Mortality_tab$Duration, Excess_tab_breast[4:5,])
formattable(Excess_tab_breast,col.names = c("Duration of disruption","Year 0","Year 1","Year 2", "Year 3", "Year 4", "Year 5"), align=c("l","c"), row.names = FALSE, digits = 3)

###############PSA###############################
#Get 95% CI for excess mortality estimates accounting for sampling uncertainty in TTI HR and stage specific incidence

#Generate random draw of ITT HR & stage distribution
#TTI_HR - using standard errors from Khorana 2019
#lung 1.032 (1.031-1.034)
#colorectal  1.005 (1.002-1.008)
#breast 1.018 (1.015-1.020) [HR for stage 1 - used for all transitions]

#RNG
set.seed(100) # seed for replication
n_rep <- 1000 #number of PSA iterations
#values of TTI for each cancer
TTI_HR_ran_lung <-  exp(rnorm(n = n_rep,mean = log(1.032),sd = 0.0008)) 
TTI_HR_ran_colorectal <-  exp(rnorm(n = n_rep,mean = log(1.005),sd = 0.0015)) 
TTI_HR_ran_breast <-  exp(rnorm(n = n_rep,mean = log(1.018),sd = 0.001)) 
#Values of number in each stage for each cancer:
v_stage_numbers_lung_ran <- rdirichlet(n_rep, Base_v_stage_numbers_lung) #values of stage distribution proportions
v_stage_numbers_lung_ran <- as.data.frame(round(v_stage_numbers_lung_ran*Base_n_incidence_lung,digits = 0)) %>%
  rename( "I" =V1, "II" = V2, "III" = V3, "IV" = V4, "unknown" = V5) #number not proportion

v_stage_numbers_colorectal_ran <- rdirichlet(n_rep, Base_v_stage_numbers_colorectal) #values of stage distribution proportions
v_stage_numbers_colorectal_ran <- as.data.frame(round(v_stage_numbers_colorectal_ran*Base_n_incidence_colorectal,digits = 0)) %>%
  rename( "I" =V1, "II" = V2, "III" = V3, "IV" = V4, "unknown" = V5) #number not proportion

v_stage_numbers_breast_ran <- rdirichlet(n_rep, Base_v_stage_numbers_breast) #values of stage distribution proportions
v_stage_numbers_breast_ran <- as.data.frame(round(v_stage_numbers_breast_ran*Base_n_incidence_breast,digits = 0)) %>%
  rename( "I" =V1, "II" = V2, "III" = V3, "IV" = V4, "unknown" = V5) #number not proportion

#Data frame to store results
PSA_output <- as.data.frame(matrix(ncol = 12,nrow = length(TTI_HR_ran_lung),data = 0)) %>%
  rename( "3_months_lung" =V1, "6_months_lung" = V2, "3_months_lung_LY" = V3, "6_months_lung_LY" = V4, "3_months_colorectal" = V5, "6_months_colorectal" = V6, "3_months_colorectal_LY" = V7, "6_months_colorectal_LY" = V8, "3_months_breast" =V9, "6_months_breast" = V10, "3_months_breast_LY" = V11, "6_months_breast_LY" = V12)

#Loop to run model at random values and keep results
j <- 1 #counter
for (x in 1:length(TTI_HR_ran_lung)){
  
iter_v_stage_numbers_lung <- c(I = v_stage_numbers_lung_ran[x,1], II = v_stage_numbers_lung_ran[x,2], III = v_stage_numbers_lung_ran[x,3], IV = v_stage_numbers_lung_ran[x,4], unknown = v_stage_numbers_lung_ran[x,5])

iter_v_stage_numbers_colorectal <- c(I = v_stage_numbers_colorectal_ran[x,1], II = v_stage_numbers_colorectal_ran[x,2], III = v_stage_numbers_colorectal_ran[x,3], IV = v_stage_numbers_colorectal_ran[x,4], unknown = v_stage_numbers_colorectal_ran[x,5])

iter_v_stage_numbers_breast <- c(I = v_stage_numbers_breast_ran[x,1], II = v_stage_numbers_breast_ran[x,2], III = v_stage_numbers_breast_ran[x,3], IV = v_stage_numbers_breast_ran[x,4], unknown = v_stage_numbers_breast_ran[x,5])


  Output_lung <- stage_shift(TTI_HR_ran_lung[x],Base_df_survival_lung,Base_n_incidence_lung,iter_v_stage_numbers_lung,0)
  Mortality_tab <- Output_lung[[1]]
  PSA_output[j,1] <- c(Mortality_tab$excess_mortality[1])
  PSA_output[j,2] <- c(Mortality_tab$excess_mortality[2])
  PSA_output[j,3] <- c(Mortality_tab$life_years_lost_10[1])
  PSA_output[j,4] <- c(Mortality_tab$life_years_lost_10[2])
  
  Output_colorectal <- stage_shift(TTI_HR_ran_colorectal[x],Base_df_survival_colorectal,Base_n_incidence_colorectal,iter_v_stage_numbers_colorectal,0)
  Mortality_tab <- Output_colorectal[[1]]
  PSA_output[j,5] <- c(Mortality_tab$excess_mortality[1]) 
  PSA_output[j,6] <- c(Mortality_tab$excess_mortality[2])
  PSA_output[j,7] <- c(Mortality_tab$life_years_lost_10[1])
  PSA_output[j,8] <- c(Mortality_tab$life_years_lost_10[2])
  
  Output_breast <- stage_shift(TTI_HR_ran_breast[x],Base_df_survival_breast,Base_n_incidence_breast,iter_v_stage_numbers_breast,0)
  Mortality_tab <- Output_breast[[1]]
  PSA_output[j,9] <- c(Mortality_tab$excess_mortality[1]) 
  PSA_output[j,10] <- c(Mortality_tab$excess_mortality[2])
  PSA_output[j,11] <- c(Mortality_tab$life_years_lost_10[1])
  PSA_output[j,12] <- c(Mortality_tab$life_years_lost_10[2])
  

j <- j + 1  
  }

#Get 95% interval for excess deaths over 5 years for 3 & 6 month delays for all 3 cancers
quantile(PSA_output[,1],probs = c(0.025,0.975))
quantile(PSA_output[,2],probs = c(0.025,0.975))
quantile(PSA_output[,5],probs = c(0.025,0.975))
quantile(PSA_output[,6],probs = c(0.025,0.975))
quantile(PSA_output[,9],probs = c(0.025,0.975))
quantile(PSA_output[,10],probs = c(0.025,0.975))

quantile(PSA_output[,3],probs = c(0.025,0.975))
quantile(PSA_output[,4],probs = c(0.025,0.975))
quantile(PSA_output[,7],probs = c(0.025,0.975))
quantile(PSA_output[,8],probs = c(0.025,0.975))
quantile(PSA_output[,11],probs = c(0.025,0.975))
quantile(PSA_output[,12],probs = c(0.025,0.975))

#####OWSA##################################################
#Lung
#Loop to generate values for line graph of excess deaths (3 & 6 month delay) increasing TTI HR range 1.001 -> 1.05

TTI_HR_list <- c(1.001,1.005,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1) #values of TTI to use
TTI_HR_excess <- as.data.frame(matrix(ncol = 3,nrow = length(TTI_HR_list),data = 0)) %>%
  rename("TTI_HR" =V1, "3_months" =V2, "6_months" = V3)

j <- 1 #counter
for (x in TTI_HR_list){
  Output_1 <- stage_shift(x,Base_df_survival_lung,Base_n_incidence_lung,Base_v_stage_numbers_lung,0)[[1]]
  TTI_HR_excess[j,] <- c(x,Output_1$excess_mortality[1],Output_1$excess_mortality[2])
  j <- j + 1
}
#Reshape to long format to use ggplot with groups
TTI_HR_excess <- TTI_HR_excess %>%
  pivot_longer(-"TTI_HR", names_to = "Delay", values_to = "Excess_deaths")

ggplot(data=TTI_HR_excess, aes(x=TTI_HR, y=Excess_deaths, group=Delay)) +
  geom_line(aes(linetype=Delay))+
  geom_point(aes(x = 1.032, y = 48.01), colour = "black") +
  geom_point(aes(x = 1.032, y = 150.19), colour = "black") +
  geom_vline(aes(xintercept = 1.032), linetype = "twodash") +
  xlab("Time to treatment initiation hazard ratio") +
  ylab("Excess deaths at 5 years in Scotland") +
  labs(linetype = "Duration of Delay") 

#Colorectal
#Loop to generate values for line graph of excess deaths (3 & 6 month delay) increasing TTI HR range 1.001 -> 1.01

TTI_HR_list <- c(1.001,1.005,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1) #values of TTI to use
TTI_HR_excess <- as.data.frame(matrix(ncol = 3,nrow = length(TTI_HR_list),data = 0)) %>%
  rename("TTI_HR" =V1, "3_months" =V2, "6_months" = V3)

j <- 1 #counter
for (x in TTI_HR_list){
  Output_1 <- stage_shift(x,Base_df_survival_colorectal,Base_n_incidence_colorectal,Base_v_stage_numbers_colorectal,0)[[1]]
  TTI_HR_excess[j,] <- c(x,Output_1$excess_mortality[1],Output_1$excess_mortality[2])
  j <- j + 1
}
#Reshape to long format to use ggplot with groups
TTI_HR_excess <- TTI_HR_excess %>%
  pivot_longer(-"TTI_HR", names_to = "Delay", values_to = "Excess_deaths")

ggplot(data=TTI_HR_excess, aes(x=TTI_HR, y=Excess_deaths, group=Delay)) +
  geom_line(aes(linetype=Delay))+
  geom_point(aes(x = 1.005, y = 14.69), colour = "black") +
  geom_point(aes(x = 1.005, y = 55.43), colour = "black") +
  geom_vline(aes(xintercept = 1.005), linetype = "twodash") +
  xlab("Time to treatment initiation hazard ratio") +
  ylab("Excess deaths at 5 years in Scotland") +
  labs(linetype = "Duration of Delay")

#Breast
#Loop to generate values for line graph of excess deaths (3 & 6 month delay) increasing TTI HR range 1.01 -> 1.1

TTI_HR_list <- c(1.001,1.005,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1) #values of TTI to use
TTI_HR_excess <- as.data.frame(matrix(ncol = 3,nrow = length(TTI_HR_list),data = 0)) %>%
  rename("TTI_HR" =V1, "3_months" =V2, "6_months" = V3)

j <- 1 #counter
for (x in TTI_HR_list){
  Output_1 <- stage_shift(x,Base_df_survival_breast,Base_n_incidence_breast,Base_v_stage_numbers_breast,0)[[1]]
  TTI_HR_excess[j,] <- c(x,Output_1$excess_mortality[1],Output_1$excess_mortality[2])
  j <- j + 1
}
#Reshape to long format to use ggplot with groups
TTI_HR_excess <- TTI_HR_excess %>%
  pivot_longer(-"TTI_HR", names_to = "Delay", values_to = "Excess_deaths")

ggplot(data=TTI_HR_excess, aes(x=TTI_HR, y=Excess_deaths, group=Delay)) +
  geom_line(aes(linetype=Delay))+
  geom_point(aes(x = 1.018, y = 38.7), colour = "black") +
  geom_point(aes(x = 1.018, y = 139.15), colour = "black") +
  geom_vline(aes(xintercept = 1.018), linetype = "twodash") +
  xlab("Time to treatment initiation hazard ratio") +
  ylab("Excess deaths at 5 years in Scotland") +
  labs(linetype = "Duration of Delay")
