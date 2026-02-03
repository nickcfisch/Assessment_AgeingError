fit_tmb<-function (obj, fn = obj$fn, gr = obj$gr, startpar = NULL, lower = -Inf, 
                   upper = Inf, getsd = TRUE, control = list(eval.max = 10000, 
                                                             iter.max = 10000, trace = 1), bias.correct = FALSE, bias.correct.control = list(sd = FALSE, 
                                                                                                                                             split = NULL, nsplit = NULL, vars_to_correct = NULL), 
                   savedir = NULL, loopnum = 2, newtonsteps = 0, n = Inf, getReportCovariance = FALSE, 
                   getJointPrecision = FALSE, getHessian = FALSE, quiet = FALSE, 
                   start_time_elapsed = as.difftime("0:0:0"), ...) 
{
  if (is.null(startpar)) 
    startpar = obj$par
  if (bias.correct == TRUE & is.null(obj$env$random)) {
    message("No random effects detected in TMB model, so overriding user input to `TMBhelper::fit_tmb` to instead specify `bias.correct=FALSE`")
    bias.correct = FALSE
  }
  if (getReportCovariance == FALSE) {
    message("Note that `getReportCovariance=FALSE` causes an error in `TMB::sdreport` when no ADREPORTed variables are present")
  }
  List = list(...)
  combine_lists = function(default, input) {
    output = default
    for (i in seq_along(input)) {
      if (names(input)[i] %in% names(default)) {
        output[[names(input)[i]]] = input[[i]]
      }
      else {
        output = c(output, input[i])
      }
    }
    return(output)
  }
  BS.control = list(sd = FALSE, split = NULL, nsplit = NULL, 
                    vars_to_correct = NULL)
  BS.control = combine_lists(default = BS.control, input = bias.correct.control)
  nlminb.control = list(eval.max = 10000, iter.max = 10000, 
                        trace = 0)
  nlminb.control = combine_lists(default = nlminb.control, 
                                 input = control)
  start_time = Sys.time()
  parameter_estimates = nlminb(start = startpar, objective = fn, 
                               gradient = gr, control = nlminb.control, lower = lower, 
                               upper = upper)
  for (i in seq(2, loopnum, length = max(0, loopnum - 1))) {
    Temp = parameter_estimates[c("iterations", "evaluations")]
    parameter_estimates = nlminb(start = parameter_estimates$par, 
                                 objective = fn, gradient = gr, control = nlminb.control, 
                                 lower = lower, upper = upper)
    parameter_estimates[["iterations"]] = parameter_estimates[["iterations"]] + 
      Temp[["iterations"]]
    parameter_estimates[["evaluations"]] = parameter_estimates[["evaluations"]] + 
      Temp[["evaluations"]]
  }
  for (i in seq_len(newtonsteps)) {
    g = as.numeric(gr(parameter_estimates$par))
    h = optimHess(parameter_estimates$par, fn = fn, gr = gr)
    parameter_estimates$par = parameter_estimates$par - solve(h, 
                                                              g)
    parameter_estimates$objective = fn(parameter_estimates$par)
  }
  parameter_estimates = parameter_estimates[c("par", 
                                              "objective", "iterations", "evaluations")]
  parameter_estimates[["time_for_MLE"]] = Sys.time() - 
    start_time
  parameter_estimates[["max_gradient"]] = max(abs(gr(parameter_estimates$par)))
  parameter_estimates[["Convergence_check"]] = ifelse(parameter_estimates[["max_gradient"]] < 
                                                        1e-04, "There is no evidence that the model is not converged", 
                                                      "The model is likely not converged")
  parameter_estimates[["number_of_coefficients"]] = c(Total = length(unlist(obj$env$parameters)), 
                                                      Fixed = length(startpar), Random = length(unlist(obj$env$parameters)) - 
                                                        length(startpar))
  parameter_estimates[["AIC"]] = TMBAIC(opt = parameter_estimates)
  if (n != Inf) {
    parameter_estimates[["AICc"]] = TMBAIC(opt = parameter_estimates, 
                                                      n = n)
    parameter_estimates[["BIC"]] = TMBAIC(opt = parameter_estimates, 
                                                     p = log(n))
  }
  parameter_estimates[["diagnostics"]] = data.frame(Param = names(startpar), 
                                                    starting_value = startpar, Lower = lower, MLE = parameter_estimates$par, 
                                                    Upper = upper, final_gradient = as.vector(gr(parameter_estimates$par)))
  if (getsd == TRUE) {
    sd_time = Sys.time()
    h = optimHess(parameter_estimates$par, fn = fn, gr = gr)
    if (is.character(try(chol(h), silent = TRUE))) {
      warning("Hessian is not positive definite, so standard errors are not available")
      if (!is.null(savedir)) {
        capture.output(parameter_estimates, file = file.path(savedir, 
                                                             "parameter_estimates.txt"))
      }
      return(list(opt = parameter_estimates, h = h))
    }
    if (bias.correct == FALSE | is.null(BS.control[["vars_to_correct"]])) {
      if (!is.null(BS.control[["nsplit"]])) {
        if (BS.control[["nsplit"]] == 1) 
          BS.control[["nsplit"]] = NULL
      }
      parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, 
                                                  par.fixed = parameter_estimates$par, hessian.fixed = h, 
                                                  bias.correct = bias.correct, bias.correct.control = BS.control[c("sd", 
                                                                                                                   "split", "nsplit")], getReportCovariance = getReportCovariance, 
                                                  getJointPrecision = getJointPrecision, ...)
    }
    else {
      if ("ADreportIndex" %in% names(obj$env)) {
        Which = as.vector(unlist(obj$env$ADreportIndex()[BS.control[["vars_to_correct"]]]))
      }
      else {
        parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, 
                                                    par.fixed = parameter_estimates$par, hessian.fixed = h, 
                                                    bias.correct = FALSE, getReportCovariance = FALSE, 
                                                    getJointPrecision = FALSE, ...)
        Which = which(rownames(summary(parameter_estimates[["SD"]], 
                                       "report")) %in% BS.control[["vars_to_correct"]])
      }
      if (!is.null(BS.control[["nsplit"]]) && BS.control[["nsplit"]] > 
          1) {
        Which = split(Which, cut(seq_along(Which), BS.control[["nsplit"]]))
      }
      Which = Which[sapply(Which, FUN = length) > 0]
      if (length(Which) == 0) 
        Which = NULL
      message(paste0("Bias correcting ", length(Which), 
                     " derived quantities"))
      parameter_estimates[["SD"]] = TMB::sdreport(obj = obj, 
                                                  par.fixed = parameter_estimates$par, hessian.fixed = h, 
                                                  bias.correct = TRUE, bias.correct.control = list(sd = BS.control[["sd"]], 
                                                                                                   split = Which, nsplit = NULL), getReportCovariance = getReportCovariance, 
                                                  getJointPrecision = getJointPrecision, ...)
    }
    parameter_estimates[["Convergence_check"]] = ifelse(parameter_estimates$SD$pdHess == 
                                                          TRUE, parameter_estimates[["Convergence_check"]], 
                                                        "The model is definitely not converged")
    parameter_estimates[["time_for_sdreport"]] = Sys.time() - 
      sd_time
    if (getHessian == TRUE) {
      parameter_estimates[["hessian"]] = h
    }
  }
  parameter_estimates[["time_for_run"]] = Sys.time() - 
    start_time + start_time_elapsed
  if (!is.null(savedir)) {
    save(parameter_estimates, file = file.path(savedir, "parameter_estimates.RData"))
    capture.output(parameter_estimates, file = file.path(savedir, 
                                                         "parameter_estimates.txt"))
  }
  if (quiet == FALSE & parameter_estimates[["Convergence_check"]] != 
      "There is no evidence that the model is not converged") {
    message("#########################")
    message(parameter_estimates[["Convergence_check"]])
    message("#########################")
  }
  return(parameter_estimates)
}


TMBAIC<-function (opt, p = 2, n = Inf) {
    k = length(opt[["par"]])
    if (all(c("par", "objective") %in% names(opt))) 
        negloglike = opt[["objective"]]
    if (all(c("par", "value") %in% names(opt))) 
        negloglike = opt[["value"]]
    Return = p * k + 2 * negloglike + 2 * k * (k + 1)/(n - k - 
        1)
    return(Return)
}