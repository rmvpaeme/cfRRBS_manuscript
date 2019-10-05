#Helper functions 

#(from cookbook for R)
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



### convert to log format, rescale, calculate statistics (sd, sem, CI), convert to original scale
log_rescaling_CI <- function(data, measurevar, groupvar=NULL, log_base=10, conf=0.95){
  require(dplyr)
  # log transform
  data_resc <- data %>% mutate(measurevar_log = log(get(measurevar), as.numeric(log_base)))
  # calculate mean for each group (log scale)
  log_mean <- data_resc %>% group_by(get(groupvar)) %>% 
    dplyr::summarise(mean_logscale = mean(measurevar_log))
  # rescale accoring to min (log scale)
  data_resc <- data_resc %>% mutate(measurevar_log_resc = measurevar_log - min(log_mean$mean_logscale))
  
  #Calculate statistics (sd, sem, 95% CI)
  tmp <- summarySE(data_resc, measurevar="measurevar_log_resc", groupvars=paste(groupvar), conf.interval=as.numeric(paste(conf)))
  tmp <- tmp %>% dplyr::rename(mean_log_resc = measurevar_log_resc) %>%
    mutate(mean_oriscale = log_base^mean_log_resc, 
           ci_lower_oriscale = log_base^(mean_log_resc-ci), 
           ci_upper_oriscale =log_base^(mean_log_resc+ci))
  full_join(tmp,data_resc, by=paste(groupvar))
}


# ggplot(test, aes(x=GraphKit, y=measurevar_oriscale, colour=RNAisolation)) + 
#   geom_point() +  
#   geom_hline(yintercept = 1, linetype="dashed",colour="grey88") +
#   geom_errorbar(aes(ymin=ci_lower_oriscale, ymax=ci_upper_oriscale), width=.1) +
#   geom_line() +
#   labs(x="", y="relative amount of genes",title="Number of genes (unfiltered)") +
#   scale_colour_manual(values=color_panel2) +
#   scale_y_continuous(limits=c(0,NA)) +
#   theme_point +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
# 
# OR
# 
# ggplot(test, aes(x=GraphKit, y=10^measurevar_log_resc, colour=RNAisolation)) + 
#   geom_point() +  
#   geom_hline(yintercept = 1, linetype="dashed",colour="grey88") +
#   geom_errorbar(aes(ymin=10^(measurevar_log_resc-ci), ymax=10^(measurevar_log_resc+ci)), width=.1) +
#   geom_line() +
#   labs(x="", y="relative amount of genes",title="Number of genes (unfiltered)") +
#   scale_colour_manual(values=color_panel2) +
#   scale_y_continuous(limits=c(0,NA)) +
#   theme_point +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))


### plot multiple graphs next to each other and share legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  require(grid)
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: http://goo.gl/K4yh

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

#p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)