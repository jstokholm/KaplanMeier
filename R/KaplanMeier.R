#' Kaplan-Meier plot
#'
#' For creating a Kaplan-Meier plot using ggplot2
#'
#' @param data data.frame
#' @param predictor Predictor of interest for plotting from @data.
#' @param event event
#' @param eventage eventage
#' @param table Display survival table below graph; default is TRUE.
#' @param timeby Cut points for the table; default is 365 (days).
#' @param years X-axis in years - else in days; default is TRUE.
#' @param reverse_table reverse the upper/lower part of table; default is FALSE.
#' @param hr Make and plot cox proportional HR; default is TRUE.
#' @param reverseHR Change case/control categories for HR; default is FALSE.
#' @param adjust_cat Categorical variables to adjust the cox HR for.
#' @param adjust_con Continues variables to adjust the cox HR for.
#' @param returns outputs the ggplot/grob object; default is TRUE.
#' @param ycut Cut y-axis at a certain percentage; default is NULL.
#' @param xlabs X-axis label
#' @param ylabs Y-axis label
#' @param main Plot title.
#' @param legend_title Legend title; default is name of predictor.
#' @param legend_names Define variable names for legend text.
#' @param mintime To get consistent min times on x axes in stratified plots.
#' @param maxtime To get consistent max times on x axes in stratified plots.
#' @param rotate Switch the predictor order; default is FALSE.
#' @param color Color plot, else black/greytone; default is TRUE.
#' @param colors Set of defined colors for predictor groups; default is NULL.
#' @param reverse_colors Flip the predictor group; default is FALSE.
#' @param annotation Paste an annotation in the corner; default is NULL.
#' @param no_legends Removes legend; default is FALSE.
#' @param no_names Removes y names; default is FALSE.
#' @param low_legend Put legends in the bottom, else in the top corner;
#' @param y_tick_n Number of y-ticks; default is NULL.
#'
#' @import ggplot2 survival gridExtra grid dplyr RColorBrewer
#' @return A ggplot / grob
#' @export

# fig1 <- KaplanMeier(data=X, predictor="fishoil", event="event_j45_5yr",eventage= "eventage_j45_5yr", timeby=365, main = NULL, xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)", legend_names=c("n3-LCPUFA","Placebo"),legend_title=NULL,returns = TRUE, hr =TRUE,table=TRUE, reverseHR=TRUE, color=TRUE)

# fig1 <- KaplanMeier(data=X, predictor=1, event="event_j45_5yr",eventage= "eventage_j45_5yr", timeby=365, main = "NULL", xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)",legend_title=NULL,returns = TRUE, hr =FALSE,table=FALSE, reverseHR=FALSE, color=FALSE)
#
#Make a pdf
#pdf('./FigsNtabs/PUFA_5YR_plot.pdf',height = 5,width = 5)
#grid.arrange(fig1)
#dev.off()
#
# Eks. stratify the data (here by prelevel PUFA tertile)
# KaplanMeier(data=X[X$tPufa =="(0.829,4.3]",], predictor="fishoil", event="event_j45_5yr",eventage= "eventage_j45_5yr", timeby=365, main = NULL, xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)", legend_names=c("n3-LCPUFA","Placebo"),legend_title=NULL,returns = FALSE, hr =TRUE,table=TRUE, reverseHR=TRUE, color=FALSE)

# Eks. adjust the cox (here by categorical: dvit and sex; and continous: solely_b)

# KaplanMeier(data=X, predictor="fishoil", event="event_j45_5yr",eventage= "eventage_j45_5yr", adjust_cat = c("dvit_all","sex"), adjust_con = "solely_b", timeby=365, main = NULL, xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)", legend_names=c("n3-LCPUFA","Placebo"),legend_title=NULL,returns = TRUE, hr =TRUE,table=TRUE, reverseHR=TRUE, color=FALSE, ycut=50)

KaplanMeier <- function(data,
                        predictor,
                        event,
                        eventage,
                        table=TRUE,
                        reverse_table=FALSE,
                        timeby = 365,
                        years = TRUE,
                        hr = TRUE,
                        reverseHR=FALSE,
                        adjust_cat = NULL,
                        adjust_con = NULL,
                        returns = TRUE,
                        ycut=NULL,
                        xlabs = "Time",
                        ylabs = "Survival probability",
                        main = "Kaplan-Meier Plot",
                        legend_title = NULL,
                        legend_names = NULL,
                        mintime = NULL,
                        maxtime = NULL,
                        rotate = TRUE,
                        color=TRUE,
                        colors=NULL,
                        reverse_colors=FALSE,
                        annotation=NULL,
                        no_legends=FALSE,
                        no_names=FALSE,
                        low_legend=FALSE,
                        y_tick_n=NULL) {

  # require(ggplot2)
  # require(survival)
  # require(gridExtra)
  # require(grid)
  # require(dplyr)
  # require(RColorBrewer)
  data <- data.frame(data)
  sfit <- survfit(Surv(data[,eventage], data[,event] == 1) ~ factor(data[,predictor]))

  if(is.null(legend_names)) legend_names <-as.character(levels(factor(data[,predictor])))
  m <- max(nchar(legend_names))
  if(is.null(legend_title)) legend_title <- "Strata"
  if(is.null(mintime)) mintime <- 0
  if(is.null(maxtime)) maxtime <- max(sfit$time)

  if(years) {
    times <- seq(mintime, round(maxtime / 365)*365, by = timeby)
    sfit$time <- sfit$time / 365
    times <- times / 365
  }
  else  times <- seq(mintime, ceiling(maxtime/timeby)*timeby, by = timeby)
  if(rotate) {
    sfit$upper <- 1 - sfit$upper
    sfit$surv <- 1 - sfit$surv
    sfit$lower <- 1 - sfit$lower
  }
  .df <- data.frame(time = sfit$time, n.risk = sfit$n.risk,
                    n.event = sfit$n.event, surv = sfit$surv, strata = summary(sfit, censored = T)$strata,
                    upper = sfit$upper, lower = sfit$lower)
  levels(.df$strata) <- legend_names
  if(rotate) {
    zeros <- data.frame(time = mintime, surv = 0, strata = factor(legend_names, levels=levels(.df$strata)),upper = 0, lower = 0)
  }
  else {
    zeros <- data.frame(time = mintime, surv = 1.00, strata = factor(legend_names, levels=levels(.df$strata)), upper = 1.00, lower = 1.00)
  }
  .df <- bind_rows(zeros, .df)
  .df[.df$time > max(times),"time"] <- max(times) #fix for curves being cut;
  d <- length(levels(.df$strata))
  if(reverse_colors) .df$strata <- factor(.df$strata, levels=rev(levels(.df$strata)))
  if(color) {
    if(is.null(colors)){
      cols  <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Dark2"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(7,"Accent"),brewer.pal(12,"Paired"),"gray")
      cols <- cols[1:d]
    }
    else cols <- colors

    p <- ggplot(.df, aes(time, surv, color = strata))+geom_step( size = 0.7) + scale_color_manual(values=cols)
  }
  else {
    p <- ggplot(.df, aes(time, surv, group = strata)) + geom_step(aes(linetype = strata), size = 0.7)
  }#ceiling(max(.df$surv)/0.05)*0.05
  p <- p +theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    scale_x_continuous(xlabs, breaks = times, limits = c(mintime, max(times)))  +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(legend.title =element_blank()) +
    theme(plot.title = element_text(size=20)) +
    theme(legend.justification=c(1,ifelse(low_legend,0,1)), legend.position=c(1,ifelse(low_legend,0,1.00))) +
    theme(legend.background = element_rect(colour = NA, fill = NA)) +
    theme(legend.key = element_rect(colour = NA, fill = NA)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(linetype = legend_title) +
    theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines")) +
    ggtitle(main)

  if(!is.null(y_tick_n)) {
    p <- p + scale_y_continuous(ylabs, limits = c(0, ifelse(is.null(ycut),1,ycut/100)), breaks = seq(0, ifelse(is.null(ycut),1,ycut/100), by = ifelse(is.null(ycut),1,ycut/100)/y_tick_n), labels = scales::percent_format(accuracy = 1))
    }
  else {
    p <- p + scale_y_continuous(ylabs, limits = c(0, ifelse(is.null(ycut),1,ycut/100)), labels = scales::percent_format(accuracy = 1))
  }


  if(!is.null(annotation)) p <- p+ annotate("text", x = -Inf, y = Inf, label = annotation,size=12,color="black", vjust=1.3, hjust=-0.3)
  if(no_legends) {
    p <- p + theme(legend.position="none")
  }

  if(hr) {
    if(is.null(adjust_cat) & is.null(adjust_con))
    {
      fit <- coxph(Surv(data[,eventage], data[,event] == 1) ~ factor(data[,predictor]))
    }
    if(!is.null(adjust_cat) & !is.null(adjust_con))
    {
      fit <- coxph(as.formula(paste0("Surv(",eventage,",",event," == 1) ~ factor(",predictor,")", paste0("+ factor(",adjust_cat,")", collapse=" + ") , paste0("+ ",adjust_con, collapse=" + " ))), data=data)
    }
    if(!is.null(adjust_cat) & is.null(adjust_con))
    {
      fit <- coxph(as.formula(paste0("Surv(",eventage,",",event," == 1) ~ factor(",predictor,")", paste0("+ factor(",adjust_cat,")", collapse=" + "))), data=data)
    }
    if(is.null(adjust_cat) & !is.null(adjust_con))
    {
      fit <- coxph(as.formula(paste0("Surv(",eventage,",",event," == 1) ~ factor(",predictor,")", paste0("+ ",adjust_con, collapse=" + " ))), data=data)
    }
    fits <- summary(fit)
    pval <- fits$coefficients[1,5]
    if(reverseHR) {
      hr <- exp(-fits$coefficients[1,1])
      hr_low <- exp(-fits$coefficients[1,1] - 1.96*fits$coefficients[1,3])
      hr_high <- exp(-fits$coefficients[1,1] + 1.96*fits$coefficients[1,3])
    }
    else {
      hr <- exp(fits$coefficients[1,1])
      hr_low <- exp(fits$coefficients[1,1] - 1.96*fits$coefficients[1,3])
      hr_high <- exp(fits$coefficients[1,1] + 1.96*fits$coefficients[1,3])
    }
    if(is.null(adjust_cat) & is.null(adjust_con)) {
      HRtext <- paste0("HR ", sprintf("%.2f",round(hr,2)), " [", sprintf("%.2f",round(hr_low,2)), "-",sprintf("%.2f",round(hr_high,2)),"], p ",ifelse(pval<0.001,"< 0.001",paste0("= ",format.pval(pval,1,0.001,nsmall=3))))
    }
    else {
      HRtext <- paste0("aHR ", sprintf("%.2f",round(hr,2)), " [", sprintf("%.2f",round(hr_low,2)), "-",sprintf("%.2f",round(hr_high,2)),"], p ",ifelse(pval<0.001,"< 0.001",paste0("= ",format.pval(pval,1,0.001,nsmall=3))))
    }
    p <-p + annotate("text", x = mean(times),y=ifelse(is.null(ycut),1,ycut/100), label=HRtext, size = 5)
  }

  if(table) {
    ## Create table graphic to include at-risk numbers
    risk.data <- data.frame(strata = summary(sfit, times = times, extend = TRUE)$strata,
                            time = summary(sfit, times = times, extend = TRUE)$time,
                            n.risk = summary(sfit, times = times, extend = TRUE)$n.risk)

    if(reverse_table) {
      data.table <- ggplot(risk.data, aes(x = time, y = rev(strata), label = format(n.risk, nsmall = 0)))+  scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = rev(legend_names))
    }
    else {
      data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0)))+  scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = legend_names)
    }
    data.table  <- data.table +
      #, color = strata)) +
      geom_text(size = 3.5) +
      theme_bw() +
      scale_x_continuous(limits = c(mintime, max(times))) +
      theme(axis.title.x = element_text(size = 10, vjust = 1), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.border = element_blank(),
            axis.text.x = element_blank(), axis.ticks = element_blank(),
            axis.text.y = element_text(face = "bold", hjust = 1))
    data.table <- data.table + theme(legend.position = "none") +
      xlab(NULL)+ ylab("")
    data.table <- data.table +
      theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines"))
    #Adjustment of y-axis for long stratanames;
    p <-p   +    theme(axis.title.y =element_text(margin=margin(r = ifelse(m < 8, 0, -6*m), unit = "pt"))) #R >3.22
    #   p <-p   +    theme(axis.title.y = element_text(vjust =ifelse(m < 8, 0, -0.4*m))) #R <=3.22
    if(no_names) {
        p <- p +  theme(axis.title.y =element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
        data.table <- data.table +  theme(axis.text.y=element_blank())
      }

    p <- ggplot_gtable(ggplot_build(p))
    data.table <- ggplot_gtable(ggplot_build(data.table))
    # maxWidth = unit.pmax(p$widths[2:3], data.table$widths[2:3])
    # p$widths[2:3] <- maxWidth
    #
    # data.table$widths[2:3] <- maxWidth
    maxWidth = unit.pmax(p$widths, data.table$widths)

    p$widths = as.list(maxWidth)
    data.table$widths = as.list(maxWidth)

    grid.arrange(p, data.table,
                 clip = FALSE, nrow = 2, ncol = 1,
                 heights = unit(c(1,d*0.06),c("null", "null"))) #fix label
    if(returns) {
      a <- arrangeGrob(p, data.table, clip = FALSE,
                       nrow = 2, ncol = 1, heights = unit(c(1,d*0.06),c("null", "null")))
      return(a)
    }
  }
  else {
    print(p)
    if(returns) return(p)
  }
}
