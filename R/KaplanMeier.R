#' Kaplan-Meier plot
#'
#' For creating a Kaplan-Meier plot using ggplot2
#'
#' @param data data.frame
#' @param predictor Predictor of interest for plotting from @data.
#' @param event event
#' @param eventage eventage
#'
#' @import ggplot2 survival gridExtra grid dplyr RColorBrewer
#' @return A ggplot / grob
#' @export

# fig1 <- KaplanMeier(data=X, predictor="fishoil", event="event_j45_5yr",eventage= "eventage_j45_5yr", timeby=365, main = NULL, xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)", ystratalabs=c("n3-LCPUFA","Placebo"),ystrataname=NULL,pval = FALSE,returns = TRUE, hr =TRUE,table=TRUE, reverseHR=TRUE, color=T)

# fig1 <- KaplanMeier(data=X, predictor=1, event="event_j45_5yr",eventage= "eventage_j45_5yr", timeby=365, main = "NULL", xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)",ystrataname=NULL,pval = FALSE,returns = TRUE, hr =FALSE,table=FALSE, reverseHR=FALSE, color=FALSE)
#
#Make a pdf
#pdf('./FigsNtabs/PUFA_5YR_plot.pdf',height = 5,width = 5)
#grid.arrange(fig1)
#dev.off()
# 
# Eks. stratify the data (here by prelevel PUFA tertile)
# KaplanMeier(data=X[X$tPufa =="(0.829,4.3]",], predictor="fishoil", event="event_j45_5yr",eventage= "eventage_j45_5yr", timeby=365, main = NULL, xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)", ystratalabs=c("n3-LCPUFA","Placebo"),ystrataname=NULL,pval = FALSE,returns = FALSE, hr =TRUE,table=TRUE, reverseHR=TRUE, color=FALSE)

# Eks. adjust the cox (here by categorical: dvit and sex; and continous: solely_b)

# KaplanMeier(data=X, predictor="fishoil", event="event_j45_5yr",eventage= "eventage_j45_5yr", adjust_cat = c("dvit_all","sex"), adjust_con = "solely_b", timeby=365, main = NULL, xlabs = "Age (years)", ylabs="Risk of persistent wheeze / asthma (%)", ystratalabs=c("n3-LCPUFA","Placebo"),ystrataname=NULL,pval = FALSE,returns = TRUE, hr =TRUE,table=TRUE, reverseHR=TRUE, color=FALSE, ycut=50)
# Changelog
# 16/3-16 JT added "maxtime" to get consistent max times on x axes in stratified plots.
# 29/3-16 JT changed anova() p values to standard parameter (z) p values


KaplanMeier <- function(data, predictor, event, eventage, 
                         adjust_cat = NULL, adjust_con = NULL,
                         table=TRUE, returns = FALSE, ycut=NULL,HR_text=NULL,
                         xlabs = "Time", ylabs = "Survival probability",
                         ystratalabs = NULL, ystrataname = NULL,
                         timeby = 365, maxtime = NULL, main = "Kaplan-Meier Plot",
                         pval = FALSE, rotate = TRUE, years = TRUE, hr = FALSE,
                         color=FALSE, colors=NULL,reverseHR=FALSE,
                         annotation=NULL,reverse_table=FALSE,
                         No_legends=FALSE, No_names=FALSE, low_legend=FALSE) {
  
  require(ggplot2)
  require(survival)
  require(gridExtra)
  require(grid)
  require(dplyr)
  require(RColorBrewer)
  
  sfit <- survfit(Surv(data[,eventage], data[,event] == 1) ~ factor(data[,predictor])) 
  
  if(is.null(ystratalabs)) {
    ystratalabs <- as.character(levels(summary(sfit)$strata))
    ystratalabs <- gsub("factor\\(data\\[\\, predictor\\]\\)\\=", "", ystratalabs)
  }
  m <- max(nchar(ystratalabs))
  if(is.null(ystrataname)) ystrataname <- "Strata"
  if(is.null(maxtime)) maxtime <- sfit$time
  
  if(years) {
    times <- seq(0, round(max(maxtime / 365))*365, by = timeby)
    sfit$time <- sfit$time / 365
    times <- times / 365
  }
  else  times <- seq(0, ceiling(maxtime/timeby)*timeby, by = timeby)
  if(rotate) {
    sfit$upper <- 1 - sfit$upper
    sfit$surv <- 1 - sfit$surv
    sfit$lower <- 1 - sfit$lower
  }
  .df <- data.frame(time = sfit$time, n.risk = sfit$n.risk,
                    n.event = sfit$n.event, surv = sfit$surv, strata = summary(sfit, censored = T)$strata,
                    upper = sfit$upper, lower = sfit$lower)
  levels(.df$strata) <- ystratalabs
  if(rotate) {
    zeros <- data.frame(time = 0, surv = 0, strata = factor(ystratalabs, levels=levels(.df$strata)),upper = 0, lower = 0)
  }  
  else {
    zeros <- data.frame(time = 0, surv = 1.00, strata = factor(ystratalabs, levels=levels(.df$strata)), upper = 1.00, lower = 1.00)
  }
  .df <- bind_rows(zeros, .df)
  d <- length(levels(.df$strata))
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
    scale_x_continuous(xlabs, breaks = times, limits = c(0, max(times))) +
    scale_y_continuous(ylabs, limits = c(0, ifelse(is.null(ycut),1,ycut/100)), labels = scales::percent_format(accuracy = 1)) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(legend.title =element_blank()) +
    theme(plot.title = element_text(size=20)) +
    theme(legend.justification=c(1,ifelse(low_legend,0,1)), legend.position=c(1,ifelse(low_legend,0,1.00))) +
    theme(legend.background = element_rect(colour = NA, fill = NA)) +
    theme(legend.key = element_rect(colour = NA, fill = NA)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(linetype = ystrataname) +
    theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines")) + 
    ggtitle(main) 
  
  
  if(!is.null(annotation)) p= p+ annotate("text", x = -Inf, y = Inf, label = annotation,size=12,color="black", vjust=1.3, hjust=-0.3)
  if(No_legends) {
    p <- p + theme(legend.position="none")
  }
  if(No_names) {
    p <- p +scale_y_continuous(NULL, limits = c(0, ifelse(is.null(ycut),1,ycut/100)), labels = scales::percent_format(accuracy = 1)) +  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  }
  if(pval) {
    sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
    pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
    p <- p + annotate("text", x = 0.6 * max(sfit$time), y = 10, label = pvaltxt)
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
    if(!is.null(HR_text)) {
      HRtext <- HR_text
    }
    p <-p + annotate("text", x = 0.5 * max(times),y=ifelse(is.null(ycut),1,ycut/100), label=HRtext, size = 5)
  }
  
  if(table) {
    ## Create table graphic to include at-risk numbers
    risk.data <- data.frame(strata = summary(sfit, times = times, extend = TRUE)$strata,
                            time = summary(sfit, times = times, extend = TRUE)$time,
                            n.risk = summary(sfit, times = times, extend = TRUE)$n.risk)
    
    if(reverse_table) {
      data.table <- ggplot(risk.data, aes(x = time, y = rev(strata), label = format(n.risk, nsmall = 0)))+  scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = rev(ystratalabs))
    }
    else {  
      data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0)))+  scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = ystratalabs)
    }  
    data.table  <- data.table +
      #, color = strata)) +
      geom_text(size = 3.5) +
      theme_bw() +
      scale_x_continuous(limits = c(0, max(times))) +
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
    if(No_names)  data.table <- data.table +  theme(axis.text.y=element_blank())
    
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