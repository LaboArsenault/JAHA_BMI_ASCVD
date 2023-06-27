library("dplyr")
library("data.table")
library("survival")
library("survminer")
library("gtsummary")
library("survivalAnalysis")
library("stringr")
library("plotly")
library("ggplot2")
library("gt")
library("gridExtra")
library("forestplot")
library(boot)
library(rms)

load("/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Data/DF_All_02.RData")

## Incidence ASCVD 1=maladie 0=pas maladie 2=maladie avant rec -----------------------------------------------------
DF_All$INC_ASCVD<-if_else(DF_All$DATE_ALL_MACE>DF_All$REC_DATE,1,2,missing = 0)
DF_All<-dplyr::filter(DF_All,!(INC_ASCVD==2))
y_end<-as.Date("2021-08-30",format = "%Y-%m-%d")
DF_All$INC_ASCVD<-if_else(DF_All$DATE_ALL_MACE>y_end&DF_All$INC_ASCVD==1,0,DF_All$INC_ASCVD)

DF_All<-DF_All%>%
  dplyr::mutate(DATE_end=do.call(pmin,c(dplyr::select(.,starts_with(c('DATE_All',"DATE_D"))), na.rm = TRUE)))
DF_All$DATE_end<-as.Date(if_else(DF_All$DATE_end>y_end|is.na(DF_All$DATE_end),as.character(y_end),as.character(DF_All$DATE_end)))
DF_All$Time_END<-difftime(DF_All$DATE_end,DF_All$REC_DATE)


## Categories BMI -----------------------------------------------------------------------------------------------------
DF_All<-DF_All%>%filter(BMI>=18.5)

DF_All$BMI_Cat<-if_else(DF_All$BMI>=18.5&DF_All$BMI<25,1,0)
DF_All$BMI_Cat<-if_else(DF_All$BMI>=25&DF_All$BMI<30,2,DF_All$BMI_Cat)
DF_All$BMI_Cat<-if_else(DF_All$BMI>=30&DF_All$BMI<35,3,DF_All$BMI_Cat)
DF_All$BMI_Cat<-if_else(DF_All$BMI>=35,4,DF_All$BMI_Cat)


## Group CVHS ---------------------------------------------------------------------------------------------------------
DF_All$CVHS<-as.factor(str_replace(DF_All$CVHS,"10","9"))

DF_All$CVHS_C<-as.character(DF_All$CVHS)
DF_All$CVHS_C<-str_replace_all(DF_All$CVHS_C,c("0|1|2"),"0")
DF_All$CVHS_C<-str_replace_all(DF_All$CVHS_C,c("3|4"),"1")
DF_All$CVHS_C<-str_replace_all(DF_All$CVHS_C,c("5|6"),"2")
DF_All$CVHS_C<-str_replace_all(DF_All$CVHS_C,c("7|8|9"),"3")
DF_All$CVHS_C<-as.factor(DF_All$CVHS_C)

## Group ref BMI/CVHS (Figure 2) --------------------------------------------------------------------------------------
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==2&DF_All$CVHS_C==0,1,0)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==3&DF_All$CVHS_C==0,2,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==4&DF_All$CVHS_C==0,3,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==1&DF_All$CVHS_C==1,4,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==2&DF_All$CVHS_C==1,5,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==3&DF_All$CVHS_C==1,6,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==4&DF_All$CVHS_C==1,7,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==1&DF_All$CVHS_C==2,8,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==2&DF_All$CVHS_C==2,9,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==3&DF_All$CVHS_C==2,10,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==4&DF_All$CVHS_C==2,11,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==1&DF_All$CVHS_C==3,12,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==2&DF_All$CVHS_C==3,13,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==3&DF_All$CVHS_C==3,14,DF_All$Ref_BMI)
DF_All$Ref_BMI<-if_else(DF_All$BMI_Cat==4&DF_All$CVHS_C==3,15,DF_All$Ref_BMI)
DF_All$Ref_BMI<-as.factor(DF_All$Ref_BMI)

## Sex-Specific WHR quartiles -----------------------------------------------------------------------------------------
DF_WHR_F<-DF_All%>%dplyr::select("IID","SEX","WHR")%>%filter(SEX==0)
q_WHR_F<-quantile(DF_WHR_F$WHR, probs = seq(0,1,1/4))
DF_WHR_F$WHR_Q4<-with(DF_WHR_F,
                      cut(WHR,
                          q_WHR_F,
                          include.lowest = T, labels = c(1:4)))

DF_WHR_M<-DF_All%>%dplyr::select("IID","SEX","WHR")%>%filter(SEX==1)
q_WHR_M<-quantile(DF_WHR_M$WHR, probs = seq(0,1,1/4))
DF_WHR_M$WHR_Q4<-with(DF_WHR_M,
                      cut(WHR,
                          q_WHR_M,
                          include.lowest = T, labels = c(1:4)))

DF_WHR_Q4<-merge(x = DF_WHR_F, y = DF_WHR_M, all = TRUE)
DF_WHR_Q4$SEX<-ifelse(is.na(DF_WHR_Q4$SEX.x),DF_WHR_Q4$SEX.y,DF_WHR_Q4$SEX.x)
DF_WHR_Q4$WHR<-ifelse(is.na(DF_WHR_Q4$WHR.x),DF_WHR_Q4$WHR.y,DF_WHR_Q4$WHR.x)
DF_WHR_Q4$WHR_Q4<-ifelse(is.na(DF_WHR_Q4$WHR_Q4.x),DF_WHR_Q4$WHR_Q4.y,DF_WHR_Q4$WHR_Q4.x)
DF_WHR_Q4<-DF_WHR_Q4%>%dplyr::select("IID","WHR","WHR_Q4")

DF_All<-merge(x = DF_All, y = DF_WHR_Q4,by = c("IID","WHR"))

## Group ref BMI/WHR (Figure 3) ---------------------------------------------------------------------------------------
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==1&DF_All$WHR_Q4==2,1,0)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==1&DF_All$WHR_Q4==3,2,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==1&DF_All$WHR_Q4==4,3,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==2&DF_All$WHR_Q4==1,4,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==2&DF_All$WHR_Q4==2,5,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==2&DF_All$WHR_Q4==3,6,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==2&DF_All$WHR_Q4==4,7,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==3&DF_All$WHR_Q4==1,8,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==3&DF_All$WHR_Q4==2,9,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==3&DF_All$WHR_Q4==3,10,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==3&DF_All$WHR_Q4==4,11,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==4&DF_All$WHR_Q4==1,12,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==4&DF_All$WHR_Q4==2,13,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==4&DF_All$WHR_Q4==3,14,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-if_else(DF_All$BMI_Cat==4&DF_All$WHR_Q4==4,15,DF_All$Ref_BMI_cat_WHR)
DF_All$Ref_BMI_cat_WHR<-as.factor(DF_All$Ref_BMI_cat_WHR)

## Group ref CVHS/WHR (Figure 4) --------------------------------------------------------------------------------------
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==2&DF_All$CVHS_C==0,1,0)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==3&DF_All$CVHS_C==0,2,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==4&DF_All$CVHS_C==0,3,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==1&DF_All$CVHS_C==1,4,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==2&DF_All$CVHS_C==1,5,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==3&DF_All$CVHS_C==1,6,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==4&DF_All$CVHS_C==1,7,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==1&DF_All$CVHS_C==2,8,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==2&DF_All$CVHS_C==2,9,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==3&DF_All$CVHS_C==2,10,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==4&DF_All$CVHS_C==2,11,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==1&DF_All$CVHS_C==3,12,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==2&DF_All$CVHS_C==3,13,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==3&DF_All$CVHS_C==3,14,DF_All$Ref_WHR)
DF_All$Ref_WHR<-if_else(DF_All$WHR_Q4==4&DF_All$CVHS_C==3,15,DF_All$Ref_WHR)
DF_All$Ref_WHR<-as.factor(DF_All$Ref_WHR)


DF_All_W<-DF_All%>%filter(SEX==0)
DF_All_M<-DF_All%>%filter(SEX==1)
# =====================================================================================================================
## Figure 1 (A) : All participants ------------------------------------------------------------------------------------
km_CVHS_M<-survfit(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS, data=DF_All)
plot_F1<-ggsurvplot(km_CVHS_M,ggtheme = theme_classic(),ylim = c(0.75,1),censor.size = 0.2,size = 0.5,legend = "right",
                    legend.title = "CVHS",legend.labs = c("10","9","8","7","6","5","4","3","2","0-1"), xlab = "Follow-up (Years)", ylab = "Event-free survival",font.legend = 12)
plot_F1$plot<-plot_F1$plot + 
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,3))+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10))

# HR Figure 1
Cox_CVHS<-coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+ETH+Townsend, data=DF_All)
Cox_CVHS%>%
  gtsummary::tbl_regression(exp = TRUE)
Cox_summary<-summary(Cox_CVHS)
HR<-Cox_summary$conf.int
colnames(HR)<-c("HR","inv","min","max")
HR<-HR[1:9,c("HR","min","max")]
ref<-c(1.00,1.00,1.00)
HR<-rbind(ref,HR)
name<-c("CVHS = 10","CVHS = 9","CVHS = 8","CVHS = 7","CVHS = 6","CVHS = 5","CVHS = 4","CVHS = 3","CVHS = 2","CVHS = 0-1")

HR<-cbind(name,HR)
HR<-as.data.frame(HR)
HR$HR<-as.numeric(HR$HR)
HR$HR<-round(HR$HR, digits = 2)
HR$min<-as.numeric(HR$min)
HR$min<-round(HR$min, digits = 2)
HR$max<-as.numeric(HR$max)
HR$max<-round(HR$max, digits = 2)

# Event rate Figure 1 
N_event_m<-as.data.frame.matrix(table(DF_All$CVHS,DF_All$INC_ASCVD))
colnames(N_event_m)<-c("No_event","Event")
N_event_m$N<-N_event_m$No_event+N_event_m$Event
N_event_m$ER<-(N_event_m$Event/N_event_m$N)*100
N_event_m$ER<-round(N_event_m$ER, digits = 2)

# Figure 1 (ER + HR)
Figure_1<-cbind(N_event_m,HR)
Figure_1$v1<-c("1.00 (ref)", paste0(sprintf("%.2f",Figure_1$HR[2:10])," (", sprintf("%.2f",Figure_1$min[2:10])," - ",sprintf("%.2f",Figure_1$max[2:10]),")"))
rownames(Figure_1)<-Figure_1$name
Figure_1$CVHS<-rev(1:10)
Figure_1[10,10]<-"0-1"

Figure_1<-Figure_1%>%dplyr::select("CVHS","N","Event","ER","v1")
colnames(Figure_1)<-c("CVHS","N total","N cases","Event rate (%)","HR (95% CI)")

table_F1<-ggtexttable(Figure_1,theme = ttheme("light"),rows = NULL)

F1_all<-ggarrange(plotlist=list(plot_F1$plot,table_F1),ncol = 1,nrow = 2,heights = c(2,1),labels = "A")

png(filename = "/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Figure_1_A.png",type = "cairo",res = 600,units = "cm", width = 28,height = 27)
F1_all
dev.off()
# =====================================================================================================================

# =====================================================================================================================
## Figure 1 (B-C) -----------------------------------------------------------------------------------------------------
Sex_loop<-c("DF_All_W","DF_All_M")
DF_loop<-data.frame(X=Sex_loop,Y=c("B","C"))

for(s in Sex_loop){
  DF<-get(s)
  x<-DF_loop[which((DF_loop[,c("X")])==s),c("Y")]
  km_CVHS_M<-survfit(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS, data=DF)
  plot_F1<-ggsurvplot(km_CVHS_M,ggtheme = theme_classic(),ylim = c(0.75,1),censor.size = 0.2,size = 0.5,legend = "right",
                      legend.title = "CVHS",legend.labs = c("10","9","8","7","6","5","4","3","2","0-1"), xlab = "Follow-up (Years)", ylab = "Event-free survival",font.legend = 12)
  plot_F1$plot<-plot_F1$plot + 
    scale_x_continuous(limits = c(0,15), breaks = seq(0,15,3))+
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  # HR Figure 1 
  Cox_CVHS<-coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+ETH+Townsend, data=DF)
  Cox_CVHS%>%
    gtsummary::tbl_regression(exp = TRUE)
  Cox_summary<-summary(Cox_CVHS)
  HR<-Cox_summary$conf.int
  colnames(HR)<-c("HR","inv","min","max")
  HR<-HR[1:9,c("HR","min","max")]
  ref<-c(1.00,1.00,1.00)
  HR<-rbind(ref,HR)
  name<-c("CVHS = 10","CVHS = 9","CVHS = 8","CVHS = 7","CVHS = 6","CVHS = 5","CVHS = 4","CVHS = 3","CVHS = 2","CVHS = 0-1")
  
  HR<-cbind(name,HR)
  HR<-as.data.frame(HR)
  HR$HR<-as.numeric(HR$HR)
  HR$HR<-round(HR$HR, digits = 2)
  HR$min<-as.numeric(HR$min)
  HR$min<-round(HR$min, digits = 2)
  HR$max<-as.numeric(HR$max)
  HR$max<-round(HR$max, digits = 2)
  
  # Event rate Figure 1 
  N_event_m<-as.data.frame.matrix(table(DF$CVHS,DF$INC_ASCVD))
  colnames(N_event_m)<-c("No_event","Event")
  N_event_m$N<-N_event_m$No_event+N_event_m$Event
  N_event_m$ER<-(N_event_m$Event/N_event_m$N)*100
  N_event_m$ER<-round(N_event_m$ER, digits = 2)
  
  # Figure 1 (ER + HR)
  Figure_1<-cbind(N_event_m,HR)
  Figure_1$v1<-c("1.00 (ref)", paste0(sprintf("%.2f",Figure_1$HR[2:10])," (", sprintf("%.2f",Figure_1$min[2:10])," - ",sprintf("%.2f",Figure_1$max[2:10]),")"))
  rownames(Figure_1)<-Figure_1$name
  Figure_1$CVHS<-rev(1:10)
  Figure_1[10,10]<-"0-1"
  
  Figure_1<-Figure_1%>%dplyr::select("CVHS","N","Event","ER","v1")
  colnames(Figure_1)<-c("CVHS","N total","N cases","Event rate (%)","HR (95% CI)")
  
  table_F1<-ggtexttable(Figure_1,theme = ttheme("light"),rows = NULL)
  
  F1_all<-ggarrange(plotlist=list(plot_F1$plot,table_F1),ncol = 1,nrow = 2,heights = c(2,1),labels = paste0(x))
  
  png(filename = paste0("/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Figure_1_",x,".png"),type = "cairo",res = 600,units = "cm", width = 28,height = 27)
  print(F1_all)
  dev.off()
}
# =====================================================================================================================

# =====================================================================================================================
## Figure 2 (A-B-C) ---------------------------------------------------------------------------------------------------
Sex_loop<-c("DF_All","DF_All_W","DF_All_M")
DF_loop<-data.frame(X=Sex_loop,Y=c("A","B","C"))

for(s in Sex_loop){
  DF<-get(s)
  x<-DF_loop[which((DF_loop[,c("X")])==s),c("Y")]
  message("** Creating Figure ", s)
  message("** Cox ",s)
  Cox<-coxph(Surv(Time_END/365.25, INC_ASCVD) ~ Ref_BMI+ETH+AGE+SEX+Townsend, data=DF)
  
  message("HR ",s)
  
  Cox_summary<-summary(Cox)
  DF_HR<-Cox_summary$conf.int
  colnames(DF_HR)<-c("HR","inv","min","max")
  DF_HR<-DF_HR[1:15,c("HR","min","max")]
  ref<-c(1.00,1.00,1.00)
  DF_HR<-rbind(ref,DF_HR)
  name<-c(0:15)
  name<-paste0("ref ",name)
  DF_HR<-cbind(name,DF_HR)
  DF_HR<-as.data.frame(DF_HR)
  DF_HR$HR<-as.numeric(DF_HR$HR)
  DF_HR$HR<-signif(DF_HR$HR, digits = 3)
  DF_HR$min<-as.numeric(DF_HR$min)
  DF_HR$min<-signif(DF_HR$min, digits = 3)
  DF_HR$max<-as.numeric(DF_HR$max)
  DF_HR$max<-signif(DF_HR$max, digits = 3)
  DF_HR$Ref<-rep(c("BMI 18.5-24.9","BMI 25-29.9","BMI 30-34.9",paste0("BMI",intToUtf8(8805),"35")),4)
  levels(DF_HR$Ref)<-c("BMI 18.5-24.9","BMI 25-29.9","BMI 30-34.9",paste0("BMI",intToUtf8(8805),"35"))
  DF_HR$CVHS<-factor(c(rep("CVHS = 8-10\n(Healthy)",4),rep("CVHS = 6-7\n(Moderately Healthy)",4),rep("CVHS = 4-5\n(Moderately Unhealthy)",4),rep("CVHS = 0-3\n(Unhealthy)",4)),
                       ordered = TRUE,levels = c("CVHS = 8-10\n(Healthy)","CVHS = 6-7\n(Moderately Healthy)","CVHS = 4-5\n(Moderately Unhealthy)","CVHS = 0-3\n(Unhealthy)"))
  
  pval<-Cox_summary$coefficients[1:15,5]
  DF_HR$pval<-c(NA,sprintf("%.2e",pval))
  
  message("**Event Rate", s)
  N_event<-as.data.frame.matrix(table(DF[,Ref_BMI],DF[,INC_ASCVD]))
  colnames(N_event)<-c("No_event","Event")
  N_event$N<-N_event$No_event+N_event$Event
  N_event$ER<-sprintf("%.2f",(N_event$Event/N_event$N)*100)
  
  ER_Table<-cbind(N_event,DF_HR)
  
  ER_Table$v1<-c("1.00 (ref)", paste0(sprintf("%.2f",ER_Table$HR[2:16])," (", sprintf("%.2f",ER_Table$min[2:16])," - ",sprintf("%.2f",ER_Table$max[2:16]),")"))
  ER_Table<-ER_Table%>%dplyr::select("N","Event","ER","v1","HR","min","max","CVHS","Ref","pval")
  colnames(ER_Table)<-c("N total","N cases","Event rate (%)","HR (95% CI)","mean","lower","upper","CVHS","Ref","P value")
  

  message("**Info Forestplot", s)
  CVHS_BMI_table<-structure(list(mean = c(NA,ER_Table[1:4,5],ER_Table[5:8,5],ER_Table[9:12,5],ER_Table[13:16,5]),
                                   lower = c(NA,ER_Table[1:4,6],ER_Table[5:8,6],ER_Table[9:12,6],ER_Table[13:16,6]),
                                   upper = c(NA,ER_Table[1:4,7],ER_Table[5:8,7],ER_Table[9:12,7],ER_Table[13:16,7])),
                              .Names = c("mean","lower","upper"),
                              row.names = c(NA,-17L),
                              class = "data.frame")
  
  tabletext_BMI<-cbind(c(NA,NA,"CVHS = 8-10\n(Healthy)",NA,NA,NA,"CVHS = 6-7\n(Moderately Healthy)",NA,NA,NA,"CVHS = 4-5\n(Moderately Unhealthy)",NA,NA,NA,"CVHS = 0-3\n(Unhealthy)",NA,NA),
                       c("BMI (kg/m²)",rep(c("BMI 18.5-24.9","BMI 25-29.9","BMI 30-34.9",paste0("BMI",intToUtf8(8805),"35")),4)),
                       c("Cases (n)",ER_Table[1:4,2],ER_Table[5:8,2],ER_Table[9:12,2],ER_Table[13:16,2]),
                       c("Total (n)",ER_Table[1:4,1],ER_Table[5:8,1],ER_Table[9:12,1],ER_Table[13:16,1]),
                       c("Event rate (%)",ER_Table[1:4,3],ER_Table[5:8,3],ER_Table[9:12,3],ER_Table[13:16,3]),
                       c("HR (95% CI)",ER_Table[1:4,4],ER_Table[5:8,4],ER_Table[9:12,4],ER_Table[13:16,4]),
                       c("P value",ER_Table[1:4,10],ER_Table[5:8,10],ER_Table[9:12,10],ER_Table[13:16,10]))
  
  message("** Creating Forest plot",s)
  
  png(filename = paste0("/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Figure_2_",x,".png"),type = "cairo",res = 600,units = "cm", width = 37,height = 15)
  plot.new()
  fn <- local({
    i = 0
    no_lines <- sum(!is.na(CVHS_BMI_table$mean))
    b_clrs = rep(c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),4)
    l_clrs = rep(c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),4)
    b_shape = rep(c(21,22,23,24),4)
    
    function(..., clr.line, clr.marker,size,pch){
      i <<- i + 1
      fpDrawPointCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i],size = unit(0.5,"snpc"), pch = b_shape[i])
    }
  })
  print(CVHS_BMI_table%>%forestplot(labeltext = tabletext_BMI,
                                      align = c("l","c","c","c","c","c","c"),
                                      fn.ci_norm = fn,
                                      is.summary = c(rep(TRUE,1),rep(FALSE,17)),
                                      clip = c(0.5,8),
                                      txt_gp = fpTxtGp(label = gpar(fontsize = 10),ticks = gpar(cex = 0.8)),
                                      graphwidth = unit(150,"mm"),
                                      colgap = unit(0.02,"npc"),
                                      vertices = TRUE,
                                      grid = structure(c(1),gp = gpar(lty = 2, col = "grey")),
                                      hrzl_lines = list("2" = gpar(lty = 1, columns = c(1:7)),
                                                        "6" = gpar(lty = 2, columns = c(1:7)),
                                                        "10" = gpar(lty = 2, columns = c(1:7)),
                                                        "14" = gpar(lty = 2, columns = c(1:7))),
                                      mar = unit(c(20,15,10,15),"mm")))
  legend("topright",legend = c("BMI 18.5-24.9","BMI 25-29.9","BMI 30-34.9",paste0("BMI",intToUtf8(8805),"35")), title = "BMI (kg/m²)",pch = c(21,22,23,24),
         col = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),pt.bg = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),pt.cex = 1.3,
         box.lwd = 1.5,box.lty = 1,box.col = "black")
  mtext(paste0(x),side = 3, adj = 0,line = 2, font = 2,cex = 1.5,at = -0.08)
  dev.off()
  message("Forestplot save")
}


# =====================================================================================================================
## Figure 3 (A-B-C) ---------------------------------------------------------------------------------------------------
Sex_loop<-c("DF_All","DF_All_W","DF_All_M")
DF_loop<-data.frame(X=Sex_loop,Y=c("A","B","C"))

for(s in Sex_loop){
  DF<-get(s)
  x<-DF_loop[which((DF_loop[,c("X")])==s),c("Y")]
  message("** Creating Figure ", s)
  message("** Cox ",s)
  Cox<-coxph(Surv(Time_END/365.25, INC_ASCVD) ~ Ref_BMI_cat_WHR+ETH+AGE+SEX+Townsend, data=DF)
  
  message("HR ",s)
  
  Cox_summary<-summary(Cox)
  DF_HR<-Cox_summary$conf.int
  colnames(DF_HR)<-c("HR","inv","min","max")
  DF_HR<-DF_HR[1:15,c("HR","min","max")]
  ref<-c(1.00,1.00,1.00)
  DF_HR<-rbind(ref,DF_HR)
  name<-c(0:15)
  name<-paste0("ref ",name)
  DF_HR<-cbind(name,DF_HR)
  DF_HR<-as.data.frame(DF_HR)
  DF_HR$HR<-as.numeric(DF_HR$HR)
  DF_HR$HR<-signif(DF_HR$HR, digits = 3)
  DF_HR$min<-as.numeric(DF_HR$min)
  DF_HR$min<-signif(DF_HR$min, digits = 3)
  DF_HR$max<-as.numeric(DF_HR$max)
  DF_HR$max<-signif(DF_HR$max, digits = 3)
  DF_HR$Ref<-rep(c("Q1","Q2","Q3","Q4"),4)
  levels(DF_HR$Ref)<-c("Q1","Q2","Q3","Q4")
  DF_HR$BMI<-factor(c(rep("BMI 18.5-24.9",4),rep("BMI 25-29.9",4),rep("BMI 30-34.9",4),rep(paste0("BMI",intToUtf8(8805),"35"),4)),
                    ordered = TRUE,levels = c("BMI 18.5-24.9","BMI 25-29.9","BMI 30-34.9",paste0("BMI",intToUtf8(8805),"35")))
  
  pval<-Cox_summary$coefficients[1:15,5]
  DF_HR$pval<-c(NA,sprintf("%.2e",pval))
  
  message("**Event Rate", s)
  N_event<-as.data.frame.matrix(table(DF[,Ref_BMI_cat_WHR],DF[,INC_ASCVD]))
  colnames(N_event)<-c("No_event","Event")
  N_event$N<-N_event$No_event+N_event$Event
  N_event$ER<-sprintf("%.2f",(N_event$Event/N_event$N)*100)
  
  ER_Table<-cbind(N_event,DF_HR)
  
  ER_Table$v1<-c("1.00 (ref)", paste0(sprintf("%.2f",ER_Table$HR[2:16])," (", sprintf("%.2f",ER_Table$min[2:16])," - ",sprintf("%.2f",ER_Table$max[2:16]),")"))
  ER_Table<-ER_Table%>%dplyr::select("N","Event","ER","v1","HR","min","max","BMI","Ref","pval")
  colnames(ER_Table)<-c("N total","N cases","Event rate (%)","HR (95% CI)","mean","lower","upper","BMI","Ref","P value")
  

  message("**Info Forestplot", s)
  CVHS_BMI_table<-structure(list(mean = c(NA,ER_Table[1:4,5],ER_Table[5:8,5],ER_Table[9:12,5],ER_Table[13:16,5]),
                                   lower = c(NA,ER_Table[1:4,6],ER_Table[5:8,6],ER_Table[9:12,6],ER_Table[13:16,6]),
                                   upper = c(NA,ER_Table[1:4,7],ER_Table[5:8,7],ER_Table[9:12,7],ER_Table[13:16,7])),
                              .Names = c("mean","lower","upper"),
                              row.names = c(NA,-17L),
                              class = "data.frame")
  
  tabletext_BMI<-cbind(c(NA,NA,"BMI 18.5-24.9",NA,NA,NA,"BMI 25-29.9",NA,NA,NA,"BMI 30-34.9",NA,NA,NA,paste0("BMI",intToUtf8(8805),"35"),NA,NA),
                       c("WHR Quartile",rep(c("Q1","Q2","Q3","Q4"),4)),
                       c("Cases (n)",ER_Table[1:4,2],ER_Table[5:8,2],ER_Table[9:12,2],ER_Table[13:16,2]),
                       c("Total (n)",ER_Table[1:4,1],ER_Table[5:8,1],ER_Table[9:12,1],ER_Table[13:16,1]),
                       c("Event rate (%)",ER_Table[1:4,3],ER_Table[5:8,3],ER_Table[9:12,3],ER_Table[13:16,3]),
                       c("HR (95% CI)",ER_Table[1:4,4],ER_Table[5:8,4],ER_Table[9:12,4],ER_Table[13:16,4]),
                       c("P value",ER_Table[1:4,10],ER_Table[5:8,10],ER_Table[9:12,10],ER_Table[13:16,10]))

  
  message("** Creating Forest plot",s)
  
  png(filename = paste0("/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Figure_3_",x,".png"),type = "cairo",res = 600,units = "cm", width = 37,height = 15)
  plot.new()
  fn <- local({
    i = 0
    no_lines <- sum(!is.na(CVHS_BMI_table$mean))
    b_clrs = rep(c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),4)
    l_clrs = rep(c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),4)
    b_shape = rep(c(21,22,23,24),4)
    
    function(..., clr.line, clr.marker,size,pch){
      i <<- i + 1
      fpDrawPointCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i],size = unit(0.5,"snpc"), pch = b_shape[i])
    }
  })
  print(CVHS_BMI_table%>%forestplot(labeltext = tabletext_BMI,
                                      align = c("l","c","c","c","c","c","c"),
                                      fn.ci_norm = fn,
                                      is.summary = c(rep(TRUE,1),rep(FALSE,17)),
                                      clip = c(0.5,8),
                                      txt_gp = fpTxtGp(label = gpar(fontsize = 10),ticks = gpar(cex = 0.8)),
                                      graphwidth = unit(150,"mm"),
                                      colgap = unit(0.02,"npc"),
                                      vertices = TRUE,
                                      grid = structure(c(1),gp = gpar(lty = 2, col = "grey")),
                                      hrzl_lines = list("2" = gpar(lty = 1, columns = c(1:7)),
                                                        "6" = gpar(lty = 2, columns = c(1:7)),
                                                        "10" = gpar(lty = 2, columns = c(1:7)),
                                                        "14" = gpar(lty = 2, columns = c(1:7))),
                                      mar = unit(c(20,15,10,15),"mm")))
  legend("topright",legend = c("Q1","Q2","Q3","Q4"), title = "WHR Quartile",pch = c(21,22,23,24),
         col = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),pt.bg = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),pt.cex = 1.3,
         box.lwd = 1.5,box.lty = 1,box.col = "black")
  mtext(paste0(x),side = 3, adj = 0,line = 2, font = 2,cex = 1.5,at = -0.08)
  dev.off()
  message("Forestplot save")
}


# =====================================================================================================================
## Figure 4 (A-B-C) ---------------------------------------------------------------------------------------------------
Sex_loop<-c("DF_All","DF_All_W","DF_All_M")
DF_loop<-data.frame(X=Sex_loop,Y=c("A","B","C"))

for(s in Sex_loop){
  DF<-get(s)
  x<-DF_loop[which((DF_loop[,c("X")])==s),c("Y")]
  message("** Creating Figure ", s)
  message("** Cox ",s)
  Cox<-coxph(Surv(Time_END/365.25, INC_ASCVD) ~ Ref_WHR+ETH+AGE+SEX+Townsend, data=DF)
  
  message("HR ",s)
  
  Cox_summary<-summary(Cox)
  DF_HR<-Cox_summary$conf.int
  colnames(DF_HR)<-c("HR","inv","min","max")
  DF_HR<-DF_HR[1:15,c("HR","min","max")]
  ref<-c(1.00,1.00,1.00)
  DF_HR<-rbind(ref,DF_HR)
  name<-c(0:15)
  name<-paste0("ref ",name)
  DF_HR<-cbind(name,DF_HR)
  DF_HR<-as.data.frame(DF_HR)
  DF_HR$HR<-as.numeric(DF_HR$HR)
  DF_HR$HR<-signif(DF_HR$HR, digits = 3)
  DF_HR$min<-as.numeric(DF_HR$min)
  DF_HR$min<-signif(DF_HR$min, digits = 3)
  DF_HR$max<-as.numeric(DF_HR$max)
  DF_HR$max<-signif(DF_HR$max, digits = 3)
  DF_HR$Ref<-rep(c("Q1","Q2","Q3","Q4"),4)
  levels(DF_HR$Ref)<-c("Q1","Q2","Q3","Q4")
  DF_HR$CVHS<-factor(c(rep("CVHS = 8-10\n(Healthy)",4),rep("CVHS = 6-7\n(Moderately Healthy)",4),rep("CVHS = 4-5\n(Moderately Unhealthy)",4),rep("CVHS = 0-3\n(Unhealthy)",4)),
                       ordered = TRUE,levels = c("CVHS = 8-10\n(Healthy)","CVHS = 6-7\n(Moderately Healthy)","CVHS = 4-5\n(Moderately Unhealthy)","CVHS = 0-3\n(Unhealthy)"))
  
  pval<-Cox_summary$coefficients[1:15,5]
  DF_HR$pval<-c(NA,sprintf("%.2e",pval))
  
  message("**Event Rate", s)
  N_event<-as.data.frame.matrix(table(DF[,Ref_WHR],DF[,INC_ASCVD]))
  colnames(N_event)<-c("No_event","Event")
  N_event$N<-N_event$No_event+N_event$Event
  N_event$ER<-sprintf("%.2f",(N_event$Event/N_event$N)*100)
  
  ER_Table<-cbind(N_event,DF_HR)
  
  ER_Table$v1<-c("1.00 (ref)", paste0(sprintf("%.2f",ER_Table$HR[2:16])," (", sprintf("%.2f",ER_Table$min[2:16])," - ",sprintf("%.2f",ER_Table$max[2:16]),")"))
  ER_Table<-ER_Table%>%dplyr::select("N","Event","ER","v1","HR","min","max","CVHS","Ref","pval")
  colnames(ER_Table)<-c("N total","N cases","Event rate (%)","HR (95% CI)","mean","lower","upper","CVHS_2","Ref","P value")
  
  message("**Info Forestplot", s)
  CVHS_2_BMI_table<-structure(list(mean = c(NA,ER_Table[1:4,5],ER_Table[5:8,5],ER_Table[9:12,5],ER_Table[13:16,5]),
                                   lower = c(NA,ER_Table[1:4,6],ER_Table[5:8,6],ER_Table[9:12,6],ER_Table[13:16,6]),
                                   upper = c(NA,ER_Table[1:4,7],ER_Table[5:8,7],ER_Table[9:12,7],ER_Table[13:16,7])),
                              .Names = c("mean","lower","upper"),
                              row.names = c(NA,-17L),
                              class = "data.frame")
  
  tabletext_BMI<-cbind(c(NA,NA,"CVHS = 8-10\n(Healthy)",NA,NA,NA,"CVHS = 6-7\n(Moderately Healthy)",NA,NA,NA,"CVHS = 4-5\n(Moderately Unhealthy)",NA,NA,NA,"CVHS = 0-3\n(Unhealthy)",NA,NA),
                       c("WHR Quartile",rep(c("Q1","Q2","Q3","Q4"),4)),
                       c("Cases (n)",ER_Table[1:4,2],ER_Table[5:8,2],ER_Table[9:12,2],ER_Table[13:16,2]),
                       c("Total (n)",ER_Table[1:4,1],ER_Table[5:8,1],ER_Table[9:12,1],ER_Table[13:16,1]),
                       c("Event rate (%)",ER_Table[1:4,3],ER_Table[5:8,3],ER_Table[9:12,3],ER_Table[13:16,3]),
                       c("HR (95% CI)",ER_Table[1:4,4],ER_Table[5:8,4],ER_Table[9:12,4],ER_Table[13:16,4]),
                       c("P value",ER_Table[1:4,10],ER_Table[5:8,10],ER_Table[9:12,10],ER_Table[13:16,10]))
  
  message("** Creating Forest plot",s)
  
  png(filename = paste0("/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Figure_4_",x,".png"),type = "cairo",res = 600,units = "cm", width = 37,height = 15)
  plot.new()
  fn <- local({
    i = 0
    no_lines <- sum(!is.na(CVHS_2_BMI_table$mean))
    b_clrs = rep(c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),4)
    l_clrs = rep(c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),4)
    b_shape = rep(c(21,22,23,24),4)
    
    function(..., clr.line, clr.marker,size,pch){
      i <<- i + 1
      fpDrawPointCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i],size = unit(0.5,"snpc"), pch = b_shape[i])
    }
  })
  print(CVHS_2_BMI_table%>%forestplot(labeltext = tabletext_BMI,
                                      align = c("l","c","c","c","c","c","c"),
                                      fn.ci_norm = fn,
                                      is.summary = c(rep(TRUE,1),rep(FALSE,17)),
                                      clip = c(0.5,8),
                                      txt_gp = fpTxtGp(label = gpar(fontsize = 10),ticks = gpar(cex = 0.8)),
                                      graphwidth = unit(150,"mm"),
                                      colgap = unit(0.02,"npc"),
                                      vertices = TRUE,
                                      grid = structure(c(1),gp = gpar(lty = 2, col = "grey")),
                                      hrzl_lines = list("2" = gpar(lty = 1, columns = c(1:7)),
                                                        "6" = gpar(lty = 2, columns = c(1:7)),
                                                        "10" = gpar(lty = 2, columns = c(1:7)),
                                                        "14" = gpar(lty = 2, columns = c(1:7))),
                                      mar = unit(c(20,15,10,15),"mm")))
  legend("topright",legend = c("BMI 18.5-24.9","BMI 25-29.9","BMI 30-34.9",paste0("BMI",intToUtf8(8805),"35")), title = "BMI (kg/m²)",pch = c(21,22,23,24),
         col = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),pt.bg = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"),pt.cex = 1.3,
         box.lwd = 1.5,box.lty = 1,box.col = "black")
  mtext(paste0(x),side = 3, adj = 0,line = 2, font = 2,cex = 1.5,at = -0.08)
  dev.off()
  message("Forestplot save")
}
# =====================================================================================================================
## Table 2 : C-stat ---------------------------------------------------------------------------------------------------
Cindex_func <- function(formula, data, indices){
  sample = as.data.frame(data[indices,])
  res.cox <- coxph(formula=as.formula(formula), data = sample, x=T)
  C1 <- as.numeric(as.character(summary(res.cox)$concordance["C"]))
  return(C1)
}

C_index_iter_func = function(data, niter, formula){
  sampling <- boot(data,Cindex_func,R=niter, formula = formula, parallel=c("multicore"), ncpus=25)
  C_index_list = as.data.frame(sampling$t)
  C_index_list = arrange(C_index_list, V1)
  C_index_inf=round(C_index_list[c(25),], digits = 3)
  C_index_sup=round(C_index_list[c(975),], digits = 3)
  result=paste0("(",C_index_inf," to ",C_index_sup,")")
  return(result)
}

c_result<-matrix(data = NA,nrow = 3,ncol = 4)
colnames(c_result)<-c("Model 1 : AGE + SEX","Model 2 : AGE + SEX + CVHS","Model 3 : AGE + SEX + CVHS + BMI","Model 4 : AGE + SEX + CVHS + WHR")
rownames(c_result)<-c("All","Female","Male")
c_result<-as.data.frame(c_result)

DF_boot<-c("DF_All","DF_All_W","DF_All_M")

# Model 1 : AGE + SEX
x<-0
for(DF in DF_boot){
  x<-x+1
  if(DF=="DF_All"){
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ AGE+SEX"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ AGE+SEX, data = DF)
  }else{
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ AGE"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ AGE, data = DF)
  }
  C1 <- round(as.numeric(as.character(summary(res.cox)$concordance["C"])),digits = 3)
  c_result[x,1]<-paste0(C1," ",c_stat)
}

# Model 2 : AGE + SEX + CVHS
x<-0
for(DF in DF_boot){
  x<-x+1
  if(DF=="DF_All"){
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX, data = DF)
  }else{
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE, data = DF) 
  }
  C1 <- round(as.numeric(as.character(summary(res.cox)$concordance["C"])),digits = 3)
  c_result[x,2]<-paste0(C1," ",c_stat)
}

# Model 3 : AGE + SEX + CVHS + BMI
x<-0
for(DF in DF_boot){
  x<-x+1
  if(DF=="DF_All"){
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+BMI"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+BMI, data = DF)
  }else{
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+BMI"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+BMI, data = DF)
  }
  C1 <- round(as.numeric(as.character(summary(res.cox)$concordance["C"])),digits = 3)
  c_result[x,3]<-paste0(C1," ",c_stat)
}

# Model 4 : AGE + SEX + CVHS + WHR
x<-0
for(DF in DF_boot){
  x<-x+1
  if(DF=="DF_All"){
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+WHR"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+WHR, data = DF)
  }else{
    DF<-get(DF)
    c_stat<-C_index_iter_func(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+WHR"))
    res.cox <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+WHR, data = DF)
  }
  C1 <- round(as.numeric(as.character(summary(res.cox)$concordance["C"])),digits = 3)
  c_result[x,4]<-paste0(C1," ",c_stat)
}


Cox_CStat_diff_CI <- function(formula,formula_noprs, data, indices){
  sample = as.data.frame(data[indices,])
  model <- coxph(as.formula(formula), data = sample,x=T)
  model_noprs <- coxph(as.formula(formula_noprs),x=T, data = sample)
  Cdiff <- as.numeric(as.character(summary(model)$concordance["C"])) - as.numeric(as.character(summary(model_noprs)$concordance["C"]))
  return(Cdiff)
}

Cox_CStat_diff_CI_iter = function(data, niter, formula, formula_noprs){
  sampling <- boot(data,Cox_CStat_diff_CI,R=niter, formula=formula, formula_noprs=formula_noprs, parallel=c("multicore"), ncpus=25)
  C2_list <- as.data.frame(sampling$t) %>% arrange(., V1)
  LOWER = round(C2_list[c(25),], digits=3)
  UPPER = round(C2_list[c(975),], digits=3)
  result = paste0(LOWER, " - ", UPPER)
  return(result)
}


c_result_diff<-matrix(data = NA,nrow = 3,ncol = 2)
colnames(c_result_diff)<-c("Model 2 vs Model 3","Model 2 vs Model 4")
rownames(c_result_diff)<-c("All","Female","Male")
c_result_diff<-as.data.frame(c_result_diff)

DF_boot<-c("DF_All","DF_All_W","DF_All_M")

# Model 2 vs Model 3
x<-0
for(DF in DF_boot){
  x<-x+1
  if(DF=="DF_All"){
    DF<-get(DF)
    c_diff<-Cox_CStat_diff_CI_iter(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+BMI"),formula_noprs = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX"))
    model <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+BMI, data = DF)
    model_noprs <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX, data = DF)
  }else{
    DF<-get(DF)
    c_diff<-Cox_CStat_diff_CI_iter(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+BMI"),formula_noprs = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE"))
    model <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+BMI, data = DF)
    model_noprs <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE, data = DF)
  }
  Cdiff <- round(as.numeric(as.character(summary(model)$concordance["C"])) - as.numeric(as.character(summary(model_noprs)$concordance["C"])),digits = 3)
  c_result_diff[x,1]<-paste0(Cdiff," (",c_diff,")")
}


# Model 2 vs Model 4
x<-0
for(DF in DF_boot){
  x<-x+1
  if(DF=="DF_All"){
    DF<-get(DF)
    c_diff<-Cox_CStat_diff_CI_iter(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+WHR"),formula_noprs = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX"))
    model <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+SEX+WHR, data = DF)
    model_noprs <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+SEX+AGE, data = DF)
  }else{
    DF<-get(DF)
    c_diff<-Cox_CStat_diff_CI_iter(data = DF,niter = 999, formula = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+WHR"),formula_noprs = paste0("Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE"))
    model <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE+WHR, data = DF)
    model_noprs <- coxph(Surv(Time_END/365.25, INC_ASCVD) ~ CVHS+AGE, data = DF)
  }
  Cdiff <- round(as.numeric(as.character(summary(model)$concordance["C"])) - as.numeric(as.character(summary(model_noprs)$concordance["C"])),digits = 3)
  c_result_diff[x,2]<-paste0(Cdiff," (",c_diff,")")
}

C_result<-c(c_result,c_result_diff)
C_result<-as.data.frame(C_result)
writexl::write_xlsx(C_result, "/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Table_2.xlsx")
# =====================================================================================================================

# =====================================================================================================================
## Figure 5 -----------------------------------------------------------------------------------------------------------

DF_CVHS<-dplyr::select(DF_All,"IID","AGE","SEX","BMI","CVHS","CVHS_C","WHR_Q4","INC_ASCVD","Time_END","WHR","ETH","Townsend")

DF_CVHS$CVHS_C<-factor(DF_CVHS$CVHS_C, labels = c("Healthy","Mod Healthy","Mod Unhealthy","Unhealthy"))

dd<<-datadist(DF_CVHS)
options(datadist='dd')


# Graph color 
f_color<-scale_color_manual(values = c("#F8766D","#7CAE00","#00BFC4","#E76BF3"))
f_yline_1<-geom_hline(yintercept = 1,linetype="dashed",color = "#4C4C4C")

label_sex<-as_labeller(c("0"= "Females","1"="Males"))
# All participants 
f_BMI_adj<-cph(Surv(Time_END/365.25, INC_ASCVD) ~ rcs(BMI,c(25,30,35))+CVHS_C+SEX+AGE+ETH+Townsend, data = DF_CVHS,x=TRUE, y=TRUE)

f_WHR_adj<-cph(Surv(Time_END/365.25, INC_ASCVD) ~ rcs(WHR,4)+CVHS_C+SEX+AGE+ETH+Townsend, data = DF_CVHS,x=TRUE, y=TRUE)


p_BMI_adj<-ggplot(Predict(f_BMI_adj,BMI,CVHS_C,fun=exp), ylab="Hazard ratio (95% CI)", xlab="BMI (kg/m²)")

p_WHR_adj<-ggplot(Predict(f_WHR_adj,WHR,CVHS_C,fun=exp), ylab="Hazard ratio (95% CI)", xlab="WHR")

plot_A_1<-p_BMI_adj+f_color+f_yline_1+labs(colour = "CVHS")+theme(legend.position = "none")+theme(plot.caption = element_blank())+labs(title = "A : All participants")+scale_y_continuous(breaks = c(0,1,2,3,4))+coord_cartesian(ylim=c(0,4))

plot_A_2<-p_WHR_adj+f_color+f_yline_1+labs(colour = "CVHS")+theme(legend.position = "none")+theme(plot.caption = element_blank())+labs(title = "  ")+coord_cartesian(ylim = c(0,4))

plot_A<-ggarrange(plot_A_1,plot_A_2,ncol = 2)

## Sex-specific 
DF_CVHS_F<-DF_CVHS%>%filter(SEX==0)
DF_CVHS_M<-DF_CVHS%>%filter(SEX==1)

# Females
dd<<-datadist(DF_CVHS_F)
options(datadist='dd')

f_BMI_adj_F<-cph(Surv(Time_END/365.25, INC_ASCVD) ~ rcs(BMI,c(25,30,35))+CVHS_C+AGE+ETH+Townsend, data = DF_CVHS_F,x=TRUE, y=TRUE)

f_WHR_adj_F<-cph(Surv(Time_END/365.25, INC_ASCVD) ~ rcs(WHR,5)+CVHS_C+AGE+ETH+Townsend, data = DF_CVHS_F,x=TRUE, y=TRUE)


p_BMI_F<-ggplot(Predict(f_BMI_adj_F,BMI,CVHS_C,fun=exp), ylab="Hazard ratio (95% CI)", xlab="BMI (kg/m²)")
p_WHR_F<-ggplot(Predict(f_WHR_adj_F,WHR,CVHS_C,fun=exp), ylab="Hazard ratio (95% CI)", xlab="WHR")

plot_B_1<-p_BMI_F+f_color+f_yline_1+labs(colour = "CVHS")+theme(legend.position = "none")+theme(plot.caption = element_blank())+labs(title = "B : Females")+coord_cartesian(ylim = c(0,12))

plot_B_2<-p_WHR_F+f_color+f_yline_1+labs(colour = "CVHS")+theme(legend.position = "none")+theme(plot.caption = element_blank())+labs(title = "  ")+coord_cartesian(ylim = c(0,12))

plot_B<-ggarrange(plot_B_1,plot_B_2,ncol = 2)

# Males
dd<<-datadist(DF_CVHS_M)
options(datadist='dd')

f_BMI_adj_M<-cph(Surv(Time_END/365.25, INC_ASCVD) ~ rcs(BMI,c(25,30,35))+CVHS_C+AGE+ETH+Townsend, data = DF_CVHS_M,x=TRUE, y=TRUE)
f_WHR_adj_M<-cph(Surv(Time_END/365.25, INC_ASCVD) ~ rcs(WHR,4)+CVHS_C+AGE+ETH+Townsend, data = DF_CVHS_M,x=TRUE, y=TRUE)

p_BMI_M<-ggplot(Predict(f_BMI_adj_M,BMI,CVHS_C,fun=exp), ylab="Hazard ratio (95% CI)", xlab="BMI (kg/m²)")
p_WHR_M<-ggplot(Predict(f_WHR_adj_M,WHR,CVHS_C,fun=exp), ylab="Hazard ratio (95% CI)", xlab="WHR")

plot_C_1<-p_BMI_M+f_color+f_yline_1+labs(colour = "CVHS")+theme(plot.caption = element_blank())+labs(title = "C : Males")+coord_cartesian(ylim = c(0,4))

plot_C_2<-p_WHR_M+f_color+f_yline_1+labs(colour = "CVHS")+theme(plot.caption = element_blank())+labs(title = "  ")+coord_cartesian(ylim = c(0,4))

plot_C<-ggarrange(plot_C_1,plot_C_2,ncol = 2, common.legend = TRUE, legend = "bottom")


# Figure 5 
Figure_5<-ggarrange(plot_A,plot_B,plot_C,nrow = 3, common.legend = TRUE, legend = "bottom", heights = c(1,1,1.1))


# Save Figure 5
png(filename = "/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Figure_5.png",type = "cairo",res = 600,units = "cm", width = 20.7,height = 35)
Figure_5
dev.off()
# =====================================================================================================================

# =====================================================================================================================
## Table 1 ------------------------------------------------------------------------------------------------------------
DF_All$N<-1
n<-sum(DF_All$N)
DF_All$INC_ASCVD<-as.numeric(DF_All$INC_ASCVD)
DF_All$Death<-as.numeric(DF_All$Death)
All<-c(n,paste0(sum(DF_All$SEX==1)," (",round(((sum(DF_All$SEX==1)/n)*100),1),"%)"),
       paste0(round(mean(DF_All$AGE),1)," (",round(sd(DF_All$AGE),1),")"),
       paste0(round(mean(DF_All$BMI),1)," (",round(sd(DF_All$BMI),1),")"),
       paste0(sprintf("%.2f",round(mean(DF_All_M$WHR),2))," (",sprintf("%.2f",round(sd(DF_All_M$WHR),2)),")"),
       paste0(sprintf("%.2f",round(mean(DF_All_W$WHR),2))," (",sprintf("%.2f",round(sd(DF_All_W$WHR),2)),")"),
       paste0(sum(DF_All$ETH==1)," (",round(((sum(DF_All$ETH==1)/n)*100),1),"%)"),
       paste0(sum(DF_All$Smoking==1)," (",round(((sum(DF_All$Smoking==1)/n)*100),1),"%)"),
       paste0(sum(DF_All$S_F_L==1)," (",round(((sum(DF_All$S_F_L==1)/n)*100),1),"%)"),
       paste0(sum(DF_All$S_AP==1)," (",round(((sum(DF_All$S_AP==1)/n)*100),1),"%)"),
       paste0(sum(DF_All$Sleep_Score_2==1)," (",round(((sum(DF_All$Sleep_Score_2==1)/n)*100),1),"%)"),
       paste0(round(mean(DF_All$Townsend,na.rm = T),1)," (",sprintf("%.1f",round(sd(DF_All$Townsend,na.rm = T),1)),")"),
       paste0(round(mean(DF_All$SBP),1)," (",round(sd(DF_All$SBP),1),")"),
       paste0(round(mean(DF_All$DBP),1)," (",round(sd(DF_All$DBP),1),")"),
       paste0(sum(DF_All$Rx_BP==1)," (",round(((sum(DF_All$Rx_BP==1)/n)*100),1),")"),
       paste0(round(mean(DF_All$CRP),1)," (",round(sd(DF_All$CRP),1),")"),
       paste0(round(mean(DF_All$TG),1)," (",sprintf("%.1f",round(sd(DF_All$TG),1)),")"),
       paste0(round(mean(DF_All$LDL_C),1)," (",round(sd(DF_All$LDL_C),1),")"),
       paste0(round(mean(DF_All$HDL_C),1)," (",round(sd(DF_All$HDL_C),1),")"),
       paste0(sum(DF_All$Rx_C==1)," (",round(((sum(DF_All$Rx_C==1)/n)*100),1),"%)"),
       paste0(round(mean(DF_All$HbA1C),1)," (",round(sd(DF_All$HbA1C),1),")"),
       paste0(sum(DF_All$Rx_D==1)," (",sprintf("%.1f",round(((sum(DF_All$Rx_D==1)/n)*100),1)),"%)"))

Table_1<-matrix(c(1:10),ncol = 10,nrow=22)
Table_1<-as.data.frame(Table_1)

n_loop<-0:9
for (x in n_loop){
  DF<-DF_All%>%filter(CVHS==x)
  DF_M<-DF%>%filter(SEX==1)
  DF_W<-DF%>%filter(SEX==0)
  n<-sum(DF$N)
  col<-c(n,paste0(sum(DF$SEX==1)," (",round(((sum(DF$SEX==1)/n)*100),1),"%)"),
         paste0(round(mean(DF$AGE),1)," (",round(sd(DF$AGE),1),")"),
         paste0(sprintf("%.1f",round(mean(DF$BMI),1))," (",sprintf("%.1f",round(sd(DF$BMI),1)),")"),
         paste0(sprintf("%.2f",round(mean(DF_M$WHR),2))," (",round(sd(DF_M$WHR),2),")"),
         paste0(sprintf("%.2f",round(mean(DF_W$WHR),2))," (",sprintf("%.2f",round(sd(DF_W$WHR),2)),")"),
         paste0(sum(DF$ETH==1)," (",sprintf("%.1f",round(((sum(DF$ETH==1)/n)*100),1)),"%)"),
         paste0(sum(DF$Smoking==1)," (",sprintf("%.1f",round(((sum(DF$Smoking==1)/n)*100),1)),"%)"),
         paste0(sum(DF$S_F_L==1)," (",round(((sum(DF$S_F_L==1)/n)*100),1),"%)"),
         paste0(sum(DF$S_AP==1)," (",round(((sum(DF$S_AP==1)/n)*100),1),"%)"),
         paste0(sum(DF$Sleep_Score_2==1)," (",round(((sum(DF$Sleep_Score_2==1)/n)*100),1),"%)"),
         paste0(sprintf("%.1f",round(mean(DF$Townsend,na.rm = T),1))," (",sprintf("%.1f",round(sd(DF$Townsend,na.rm = T),1)),")"),
         paste0(sprintf("%.1f",round(mean(DF$SBP),1))," (",sprintf("%.1f",round(sd(DF$SBP),1)),")"),
         paste0(sprintf("%.1f",round(mean(DF$DBP),1))," (",sprintf("%.1f",round(sd(DF$DBP),1)),")"),
         paste0(sum(DF$Rx_BP==1)," (",round(((sum(DF$Rx_BP==1)/n)*100),1),")"),
         paste0(round(mean(DF$CRP),1)," (",round(sd(DF$CRP),1),")"),
         paste0(round(mean(DF$TG),1)," (",round(sd(DF$TG),1),")"),
         paste0(round(mean(DF$LDL_C),1)," (",round(sd(DF$LDL_C),1),")"),
         paste0(sprintf("%.1f",round(mean(DF$HDL_C),1))," (",sprintf("%.1f",round(sd(DF$HDL_C),1)),")"),
         paste0(sum(DF$Rx_C==1)," (",sprintf("%.1f",round(((sum(DF$Rx_C==1)/n)*100),1)),"%)"),
         paste0(sprintf("%.1f",round(mean(DF$HbA1C),1))," (",sprintf("%.1f",round(sd(DF$HbA1C),1)),")"),
         paste0(sum(DF$Rx_D==1)," (",sprintf("%.1f",round(((sum(DF$Rx_D==1)/n)*100),1)),"%)"))
  y<-x+1
  Table_1[,y]<-col
}

rows<-c("N","Male*","Age, years, mean (SD)","BMI, kg/m², mean (SD)","Waist-to-hip ratio (Male), cm, mean (SD)","Waist-to-hip ratio (Female), cm, mean (SD)","White ethnicity*","Current smoking*","Low fruit and vegetable intake*","Low physical activity*","Poor sleep quality*",
        "Deprivation index, mean (SD)","Systolic BP, mm Hg, mean (SD)","Diastolic BP, mm Hg, mean (SD)","Antihypertensive medications*","CRP, mg/L, mean (SD)","Triglycerides, mmol/L, mean (SD)","LDL-C, mmol/L, mean (SD)","HDL-C, mmol/L, mean (SD)","Cholesterol lowering medications*","HbA1c, mmol/mol, mean (SD)","On insulin*")

Table1<-data.frame(rows,All,Table_1)
colnames(Table1)<-c("Factor","Total",paste0("CVHS = ",rev(c(2:10))),"CVHS = 0-1")

writexl::write_xlsx(Table1, "/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Plots/Table_1.xlsx")
# =====================================================================================================================






