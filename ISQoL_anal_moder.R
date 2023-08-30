# Prelude- Clean data and Load packages
#first run Decipher_PEC.R script
close.screen(all=TRUE)
rm(list=ls())
header=1;

library(sl3)
library(SuperLearner)
library(tmle)
library(earth)
library(glmnet)
library(xgboost)
library(dplyr)
library(sjPlot)

# load data, remove missing data
load(file="ISQoL.RData") 

#Not Dummy for table1
# mat_A_W_Table1 <- df_all %>% 
#   dplyr::select(c(Study.Subject.ID,A,
#     Age,Sex,Education,CENTRE,#Employment,Adresseur,
#     RQTH,Fam,Marginalisation,Forensic,
#     Dx,Dx2,Dx_SOMA,Addictions,
#     GAF,CGI,
#     First_Contact,N_Admissions,
#     
#     SQoL18_SEL:SQoL18_TOT,
#     WEMWBS_TOT,SERS_TOT,
#     IS_TOT,IS_Sympt,IS_Disease,IS_Treatment,
#     ISMI_TOT,#ISMI_Alien,ISMI_Approb,ISMI_Discrim,ISMI_Isol,ISMI_Intern,ISMI_Resist,
#     STORI)) %>%
#   #filter(Dx=="SCZ") %>% # & !is.na(A))%>%
#   filter(Dx=="SCZ"  & !is.na(A))%>%
#   as.data.frame(stringsAsFactors = TRUE) 

df_all1 <- filter(df_all,Dx=="SCZ"  & !is.na(A)) %>%
  select(-c(Dx,STORI_REBUILDING,STORI_GROWTH,WEMWBS_TOT,SERS_TOT, starts_with("ISMI"))) %>% #ISMI_Alien:ISMI_Resist,
  mutate(A=case_when(IS_TOT>=9 ~ 1,
                     IS_TOT<9 ~0))

##############################
# Missing data and outliers: evaluation
##############################

#Percentage of missing data per observation
missing_per_row = data.frame(Study.Subject.ID=df_all$Study.Subject.ID,
                             pct_miss=rowSums(is.na(df_all)/ncol(df_all))
)

hist(missing_per_row$pct_miss)

# take out subj who miss a lot of data
who_miss_10=filter(missing_per_row,round(pct_miss,1)>0.2)#
df_all<-anti_join(df_all_SAI_SAII,who_miss_10,by="Study.Subject.ID")

#check outliers
for (i in (3:ncol(df_all))) {
  
 boxplot( df_all[,i])

}

#distribution of outcomes
hist(df_all$SQoL18_TOT)


##############################
# Univariate and Bivariate analysis aka Table 1
##############################

#remove missing values of outcome only
df_all=df_all1 %>%
  filter_at(vars(starts_with("SQoL")), all_vars(!is.na(.)))
df_all=df_all1[complete.cases(df_all1),]

#Table1_data <- mat_A_W_Table1[complete.cases(mat_A_W_Table1),]
Table1_data <- dplyr::select(df_all,-Study.Subject.ID) %>%
   mutate(across(c(Sex_Male:Addictions_Addictions), as.factor)
          )

df_final <- df_final %>%
   mutate(across(c(Sex_Male:Addictions_Addictions), as.factor)
          )

library(table1)
my.render.cont <- function(x) {
    with(stats.default(x), 
         sprintf("%0.2f (%0.1f)", MEAN, SD))
}


pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For numeric variables, perform a standard 2-sample t-test
        p <- t.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits=2, eps=0.0001)))
}

table1(~ .
       | A,
       data=df_final,  
       rowlabelhead = "Variables",
       overall=F, extra.col=list(`P-value`=pvalue)
       #,render.continuous=my.render.cont
       ) -> Table1_final

##############################
# Missing data and outliers: correction
##############################
#median and missingness indicator variable
 task <- make_sl3_Task(
    data=df_all,
    covariates = colnames(dplyr::select(df_all,-c(starts_with("SQoL18_"),
                                               starts_with("ISMI")
                                               )
                                        )
                          ),
    outcome = colnames(dplyr::select(df_all,c(starts_with("SQoL18_"),
                                           starts_with("ISMI")
                                           )
                                     )
                       )
  )
 
### same strategy without sl3
no_SL3=0
if (no_SL3) {
dummy_cols(
  .data=  dplyr::select(mat,-c(starts_with("SQoL"),
                               starts_with("ISMI"))),
  select_columns = NULL,
  remove_first_dummy = TRUE,#to avoid multicol
  remove_most_frequent_dummy = FALSE,
  ignore_na = TRUE,
  split = NULL,
  remove_selected_columns = TRUE#removes the columns used to generate the dummy columns
  ) %>%
  mutate(across(everything(.),
                ~ replace(., is.na(.), median(.))
                )
         )
 #select everything(ends_with("_NA"))
   
   
df_all_NA=select(df_all,ends_with("_NA"))

#function to find the mode
find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

#replace NA by median and mode, then cbind NA_dum 
tmp<-df_all %>%

  group_by() %>%
  mutate(across(where(is.numeric & !starts_with("SQoL")),
                ~case_when(!is.na(.x)~as.numeric(.x),
                           is.na(.x)~as.numeric(median(.,na.rm=TRUE)))
                ),
         across(where(is.factor),
                ~case_when(!is.na(.x)~as.factor(.x),
                           is.na(.x)~factor(find_mode(.))
                           )
                )
  ) #%>% cbind(NA_dum)

}

##############################
## SAI - ATE
##############################
sensitivity=1# if wants to do a complete cases analysis
sensitivity_weight=0

#ML algorithms 
SL.library.A=c("SL.glmnet","SL.earth","SL.ranger")
SL.library.Delta=c("SL.glmnet","SL.earth","SL.ranger")
SL.library.Y=list("SL.glmnet","SL.ranger","SL.earth")

#### for sensitivity analysis use only
if (sensitivity) {
  df_all = df_all %>% select(-starts_with("ISMI"))
  df_final=df_all[complete.cases(df_all),]# if wants to do a complete cases analysis
} else {df_final=task$data}
####
df_final=select(df_final,-c(Study.Subject.ID,CENTRE,IS_TOT))


A=df_final$A

psycho_var=c("WEMWBS","SQoL","SERS","ISMI","IS","A","STORI")
W=select(df_final,-c(starts_with("SQoL"),
                     starts_with("IS"),
                     A
                     )
         )#WEMWBS_TOT,ISMI_TOT,SERS_TOT,STORI))

#SL estimation of p(A|W)
SL_onlyA=SuperLearner(Y = A,X = W,family = binomial(),SL.library = SL.library.A, 
                      method="method.NNLS")

#summary(SL_onlyA$SL.predict/sum(SL_onlyA$SL.predict))
#hist(1/SL_onlyA$SL.predict)

##### for sensitivity analysis (on weights) use only
if (sensitivity_weight) {
threshold=100 #threshold for removing observations based on weights
df_final=df_final[1/SL_onlyA$SL.predict<threshold,]
}
#####

#function to calculate effect of A (insight) on QoL measures based on tmle
run_tmle <- function(Y, g1W) {

SL_onlyY1=SuperLearner(Y = Y,X = cbind(A=1,W),family = gaussian(),SL.library = SL.library.Y, 
                      method="method.NNLS")
SL_onlyY0=SuperLearner(Y = Y,X = cbind(A=0,W),family = gaussian(),SL.library = SL.library.Y, 
                      method="method.NNLS")
Q=cbind(SL_onlyY0$SL.predict,SL_onlyY1$SL.predict)

result <- tmle(Y,A,W,
     # Q.SL.library=SL.library.Y ,
     # g.SL.library=SL.library.A ,
     Q=Q,
     g1W=g1W,
     g.Delta.SL.library = SL.library.Delta,
     Delta=ifelse(is.na(Y),yes = 0, no=1),
     V=10,
     family="gaussian", fluctuation="logistic", prescreenW.g = FALSE)
return(list(Q,result,SL_onlyY0,SL_onlyY1))
}
(results_TOT <- run_tmle(Y=df_final$SQoL18_TOT,
                         g1W=SL_onlyA$SL.predict
                         ))
(AUT <- run_tmle(Y=df_final$SQoL18_AUT,
                 g1W=SL_onlyA$SL.predict
                         ))
(RES <- run_tmle(Y=df_final$SQoL18_RES,
                 g1W=SL_onlyA$SL.predict
                         ))
(FAM <- run_tmle(Y=df_final$SQoL18_FAM,
                 g1W=SL_onlyA$SL.predict
                         ))
(FRI <- run_tmle(Y=df_final$SQoL18_FRI,
                 g1W=SL_onlyA$SL.predict
                         ))
(SEL <- run_tmle(Y=df_final$SQoL18_SEL,
                 g1W=SL_onlyA$SL.predict
                         ))
(PSY <- run_tmle(Y=df_final$SQoL18_PSY,
                 g1W=SL_onlyA$SL.predict
                         ))
(PHY <- run_tmle(Y=df_final$SQoL18_PHY,
                 g1W=SL_onlyA$SL.predict
                         ))
(ROM <- run_tmle(Y=df_final$SQoL18_ROM,
                 g1W=SL_onlyA$SL.predict
                         ))

## SAI - model checks and print results
summary(1/SL_onlyA$SL.predict)
summary(SL_onlyA$SL.predict)

SL_onlyA
cvAUC::AUC(SL_onlyA$SL.predict,A)

all_results=c("SQoL18_TOT","SQoL18_AUT","SQoL18_FAM","SQoL18_FRI","SQoL18_PHY","SQoL18_PSY","SQoL18_RES","SQoL18_ROM","SQoL18_SEL")
all_outcomes=c("results_TOT","AUT","FAM","FRI","PHY","PSY","RES","ROM","SEL")

Table2=NULL

for (i in 1:length(all_results)) {
  SL_onlyY=SuperLearner(Y = df_final[[all_results[i]]],X = cbind(W),family =gaussian(),SL.library = SL.library.Y, method="method.NNLS")
  print(all_results[i])

print(caret::R2(pred=get(all_outcomes[[i]])[1][[1]][,1],df_final[[all_results[i]]]))
print(caret::R2(pred=get(all_outcomes[[i]])[1][[1]][,2],df_final[[all_results[i]]]))
print(get(all_outcomes[[i]])[3][[1]])
print(get(all_outcomes[[i]])[4][[1]])

#print(get(all_outcomes[[i]])[2][[1]])

tmp=get(paste0(all_outcomes[i]))[[2]]$estimates$ATE    
tmp$lb=tmp$CI[1]
tmp$ub=tmp$CI[2]
tmp$CI=NULL
tmp$se.psi=sqrt(tmp$var.psi)

Table2=rbind(Table2,as.data.frame(tmp))

}
Table2[,c(1,6,4,5,3)]
tab_df(Table2[,c(1,6)])
tab_df(Table2[,c(1,3)], digits=4)


###################################
## SAII - MSM
###################################

SL.library.A=c("SL.glmnet","SL.earth","SL.ranger")
SL.library.Delta=c("SL.glmnet","SL.earth","SL.ranger")
SL.library.Y=list("SL.glmnet","SL.earth","SL.ranger")

#### if wants to do a complete cases analysis
if (sensitivity) {
  df_all = df_all %>% select(-starts_with("ISMI"))
df_final=df_all[complete.cases(df_all),]
} else {df_final=task$data}
####
df_final=select(df_final,-c(Study.Subject.ID,CENTRE,IS_TOT))

A=df_final$A

psycho_var=c("WEMWBS","SQoL","SERS","ISMI","IS","A","STORI")
W=select(df_final,-c(starts_with("SQoL"),
                     starts_with("IS"),
                     starts_with("First_"),
                     A
                     )
         )#WEMWBS_TOT,ISMI_TOT,SERS_TOT,STORI))


potential_moderators=W %>%
  select(-c(starts_with("delta"))) %>%
  colnames()
potential_moderators=c("Age","GAF","CGI","Illness_Duration","Sex_Male","Education__12_years")

potential_outcomes=c("SQoL18_SEL", "SQoL18_AUT", "SQoL18_PSY","SQoL18_PHY","SQoL18_TOT")


#function to calculate HTE based on tmle
moder_anal<-function(which_V, Y) {
V<-W %>%
  select(all_of(which_V)) %>%
  rename(V=all_of(which_V)) #%>%
  # mutate(across(starts_with("SQoL18"),~cut(.x, breaks=2,include.lowest=TRUE, ordered_result = TRUE, labels=c(1,0)))) %>%
  # mutate(across(starts_with("SQoL18"),~as.numeric(.)))
  #mutate(which_V=cut(which_V, breaks=2,include.lowest=TRUE, ordered_result = TRUE, labels=c(0,1))) 

W=select(W,-all_of(which_V))

results_MSM <- tmleMSM(Y, A, W, V, 
        MSM="A*V",
        T = rep(1,length(Y)), 
        Delta=ifelse(is.na(Y),yes = 0, no=1),
        v = NULL, Q = NULL, Qform = NULL, Qbounds = c(-Inf, Inf), 
        Q.SL.library = SL.library.Y, g.SL.library = SL.library.A,#same as Delta.SL.library
        cvQinit = FALSE, hAV = NULL, hAVform = NULL, g1W = NULL, 
        gform = NULL, pDelta1 = NULL, g.Deltaform = NULL, 
        ub = 1/0.025, family = "gaussian", fluctuation = "logistic", 
        alpha  = 0.995, id = 1:length(Y), V_SL = 10, inference = TRUE, 
        verbose = TRUE, Q.discreteSL = FALSE, g.discreteSL = FALSE) 
}


for (i in (1:length(potential_outcomes))) {
  
  for (j in (1:length(potential_moderators))) {
      print(potential_outcomes[i])
      print(potential_moderators[j])
      
      assign(paste("results_MSM",potential_outcomes[i],potential_moderators[j],sep="_"),
                   moder_anal(which_V=potential_moderators[j],
                              Y=df_final[[potential_outcomes[i]]]#because data.table
                              )
                   )
      
  }
  
}


## SAII - check positivity and display results

#positivity assumption
all_weights=NULL
Table3=NULL

summary(1/results_MSM_SQoL18_TOT_Age$g$g1W)
summary(results_MSM_SQoL18_TOT_CGI$g$g1W)
summary(1/results_MSM_SQoL18_TOT_CGI$g$g1W)

potential_moderators=c("Age","Sex_Male","Education__12_years","GAF","CGI","Illness_Duration")

potential_outcomes=c("SQoL18_TOT", "SQoL18_AUT","SQoL18_PHY", "SQoL18_PSY","SQoL18_SEL")

#gather variables starting with results_MSM_SQoL18 and display results
for (i in (1:length(potential_outcomes))) {
Table3=NULL
Supp_Table3=NULL

  for (j in (1:length(potential_moderators))) {
      print(potential_outcomes[i])
      print(potential_moderators[j])
      
      print(summary(get(paste0(
        paste("results_MSM",potential_outcomes[i],potential_moderators[j],sep="_"),
        "")
      )
      )
      )

      print(get(paste("results_MSM",potential_outcomes[i],potential_moderators[j],sep="_")))
      
      tmp=get(paste("results_MSM",potential_outcomes[i],potential_moderators[j],sep="_"))
      all_weights=rbind(all_weights,summary(tmp$g$g1W))
      
      T3_tmp=NULL
      T3_tmp$psi=tmp$psi[4]
      T3_tmp$se=tmp$se[4]
      T3_tmp$lb=tmp$lb[4]
      T3_tmp$ub=tmp$ub[4]
      T3_tmp$pvalue=tmp$pvalue[4]
      
      Sup_tmp=data.frame(Parameter=c("Intercept","A","V","A:V"),Mean=tmp$psi,SE=tmp$se,"p-value"=tmp$pvalue,lb=tmp$lb,ub=tmp$ub)
      rownames(Sup_tmp)=NULL
      Sup_tmp[,1]=gsub("A","Insight",Sup_tmp[,1])
      Sup_tmp[,1]=gsub("V",potential_moderators[j],Sup_tmp[,1])


Table3=rbind(Table3,as.data.frame(T3_tmp))
Supp_Table3=rbind(Supp_Table3, as.data.frame(Sup_tmp))

  }
assign(paste0("Table3","_",i),Table3)
assign(paste0("Supp_Table3","_",i),Supp_Table3)

}
  tab_df(Supp_Table3_4, show.rownames=TRUE, digits = 2)

  tab_df(as.data.frame(Supp_Table3_4$`p-value`),col.header = NULL,digits = 4)

  tab_df(Table3_1)