####### Table 1. Clinical characteristics of patients admitted to the hospital with infections. 
########### demographic summary statistics ############
## AGE
quantile(dt$age_admission, probs = seq(0, 1, 0.25))
quantile(subset(dt,risk==0)$age_admission, probs = seq(0, 1, 0.25))
quantile(subset(dt,risk==1)$age_admission, probs = seq(0, 1, 0.25))

## Freqs counts
dt_f = select(dt, GRID,GENDER,risk,sepsis,septic_shock,cardio,respiratory,renal,hepatic,blood,DEATH,
        chf,cpd,cvd,demen,diabwith,diabwithout,hemi,hiv,mi,malig,mild_liver,severe_liver,pud,pvd,
        renal_comor,rheum,solid)
dx = names(dt_f)[-1]
dt_f = dt_f %>% mutate_if(is.numeric, as.character)

## Freqs with dynamtic column names
tb_all = data.frame()
for(i in 1:length(dx)){
  out = dt %>%
    group_by(!!sym(dx[i])) %>%
    summarise(n = n()) %>%
    mutate(Freq = n/sum(n))
    out$dx = dx[i]
    names(out)[1] = "level"
    tb_all = rbind(tb_all,out)
}
tb_all = mutate(tb_all,count=paste(n,"(",round(Freq,3)*100,")",sep=""))
tb_all = select(tb_all,dx,level,count)
write.table(tb_all,"demo_freq_risk.txt",row.names=F,quote=F)

#### Freqs risk=1
dt_f_risk= subset(dt_f, risk==1)
tb_all = data.frame()
for(i in 1:length(dx)){
   out =  dt_f_risk %>%
    group_by(!!sym(dx[i])) %>%
    summarise(n = n()) %>%
    mutate(Freq = n/sum(n))
    out$dx = dx[i]
    names(out)[1] = "level"
    tb_all = rbind(tb_all,out)
}
tb_all = mutate(tb_all,count=paste(n,"(",round(Freq,3)*100,")",sep=""))
tb_all = select(tb_all,dx,level,count)
write.table(tb_all,"demo_freq_risk1.txt",row.names=F,quote=F)

#### Freqs risk=0
dt_f_risk= subset(dt_f, risk==0)
tb_all = data.frame()
for(i in 1:length(dx)){
   out =  dt_f_risk %>%
    group_by(!!sym(dx[i])) %>%
    summarise(n = n()) %>%
    mutate(Freq = n/sum(n))
    out$dx = dx[i]
    names(out)[1] = "level"
    tb_all = rbind(tb_all,out)
}
tb_all = mutate(tb_all,count=paste(n,"(",round(Freq,3)*100,")",sep=""))
tb_all = select(tb_all,dx,level,count)
write.table(tb_all,"demo_freq_risk0.txt",row.names=F,quote=F)

#########  Figure 1. Association between APOL1 high-risk genotypes and the risk of sepsis, 
######################and septic shock, and individual sepsis organ dysfunction criteria.
############ logistic regression: pheno ~ risk + covariates ############# figure 1.
phenotypes <- c("sepsis","septic_shock","renal","respiratory","blood","cardio","hepatic","DEATH")
tb_all = data.frame()
for(i in 1:length(phenotypes)){
  y = phenotypes[i]
  vars = c("risk","GENDER","age_admission","PC1","PC2","PC3")
  f = as.formula(paste(y, paste(vars, collapse=" + "), sep=" ~ "))
  reg = glm(f,family=binomial,data=dt)
  out = summary(reg)
  beta = reg$coefficients
  odds = exp(beta)
  ci = exp(confint.default(reg))
  tb = cbind(out$coefficients,odds, ci)
  colnames(tb) = c("beta","ste","z","p","odds","ci_lower","ci_upper")
  tb = data.frame(pheno=phenotypes[i],tb)
  tb_all = rbind(tb_all,tb)
}
tb_all$var = row.names(tb_all)
tb_all_risk = filter(tb_all, grepl("risk",var))
write.table(tb_all,file="pheno_risk_covariates_details.txt",row.names=F,quote=F)
write.table(tb_all_risk, file="pheno_risk_covariates.txt",row.names=F,quote=F)

######### Figure 2. Associations between APOL1 high-risk genotypes and the risk of sepsis 
############and renal dysfunction before and after adjustment for renal disease, 
############and exclusion of patients with pre-existing severe renal disease.
############## remove kidney failure individuals ################
kidney_failure_ind = read.table("kidney_failure_ind.txt",header=T)
dt_new = subset(dt,! GRID %in% kidney_failure_ind$GRID)

############# 1. APOL1 risk vs sepsis, adjust by age, gender and PCs (done)
reg = glm(sepsis~risk+GENDER+age_admission+PC1+PC2+PC3,family=binomial,data=dt)
    out = summary(reg)
    beta = reg$coefficients
    odds = exp(beta)
    ci = exp(confint.default(reg))
    tb = cbind(out$coefficients,odds, ci)
    colnames(tb) = c("beta","ste","z","p","odds","ci_lower","ci_upper")
    tb = data.frame(pheno="sepsis",tb)
write.table(tb, file="sepsis_risk_covariates.txt",quote=F,row.names=T)

############# renal ~ risk + sever_renal + covariates ###################
reg = glm(renal~risk+sever_renal+GENDER+age_admission+PC1+PC2+PC3,family=binomial,data=dt)
out = summary(reg)
beta = reg$coefficients
odds = exp(beta)
ci = exp(confint.default(reg))
tb = cbind(out$coefficients,odds, ci)
colnames(tb) = c("beta","ste","z","p","odds","ci_lower","ci_upper")
tb = data.frame(pheno="renal",tb)
write.table(tb, file="renal_risk_sever_renal_covariates.txt",quote=F)

############ APOL1 risk vs sepsis, exclude 458 patients with severe renal disease, adj. by age, gendr and PCs.
reg = glm(sepsis~risk+GENDER+age_admission+PC1+PC2+PC3,family=binomial,data=dt_new)
    out = summary(reg)
    beta = reg$coefficients
    odds = exp(beta)
    ci = exp(confint.default(reg))
    tb = cbind(out$coefficients,odds, ci)
    colnames(tb) = c("beta","ste","z","p","odds","ci_lower","ci_upper")
    tb = data.frame(pheno="sepsis",tb)
write.table(tb, file="sepsis_risk_covariates_rmKidneyfailure.txt",quote=F,row.names=T)

############ sepsis ~ risk + sever_renal + covaiates ###################
reg = glm(sepsis~risk+sever_renal+GENDER+age_admission+PC1+PC2+PC3,family=binomial,data=dt)
    out = summary(reg)
    beta = reg$coefficients
    odds = exp(beta)
    ci = exp(confint.default(reg))
    tb = cbind(out$coefficients,odds, ci)
    colnames(tb) = c("beta","ste","z","p","odds","ci_lower","ci_upper")
    tb = data.frame(pheno="sepsis",tb)
write.table(tb, file="sepsis_risk_sever_renal_covariates.txt",quote=F,row.names=T)

############## renal ~ risk + covariates after removing individuals with kidney failure ###################
reg = glm(renal~risk+GENDER+age_admission+PC1+PC2+PC3,family=binomial,data=dt_new)
out = summary(reg)
beta = reg$coefficients
odds = exp(beta)
ci = exp(confint.default(reg))
tb1 = cbind(out$coefficients,odds, ci)
colnames(tb1) = c("beta","ste","z","p","odds","ci_lower","ci_upper")
tb1 = data.frame(pheno="renal",tb1)
write.table(tb1, file="renal_risk_covariates_rmkidneyfailure_short2.txt",quote=F)

table(dt$sepsis)
table(dt$risk)
table(dt$sepsis,dt$risk)

table(dt_new$sepsis)
table(dt_new$risk)
table(dt_new$sepsis,dt_new$risk)

#####4.	Figure 3. Restricted PheWAS associations between APOL1 high-risk genotype and sepsis-related phenotypes.
############ phewas #############

phe_table_b<-createPhewasTable(phecode_03252022_mega_pc_b_apol1, min.code.count = 2, add.exclusions=T, translate=F)

phe_table_b

#14,713 rows 

all<-merge(phe_table_b,mega_demo_pc_b_apol1 )

all

#df [14,713 × 1,912]

outdir="APOL1/”

phenotypes<-names(phe_table_b[2:1876])

predictors<-c("risk")

covariates<-c("GENDER_SOURCE_VALUE","AGE_LAST_VISIT","EHR_LENGTH","PC1","PC2","PC3")

output<-phewas(phenotypes, predictors, data=all, covariates,additive.genotypes = F)

# now lets annotate our phewas

annotated<-addPhecodeInfo(output)

head(annotated) 

tabletitle=("APOL1_aa_phewas_adj_gender_agelastvisit_ehr_pc123.txt")

write.table(output, file=paste(outdir,tabletitle),quote=F,row.names=F, sep="\t")

annotatedtitle=("APOL1_aa_phewas_adj_gender_agelastvisit_ehr_pc123_annotated.txt")

write.table(annotated, file=paste(outdir,annotatedtitle),quote=F,row.names=F, sep="\t")

# now lets try to plot our data

# jpeg(jpgoutdir)

jpgtitle = ("APOL1 in AA")

phewasManhattan(output, annotate.level=1e-3, title=jpgtitle,OR.direction = T)

ggsave(jpgoutdir)

dev.off()

##########phewas after excluding patients with severe renal disease ###############

all_exkfc <- all[!all$ID%in%kfc_aa_grid$MEGA_GRID,]

all_exkfc

#12,547

outdir="APOL1/"

phenotypes<-names(phe_table_b[2:1876])

predictors<-c("risk")

covariates<-c("GENDER_SOURCE_VALUE","AGE_LAST_VISIT","EHR_LENGTH","PC1","PC2","PC3")

output<-phewas(phenotypes, predictors, data=all_exkfc, covariates,additive.genotypes = F)

# now lets annotate our phewas

annotated<-addPhecodeInfo(output)

head(annotated)

tabletitle=("APOL1_aa_exkidneyfailure_phewas_adj_gender_agelastvisit_ehr_pc123.txt")

write.table(output, file=paste(outdir,tabletitle),quote=F,row.names=F, sep="\t")


annotatedtitle=("APOL1_aa_exkidneyfailure_phewas_adj_gender_agelastvisit_ehr_pc123_annotated.txt")

write.table(annotated, file=paste(outdir,annotatedtitle),quote=F,row.names=F, sep="\t")

 

# now lets try to plot our data

# Need find a better way to save plot.

jpgoutdir=("APOL1/apol1_aa_exkidneyfailure.jpg")

jpgoutdir

# jpeg(jpgoutdir)

jpgtitle = ("APOL1 in AA (exclude Kidney Failure)")

phewasManhattan(output, annotate.level=1e-3, title=jpgtitle,OR.direction = T)

ggsave(jpgoutdir)

dev.off()

```
