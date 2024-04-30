###parameter
arg=commandArgs(T)
working_dictory=arg[1]###your working path
span_size_length=arg[2]###the span size length of TD, the default is 1000
breaks_seq=arg[3]###the interval size of TD, the default is seq(-1,6,0.1)
fold_chang=arg[4]###the fold change cutoff of upregulated or downregulated genes 
adj_p=arg[5]###the adjust p value cutoff of upregulated or downregulated genes 

##figure 1
#############recognize TD
options(stringsAsFactors=F)
setwd(working_dictory)
all_cnv<-read.table(file="TCGA-LUAD.masked_cnv.txt",sep="\t",header = T)
sample_name<-unique(all_cnv$sample)
sample_cnv<-list()##每个样本的CNV信息
for(i in 1:length(sample_name)){
sample_cnv[[i]]<-all_cnv[all_cnv$sample==sample_name[i],]
}
names(sample_cnv)<-sample_name

chr_name<-unique(sample_cnv[[1]]$Chrom)###23 chromosomes
##conditions：0.length of amplicon larger than 100b 1. The log rd of amplicon; 2 is greater than 0.3; 3. The difference between two adjacent fragments is less than 0.3;
all_sample_chr_cv<-list()###Record the CNV information of each chromosome of each patient, whether there is TD
for(w in 1:length(sample_cnv)){
sample_chr_cnv<-list()##CNV information for each chromosome of each patient
for(j in 1:length(chr_name)){
sample_chr_cnv[[j]]<-sample_cnv[[w]][sample_cnv[[w]]$Chrom==chr_name[j],]
}
names(sample_chr_cnv)<-chr_name

for(z in 1:length(sample_chr_cnv)){
if(dim(sample_chr_cnv[[z]])[1]>1){
TD_likegroup<-c()
for(i in 1:(dim(sample_chr_cnv[[z]])[1]-1)){
if(all(sample_chr_cnv[[z]]$End[i]-sample_chr_cnv[[z]]$Start[i]>=100,sample_chr_cnv[[z]]$value[i]>0,sample_chr_cnv[[z]]$value[i]-max(sample_chr_cnv[[z]]$value[i-1],sample_chr_cnv[[z]]$value[i+1])>0.3,abs(sample_chr_cnv[[z]]$value[i-1]-sample_chr_cnv[[z]]$value[i+1])<=0.3)){
TD_likegroup[i]<-1
}else{TD_likegroup[i]<-0}
}
if(all(sample_chr_cnv[[z]]$End[dim(sample_chr_cnv[[z]])[1]]-sample_chr_cnv[[z]]$Start[dim(sample_chr_cnv[[z]])[1]]>=100,sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]]>0,sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]]-sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]-1]>0.3,abs(sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]-1])<=0.3)){
TD_likegroup[dim(sample_chr_cnv[[z]])[1]]<-1
}else{TD_likegroup[dim(sample_chr_cnv[[z]])[1]]<-0}
sample_chr_cnv[[z]]$TD_likegroup<-TD_likegroup
}else{sample_chr_cnv[[z]]$TD_likegroup<-0}
}
all_sample_chr_cv[[w]]<-sample_chr_cnv
}
names(all_sample_chr_cv)<-sample_name


#######Calculate the total of TDS for each sample
all_TD_frame<-list()
for(i in 1:length(all_sample_chr_cv)){
TD_frame<-data.frame()
for(j in 1:length(all_sample_chr_cv[[i]])){

TD_frame<-rbind(TD_frame,all_sample_chr_cv[[i]][[j]],stringsAsFactors=F)
}
all_TD_frame[[i]]<-TD_frame
}
names(all_TD_frame)<-sample_name
save(all_TD_frame,file="all_TD_frame.Rdata")####Record each patient's TD segment
TD_num<-c()
for(i in 1:length(all_TD_frame)){
TD_num[i]<-sum(all_TD_frame[[i]]$TD_likegroup)
}
names(TD_num)<-sample_name


####Calculate the mean TD of each chromosome in all samples
chr_mean<-c()
for(i in 1:length(chr_name)){
chr_num<-c()
for(j in 1:length(all_sample_chr_cv)){
chr_num[j]<-sum(all_sample_chr_cv[[j]][[i]]$TD_likegroup)
}
chr_mean[i]<-sum(chr_num)/length(all_sample_chr_cv)
}
names(chr_mean)<-chr_name


##################Calculate the TDP score for each sample
TDP_score<-c()
for(i in 1:length(all_sample_chr_cv)){
chr_abs<-c()
for(j in 1:length(chr_name)){
chr_abs[j]<-abs(sum(all_sample_chr_cv[[i]][[j]]$TD_likegroup)-chr_mean[j])
}
TDP_score[i]<--sum(chr_abs)/TD_num[i]
}
less20_name<-names(TD_num[TD_num<=20])
names(TDP_score)<-sample_name
save(TDP_score,less20_name,file="TDP_score.Rdata")


##########The relationship between TDP score and OS or RFS was simulated and the survival of TDP group was found to be better
setwd("E:/NSCLC_TDP/LUAD/CNV")
options(stringsAsFactors=F)
load("TDP_score.Rdata")
library(mclust)
set.seed(12345678)
TDP_score<-TDP_score[setdiff(names(TDP_score),less20_name)]
re_TDP_score<-Mclust(TDP_score)
group<-re_TDP_score$classification##4个类别
dat<-data.frame(score=TDP_score,group=group)
value_a<-re_TDP_score[["parameters"]][["mean"]][["1"]]
value_b<-re_TDP_score[["parameters"]][["mean"]][["2"]]

TDP_group<-rownames(dat[dat[,1]> value_b,])
NonTDP_group<-rownames(dat[dat[,1]<=value_a,])

hist(TDP_score,breaks = seq(-2,0,0.05),freq = F)
library(mixtools)
mixmdl=normalmixEM(TDP_score,mu=re_TDP_score$parameters$mean,sigma=re_TDP_score$parameters$variance$sigmasq,epsilon=1e-3)
plot(mixmdl,which=2)
abline(v=value_a,lwd=2,col="blue")
abline(v=value_b,lwd=2,col="red")
save(TDP_group,NonTDP_group,file="patient_TDP_group.Rdata")

##########TDP patients were divided according to amplification length
setwd("E:/NSCLC_TDP/LUAD/CNV")
load("patient_TDP_group.Rdata")
load("all_TD_frame.Rdata")
options(stringsAsFactors=F)
TDP_frame<-all_TD_frame[TDP_group]##TDP patient
span_size_<-list()##The length distribution of TDS in each sample
for(i in 1:length(TDP_frame)){
test_vec<-TDP_frame[[i]][TDP_frame[[i]]$TD_likegroup==1,]
span_size<-test_vec$End-test_vec$Start
span_size_[[i]]<-10(span_size/span_size_length)
}
names(span_size_)<-names(TDP_frame)

hist(unlist(span_size_),fre=F,breaks = breaks_seq)

library(mclust)
set.seed(12345678)
all_TD_size<-unlist(span_size_)
span_size_distri<-Mclust(all_TD_size)

library(mixtools)
mixmdl=normalmixEM(all_TD_size,mu=span_size_distri$parameters$mean,sigma=span_size_distri$parameters$variance$sigmasq,epsilon=1e-3)
plot(mixmdl,which=2)


####Range of fragments in a category
class_0<-range(all_TD_size[span_size_distri$classification==1])##-1~0.414
class_1<-range(all_TD_size[span_size_distri$classification==1])##0.416~0.846
class_2<-range(all_TD_size[span_size_distri$classification==2])##1.85~2.71
class_3<-range(all_TD_size[span_size_distri$classification==3])##2.71~3.85

group_list<-list()
for(i in 1:length(span_size_)){
group_list[[i]]<-c(length(which(class_1[1]<span_size_[[i]]&span_size_[[i]]<class_1[2])),length(which(class_2[1]<span_size_[[i]]&span_size_[[i]]<class_2[2])),length(which(class_3[1]<span_size_[[i]]&span_size_[[i]]<class_3[2])))
}
names(group_list)<-names(span_size_)
group_set<-as.data.frame(group_list)
rownames(group_set)<-c("small_span","class_2","large_span")
colnames(group_set)<-names(group_list)

####Determine which group the patient belongs to
sample_sum<-apply(group_set,2,sum)#####The sum of 1,2, and 3 classes per patient
names(sample_sum)<-colnames(group_set)

###A proportional matrix, using 1/4 to determine which group you belong to
new_group_list<-matrix(0,nrow=dim(group_set)[1],ncol=dim(group_set)[2])
for(i in 1:dim(group_set)[2]){
new_group_list[,i]<-group_set[,i]/sample_sum[i]
}
new_group_list<-as.data.frame(new_group_list)
rownames(new_group_list)<-rownames(group_set)
colnames(new_group_list)<-colnames(group_set)

group_1<-colnames(new_group_list[,new_group_list[1,]>=1/4])
group_2<-colnames(new_group_list[,new_group_list[2,]>=1/4])
group_3<-colnames(new_group_list[,new_group_list[3,]>=1/4])
group_1_2<-colnames(new_group_list[,new_group_list[1,]>=1/4&new_group_list[2,]>=1/4])
group_1_3<-colnames(new_group_list[,new_group_list[1,]>=1/4&new_group_list[3,]>=1/4])
group_2_3<-colnames(new_group_list[,new_group_list[2,]>=1/4&new_group_list[3,]>=1/4])

only_1<-setdiff(group_1,unique(union(group_1_2,group_1_3)))
only_2<-setdiff(group_2,unique(union(group_1_2,group_2_3)))
only_3<-setdiff(group_3,unique(union(group_1_3,group_2_3)))
only_1_2<-setdiff(group_1_2,intersect(intersect(group_1,group_2),group_3))
only_1_3<-setdiff(group_1_3,intersect(intersect(group_1,group_2),group_3))
only_2_3<-group_2_3
######It's divided into six different groups
new_group_list[4,only_1]<-1
new_group_list[4,only_2]<-2
new_group_list[4,only_3]<-3
new_group_list[4,only_1_2]<-4
new_group_list[4,only_1_3]<-5
new_group_list[4,only_2_3]<-6
final_group<-as.numeric(new_group_list[4,])
names(final_group)<-colnames(new_group_list)
save(NonTDP_group,only_1,only_2,only_3,only_1_2,only_1_3,only_2_3,file="patient_group.Rdata")
##########The genome complexity was calculated for each group
complexity<-function(TDP_group){
TDP_group_genome<-lapply(all_TD_frame[TDP_group],dim)
num_amp_del<-c()
for(i in 1:length(TDP_group_genome)){
num_amp_del[i]<-TDP_group_genome[[i]][1]
}
names(num_amp_del)<-names(TDP_group_genome)
return(num_amp_del)
}

NonTDP_group_complexity<-complexity(NonTDP_group)
only_1_2_complexity<-complexity(only_1_2)
only_1_3_complexity<-complexity(only_1_3)
only_2_complexity<-complexity(only_2)
only_3_complexity<-complexity(only_3)
only_2_3_complexity<-complexity(only_2_3)
wilcox.test(only_2_complexity,NonTDP_group_complexity,alternative = "greater")

df<-data.frame(c(NonTDP_group_complexity,only_1_2_complexity,only_1_3_complexity,only_2_complexity,
only_3_complexity,only_2_3_complexity),c(rep("NonTDP_group_complexity",length(NonTDP_group_complexity)),
rep("only_1_2_complexity",length(only_1_2_complexity)),rep("only_1_3_complexity",length(only_1_3_complexity)),
rep("only_2_complexity",length(only_2_complexity)),rep("only_3_complexity",length(only_3_complexity)),
rep("only_2_3_complexity",length(only_2_3_complexity))))
colnames(df)<-c("num","group")


library(ggpubr)

library(digest)
p<-ggboxplot(df, "group", "num",color = "group", size=2,palette =c("#00AFBB", "#E7B800", "#FC4E07","#0000CD",'#8B658B','#D02090'),
add = "jitter", shape = "group",ylim =c(0,1000))


##########figure 6
##clinical informations
options(stringsAsFactors=F)
setwd("E:/NSCLC_TDP/LUAD/clin")
LUAD_clin<-read.table(file="TCGA-LUAD.GDC_phenotype.txt",sep="\t",header=T,na.strings=c("","NA"),quote='')
rownames(LUAD_clin)<-LUAD_clin[,1]

new_LUAD_clin<-LUAD_clin[intersect(TDP_group,rownames(LUAD_clin)),]
new_final_group<-final_group[intersect(rownames(new_LUAD_clin),names(final_group))]
new_LUAD_clin<-cbind(new_LUAD_clin,new_final_group)
##Group and TDP group analysis
attach(new_LUAD_clin)
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage i"]<-"Stage I"
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage ia"]<-"Stage I"
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage ib"]<-"Stage I"
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage iia"]<-"Stage II"
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage iib"]<-"Stage II"
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage iiia"]<-"Stage III"
tumor_stage.diagnoses[tumor_stage.diagnoses=="stage iiib"]<-"Stage III"
tumor_stage.diagnoses[tumor_stage.diagnoses=="not reported"]<-NA
barplot(table(new_final_group,tumor_stage.diagnoses),col=2:6)
barplot(table(new_final_group,primary_therapy_outcome_success),col=2:6)
barplot(table(new_final_group,pathoic_M),col=2:6)
boxplot(age_at_initial_pathoic_diagnosis~new_final_group,col = 2:6)


########prognosis
new_LUAD_clin$RFS<-new_LUAD_clin$days_to_new_tumor_event_after_initial_treatment
new_LUAD_clin$RFS_event<-new_LUAD_clin$new_tumor_event_after_initial_treatment
new_LUAD_clin<-new_LUAD_clin[complete.cases(new_LUAD_clin$RFS_event),]
for(i in 1:dim(new_LUAD_clin)[1]){
if(new_LUAD_clin$RFS_event[i]=='NO'){
new_LUAD_clin$RFS[i]<-new_LUAD_clin$days_to_last_follow_up.diagnoses[i]
}
}
new_LUAD_clin$RFS_event[new_LUAD_clin$RFS_event=='NO']<-0
new_LUAD_clin$RFS_event[new_LUAD_clin$RFS_event=='YES']<-1

library(survival)
attach(new_LUAD_clin)
ovsur<-Surv(RFS,as.numeric(RFS_event))                                            #################建立生存对象
sur_stat<-summary(coxph(ovsur~new_LUAD_clin$new_final_group))
srisk<-cbind(as.numeric(RFS_event),as.numeric(RFS),new_LUAD_clin$new_final_group)
rownames(srisk)<-rownames(new_LUAD_clin)
colnames(srisk)<-c("status","days","group")
srisk<-data.frame(srisk)

##srisk<-srisk[srisk$group==5|srisk$group==6,]

risk.fit<-survfit(Surv(as.numeric(srisk$days/30),as.numeric(srisk$status))~srisk$group)
t_results<-survdiff(Surv(as.numeric(srisk$days/30),as.numeric(srisk$status))~srisk$group);
p<-pchisq(t_results$chisq, 1 , lower.tail = F)
print(p)          #########0.01830103
plot(risk.fit,col = c(2,4,5,6),main=paste("Survival analysis based on TDP signature, ","p value=",p,sep=""),xlab = 'Time (months)',ylab = 'Survival Probability')


##########figure 2
options(stringsAsFactors=F)
setwd("E:/NSCLC_TDP/LUAD/CNV")
all_cnv<-read.table(file="TCGA-LUAD.masked_cnv.txt",sep="\t",header = T)
sample_name<-unique(all_cnv$sample)
sample_cnv<-list()##CNV information for each sample
for(i in 1:length(sample_name)){
sample_cnv[[i]]<-all_cnv[all_cnv$sample==sample_name[i],]
}
names(sample_cnv)<-sample_name
chr_name<-unique(sample_cnv[[1]]$Chrom)###23 chromosomes
##conditions：0.length of amplicon larger than 100b 1. The log rd of amplicon; 2 is greater than 0.3; 3. The difference between two adjacent fragments is less than 0.3;
all_sample_chr_cv<-list()
for(w in 1:length(sample_cnv)){
sample_chr_cnv<-list()
for(j in 1:length(chr_name)){
sample_chr_cnv[[j]]<-sample_cnv[[w]][sample_cnv[[w]]$Chrom==chr_name[j],]
}
names(sample_chr_cnv)<-chr_name

for(z in 1:length(sample_chr_cnv)){
if(dim(sample_chr_cnv[[z]])[1]>1){
TD_likegroup<-c()
for(i in 1:(dim(sample_chr_cnv[[z]])[1]-1)){
if(all(sample_chr_cnv[[z]]$End[i]-sample_chr_cnv[[z]]$Start[i]>=100,sample_chr_cnv[[z]]$value[i]>0,sample_chr_cnv[[z]]$value[i]-max(sample_chr_cnv[[z]]$value[i-1],sample_chr_cnv[[z]]$value[i+1])>0.3,abs(sample_chr_cnv[[z]]$value[i-1]-sample_chr_cnv[[z]]$value[i+1])<=0.3)){
TD_likegroup[i]<-1
}else{TD_likegroup[i]<-0}
}
if(all(sample_chr_cnv[[z]]$End[dim(sample_chr_cnv[[z]])[1]]-sample_chr_cnv[[z]]$Start[dim(sample_chr_cnv[[z]])[1]]>=100,sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]]>0,sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]]-sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]-1]>0.3,abs(sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]-1])<=0.3)){
TD_likegroup[dim(sample_chr_cnv[[z]])[1]]<-1
}else{TD_likegroup[dim(sample_chr_cnv[[z]])[1]]<-0}
sample_chr_cnv[[z]]$TD_likegroup<-TD_likegroup
}else{sample_chr_cnv[[z]]$TD_likegroup<-0}
}
all_sample_chr_cv[[w]]<-sample_chr_cnv
}
names(all_sample_chr_cv)<-sample_name


load("E:/NSCLC_TDP/LUAD/CNV/patient_group.Rdata")
group_1<-unique(c(only_1,only_1_2,only_1_3))


#######UCSC selects the TDP region that is effectively amplified
UCSC_paragraph<-function(chr_cnv,GENE,chrname){
invest_cnv<-list()#######Stores the data frames involved in each OG/TSG
for(j in 1:length(chr_cnv)){
start_row<-c()
end_row<-c()
for(i in 1:dim(chr_cnv[[j]][[chrname]])[1]){
if(chr_cnv[[j]][[chrname]][i,]$Start<=GENE[1]&GENE[1]<=chr_cnv[[j]][[chrname]][i,]$End){
start_row<-i
}
if(chr_cnv[[j]][[chrname]][i,]$Start<=GENE[2]&GENE[2]<=chr_cnv[[j]][[chrname]][i,]$End){
end_row<-i
}
}
if(length(start_row)!=0&length(end_row)!=0){
invest_cnv[[j]]<-chr_cnv[[j]][[chrname]][start_row:end_row,]
}
}
names(invest_cnv)<-names(chr_cnv)
final_frame<-data.frame()########Merge the matrices involved in OG/TSG
for(z in 1:length(invest_cnv)){
test_frame<-as.data.frame(invest_cnv[[z]])
final_frame<-rbind(final_frame,test_frame)
}
TDPamp_region<-final_frame[final_frame$TD_likegroup==1,]

for(q in 1:dim(TDPamp_region)[1]){
if(all(TDPamp_region[q,]$Start<=GENE[1],GENE[2]<=TDPamp_region[q,]$End)){
TDPamp_region$effect[q]<-"DUP"
}
if(all(TDPamp_region[q,]$Start>=GENE[1],GENE[2]>=TDPamp_region[q,]$End)){
TDPamp_region$effect[q]<-"DT"
}
if(all(TDPamp_region[q,]$Start>=GENE[1],GENE[2]<=TDPamp_region[q,]$End)){
TDPamp_region$effect[q]<-"ST"
}
if(all(TDPamp_region[q,]$Start<=GENE[1],GENE[2]>=TDPamp_region[q,]$End)){
TDPamp_region$effect[q]<-"ST"
}
}

return(TDPamp_region)
}
######
########Study OG and TSG in only_2
only_2_chr_cnv<-all_sample_chr_cv[only_2]
##OG
MAPK1<-c(22108789,22221970)
MAPK1_chr<-"22"
MAPK1_group_1<-UCSC_paragraph(only_2_chr_cnv,MAPK1,MAPK1_chr)

MLLT11<-c(151030234,151040970)
MLLT11_chr<-"1"
MLLT11_group_1<-UCSC_paragraph(only_2_chr_cnv,MLLT11,MLLT11_chr)


##TSG
TMCO2<-c(40711619,40717363)
TMCO2_chr<-"1"
TMCO2_group_1<-UCSC_paragraph(only_2_chr_cnv,TMCO2,TMCO2_chr)

EXO5<-c(40974413,40982228)
EXO5_chr<-"1"
EXO5_group_1<-UCSC_paragraph(only_2_chr_cnv,EXO5,EXO5_chr)


########The sample heat map influenced by TSG/OG was constructed
library(pheatmap)
hotmap<-function(x){
  pheatmap(x,cluster_row = FALSE,
           color = c("#FFFFFF","#C00000","#64A36F","#E36488","#2980B9"),
           legend_breaks = 0:4,
           legend_labels=c("NO_TD","DUP in OG","ST","DT","DUP in TSG"),
           cutree_rows=dim(x)[1],show_colnames=F,treeheight_row=0, treeheight_col=0,cellheight=20)
}

####group_1
effct_matrix<- matrix(0,nrow = 4,ncol =length(only_2_chr_cnv))
rownames(effct_matrix)<-c("MAPK1","MLLT11","TMCO2","EXO5")
colnames(effct_matrix)<-names(only_2_chr_cnv)
test_group<-list(MAPK1_group_1,MLLT11_group_1,TMCO2_group_1,EXO5_group_1)
names(test_group)<-rownames(effct_matrix)
#OG
for(i in 1:2){
for(j in 1:dim(test_group[[i]])[1]){
if(test_group[[i]][j,]$effect=="DUP"){
effct_matrix[i,test_group[[i]][j,]$sample]<-1
}
if(test_group[[i]][j,]$effect=="ST"){
effct_matrix[i,test_group[[i]][j,]$sample]<-2
}
if(test_group[[i]][j,]$effect=="DT"){
effct_matrix[i,test_group[[i]][j,]$sample]<-3
}
}
}
#TSG
for(i in 3:4){
for(j in 1:dim(test_group[[i]])[1]){
if(test_group[[i]][j,]$effect=="DUP"){
effct_matrix[i,test_group[[i]][j,]$sample]<-4
}
if(test_group[[i]][j,]$effect=="ST"){
effct_matrix[i,test_group[[i]][j,]$sample]<-2
}
if(test_group[[i]][j,]$effect=="DT"){
effct_matrix[i,test_group[[i]][j,]$sample]<-3
}
}
}
setwd("E:/NSCLC_TDP/LUAD/amp_effect")
save(effct_matrix,file="only_2_effctmatrix.Rdata")
phat <- hotmap(effct_matrix)

final_effect<-c()
for(i in 1:dim(effct_matrix)[2]){
if(any(effct_matrix[,i]%in%1)&any(effct_matrix[,i]%in%4)){
final_effect[i]<-5
}
if(any(effct_matrix[,i]%in%1)&!any(effct_matrix[,i]%in%4)){
final_effect[i]<-1
}
if(!any(effct_matrix[,i]%in%1)&any(effct_matrix[,i]%in%4)){
final_effect[i]<-4
}
if(all(!any(effct_matrix[,i]%in%1),!any(effct_matrix[,i]%in%4),any(effct_matrix[,i]%in%2))){
final_effect[i]<-2
}
if(all(!any(effct_matrix[,i]%in%1),!any(effct_matrix[,i]%in%4),any(effct_matrix[,i]%in%3))){
final_effect[i]<-3
}
if(all(!any(effct_matrix[,i]%in%1),!any(effct_matrix[,i]%in%4),!any(effct_matrix[,i]%in%2),!any(effct_matrix[,i]%in%3))){
final_effect[i]<-0
}
}
names(final_effect)<-colnames(effct_matrix)
table(final_effect)
final_effect<-final_effect[phat$tree_col$order]
pheatmap(t(final_effect),cluster_row = FALSE,cluster_col = FALSE,
           color = c("#FFFFFF","#C00000","#64A36F","#E36488","#2980B9","#000000"),
           legend_breaks = 0:5,
           legend_labels=c("NO_TD","DUP in OG","ST","DT","DUP in TSG","DUP in TSG&OG"),
           cutree_rows=dim(final_effect)[1],show_colnames=F,treeheight_row=0, treeheight_col=0,cellheight=20)


########OG and TSG in only_2_3 were studied
only_2_3_chr_cnv<-all_sample_chr_cv[only_2_3]
##OG
KRAS<-c(25357723,25403870)
KRAS_chr<-"12"
KRAS_group_1<-UCSC_paragraph(only_2_3_chr_cnv,KRAS,KRAS_chr)

EGFR<-c(55086714,55324313)
EGFR_chr<-"7"
EGFR_group_1<-UCSC_paragraph(only_2_3_chr_cnv,EGFR,EGFR_chr)

CAPRIN2<-c(30862486,30907885)
CAPRIN2_chr<-"12"
CAPRIN2_group_1<-UCSC_paragraph(only_2_3_chr_cnv,CAPRIN2,CAPRIN2_chr)

##TSG
KMT2B<-c(36208921,36229779)
KMT2B_chr<-"19"
KMT2B_group_1<-UCSC_paragraph(only_2_3_chr_cnv,KMT2B,KMT2B_chr)



library(pheatmap)
hotmap<-function(x){
  pheatmap(x,cluster_row = FALSE,
           color = c("#FFFFFF","#C00000","#64A36F","#E36488","#2980B9"),
           legend_breaks = 0:4,
           legend_labels=c("NO_TD","DUP in OG","ST","DT","DUP in TSG"),
           cutree_rows=dim(x)[1],show_colnames=F,treeheight_row=0, treeheight_col=0,cellheight=20)
}

####group_1
effct_matrix<- matrix(0,nrow = 4,ncol =length(only_2_3_chr_cnv))
rownames(effct_matrix)<-c("KRAS","EGFR","CAPRIN2","KMT2B")
colnames(effct_matrix)<-names(only_2_3_chr_cnv)
test_group<-list(KRAS_group_1,EGFR_group_1,CAPRIN2_group_1,KMT2B_group_1)
names(test_group)<-rownames(effct_matrix)
#OG
for(i in 1:3){
for(j in 1:dim(test_group[[i]])[1]){
if(test_group[[i]][j,]$effect=="DUP"){
effct_matrix[i,test_group[[i]][j,]$sample]<-1
}
if(test_group[[i]][j,]$effect=="ST"){
effct_matrix[i,test_group[[i]][j,]$sample]<-2
}
if(test_group[[i]][j,]$effect=="DT"){
effct_matrix[i,test_group[[i]][j,]$sample]<-3
}
}
}
#TSG
for(i in 4){
for(j in 1:dim(test_group[[i]])[1]){
if(test_group[[i]][j,]$effect=="DUP"){
effct_matrix[i,test_group[[i]][j,]$sample]<-4
}
if(test_group[[i]][j,]$effect=="ST"){
effct_matrix[i,test_group[[i]][j,]$sample]<-2
}
if(test_group[[i]][j,]$effect=="DT"){
effct_matrix[i,test_group[[i]][j,]$sample]<-3
}
}
}
setwd("E:/NSCLC_TDP/LUAD/amp_effect")
save(effct_matrix,file="only_2_3_effctmatrix.Rdata")
phat <- hotmap(effct_matrix)

final_effect<-c()
for(i in 1:dim(effct_matrix)[2]){
if(any(effct_matrix[,i]%in%1)&any(effct_matrix[,i]%in%4)){
final_effect[i]<-5
}
if(any(effct_matrix[,i]%in%1)&!any(effct_matrix[,i]%in%4)){
final_effect[i]<-1
}
if(!any(effct_matrix[,i]%in%1)&any(effct_matrix[,i]%in%4)){
final_effect[i]<-4
}
if(all(!any(effct_matrix[,i]%in%1),!any(effct_matrix[,i]%in%4),any(effct_matrix[,i]%in%2))){
final_effect[i]<-2
}
if(all(!any(effct_matrix[,i]%in%1),!any(effct_matrix[,i]%in%4),any(effct_matrix[,i]%in%3))){
final_effect[i]<-3
}
if(all(!any(effct_matrix[,i]%in%1),!any(effct_matrix[,i]%in%4),!any(effct_matrix[,i]%in%2),!any(effct_matrix[,i]%in%3))){
final_effect[i]<-0
}
}
names(final_effect)<-colnames(effct_matrix)
table(final_effect)
final_effect<-final_effect[phat$tree_col$order]
pheatmap(t(final_effect),cluster_row = FALSE,cluster_col = FALSE,
           color = c("#FFFFFF","#C00000","#64A36F","#E36488","#2980B9","#000000"),
           legend_breaks = 0:5,
           legend_labels=c("NO_TD","DUP in OG","ST","DT","DUP in TSG","DUP in TSG&OG"),
           cutree_rows=dim(final_effect)[1],show_colnames=F,treeheight_row=0, treeheight_col=0,cellheight=20)
		   
##########figure 3
###############Extract the MAF file from the interim sample
options(stringsAsFactors = F)
setwd("E:/NSCLC_TDP/LUAD/mutation")
load("E:/NSCLC_TDP/LUAD/CNV/patient_group.Rdata")#######Each group of patients
load("E:/NSCLC_TDP/LUAD/CNV/patient_TDP_group.Rdata")
group_1<-unique(c(only_1,only_1_2,only_1_3))
tumor_sample<-read.table(file="tumor_sample.txt",header = T,sep = "\t")####all samples
tumor_sample_1<-strtrim(tumor_sample[,1],16)
only_2_mut<-intersect(unique(tumor_sample_1),only_2)
only_2_3_mut<-intersect(unique(tumor_sample_1),only_2_3)


only_2index<-c()#######Indicates which ordinals are classified samples
for(i in 1:length(only_2_mut)){
index1<-grep(only_2_mut[i],tumor_sample[,1])
only_2index<-c(only_2index,index1)
}
only_2index_frame<-matrix(0,nrow = dim(tumor_sample)[1],ncol = 2)
only_2index_frame[only_2index,]<-1
only_2index_frame[,2]<-tumor_sample[,1]#########1Representation classification sample
write.table(only_2index_frame,file="only_2index_frame.txt",col.names = F,row.names = F,sep="\t",quote=F)

only_2_3index<-c()######Indicates which ordinals are classified samples
for(i in 1:length(only_2_3_mut)){
index1<-grep(only_2_3_mut[i],tumor_sample[,1])
only_2_3index<-c(only_2_3index,index1)
}
only_2_3index_frame<-matrix(0,nrow = dim(tumor_sample)[1],ncol = 2)
only_2_3index_frame[only_2_3index,]<-1
only_2_3index_frame[,2]<-tumor_sample[,1]#########1Representation classification sample
write.table(only_2_3index_frame,file="only_2_3index_frame.txt",col.names = F,row.names = F,sep="\t",quote=F)

###Mutsig
#####
#########Co-occurring mutated genes
options(stringsAsFactors=F)
setwd("E:/NSCLC_TDP/LUAD/mutation")
LUADonly_2<-read.table(file="only_2_mutation_gene.txt",sep="\t",header = T)
LUADonly_2_3<-read.table(file="only_2_3_mutation_gene.txt",sep="\t",header = T)


OGgene<-read.table(file="E:/TDP/oncogene.txt",sep ="\t",header = T)
TSGgene<-read.table(file="E:/TDP/TSG.txt",sep="\t",header = T)
only_2_mutation<-c(intersect(LUADonly_2[,1],OGgene[,1]),intersect(LUADonly_2[,1],TSGgene[,1]))
only_2_3_mutation<-c(intersect(LUADonly_2_3[,1],OGgene[,1]),intersect(LUADonly_2_3[,1],TSGgene[,1]))


######clustprofile
library(clusterProfiler)
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list_file=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
final_gene<-list(LUADonly_2[,1],LUADonly_2_3[,1])
names(final_gene)<-c("LUADonly_2","LUADonly_2_3")

KEGG_result<-list()
GO_result<-list()
for(i in 1:length(final_gene)){
final_gene_entriz=na.omit(list_file[match(final_gene[[i]],list_file[,"SYMBOL"]),][,2])

KEGG_gene<-enrichKEGG(final_gene_entriz,organism = "human", pvalueCutoff = 0.05)
GO_gene<-enrichGO(final_gene_entriz,'org.Hs.eg.db',ont = "BP", pvalueCutoff = 0.05)
KEGG_result[[i]]<-KEGG_gene@result
GO_result[[i]]<-GO_gene@result
}
names(KEGG_result)<-c("LUADonly_2","LUADonly_2_3")
names(GO_result)<-c("LUADonly_2","LUADonly_2_3")
save(KEGG_result,GO_result,file="kegg&go.Rdata")

###Common path
options(stringsAsFactors=F)
setwd("E:/NSCLC_TDP/LUAD/mutation")
load("kegg&go.Rdata")
observed_path<-c('Natural killer cell mediated cytotoxicity','Mitophagy - animal','Toll-like receptor signaling pathway',
'Autophagy - animal','Antigen processing and presentation','PD-L1 expression and PD-1 checkpoint pathway in cancer')
i=6
all_path<-KEGG_result[[i]]

path_index<-match(observed_path,all_path[,2])
new_kegg<-data.frame(all_path[path_index,c(2,5,8,9)],group=rep(i,length(observed_path)))
#all_kegg<-new_kegg
all_kegg<-rbind(all_kegg,new_kegg)
save(all_kegg,file="all_kegg.Rdata")
final_frame<-all_kegg
final_frame<-na.omit(final_frame)
final_frame$pvalue[final_frame$pvalue> 0.1]<-0.1
final_frame$Description<-factor(final_frame$Description,levels=rev(observed_path))

# library(pheatmap)
# pheatmap(final_frame,cluster_row = FALSE,cluster_col = FALSE,color = colorRampPalette(c("#AF1E23","white","#113F8C"))(50))
library(ggplot2)
p = ggplot(final_frame,aes(group,Description))
pbubble = p+ geom_point(aes(size=Count,color=pvalue))
# 设置渐变色
pr = pbubble+scale_color_gradient(low="red",high = "blue")
# 绘制p气泡图
pr = pr+labs(color=expression(pvalue),size="Count",  
                           x="group",y="Pathway name",title="Pathway enrichment")
pr + theme_bw()

########The immune pathway of mutation map BRCA_only_2_3
####only_2_3
setwd("E:/NSCLC_TDP/LUAD/mutation/only_2_3")
library(dplyr)
new.mutations <- read.table("only_2_3.mutations.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F)
genes <- only_2_3_mutation
new.mutations_1 <- new.mutations[which(new.mutations$Hugo_Symbol %in% genes),]
all_tumor <- unique(new.mutations$Tumor_Sample_Barcode)

#new.mutations_2 <- new.mutations_1[1:4,c(1,2,6)] %>% spread(Hugo_Symbol,Variant_Classification)
new.mutations_1[which(new.mutations_1[,2] == "Missense_Mutation"),2] <- "Missense"
new.mutations_1[which(new.mutations_1[,2] == "Nonsense_Mutation"),2] <- "Nonsense"
new.mutations_1[which(new.mutations_1[,2] == "Splice_Region"),2] <- "Splice Site"
new.mutations_1[which(new.mutations_1[,2] == "Splice_Site"),2] <- "Splice Site"
new.mutations_1[which(new.mutations_1[,2] == "Frame_Shift_Ins"),2] <- "Frame Shift"
new.mutations_1[which(new.mutations_1[,2] == "Frame_Shift_Del"),2] <- "Frame Shift"
new.mutations_1[which(new.mutations_1[,2] == "In_Frame_Del"),2] <- "In frame indel"
new.mutations_1[which(new.mutations_1[,2] == "In_Frame_Ins"),2] <- "In frame indel"
new.mutations_1[which(new.mutations_1[,2] == "3'UTR"),2] <- "Syn"
new.mutations_1[which(new.mutations_1[,2] == "5'UTR"),2] <- "Syn"
new.mutations_1[which(new.mutations_1[,2] == "3'Flank"),2] <- "Syn"
new.mutations_1[which(new.mutations_1[,2] == "5'Flank"),2] <- "Syn"
new.mutations_1[which(new.mutations_1[,2] == "Intron"),2] <- "Syn"
new.mutations_1[which(new.mutations_1[,2] == "Translation_Start_Site"),2] <- "Other non syn"
new.mutations_1[which(new.mutations_1[,2] == "Nonstop_Mutation"),2] <- "Other non syn"
new.mutations_1[which(new.mutations_1[,2] == "Silent"),2] <- "Other non syn"

genes_tumor_mutation <- matrix(0,nrow = length(genes),ncol = length(all_tumor))
rownames(genes_tumor_mutation) <- genes
colnames(genes_tumor_mutation) <- all_tumor
for (gene in genes) {
  for (tumor in all_tumor) {
    new.mutations_2 <- filter(new.mutations_1,Hugo_Symbol == gene & Tumor_Sample_Barcode == tumor)
    if(dim(new.mutations_2)[1] >= 1){
      for(i in 1:dim(new.mutations_2)[1]){
      genes_tumor_mutation[new.mutations_2$Hugo_Symbol[i],tumor] <- unique(new.mutations_2$Variant_Classification[i])
    }
	}
  }
}
write.table(genes_tumor_mutation,file = "genes_tumor_mutation.txt",quote = F,sep = "\t")

####################################
new_genes_tumor_mutation <- genes_tumor_mutation
colnames(new_genes_tumor_mutation) <- all_tumor
Variant_Classification <- unique(new.mutations_1$Variant_Classification)
new_genes_tumor_mutation_num <- matrix(0,nrow = length(genes),ncol = length(all_tumor))
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "Syn")] = 1
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "Missense")] <- 2
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "Splice Site")] <- 3
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "Nonsense")] <- 4
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "Frame Shift")] <- 5
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "In frame indel")] <- 6
new_genes_tumor_mutation_num[which(new_genes_tumor_mutation == "Other non syn")] <- 7
rownames(new_genes_tumor_mutation_num) <- rownames(new_genes_tumor_mutation)
colnames(new_genes_tumor_mutation_num) <- colnames(new_genes_tumor_mutation)

library(pheatmap)
hotmap<-function(x){
  pheatmap(x,cluster_row = FALSE,
           color = c("#EFF0F4","#64A36F","#2980B9","#E36488","#C00000","#F87829","#F9FA9B","#4F323B"),
           legend_breaks = 0:7,
           legend_labels=c("NON","Syn","Missense","Splice Site","Nonsense","Frame Shift","In frame indel","Other non syn"),
           cutree_rows=dim(x)[1],show_colnames=F,treeheight_row=0, treeheight_col=0)
}
load('all_kegg.Rdata')
BRCA_only_2_3_kegg<-all_kegg[all_kegg$group==4,]
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list_file=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
all_immune_gene<-unlist(strsplit(BRCA_only_2_3_kegg$geneID,'/'))
all_immune_symbol=na.omit(list_file[match(all_immune_gene,list_file[,"ENTREZID"]),][,3])

phat <- hotmap(new_genes_tumor_mutation_num)

##############figure 4
#########################FPKM&CNV
############LUAD FPKM
options(stringsAsFactors = F)
setwd("/Share2/home/lanxun3/jacklee/NSCLC_TDP/LUAD")
#setwd("E:/TDP/BRCA/TCGA/protein")
all_exp<-read.table(file="TCGA-LUAD.htseq_fpkm.txt",sep = "\t",header = T,na.strings = c("NA"," "),quote = "")
rownames(all_exp)<-strtrim(all_exp[,1],15)
all_exp<-all_exp[,-1]
procoding_gene<-read.table(file="/Share2/home/lanxun3/jacklee/NSCLC_TDP/hg_19_pro.txt",sep = "\t",header = T,na.strings = c("NA"," "),quote = "")
pro_exp<-all_exp[intersect(procoding_gene[,1],rownames(all_exp)),]
save(pro_exp,file="pro_exp.Rdata")

########group_1
options(stringsAsFactors = F)
setwd("/Share2/home/lanxun3/jacklee/NSCLC_TDP/LUAD")
load("patient_group.Rdata")
load('pro_exp.Rdata')
SSG<-c(only_1_2,only_1_3)
colnames(pro_exp)<-gsub('.','-',colnames(pro_exp),fixed=T)
deg_matrix<-pro_exp[,c(SSG,intersect(only_2_3,colnames(pro_exp)))]

library(limma)
group_list=c(rep(0,2),rep(1,54))
exprSet<-deg_matrix
group_info<-data.frame(colnames(exprSet),group_list)

design <- model.matrix(~factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
v <- voom(exprSet,design,normalize="quantile")
fit <- lmFit(v,design)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom = na.omit(tempOutput)
DEG_voom_can <- DEG_voom[with(DEG_voom, (adj.P.Val<0.01&abs(FC)>=1)), ]
proper_gene<-rownames(DEG_voom_can[order(DEG_voom_can$logFC,decreasing = T),])
final_frame<-deg_matrix[proper_gene,]
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
final_symbol=na.omit(list[match(rownames(final_frame),list[,"ENSEMBL"]),][,3])
rownames(final_frame)<-final_symbol

library(pheatmap)
new_frame<-scale(t(final_frame))
new_frame[new_frame< -2]<- -2
new_frame[new_frame> 2]<- 2
pdf('DEG_SSG&only_2_3.pdf')
pheatmap(t(new_frame),cluster_cols = F,cluster_rows = F,show_rownames = T,
color = colorRampPalette(c("#113F8C","white","#AF1E23"))(50),show_colnames = F,main="DEG_SSG&only_2_3")
dev.off()

###Functional enrichment
up_genes<-rownames(DEG_voom[with(DEG_voom, (adj.P.Val<adj_p & logFC>=fold_chang)), ])##Upregulated gene set
down_genes<-rownames(DEG_voom[with(DEG_voom, (adj.P.Val<adj_p & logFC<= -fold_chang)), ])##Down-regulated gene set

up_genes_entriz=na.omit(list[match(up_genes,list[,"ENSEMBL"]),][,2])
down_genes_entriz=na.omit(list[match(down_genes,list[,"ENSEMBL"]),][,2])

##up
library(clusterProfiler)
KEGG_gene<-enrichKEGG(down_genes_entriz,organism = "hsa", pvalueCutoff = 0.05)
GO_gene<-enrichGO(down_genes_entriz,'org.Hs.eg.db',ont = "BP", pvalueCutoff = 0.05)
pdf('down_kegg.pdf')
library(ggplot2)
ggplot(data=KEGG_gene@result[1:20,],aes(x=Description,y=-log10(p.adjust)))+
geom_bar(position=position_dodge(), stat="identity")+coord_flip()+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")
dev.off()

pdf('down_GO.pdf')
ggplot(data=GO_gene@result[1:20,],aes(x=Description,y=-log10(p.adjust)))+
geom_bar(position=position_dodge(), stat="identity")+coord_flip()+theme_bw()+theme(axis.line = element_line(size=1, colour = "black"),panel.border = element_blank()
,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")
dev.off()


############figure 5
##############CCLE drug
##########LUAD patient copy number data was identified in all sample data
setwd("E:/NSCLC_TDP/CCLE/LUAD")
options(stringsAsFactors=F)
LUAD_sample<-read.table(file="E:/NSCLC_TDP/CCLE/clin/LUAD_sample.txt",sep="\t",header = T)
all_sample_cnv<-read.table(file="E:/NSCLC_TDP/CCLE/CNV/CCLE_copynumber_2013-12-03.seg.txt",sep="\t",header = T)
LUAD_cnv<-all_sample_cnv[all_sample_cnv$CCLE_name%in%LUAD_sample$CCLE.name,]
colnames(LUAD_cnv)<-c("sample","chr","start","end","num_probes","value")


all_cnv<-LUAD_cnv
sample_name<-unique(all_cnv$sample)
sample_cnv<-list()
for(i in 1:length(sample_name)){
sample_cnv[[i]]<-all_cnv[all_cnv$sample==sample_name[i],]
}
names(sample_cnv)<-sample_name

chr_name<-unique(sample_cnv[[1]]$chr)

all_sample_chr_cv<-list()
for(w in 1:length(sample_cnv)){
sample_chr_cnv<-list()
for(j in 1:length(chr_name)){
sample_chr_cnv[[j]]<-sample_cnv[[w]][sample_cnv[[w]]$chr==chr_name[j],]
}
names(sample_chr_cnv)<-chr_name

for(z in 1:length(sample_chr_cnv)){
if(dim(sample_chr_cnv[[z]])[1]>1){
TD_likegroup<-c()
for(i in 1:(dim(sample_chr_cnv[[z]])[1]-1)){
if(all(sample_chr_cnv[[z]]$end[i]-sample_chr_cnv[[z]]$start[i]>=100,sample_chr_cnv[[z]]$value[i]>0,sample_chr_cnv[[z]]$value[i]-max(sample_chr_cnv[[z]]$value[i-1],sample_chr_cnv[[z]]$value[i+1])>0.3,abs(sample_chr_cnv[[z]]$value[i-1]-sample_chr_cnv[[z]]$value[i+1])<=0.3)){
TD_likegroup[i]<-1
}else{TD_likegroup[i]<-0}
}
if(all(sample_chr_cnv[[z]]$end[dim(sample_chr_cnv[[z]])[1]]-sample_chr_cnv[[z]]$start[dim(sample_chr_cnv[[z]])[1]]>=100,sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]]>0,sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]]-sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]-1]>0.3,abs(sample_chr_cnv[[z]]$value[dim(sample_chr_cnv[[z]])[1]-1])<=0.3)){
TD_likegroup[dim(sample_chr_cnv[[z]])[1]]<-1
}else{TD_likegroup[dim(sample_chr_cnv[[z]])[1]]<-0}
sample_chr_cnv[[z]]$TD_likegroup<-TD_likegroup
}else{sample_chr_cnv[[z]]$TD_likegroup<-0}
}
all_sample_chr_cv[[w]]<-sample_chr_cnv
}
names(all_sample_chr_cv)<-sample_name


#######Calculate the total of TDS for each sample
all_TD_frame<-list()
for(i in 1:length(all_sample_chr_cv)){
TD_frame<-data.frame()
for(j in 1:length(all_sample_chr_cv[[i]])){

TD_frame<-rbind(TD_frame,all_sample_chr_cv[[i]][[j]],stringsAsFactors=F)
}
all_TD_frame[[i]]<-TD_frame
}
names(all_TD_frame)<-sample_name
save(all_TD_frame,file="all_TD_frame.Rdata")####Record each patient's TD segment
TD_num<-c()
for(i in 1:length(all_TD_frame)){
TD_num[i]<-sum(all_TD_frame[[i]]$TD_likegroup)
}
names(TD_num)<-sample_name


####Calculate the mean TD of each chromosome in all samples
chr_mean<-c()
for(i in 1:length(chr_name)){
chr_num<-c()
for(j in 1:length(all_sample_chr_cv)){
chr_num[j]<-sum(all_sample_chr_cv[[j]][[i]]$TD_likegroup)
}
chr_mean[i]<-sum(chr_num)/length(all_sample_chr_cv)
}
names(chr_mean)<-chr_name


##################Calculate the TDP score for each sample
TDP_score<-c()
for(i in 1:length(all_sample_chr_cv)){
chr_abs<-c()
for(j in 1:length(chr_name)){
chr_abs[j]<-abs(sum(all_sample_chr_cv[[i]][[j]]$TD_likegroup)-chr_mean[j])
}
TDP_score[i]<--sum(chr_abs)/TD_num[i]
}
less20_name<-names(TD_num[TD_num<=20])
names(TDP_score)<-sample_name
save(TDP_score,less20_name,file="TDP_score.Rdata")


library(mclust)
#TDP_score<-TDP_score[setdiff(names(TDP_score),less20_name)]
set.seed(12345678)
re_TDP_score<-Mclust(TDP_score)
group<-re_TDP_score$classification##4个类别
dat<-data.frame(score=TDP_score,group=group)
value_a<-re_TDP_score[["parameters"]][["mean"]][["1"]]
value_b<-re_TDP_score[["parameters"]][["mean"]][["2"]]

TDP_group<-rownames(dat[dat[,1]> value_b,])
NonTDP_group<-rownames(dat[dat[,1]<=value_a,])

hist(TDP_score,breaks = seq(-2.5,0,0.05),freq = F)
library(mixtools)
mixmdl=normalmixEM(TDP_score,mu=re_TDP_score$parameters$mean,sigma=re_TDP_score$parameters$variance$sigmasq,epsilon=1e-3)
plot(mixmdl,which=2)
abline(v=value_a,lwd=2,col="blue")
abline(v=value_b,lwd=2,col="red")
save(TDP_group,NonTDP_group,file="patient_TDP_group.Rdata")


##########TDP patients were divided according to amplification length
TDP_frame<-all_TD_frame[TDP_group]
span_size_log<-list()
for(i in 1:length(TDP_frame)){
test_vec<-TDP_frame[[i]][TDP_frame[[i]]$TD_likegroup==1,]
span_size<-test_vec$end-test_vec$start
span_size_log[[i]]<-log10(span_size/1000)
}
names(span_size_log)<-names(TDP_frame)

hist(unlist(span_size_log),fre=F,breaks = seq(-1,6,0.1))

library(mclust)
set.seed(12345678)
all_TD_size<-unlist(span_size_log)
span_size_distri<-Mclust(all_TD_size)
library(mixtools)
mixmdl=normalmixEM(all_TD_size,mu=span_size_distri$parameters$mean,sigma=span_size_distri$parameters$variance$sigmasq,epsilon=1e-3)
plot(mixmdl,which=2)

####Range of fragments in a category
class_0<-range(all_TD_size[span_size_distri$classification==1])##-1~0.414
class_1<-range(all_TD_size[span_size_distri$classification==c(1,2)])##0.416~0.846
class_2<-range(all_TD_size[span_size_distri$classification==3])##1.85~2.71
class_3<-range(all_TD_size[span_size_distri$classification==4])##2.71~3.85

group_list<-list()
for(i in 1:length(span_size_log)){
group_list[[i]]<-c(length(which(class_1[1]<span_size_log[[i]]&span_size_log[[i]]<class_1[2])),length(which(class_2[1]<span_size_log[[i]]&span_size_log[[i]]<class_2[2])),length(which(class_3[1]<span_size_log[[i]]&span_size_log[[i]]<class_3[2])))
}
names(group_list)<-names(span_size_log)
group_set<-as.data.frame(group_list)
rownames(group_set)<-c("small_span","class_2","large_span")
colnames(group_set)<-names(group_list)


sample_sum<-apply(group_set,2,sum)
names(sample_sum)<-colnames(group_set)


new_group_list<-matrix(0,nrow=dim(group_set)[1],ncol=dim(group_set)[2])
for(i in 1:dim(group_set)[2]){
new_group_list[,i]<-group_set[,i]/sample_sum[i]
}
new_group_list<-as.data.frame(new_group_list)
rownames(new_group_list)<-rownames(group_set)
colnames(new_group_list)<-colnames(group_set)

group_1<-colnames(new_group_list[,new_group_list[1,]>=1/4])
group_2<-colnames(new_group_list[,new_group_list[2,]>=1/4])
group_3<-colnames(new_group_list[,new_group_list[3,]>=1/4])
group_1_2<-colnames(new_group_list[,new_group_list[1,]>=1/4&new_group_list[2,]>=1/4])
group_1_3<-colnames(new_group_list[,new_group_list[1,]>=1/4&new_group_list[3,]>=1/4])
group_2_3<-colnames(new_group_list[,new_group_list[2,]>=1/4&new_group_list[3,]>=1/4])

only_1<-setdiff(group_1,unique(union(group_1_2,group_1_3)))
only_2<-setdiff(group_2,unique(union(group_1_2,group_2_3)))
only_3<-setdiff(group_3,unique(union(group_1_3,group_2_3)))
only_1_2<-setdiff(group_1_2,intersect(intersect(group_1,group_2),group_3))
only_1_3<-setdiff(group_1_3,intersect(intersect(group_1,group_2),group_3))
only_2_3<-group_2_3

new_group_list[4,only_1]<-1
new_group_list[4,only_2]<-2
new_group_list[4,only_3]<-3
new_group_list[4,only_1_2]<-4
new_group_list[4,only_1_3]<-5
new_group_list[4,only_2_3]<-6
final_group<-as.numeric(new_group_list[4,])
names(final_group)<-colnames(new_group_list)
save(NonTDP_group,only_1,only_2,only_3,only_1_2,only_1_3,only_2_3,file="patient_group.Rdata")



##########Drug tolerance in different groups of patients
setwd("E:/NSCLC_TDP/CCLE/LUAD")
options(stringsAsFactors=F)
LUAD_drug<-read.table(file="E:/NSCLC_TDP/CCLE/clin/LUNG_drug.txt",sep="\t",header = T)
load("patient_group.Rdata")
drug_list<-list()
for(i in 1:length(unique(LUAD_drug$Compound))){
drug_list[[i]]<-LUAD_drug[LUAD_drug$Compound==unique(LUAD_drug$Compound)[i],]
}
names(drug_list)<-unique(LUAD_drug$Compound)

three_p<-data.frame()
for(i in c(1:24)){
only_1_IC<-drug_list[[i]][drug_list[[i]]$CCLE.Cell.Line.Name%in%only_1,]$IC50..uM.
only_2_3_IC<-drug_list[[i]][drug_list[[i]]$CCLE.Cell.Line.Name%in%only_2_3,]$IC50..uM.
NonTDP_group_IC<-drug_list[[i]][drug_list[[i]]$CCLE.Cell.Line.Name%in%NonTDP_group,]$IC50..uM.
only_1_2_IC<-drug_list[[i]][drug_list[[i]]$CCLE.Cell.Line.Name%in%only_1_2,]$IC50..uM.
only_1_3_IC<-drug_list[[i]][drug_list[[i]]$CCLE.Cell.Line.Name%in%only_1_3,]$IC50..uM.
group_1_IC<-c(only_1_IC,only_1_2_IC,only_1_3_IC)
c1<-wilcox.test(only_1_2_IC,only_1_IC,alternative = "greater")
three_p<-c(three_p,c1$p.value)

}######i=5,11



df<-data.frame(c(only_1_2_IC,only_1_IC),c(rep("only_1_2_IC",length(only_1_2_IC)),rep("only_1_IC",length(only_1_IC))))
colnames(df)<-c("num","group")

library(ggpubr)

library(digest)
compar<-list(c("only_2_3_IC","only_1_IC"),c("only_2_3_IC","only_1_2_IC"),c("only_2_3_IC","group_1_IC"))
p<-ggboxplot(df, "group", "num",color = "group", size=2,palette =c("#00AFBB", "#E7B800", "#FC4E07","#0000CD"),add = "jitter", shape = "group")
p1<-ggviolin(df, "group", "num",color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07","#0000CD"),add = "jitter", shape = "group")
p






























