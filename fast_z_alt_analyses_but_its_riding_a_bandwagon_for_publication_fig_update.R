#sanity check: what does coverage claim to be z-linked?
setwd("/Users/Andrew/Documents/Manduca_Demography")
library("data.table")

msamp<-fread("S35_median_cov.txt",header=F,stringsAsFactors=F)
colnames(msamp)<-c("Scaffold","Cov")
m_norm<-msamp$Cov/mean(msamp$Cov)
fsamp<-fread("LK5_median_cov.txt",header=F,stringsAsFactors=F)
colnames(fsamp)<-c("Scaffold","Cov")
f_norm<-fsamp$Cov/mean(fsamp$Cov)


cov_trans<-log((m_norm/f_norm),base=2)

cov_check<-as.data.frame(cbind(msamp$Scaffold,cov_trans),stringsAsFactors=F)
colnames(cov_check)<-c("Scaffold","Log2Cov")
zkmanmaster<-read.csv("ZENKOKU_m_sexta_genome_mat.csv",header=T, stringsAsFactors=F)
zk_hack<-as.data.frame(cbind(zkmanmaster$Scaffold,zkmanmaster$ZorA,zkmanmaster$Chromosome),stringsAsFactors=F)
colnames(zk_hack)<-c("Scaffold","ZorA","Chromosome")
scaf_pred <- zk_hack[match(unique(zk_hack$Scaffold), zk_hack$Scaffold),]


zcheck<-as.data.frame(merge(cov_check,scaf_pred,by="Scaffold"))
zcheck$Log2Cov<-as.numeric(as.character(zcheck$Log2Cov))

scaf_length<-fread("Manduca_OGS2.bed",header=F,stringsAsFactors=F)
scaf_length<-cbind(scaf_length[,1],scaf_length[,3])
colnames(scaf_length)<-c("Scaffold","Length")
zcheck<-merge(zcheck,scaf_length,by="Scaffold")

head(zcheck,1)
#per the Neo-Z, we'll only consider scaffolds over the N90 length
zcheck2<-zcheck[which(zcheck$Length>46400),]
put_z<-zcheck2[which(zcheck2$Log2Cov > 0.7),]
put_z <- put_z[is.finite(put_z$Log2Cov),]
put_no<-zcheck2[which(zcheck2$Log2Cov < 0.7),]
put_no <- put_no[is.finite(put_no$Log2Cov),]
put_no[which(put_no$Log2Cov==min(put_no$Log2Cov)),]
#incidentally, there's a few strongly female biased scaffolds
#        Scaffold          Log2Cov ZorA Chromosome Length
#419 scaffold00419 -3.13415    A          0                227848
#502 scaffold00505    -Inf         A          0                 164272
#815 scaffold00850    -Inf         A          0                   62167
#939 scaffold01000    -Inf         A          0                   46784

table(put_z$ZorA,put_z$Chromosome)
#     0  1 13 25
#  A  7  0  1  1
#  Z  0 27  0  0
#so we're recovering 9 additional Z linked scaffolds, 7 of which were previously unplaced...feels good. chr13 and 25 might merit investigation


table(put_no$ZorA,put_no$Chromosome)
#      0  10  11  12  13  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27  28   3   4   5   6   7   8   9
#  A 157  26  61  34  30  17  35  23  23  26  31   6  24  19  42  22  47  27  24  16  22  20  30  28  17  23  19  24
#furthermore, we've got no erroneously Z-linked things


hist(zcheck2$Log2Cov,breaks=30)
abline(v=0.7)

#okay what if we update Z-linkage based on coverage
new_z_scafs<-unique(put_z$Scaffold)

zk_ishin<-zkmanmaster

zk_ishin[which(zk_ishin$Scaffold%in%new_z_scafs),14]<-"Z"
table(zk_ishin$ZorA)
table(zkmanmaster$ZorA)
#from this we pick up 43 genes
#the sum of newly recovered Z-linked sequence is 2,128,494

zkmanmaster<-zk_ishin

library("PMCMR")

expressmaster<-read.csv("Msexta_pan_RNA_specificity_update.csv",header=T, stringsAsFactors=F)

NItgCalc<-function(dn,ds,pn,ps)
{
	unbiasedNI<-sum((ds*pn)/(ps+ds),na.rm=T)/sum((ps*dn)/(ps+ds),na.rm=T)
	return(unbiasedNI)
}


omni<-merge(zkmanmaster,expressmaster,by="Gene",all.x=T,all.y=T)

Sex.Bias<-rep(0,nrow(omni))
Dn<-omni$dN/omni$Non.sites
Ds<-omni$dS/omni$Syn.sites
Dn.Ds<-Dn/Ds
Pn<-omni$pN/omni$Non.sites
Ps<-omni$pS/omni$Syn.sites
Pn.Ps<-Pn/Ps
dos<-(omni$dN/(omni$dN+omni$dS))-(omni$pN/(omni$pN+omni$pS))
alph<-1-(Pn.Ps/Dn.Ds)
omni2<-cbind(omni,Sex.Bias,Dn.Ds,Pn.Ps,alph,dos)
#in this version, we'll throw a switch for SPMs
#first we want to follow the field convention (read: from the other 3 papers) of counting 1.5x sex difference as sex biased. By my maths, that's a 70-30 split
#second we want to ask the data what the best cutoffs are, then use that
#looking at the distributions, 80-20 or even 85-15 is appropriate

omni2[which(omni2$SPM.Female < 0.30),35]<-"Male"
omni2[which(omni2$SPM.Female > 0.70),35]<-"Female"
omni2[which(omni2$SPM.Female < 0.70 & omni2$SPM.Female > 0.30),35]<-"Unbiased"
#here's the switch, we'll use 85-15 for strict
#omni2[which(omni2$SPM.Female < 0.15),35]<-"Male"
#omni2[which(omni2$SPM.Female > 0.85),35]<-"Female"
#omni2[which(omni2$SPM.Female < 0.85 & omni2$SPM.Female > 0.15),35]<-"Unbiased"

om3<-omni2[is.finite(omni2$Dn.Ds),]
#apparently we don't have bias data for 80 genes?
om3<-om3[which(om3$Sex.Bias != 0),]
om3$Sex.Bias<-as.factor(om3$Sex.Bias)
om3$ZorA<-as.factor(om3$ZorA)
om3a<-om3[which(om3$ZorA=="A"),]
om3z<-om3[which(om3$ZorA=="Z"),]

#is the Z biased in composition?
z_bias<-table(om3$ZorA,om3$Sex.Bias)
z_chi<-chisq.test(z_bias)
z_chi
round(z_chi$residuals,3)
#70-30, classic:
#X-squared = 47.366, df = 2, p-value = 5.183e-11
#yes

#    Female    Male   Unbiased
#  A  0.664   -1.215    0.387
#  Z -3.108    5.688   -1.813

#85-15, strict:
#X-squared = 36.762, df = 2, p-value = 1.04e-08
#    Female   Male Unbiased
# A  0.271 -1.156    0.440
# Z -1.269  5.414   -2.058

#Z is enriched for male bias, deplete for female bias in both versions


#start by showing a fast-Z
#this obvs doesn't change
kruskal.test(Dn.Ds ~ ZorA, data=om3)
#data:  Dn.Ds by ZorA
#Kruskal-Wallis chi-squared = 6.8808, df = 1, p-value = 0.008713


median(om3a$Dn.Ds)
median(om3z$Dn.Ds)

#is it Dn, Ds, or both that give us that?
wilcox.test(om3a$dN/om3a$Non.sites,om3z$dN/om3z$Non.sites)
#W = 2624000, p-value = 7.873e-08
median(om3a$dN/om3a$Non.sites)
#0.003951175
median(om3z$dN/om3z$Non.sites)
#0.005244755

#is it Dn, Ds, or both that give us that?
wilcox.test(om3a$dS/om3a$Syn.sites,om3z$dS/om3z$Syn.sites)
#W = 2774500, p-value = 0.0005777
median(om3a$dS/om3a$Syn.sites)
#0.01588916
median(om3z$dS/om3z$Syn.sites)
#0.0186722

#both Dn and Ds are higher on the Z, but Dn moreso = higher dN/dS




kruskal.test(Dn.Ds ~ Sex.Bias, data=om3a)
#70-30
#Kruskal-Wallis chi-squared = 26.261, df = 2, p-value = 1.984e-06
#85-15
#Kruskal-Wallis chi-squared = 27.37, df = 2, p-value = 1.14e-06
posthoc.kruskal.nemenyi.test(Dn.Ds ~ Sex.Bias, data=om3a)
#70-30
#data:  Dn.Ds by Sex.Bias 

 #              Female    Male   
#Male           4.2e-05   -      
#Unbiased       0.63      8.1e-06

#85-15
#               Female    Male   
#Male           1.7e-05   -      
#Unbiased       0.078     2.0e-05

#there is a sex-specific expresion effect on the autos M > UB =F
#note that UB and F trended different from each other in the more strict set, but not strongly

table(om3a$Sex.Bias)
#70-30
# Female       Male   Unbiased 
 #    1856     2477     7219 
 
 #85-15
#  Female     Male Unbiased 
#       713     1748     9091 

kruskal.test(Dn.Ds ~ Sex.Bias, data=om3z)
#70-30
#Kruskal-Wallis chi-squared = 0.46808, df = 2, p-value = 0.7913

#85-15
#Kruskal-Wallis chi-squared = 0.97696, df = 2, p-value = 0.6136

#either way, we still don't see an effect of sex-biased expression in Dn/Ds on the Z


table(om3z$Sex.Bias)

#70-30 bias before updated Z-linkage:
# Female     Male   Unbiased 
# 48         168    284 

#with the updated linkage
# Female     Male   Unbiased 
# 55         177    295 


#but under 85-15, things get rough:
#  Female     Male    Unbiased 
#  25         131     371 

#let's visualize
#split out Z and A
omnia<-omni[which(omni$ZorA=="A"),]
omniz<-omni[which(omni$ZorA=="Z"),]

#again we'll have to throw a switch here, 70-30 or 85-15

#here's 70-30
mhia<-omnia[which(omnia$SPM.Female < 0.3),]
#3219  male-biased genes
fhia<-omnia[which(omnia$SPM.Female > 0.7),]
#2450 female
bofa<-omnia[which(omnia$SPM.Female < 0.7 & omnia$SPM.Female > 0.3),]
#9099 in the middle

mhiz<-omniz[which(omniz$SPM.Female < 0.3),]
#207 male-biased genes
fhiz<-omniz[which(omniz$SPM.Female > 0.7),]
#65 female
bofz<-omniz[which(omniz$SPM.Female < 0.7 & omniz$SPM.Female > 0.3),]
#341

#85-15
#mhia<-omnia[which(omnia$SPM.Female < 0.15),]
#2242  male-biased genes
#fhia<-omnia[which(omnia$SPM.Female > 0.85),]
#982 female
#bofa<-omnia[which(omnia$SPM.Female < 0.85 & omnia$SPM.Female > 0.15),]
#11544 in the middle

#mhiz<-omniz[which(omniz$SPM.Female < 0.15),]
#151 male-biased genes
#fhiz<-omniz[which(omniz$SPM.Female > 0.85),]
#29 female
#bofz<-omniz[which(omniz$SPM.Female < 0.85 & omniz$SPM.Female > 0.15),]
#433


#FIGURE 1
#now let's look at Z and A side-by-side
dummy_male<-rep("dodgerblue",15)
dummy_ub<-rep("dimgray",20)
dummy_fem<-rep("red3",15)
#dummy_male<-rep("dodgerblue",8)
#dummy_ub<-rep("dimgray",35)
#dummy_fem<-rep("red3",7)
dum_cols<-c(dummy_male,dummy_ub,dummy_fem)
layout(matrix(c(1,2,3,1,4,4,1,5,5),3,3,byrow=T))
#B, L, T, R
par(mai=c(0.4,0.65,0.5,0.01))
#main =  expression(paste("Faster-Z in ", italic("Manduca sexta")))
boxplot(om3z$Dn.Ds,notch=T,outline=F,xlim=c(0.7,1.8),las=1,ylab="",cex.lab=2,main="",ylim=c(0.0,1.52),cex.axis=1.4,col="gold3",cex.main=1.4)
boxplot(om3a$Dn.Ds, add=T, at =1.5,notch=T,outline=F, las=1,cex.axis=1.4, col="gray27")
mtext(side=2, line = 2, "dN/dS", cex=1.7)
text(1,1.45,"*",cex=3.2)
text(1,1.49,"*",cex=3.2)
mtext(text=c("Z", "Autos"),side=1, at=c(1,1.5),line=1.1,cex=1.35)
mtext(side =2, line =2.65,"A", at=1.65, cex=2, las=1)
mtext(side = 3, line = 2, at =1.25, "Faster-Z in", cex=1.5)
mtext(side = 3, line = 0.45, at =1.25, expression(italic("Manduca sexta")), cex=1.5)
#B, L, T, R
par(mai=c(0.55,0.7,0.5,0.01))
#1
hist(omniz$SPM.Female, breaks=50, col=dum_cols, main="Z Chromosome", xlab="SPM in Females",las=1,cex.lab=1.6, cex.main=1.7,cex.axis=1.1)
abline(v=0.3, lty=3)
abline(v=0.7, lty=3)
mtext(side =2, las=1, line =2.65,"B", at=144, cex=2)
#2
par(mai=c(0.55,0.3,0.5,0.01))
hist(omnia$SPM.Female, breaks=50, col=dum_cols, main="Autosomes", xlab="SPM in Females",las=1,cex.lab=1.6,cex.main=1.7,ylab="",cex.axis=1.1)
abline(v=0.3, lty=3)
abline(v=0.7, lty=3)
#3
#Zs
par(mai=c(0.04,0.7,0.06,0.1))
boxplot(((mhiz$dN)/(mhiz$Non.sites))/((mhiz$dS)/(mhiz$Syn.sites)),outline=F,notch=T,xlim=c(0,5.5),col="dodgerblue",
ylab="dN/dS",main="",at=0.5,ylim=c(0,1.95),las=1,cex.lab=2,cex.main=1.5,cex.axis=1.2)
#abline(v=0.68, lty=3)
boxplot((fhiz$dN/(fhiz$Non.sites))/(fhiz$dS/(fhiz$Syn.sites)),outline=F,notch=T,add=T,at=1.7,col="red3",las=1,cex.axis=1.2)
#abline(v=1.5,lty=3)
boxplot((bofz$dN/bofz$Non.sites)/(bofz$dS/bofz$Syn.sites),outline=F,notch=T,add=T,at=1.1,col="dimgray",las=1,cex.axis=1.2)
#mtext(text=c("M","UB","F"),side=1, at=c(0.2,1.1,2),line=0.45)
segments(x0=-0.1,y0=median(om3z$Dn.Ds),x1=2.7, y1=median(om3z$Dn.Ds),lty=3)
#autos
boxplot(((mhia$dN)/(mhia$Non.sites))/((mhia$dS)/(mhia$Syn.sites)),outline=F,notch=T,add=T,at=3.8,col="dodgerblue",las=1,cex.axis=1.2)
#abline(v=4.85,lty=3)
boxplot((fhia$dN/(fhia$Non.sites))/(fhia$dS/(fhia$Syn.sites)),outline=F,notch=T,add=T,at=4.9,col="red3",las=1,cex.axis=1.2)
#abline(v=3.925, lty=3)
boxplot((bofa$dN/bofa$Non.sites)/(bofa$dS/bofa$Syn.sites),outline=F,notch=T,add=T,at=4.35,col="dimgray",las=1,cex.axis=1.2)
segments(x0=3,y0=median(om3a$Dn.Ds),x1=6, y1=median(om3a$Dn.Ds),lty=3)
#text(4.8,1.08,"p = 0.09",cex=1.2)
#segments(3.5,1.76,5.2,1.76)
#segments(3.5,1.34,4.35,1.34)
#segments(5.2,1.0,4.35,1.0)
#text(4.8,1.08,"p = 0.09",cex=1.2)
#text(3.87,1.41,"p < 0.001",cex=1.2)
#text(4.35,1.86,"p < 0.0001",cex=1.2)
#mtext(text=c("M","UB","F"),side=1, at=c(3.5,4.35,5.2), line=0.45)
#mtext(text=c("Z Chromosome"),side=1, at=1.5,line=2.1,cex=1.5)
#mtext(text=c("Autosomes"),side=1, at=4.5,line=2.1,cex=1.5)
text(3.8,1.73,"*",cex=3.0)
text(3.8,1.853,"*",cex=3.0)
mtext(side =2, line =2.65,"C", at=2.3, cex=2, las=1)

#4 Polymorphism
par(mai=c(0.33,0.7,0.02,0.1))
boxplot(((mhiz$pN)/(mhiz$Non.sites))/((mhiz$pS)/(mhiz$Syn.sites)),outline=F,notch=T,xlim=c(0,5.5),col="dodgerblue",
ylab="pN/pS",main="",at=0.5,ylim=c(0,1.79),las=1,cex.lab=2,cex.main=1.5, cex.axis=1.2)
#abline(v=0.68, lty=3)
boxplot((fhiz$pN/(fhiz$Non.sites))/(fhiz$pS/(fhiz$Syn.sites)),outline=F,notch=T,add=T,at=1.7,col="red3",las=1, cex.axis=1.2)
#abline(v=1.5,lty=3)
boxplot((bofz$pN/bofz$Non.sites)/(bofz$pS/bofz$Syn.sites),outline=F,notch=T,add=T,at=1.1,col="dimgray",las=1, cex.axis=1.2)
mtext(text=c("M","UB","F"),side=1, at=c(0.2,1.1,2),line=0.75, cex=1.3)
segments(x0=-0.1,y0=median(om3z$Pn.Ps,na.rm=T),x1=2.7, y1=median(om3z$Pn.Ps,na.rm=T),lty=3)
#autos
boxplot(((mhia$pN)/(mhia$Non.sites))/((mhia$pS)/(mhia$Syn.sites)),outline=F,notch=T,add=T,at=3.8,col="dodgerblue",las=1, cex.axis=1.2)
#abline(v=4.85,lty=3)
boxplot((fhia$pN/(fhia$Non.sites))/(fhia$pS/(fhia$Syn.sites)),outline=F,notch=T,add=T,at=4.9,col="red3",las=1, cex.axis=1.2)
#abline(v=3.925, lty=3)
boxplot((bofa$pN/bofa$Non.sites)/(bofa$pS/bofa$Syn.sites),outline=F,notch=T,add=T,at=4.35,col="dimgray",las=1,cex.axis=1.2)
mtext(text=c("M","UB","F"),side=1, at=c(3.5,4.35,5.2), line=0.75, cex = 1.3)
#mtext(text=c("Z Chromosome"),side=1, at=1.5,line=2.1,cex=1.5)
#mtext(text=c("Autosomes"),side=1, at=4.5,line=2.1,cex=1.5)
segments(x0=3,y0=median(om3a$Pn.Ps,na.rm=T),x1=6, y1=median(om3a$Pn.Ps,na.rm=T),lty=3)
mtext(side =2, line =2.65,"D", at=1.8, cex=2, las=1)
text(3.8,1.53,"*",cex=3.0)
text(3.8,1.653,"*",cex=3.0)

###########End Figure 1

#Now we'll consider polymorphism

kruskal.test(Pn.Ps ~ ZorA, data=om3)
#with the new scaffolds:
#Kruskal-Wallis chi-squared = 2.5695, df = 1, p-value = 0.1089
#no Z effect here


kruskal.test(Pn.Ps ~ Sex.Bias, data=om3)
#70-30
#Kruskal-Wallis chi-squared = 43.453, df = 2, p-value = 3.667e-10

#85-15
#Kruskal-Wallis chi-squared = 41.853, df = 2, p-value = 8.161e-10
#Sex bias DOES however matter

posthoc.kruskal.nemenyi.test(Pn.Ps ~ Sex.Bias, data=om3)

#data:  Pn.Ps by Sex.Bias 
#70-30
#              Female    Male   
#Male          0.0028     -      
#Unbiased      0.1353    1.4e-10
#M > UB
#M >F
#M >F= UB

#85-15
#              Female       Male   
#Male          0.054        -      
#Unbiased      0.260        4.5e-10
#M > UB
#M >F (more weakly)
#M >F= UB
#diff numbers, same inference from the strict analysis



#so with different patterns of divergence and polymorphism this raises the question of whether alpha is different
#I don't trust simple, per-gene alpha calcs, we'll used are per-class summed alphas, just like we used to in the old days
#question 1: does evolution differ for sex-biased genes between the Z and autos
#
maa<-1-NItgCalc(mhia$dN,mhia$dS,mhia$pN,mhia$pS)
faa <-1-NItgCalc(fhia$dN,fhia$dS,fhia$pN,fhia$pS)
baa <-1-NItgCalc(bofa$dN,bofa$dS,bofa$pN,bofa$pS)

mza<-1-NItgCalc(mhiz$dN,mhiz$dS,mhiz$pN,mhiz$pS)
fza <-1-NItgCalc(fhiz$dN,fhiz$dS,fhiz$pN,fhiz$pS)
bza <-1-NItgCalc(bofz$dN,bofz$dS,bofz$pN,bofz$pS)

#and for just autos vs z
#omnia
auto_a<-1-NItgCalc(omnia$dN,omnia$dS,omnia$pN,omnia$pS)
z_a<-1-NItgCalc(omniz$dN,omniz$dS,omniz$pN,omniz$pS)






#Alpha point-ests. 
#all autos: -0.04579666, all Z: 0.04166761


bsstata<-rep(0,size)
mwhole<-omnia
for(i in 1:size)
	{
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#all autos CIs
#(-0.06382467, -0.02816349)


bsstata<-rep(0,size)
mwhole<-omniz
for(i in 1:size)
	{
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#CIs
#(-0.0331725,0.1141777)

#70-30
#                  M-biased              Unbiased           F-biased
#Autos   -0.04718937          -0.04115905     -0.06242963
#Z            -0.1180537            0.1270535         -0.1218596

#85-15
#                  M-biased              Unbiased           F-biased
#Autos    -0.03973161          -0.04666306        -0.05242167
#Z            -0.1566462            0.09248431          0.01032131


#will we get an auto vs z effect?
testor<-abs((auto_a-z_a))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(omnia,omniz)
pa<-nrow(omnia)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}

hist(bsstata)
abline(v=testor)
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.0390
#boom adaptive Z!

#so do our results per partition hold up?
# Z effect for unbiased being more adaptive, maaaybeee less for m-biased on the Z 
#let's permute
#Male-biased auto vs unbiased auto
testor<-abs((maa-mza))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(mhia,mhiz)
pa<-nrow(mhia)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}

hist(bsstata)
abline(v=testor)
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#70-30
# 0.3704
#85-15
#0.1976
#nope, m-biased things look the sameish

testor<-abs((baa-bza))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(bofa,bofz)
pa<-nrow(bofa)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}

hist(bsstata)
abline(v=testor)
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#70-30
#0.0068
#85-15
#0.0111
unbiased Z > unbiased A

testor<-abs((faa-fza))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(fhia,fhiz)
pa<-nrow(fhia)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}

hist(bsstata)
abline(v=testor)
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#70-30
#0.7004
#85-15
#0.7956
#F is F...


#to visualize, we need to get bootstrap CIs

size<-10000
bsstata<-rep(0,size)
mwhole<-mhiz
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#M Z
#(-0.3018934, 0.03571681)


bsstata<-rep(0,size)
mwhole<-mhia
for(i in 1:size)
	{
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#M auto
#(-0.08588784,  -0.0102056)


bsstata<-rep(0,size)
mwhole<-bofz
for(i in 1:size)
	{
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#unbiased Z
#(0.03924787, 0.2166495)


bsstata<-rep(0,size)
mwhole<-bofa
for(i in 1:size)
	{
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#unbiased autos
#(-0.06428019,  -0.01926753)

size<-10000
bsstata<-rep(0,size)
mwhole<-fhiz
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#F Z
#(-0.5267663, 0.1532903)


size<-10000
bsstata<-rep(0,size)
mwhole<-fhia
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN,bsap$pS)
	}
bsstata<-sort(bsstata)
bsstata[250]
bsstata[9750]
#F auto
#(-0.1131685, -0.0163059)

#FIGURE 2A
alphvec<-c(z_a, auto_a, mza, maa, bza, baa, fza, faa)
dumx<-c(0.5, 0.75, 1.5,1.75, 2.25, 2.5, 3, 3.25)
plot(dumx,alphvec, pch=c(16,15), ylim=c(-0.55,0.3), xlim=c(0.2, 3.5), xaxt='n', xlab="", ylab="", las=1, cex=c(2.5,2.5,2,2,2,2,2,2),
main=expression(paste("Adaptive evolution across the ", italic("Manduca sexta"), " genome")),cex.main=1.4, col=c("gold3","black"))
legend(0.3, -0.33, c("Z","Autosomes"),pch=c(16,15),cex=1.2, col=c("gold3","black"))
mtext(side = 3, line = 1, "A", cex=2, at = 0.07)
#alpha
mtext(side = 2, line = 3, expression(alpha), las =1, cex=2)
abline(v=1.125, lty=3)
mtext(side = 1, line = 3, "Whole Genome", at=0.48, cex=1.9)
#mtext(side = 1, line = 2, "Autos", at=0.35, cex=1.4)
#mtext(side = 1, line = 2, "Z", at=0.85, cex=1.4)

mtext(side = 1, line = 1, "Male-biased", at=1.52, cex=1.4)
mtext(side = 1, line = 1, "Unbiased", at=2.24, cex=1.4)
mtext(side = 1, line = 1, "Female-biased", at=3.2, cex=1.4)
mtext(side = 1, line = 3, "Sex-bias class", at=2.5, cex=1.9)

#CIs 
#all autos CIs
#(-0.06382467, -0.02816349)
segments(x0=.75 , y0 =-0.06382467 , x1=.75 , y1=-0.02816349)
segments(x0=.65 , y0 = -0.06382467, x1= .85, y1=-0.06382467)
segments(x0=.65 , y0 =-0.02816349 , x1= .85, y1=-0.02816349)

#all Z CIs
#(-0.0331725,0.1141777)
segments(x0=.5 , y0 =-0.0331725 , x1=.5 , y1=0.1141777,col="gold3")
segments(x0=.4 , y0 = -0.0331725, x1=.6 , y1=-0.0331725,col="gold3")
segments(x0=.4 , y0 =0.1141777 , x1=.6 , y1=0.1141777,col="gold3")

#M auto
#(-0.08588784,  -0.0102056)
segments(x0= 1.75, y0 =-0.08588784 , x1= 1.75, y1= -0.0102056)
segments(x0=1.7 , y0 =-0.08588784 , x1=1.8 , y1=-0.08588784)
segments(x0=1.7 , y0 = -0.0102056 , x1=1.8 , y1= -0.0102056)

#M Z
#(-0.3018934, 0.03571681)
segments(x0=1.5 , y0 =-0.3018934 , x1=1.5 , y1=0.03571681,col="gold3")
segments(x0=1.45, y0 =0.03571681 , x1=1.55 , y1=0.03571681,col="gold3")
segments(x0= 1.45, y0 = -0.3018934, x1= 1.55, y1=-0.3018934,col="gold3")

#unbiased autos
#(-0.06428019,  -0.01926753)
segments(x0=2.5 , y0 =-0.06428019 , x1=2.5 , y1=-0.01926753)
segments(x0=2.45 , y0 = -0.06428019, x1=2.55 , y1=-0.06428019)
segments(x0=2.45 , y0 =-0.01926753 , x1=2.55 , y1=-0.01926753)

#unbiased Z
#(0.03924787, 0.2166495)
segments(x0=2.25 , y0 =0.03924787 , x1=2.25 , y1=0.2166495,col="gold3")
segments(x0=2.2 , y0 =0.03924787 , x1=2.3 , y1=0.03924787,col="gold3")
segments(x0=2.2 , y0 = 0.2166495, x1=2.3 , y1=0.2166495,col="gold3")

#F A
segments(x0=3.25 , y0 =-0.1131685 , x1=3.25 , y1=-0.0163059)
segments(x0=3.2 , y0 =-0.1131685 , x1=3.3 , y1=-0.1131685)
segments(x0=3.2 , y0 =-0.0163059 , x1=3.3, y1=-0.0163059)
#F Z
#(-0.5267663, 0.1532903)
segments(x0=3 , y0 =-0.5267663 , x1=3 , y1=0.1532903,col="gold3")
segments(x0=2.95 , y0 =-0.5267663 , x1=3.05 , y1=-0.5267663,col="gold3")
segments(x0=2.95 , y0 =0.1532903 , x1=3.05, y1=0.1532903,col="gold3")

#####end figure 2

#and the final question: what's the effect ratio of Zs to As?
msz4x<-fread("/Users/Andrew/Documents/Manduca_Demography/Msexta_4x_z_thetas.thetas",header=T,stringsAsFactors=F)
zth<-exp(msz4x$Watterson)
msa4x<-fread("/Users/Andrew/Documents/Manduca_Demography/Msexta_4x_auto_thetas.thetas",header=T,stringsAsFactors=F)
ath<-exp(msa4x$Watterson)

median(zth)/median(ath)
#0.3747504

mean(zth)/mean(ath)
#0.4432476

#so there's definitely an extra reduction of population size going on.




#hacking out a datatable to be published with the paper
ompaper<-cbind(omni2$Gene, omni2$Scaffold, omni2$Chromosome, omni2$ZorA, omni2$Non.sites, omni2$Syn.sites,omni2$pN,omni2$pS,omni2$dN,omni2$dS,omni2$SPM.Female)
colnames(ompaper)<-c("Gene","Scaffold","Chromosome","ZorA","Non.sites","Syn.sites","pN","pS","dN","dS","SPM.Female")

#write.csv(ompaper,"Manduca_genome_matrix_for_pub.csv",quote=F,row.names=F)
