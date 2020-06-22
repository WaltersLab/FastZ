#can we double-down and get a monarch fast Z?
#first we need to turn Megan's fpkms into spms
library("PMCMR")
dp_raw<-read.table("/Users/Andrew/Documents/ExpressionWork/Danaus_Tissue_Specificity/Input_Data/danaus.fpkm.tmm.matrix.txt",header=T)
#first let's average across reps, then we can parse by sex or tissue afterwards
hd_f<-rowSums(dp_raw[,4:6])/3
hd_m<-rowSums(dp_raw[,7:9])/3
mg_f<-rowSums(dp_raw[,10:12])/3
mg_m<-rowSums(dp_raw[,13:15])/3
th_f<-rowSums(dp_raw[,16:18])/3
th_m<-rowSums(dp_raw[,19:21])/3
ov<-rowSums(dp_raw[,22:24])/3
te<-rowSums(dp_raw[,25:27])/3
ag<-rowSums(dp_raw[,28:30])/3

dp_avg<-cbind(dp_raw[,1:3],hd_f, hd_m, mg_f, mg_m, th_f, th_m, ov, te, ag)


#we may want to come back for the "anatomy of a genome" work, but for now we'll just calculate SPM for sex-bias and ovaries
ov_spm<-ov^2/(ov^2)
###

#sex specificity here
f_all<-(hd_f+ mg_f+th_f+ov)
m_all<-(hd_m+ mg_m+th_m+te+ag)
spm_f<-(f_all^2/(f_all^2+m_all^2))
hist(spm_f,breaks=50)
abline(v=0.5)

#it's a little skewed but it'll do
dp_spm_sex<-cbind(dp_raw[,1:3], spm_f)
colnames(dp_spm_sex)<-c("Gene", "bad","Chromosome","SPM_F")
dp_spm_sex$Chromosome[is.na(dp_spm_sex$Chromosome)] <- 0

#now we load in our Danaus dataframes
danmaster<-read.csv(file="/Users/Andrew/Downloads/evoconfiles/NewKanzen.csv",header=T,stringsAsFactors=F)
dan_append<-read.csv(file="/Users/Andrew/Documents/danausorthologsformanduca.csv",header=T)
colnames(dan_append)<-c("Gene","ZorA","Chromosome","Manduca.Sperm.Ortho")
danmaster<-as.data.frame(danmaster,stringsAsFactors=F)
dantranscend<-merge(danmaster,dan_append,by=c("Gene","Chromosome","ZorA"),all.x=T, all.y=T)

dantran2 <-read.csv(file="/Users/Andrew/Downloads/Z_bois_dplex_genome_mat.csv",header=T,stringsAsFactors=F)
dantran2<-merge(dantran2,dp_spm_sex,by=c("Gene","Chromosome"),all.x=T, all.y=T)
#except I fucked up pn/ps and dn/ds in the original calcs
dantran2$pN.pS<-(dantran2$pN.M/dantran2$Non.sites)/(dantran2$pS.M/dantran2$Syn.sites)
dantran2$dN.dS<-(dantran2$dN/dantran2$Non.sites)/(dantran2$dS/dantran2$Syn.sites)
#it is rectified 
dantran2$ZorA<-as.factor(dantran2$ZorA)


#1 IS THERE A FAST Z?

#Yes
kruskal.test(dN.dS ~ ZorA, data=dantran2)
#data:  dN.dS by ZorA
#Kruskal-Wallis chi-squared = 9.7227, df = 2, p-value = 0.00774

#1.5 WHAT ABOUT THE NEO Z...WHERE IS THE FAST Z?
posthoc.kruskal.nemenyi.test(dN.dS ~ ZorA, data=dantran2)
#         A           Z_anc 
#Z_anc    0.9957      -     
#Z_neo    0.0055      0.0533


boxplot(dN.dS ~ ZorA, data=dantran2,outline=F, notch=T)

#SO THERE'S A FAST Z, BUT IT'S DRIVEN ENTIRELY BY THE NEO

#2 HOW DOES SEX BIAS AFFECT EVOLUTION?

Sex.Bias<-rep(0,nrow(dantran2))
dantran3<-cbind(dantran2,Sex.Bias)
dantran3$ZorA<-as.factor(dantran3$ZorA)
#installing a cutoff switch, 70-30
dantran3[which(dantran3$SPM_F< 0.30),23]<-"Male"
dantran3[which(dantran3$SPM_F>0.70),23]<-"Female"
dantran3[which(dantran3$SPM_F< 0.70&dantran3$SPM_F> 0.30 ),23]<-"Unbiased"
#for our strict analysis, this is a weird asymetric one, but I think 0-10 M, 10-85 UB, 85-100 F is reasonable from the data
#still, we'll just do 85-15s for now
#dantran3[which(dantran3$SPM_F< 0.15),23]<-"Male"
#dantran3[which(dantran3$SPM_F>0.85),23]<-"Female"
#dantran3[which(dantran3$SPM_F< 0.85&dantran3$SPM_F> 0.15 ),23]<-"Unbiased"
#exclude things for which we don't have sex bias data
dantran3<-dantran3[which(dantran3$Sex.Bias !=0),]
dantran3$Sex.Bias<-as.factor(dantran3$Sex.Bias)

kruskal.test(dN.dS ~ Sex.Bias, data=dantran3)
#data:  dN.dS by Sex.Bias
#lax: Kruskal-Wallis chi-squared = 266.46, df = 2, p-value < 2.2e-16
#strict: Kruskal-Wallis chi-squared = 382.07, df = 2, p-value < 2.2e-16
#YES THERE'S A VERY STRONG DIFFERENCE

posthoc.kruskal.nemenyi.test(dN.dS ~ Sex.Bias, data=dantran3)
#70-30
#             Female     Male   
#Male         0.59       -      
#Unbiased     2.2e-14    < 2e-16

#strict
#             Female    Male   
#Male         0.3       -      
#Unbiased     2.9e-14   < 2e-16

#AND IT COMES FROM SEX-LIMITED/BIASED GENES EVOLVING DIFFERENTLY FROM UNBIASED, BUT NOT EACH OTHER


#in this case, both male and female are different from unbiased, but not from each other :o
boxplot(dN.dS ~ Sex.Bias, data=dantran3,outline=F, notch=T)
#both are elevated 


#2.5 HOW DOES SEX-BIAS INTERACT WITH THE (NEO) Z?
#now let's split Z and autos
dan_za<-dantran3[which(dantran3$ZorA=="Z_anc"),]
dan_au<-dantran3[which(dantran3$ZorA=="A"),]
dan_zn<-dantran3[which(dantran3$ZorA=="Z_neo"),]

kruskal.test(dN.dS ~ Sex.Bias, data=dan_zn)
#data:  dN.dS by Sex.Bias
#70-30 Kruskal-Wallis chi-squared = 11.858, df = 2, p-value = 0.002662
#strict Kruskal-Wallis chi-squared =21.137, df = 2, p-value = 2.571e-05
posthoc.kruskal.nemenyi.test(dN.dS ~ Sex.Bias, data=dan_zn)
#70-30
#            Female     Male  
#Male        0.0383     -     
#Unbiased    0.0022     0.2721

#strict
#           Female     Male   
#Male       0.96       -      
#Unbiased   0.08       7.2e-05
#well this is messy..M = F, F = UB, M > UB

#ON THE NEO, WE RECOVER THE CHARLESWORTHIAN PREDICTION (under 70-30, )
#F > UB = M

#IN OTHER WORDS, HEMIZYGOSITY OF F-BIASED THINGS DRIVES THEIR EVOLUTION

#female-biased things drive higher Dn/Ds on the neo!

boxplot(dN.dS ~ Sex.Bias, data=dan_zn, outline=F, notch=T)


kruskal.test(dN.dS ~ Sex.Bias, data=dan_za)
#data:  dN.dS by Sex.Bias
#70-30 Kruskal-Wallis chi-squared = 9.9899, df = 2, p-value = 0.006772
#strict: Kruskal-Wallis chi-squared = 24.706, df = 2, p-value = 4.316e-06
posthoc.kruskal.nemenyi.test(dN.dS ~ Sex.Bias, data=dan_za)
#70-30
#                  Female    Male  
#Male              0.6144    -      
#Unbiased          0.8151    0.0047

#                  Female     Male   
#Male              0.95      -      
#Unbiased          0.29   3.3e-06

boxplot(dN.dS ~ Sex.Bias, data=dan_za, outline=F, notch=T)
#whereas male-biased is elevated on the anc_Z....
#M > UB 
#M = F
#F = UB
# ON THE ANCESTRAL Z HOWEVER, MALE-BIASED THINGS HAVE HIGHER DN/DS


kruskal.test(dN.dS ~ Sex.Bias, data=dan_au)
#data:  dN.dS by Sex.Bias
#70-30: Kruskal-Wallis chi-squared = 249.25, df = 2, p-value < 2.2e-16
#strict: Kruskal-Wallis chi-squared = 341.27, df = 2, p-value < 2.2e-16
posthoc.kruskal.nemenyi.test(dN.dS ~ Sex.Bias, data=dan_au)
#                  Female     Male   
#Male              0.75       -      
#Unbiased          3.5e-14    < 2e-16

#strict
#                     Female  Male   
#Male                 0.28    -      
#Unbiased             3.0e-14   < 2e-16
boxplot(dN.dS ~ Sex.Bias, data=dan_au, outline=F, notch=T)
#and on autos M = F > UB
#AND ON THE AUTOS, WE SEE SEX-LIMITED EFFECT, BUT NOT A DIFF BETWEEN M AND F




#now how about polys?
#
#POLYS
#
#
#3 HOW DOES POLYMORPHISM SORT ACROSS GENE CLASS AND GENOMIC LOCATION?
kruskal.test(pN.pS ~ ZorA, data=dantran2)
#data:  pN.pS by ZorA
#Kruskal-Wallis chi-squared = 34.184, df = 2, p-value = 3.776e-08
posthoc.kruskal.nemenyi.test(pN.pS ~ ZorA, data=dantran2)
#           A           Z_anc 
#Z_anc      3.9e-07     -    
#Z_neo      0.02        0.27 

boxplot(pN.pS ~ ZorA, data=dantran2,outline=F, notch=T)
#SO SCALED POLYMORPHISM IS LOWER ON THE Z, NEO OR ANC



kruskal.test(pN.pS ~ Sex.Bias, data=dantran3)
#70-30: Kruskal-Wallis chi-squared = 319.42, df = 2, p-value < 2.2e-16
#strict: Kruskal-Wallis chi-squared = 348.41, df = 2, p-value < 2.2e-16
posthoc.kruskal.nemenyi.test(pN.pS ~ Sex.Bias, data=dantran3)
#70-30
#                  Female    Male   
#Male         2.3e-07      -      
#Unbiased < 2e-16      < 2e-16

#strict
#               Female     Male   
#Male           1.1e-07    -      
#Unbiased      < 2e-16    < 2e-16
#a sex-limited effect + a m vs f effect

kruskal.test(pN.pS ~ Sex.Bias, data=dan_zn)
#data:  pN.pS by Sex.Bias
#70-30: Kruskal-Wallis chi-squared = 5.7469, df = 2, p-value = 0.0565
#it's a bit closer
#strict: Kruskal-Wallis chi-squared = 7.2997, df = 2, p-value = 0.02599
#only relevant for strict
#posthoc.kruskal.nemenyi.test(pN.pS ~ Sex.Bias, data=dan_zn)
#            Female   Male
#Male        0.25     -   
#Unbiased    0.91     0.024

kruskal.test(pN.pS ~ Sex.Bias, data=dan_za)
#data:  pN.pS by Sex.Bias
#70-30: Kruskal-Wallis chi-squared = 2.7023, df = 2, p-value = 0.2589
#strict: Kruskal-Wallis chi-squared = 4.8797, df = 2, p-value = 0.08717

kruskal.test(pN.pS ~ Sex.Bias, data=dan_au)
#data:  pN.pS by Sex.Bias
#70-30: Kruskal-Wallis chi-squared = 323.92, df = 2, p-value < 2.2e-16
#strict: Kruskal-Wallis chi-squared = 349.61, df = 2, p-value < 2.2e-16

posthoc.kruskal.nemenyi.test(pN.pS ~ Sex.Bias, data=dan_au)
#70-30
#               Female     Male  
#Male           1.8e-06    -     
#Unbiased      <2e-16     <2e-16

#strict
#                 Female     Male  
#Male            8.8e-06     -     
#Unbiased       <2e-16      <2e-16

#the previous pattern is driven by the autos exclusively

boxplot(pN.pS ~ Sex.Bias, data=dan_au,outline=F, notch=T)

#and finally alpha
#ALPHA
#again, we'll make some switches
#
mhia<-dan_au[which(dan_au$SPM_F < 0.3),]
#4721 male-biased genes
fhia<-dan_au[which(dan_au$SPM_F > 0.7),]
#1248 female
bofa<-dan_au[which(dan_au$SPM_F < 0.7 & dan_au$SPM_F > 0.3),]
#7529
mhiza<-dan_za[which(dan_za$SPM_F < 0.3),]
#189 male-biased genes
fhiza<-dan_za[which(dan_za$SPM_F > 0.7),]
#21 female
bofza<-dan_za[which(dan_za$SPM_F < 0.7 & dan_za$SPM_F > 0.3),]
#391 


mhizn<-dan_zn[which(dan_zn$SPM_F < 0.3),]
#184 male-biased genes
fhizn<-dan_zn[which(dan_zn$SPM_F > 0.7),]
#41 female
bofzn<-dan_zn[which(dan_zn$SPM_F < 0.7 & dan_zn$SPM_F > 0.3),]
#243

#strict
mhizn<-dan_zn[which(dan_zn$SPM_F < 0.15),]
#84 male-biased genes
fhizn<-dan_zn[which(dan_zn$SPM_F > 0.85),]
#15 female yiiiiiiikessss
bofzn<-dan_zn[which(dan_zn$SPM_F < 0.85 & dan_zn$SPM_F > 0.15),]
#369

#strict
mhia<-dan_au[which(dan_au$SPM_F < 0.15),]
#2794 male-biased genes
fhia<-dan_au[which(dan_au$SPM_F > 0.85),]
#655 female
bofa<-dan_au[which(dan_au$SPM_F < 0.85 & dan_au$SPM_F > 0.15),]
#10049

#strict
mhiza<-dan_za[which(dan_za$SPM_F < 0.15),]
#189 male-biased genes
fhiza<-dan_za[which(dan_za$SPM_F > 0.85),]
#21 female
bofza<-dan_za[which(dan_za$SPM_F < 0.85 & dan_za$SPM_F > 0.15),]
#391 



#fun fact, each half of the Z is better annotated for sex biased than the Manduca Z...prolly the better RNAseq resolution (not if strict, lol)

#first let's ask the 10,000ft question: is (some part of) the Z more adaptively evolving than the autos? (this uses Zbro data)
tz<-rbind(dan_zn,dan_za)
tza<-1-NItgCalc(tz$dN,tz$dS,tz$pN.M,tz$pS.M)
aaa<-1-NItgCalc(dan_au$dN,dan_au$dS,dan_au$pN.M,dan_au$pS.M)
#           Auto                       Z
#-0.3208786        -0.1294554

testor<-abs((aaa-tza))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(dan_au,tz)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(dan_au)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#1e-04
#Adaptive Z!


bsstata<-rep(0,size)
mwhole<-dan_au
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#auto CIs
#(-0.3466682, -0.2965795)

bsstata<-rep(0,size)
mwhole<-tz
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#wholeZ CIs
#(-0.1973246,-0.05952206)

#okay, 5,000 ft now
azn<-1-NItgCalc(dan_zn$dN,dan_zn$dS,dan_zn$pN.M,dan_zn$pS.M)
aza<-1-NItgCalc(dan_za$dN,dan_za$dS,dan_za$pN.M,dan_za$pS.M)
aaa<-1-NItgCalc(dan_au$dN,dan_au$dS,dan_au$pN.M,dan_au$pS.M)
#            Auto                Anc                       Neo
#-0.3208786     -0.1920317           -0.07206235
#that looks promising...a quick permute?

bsstata<-rep(0,size)
mwhole<-dan_zn
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#neo CI
#(-0.1707368, 0.01237971)

bsstata<-rep(0,size)
mwhole<-dan_za
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#anc CI
#(-0.277029, -0.08532541)

testor<-abs((aaa-aza))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(dan_au,dan_za)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(dan_au)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.0338
#anc more adaptive than autos...


testor<-abs((aaa-azn))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(dan_au,dan_zn)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(dan_au)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#5e-04
#you better bet the Neo Z is more adaptive than the autos

testor<-abs((aza-azn))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(dan_za,dan_zn)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(dan_za)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.079
#neo and anc not differentiable in this way




maa<-1-NItgCalc(mhia$dN,mhia$dS,mhia$pN.M,mhia$pS.M)
faa <-1-NItgCalc(fhia$dN,fhia$dS,fhia$pN.M,fhia$pS.M)
baa <-1-NItgCalc(bofa$dN,bofa$dS,bofa$pN.M,bofa$pS.M)

mzaa<-1-NItgCalc(mhiza$dN,mhiza$dS,mhiza$pN.M,mhiza$pS.M)
fzaa <-1-NItgCalc(fhiza$dN,fhiza$dS,fhiza$pN.M,fhiza$pS.M)
bzaa <-1-NItgCalc(bofza$dN,bofza$dS,bofza$pN.M,bofza$pS.M)

mzna<-1-NItgCalc(mhizn$dN,mhizn$dS,mhizn$pN.M,mhizn$pS.M)
fzna <-1-NItgCalc(fhizn$dN,fhizn$dS,fhizn$pN.M,fhizn$pS.M)
bzna <-1-NItgCalc(bofzn$dN,bofzn$dS,bofzn$pN.M,bofzn$pS.M)


#just alpha me up fam

#alpha according to the 70-30 split
#             M-biased           Unbiased          F-biased
#Autos        -0.2563139         -0.3561465        -0.3527995
#Z_anc        -0.04415237        -0.2941125        -0.4299554
#Z_neo        -0.1132225         -0.09008383        0.1362487

#alpha, but strict bias
#             M-biased            Unbiased           F-biased
#Autos        -0.1633816          -0.3720162         -0.3106446
#Z_anc         0.0128717          -0.273206          -0.9047887
#Z_neo        -0.1206785          -0.08237416         0.2616075



#hell ya, so the neo-Z follows our charlesworth expectations
#and the Z bro set shows the same relationships..
#well, as long as it's signficant...let's permute 

#Neo vs autos first
#male-biased
testor<-abs((maa-mzna))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(mhia,mhizn)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(mhia)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val

#new
#0.1912
#m-auto = m-neo
#strict 
#0.8916
 

testor<-abs((faa-fzna))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(fhia,fhizn)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(fhia)
pb<-nrow(mwhole)
for(i in 1:size)
	{
		#next we'll shuffle our combined geneset
		bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  	bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
	}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#new: 0.0442
#we got it either way!
#strict: 0.1419....not enough power D:


testor<-abs((baa-bzna))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(bofa,bofzn)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
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
#new: 0.004
#new strict: 2e-04
#and unbiased neo > unbiased autos

#now anc vs autos
#male-biased
testor<-abs((maa-mzaa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(mhia,mhiza)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
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

#new: 0.0177
#m-anc > m-auto
#strict: 0.083
#m-auto = m-anc..but it trends


testor<-abs((faa-fzaa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(fhia,fhiza)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
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
#new:0.7401
#strict: 0.0856
#F anc = F auto


testor<-abs((baa-bzaa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(bofa,bofza)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
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
#new: 0.482
#strict: 0.1828
#ub z_anc = ub auto


#and finally anc vs neo
#male-biased
testor<-abs((mzna-mzaa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(mhizn,mhiza)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(mhiza)
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
#new: 0.4996
#strict: 0.3266
#m anc = m neo


testor<-abs((fzna-fzaa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(fhizn,fhiza)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(fhiza)
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
#new: 0.0081
#strict: 0023
#neo and anc different for f-biased

#unbiased
testor<-abs((bzna-bzaa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(bofzn,bofza)
mwhole$pN<-mwhole$pN.M
mwhole$pS<-mwhole$pS.M
pa<-nrow(bofza)
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
#new: 0.0480
#strict: 0.0211
#and unbiased is more adaptive on the neo than the anc

#in conclusion...the neo Z is where the adaptation is occuring...


#now we need to get some bootstrapped CIs for plotting

size<-1000
bsstata<-rep(0,size)
mwhole<-mhiza
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#M Z A
#(-0.104029, 0.2109487)
#( -0.1862542, 0.09125801)

bsstata<-rep(0,size)
mwhole<-mhizn
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#M Z N
#(-0.2579251, 0.01972085)


bsstata<-rep(0,size)
mwhole<-mhia
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#M A 
(-0.2973143, -0.2128369)



bsstata<-rep(0,size)
mwhole<-bofza
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
# B Z A
(-0.3905101, -0.1743549)

bsstata<-rep(0,size)
mwhole<-bofzn
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#B Z N 
#(-0.2335471, 0.0500527)

bsstata<-rep(0,size)
mwhole<-bofa
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#B A
#(-0.3880204, -0.3241477)

bsstata<-rep(0,size)
mwhole<-fhiza
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#F Z A
#(-0.9687395, -0.1006715)

bsstata<-rep(0,size)
mwhole<-fhizn
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#F Z N
#(-0.09668151, 0.3056676)


bsstata<-rep(0,size)
mwhole<-fhia
for(i in 1:size)
	{
		#next sample with replacement to bootstrap
		bsap<-mwhole[sample(nrow(mwhole),replace=T),]
  	bsstata[i]<-1-NItgCalc(bsap$dN,bsap$dS,bsap$pN.M,bsap$pS.M)
	}
bsstata<-sort(bsstata)
bsstata[25]
bsstata[975]
#F A 
#(-0.4368034, -0.2739353)



#####working line...for number and playing with vizualization 
median(dan_au$dN.dS,na.rm=T)
median(dan_zn$dN.dS,na.rm=T)
median(dan_za$dN.dS,na.rm=T)

median(dan_au$pN.pS,na.rm=T)
median(dan_zn$pN.pS,na.rm=T)
median(dan_za$pN.pS,na.rm=T)


dan_props<-table(dantran3$ZorA,dantran3$Sex.Bias)/rowSums(table(dantran3$ZorA,dantran3$Sex.Bias))

danp2<-t(dan_props)
barplot(danp2,col=c("red3","dodgerblue","dimgray"))

man_props<-table(om3$ZorA,om3$Sex.Bias)/rowSums(table(om3$ZorA,om3$Sex.Bias))

manp2<-t(man_props)
barplot(manp2,col=c("red3","dodgerblue","dimgray"))


#BOTTOM OF FIGURE 1 
#1 2 3 4
#1 5 5 5
#1 6 6 6
layout(matrix(c(1,2,3,4,1,5,5,5,1,6,6,6),3,4,byrow=T))
dummy_male<-rep("dodgerblue",15)
dummy_ub<-rep("dimgray",20)
dummy_fem<-rep("red3",15)
dum_cols<-c(dummy_male,dummy_ub,dummy_fem)
#1
par(mai=c(0.4,0.65,0.5,0.05))
boxplot(dan_za$dN.dS,notch=T,outline=F,xlim=c(0.65,2.35),las=1,ylab="",cex.lab=2,main="",ylim=c(0,0.52),cex.axis=1.4,col="gold3")
boxplot(dan_zn$dN.dS, add=T, at =1.5,notch=T,outline=F, las=1,cex.axis=1.4,col="yellowgreen")
boxplot(dan_au$dN.dS, add=T, at =2,notch=T,outline=F, las=1,cex.axis=1.4,col="gray27")
mtext(side=2, line = 2.45, "dN/dS", cex=1.7)
text(1.5,0.36,"*",cex=3.2)
mtext(side =2, line =2.65,"E", at=0.553, cex=2, las=1)
mtext(side = 3, line = 1.98, at =1.5, "Faster-Z in", cex=1.3)
mtext(side = 3, line = 0.3, at =1.5, expression(italic("Danaus plexippus")), cex=1.3)
#
mtext(text=c("Anc-Z", "Neo-Z","Autos"),side=1, at=c(0.85,1.5,2.15),line=0.45,cex=0.93)
#B, L, T, R
par(mai=c(0.6,0.3,0.5,0.08))
#2
hist(dan_za$SPM_F, breaks=50, col=dum_cols, main="Ancestral Z", xlab="SPM in Females",las=1,cex.lab=1.6,cex.main=1.9,cex.axis=1.2)
abline(v=0.3, lty=3)
abline(v=0.7, lty=3)
mtext(side =2, line =0,"F", at=150, cex=2, las=1)
#3
par(mai=c(0.6,0.1,0.5,0.08))
hist(dan_zn$SPM_F, breaks=50, col=dum_cols, main="Neo Z", xlab="SPM in Females",las=1,cex.lab=1.6,ylab="",cex.main=1.9,cex.axis=1.2)
abline(v=0.3, lty=3)
abline(v=0.7, lty=3)
#4
par(mai=c(0.6,0.1,0.35,0.08))
hist(dan_au$SPM_F, breaks=50, col=dum_cols, main="Autosomes", xlab="SPM in Females",las=1,cex.lab=1.6,ylab="",cex.main=1.9,cex.axis=1.2)
abline(v=0.3, lty=3)
abline(v=0.7, lty=3)
#5
#Zs
par(mai=c(0.05,0.65,0.05,0.01))
boxplot(((mhiza$dN)/(mhiza$Non.sites))/((mhiza$dS)/(mhiza$Syn.sites)),outline=F,notch=T,xlim=c(0,8),col="dodgerblue",
ylab="dN/dS",main="",at=0.2,ylim=c(0,0.54),las=1,cex.lab=2,cex.main=1.5,cex.axis=1.3)
boxplot((fhiza$dN/(fhiza$Non.sites))/(fhiza$dS/(fhiza$Syn.sites)),outline=F,notch=T,add=T,at=1.4,col="red3",las=1,cex.axis=1.3)
boxplot((bofza$dN/bofza$Non.sites)/(bofza$dS/bofza$Syn.sites),outline=F,notch=T,add=T,at=0.8,col="dimgray",las=1,cex.axis=1.3)
#mtext(text=c("M","UB","F"),side=1, at=c(0.2,0.8,1.4),line=0.45)
segments(x0=-0.1,y0=median(dan_za$dN.dS,na.rm=T),x1=1.8, y1=median(dan_za$dN.dS,na.rm=T),lty=3)
segments(0.2,0.315,0.8,0.315)
text(0.5,0.333,"p < 0.01",cex=1.2)


#neo-Z
boxplot((mhizn$dN/(mhizn$Non.sites))/(mhizn$dS/(mhizn$Syn.sites)),outline=F,notch=T,add=T,at=3.5,col="dodgerblue",las=1,cex.axis=1.3)
boxplot((fhizn$dN/(fhizn$Non.sites))/(fhizn$dS/(fhizn$Syn.sites)),outline=F,notch=T,add=T,at=4.7,col="red3",las=1,cex.axis=1.3)
boxplot((bofzn$dN/bofzn$Non.sites)/(bofzn$dS/bofzn$Syn.sites),outline=F,notch=T,add=T,at=4.1,col="dimgray",las=1,cex.axis=1.3)
#mtext(text=c("M","UB","F"),side=1, at=c(3.5,4.1,4.7),line=0.45)
segments(x0=2.8,y0=median(dan_zn$dN.dS,na.rm=T),x1=5, y1=median(dan_zn$dN.dS,na.rm=T),lty=3)
text(4.7,0.515,"*",cex=3.2)

#autos
boxplot(((mhia$dN)/(mhia$Non.sites))/((mhia$dS)/(mhia$Syn.sites)),outline=F,notch=T,add=T,at=6.5,col="dodgerblue",las=1,cex.axis=1.3)
boxplot((fhia$dN/(fhia$Non.sites))/(fhia$dS/(fhia$Syn.sites)),outline=F,notch=T,add=T,at=7.7,col="red3",las=1,cex.axis=1.3)
boxplot((bofa$dN/bofa$Non.sites)/(bofa$dS/bofa$Syn.sites),outline=F,notch=T,add=T,at=7.1,col="dimgray",las=1,cex.axis=1.3)
segments(x0=6,y0=median(dan_au$dN.dS,na.rm=T),x1=8, y1=median(dan_au$dN.dS,na.rm=T),lty=3)
text(7.1,0.3,"*",cex=3.2)
text(7.1,0.36,"*",cex=3.2)
text(7.1,0.42,"*",cex=3.2)
text(7.1,0.48,"*",cex=3.2)
#mtext(text=c("M","UB","F"),side=1, at=c(6.5,7.1,7.7),line=0.45)
#mtext(text=c("Ancestral Z"),side=1, at=1.1,line=2.1,cex=1.5)
#mtext(text=c("Neo Z"),side=1, at=4.1,line=2.1,cex=1.5)
#mtext(text=c("Autosomes"),side=1, at=7.1,line=2.1,cex=1.5)
mtext(side =2, line =2.45,"G", at=0.55, cex=2, las=1)

#6 Polys
par(mai=c(0.35,0.65,0.05,0.01))
boxplot(((mhiza$pN.M)/(mhiza$Non.sites))/((mhiza$pS.M)/(mhiza$Syn.sites)),outline=F,notch=T,xlim=c(0,8),col="dodgerblue",
ylab="pN/pS",main="",at=0.2,ylim=c(0,0.8),las=1,cex.lab=2,cex.main=1.5,cex.axis=1.3)
boxplot((fhiza$pN.M/(fhiza$Non.sites))/(fhiza$pS.M/(fhiza$Syn.sites)),outline=F,notch=T,add=T,at=1.4,col="red3",las=1,cex.axis=1.3)
boxplot((bofza$pN.M/bofza$Non.sites)/(bofza$pS.M/bofza$Syn.sites),outline=F,notch=T,add=T,at=0.8,col="dimgray",las=1,cex.axis=1.3)
mtext(text=c("M","UB","F"),side=1, at=c(0.2,0.8,1.4),line=0.45)
segments(x0=-0.1,y0=median(dan_za$pN.pS,na.rm=T),x1=1.8, y1=median(dan_za$pN.pS,na.rm=T),lty=3)



#neo-Z
boxplot((mhizn$pN.M/(mhizn$Non.sites))/(mhizn$pS.M/(mhizn$Syn.sites)),outline=F,notch=T,add=T,at=3.5,col="dodgerblue",las=1,cex.axis=1.3)
boxplot((fhizn$pN.M/(fhizn$Non.sites))/(fhizn$pS.M/(fhizn$Syn.sites)),outline=F,notch=T,add=T,at=4.7,col="red3",las=1,cex.axis=1.3)
boxplot((bofzn$pN.M/bofzn$Non.sites)/(bofzn$pS.M/bofzn$Syn.sites),outline=F,notch=T,add=T,at=4.1,col="dimgray",las=1,cex.axis=1.3)
mtext(text=c("M","UB","F"),side=1, at=c(3.5,4.1,4.7),line=0.45)
segments(x0=2.8,y0=median(dan_zn$pN.pS,na.rm=T),x1=5, y1=median(dan_zn$pN.pS,na.rm=T),lty=3)


#autos
boxplot(((mhia$pN.M)/(mhia$Non.sites))/((mhia$pS.M)/(mhia$Syn.sites)),outline=F,notch=T,add=T,at=6.5,col="dodgerblue",las=1,cex.axis=1.3)
boxplot((fhia$pN.M/(fhia$Non.sites))/(fhia$pS.M/(fhia$Syn.sites)),outline=F,notch=T,add=T,at=7.7,col="red3",las=1,cex.axis=1.3)
boxplot((bofa$pN.M/bofa$Non.sites)/(bofa$pS.M/bofa$Syn.sites),outline=F,notch=T,add=T,at=7.1,col="dimgray",las=1,cex.axis=1.3)
segments(x0=6,y0=median(dan_au$pN.pS,na.rm=T),x1=8, y1=median(dan_au$pN.pS,na.rm=T),lty=3)
text(6.42,0.3,"b",cex=1.6)
text(7.02,0.25,"a",cex=1.6)
text(7.62,0.4,"c",cex=1.6)
mtext(text=c("M","UB","F"),side=1, at=c(6.5,7.1,7.7),line=0.45)
#mtext(text=c("Ancestral Z"),side=1, at=1.1,line=2.1,cex=1.5)
#mtext(text=c("Neo Z"),side=1, at=4.1,line=2.1,cex=1.5)
#mtext(text=c("Autosomes"),side=1, at=7.1,line=2.1,cex=1.5)
mtext(side =2, line =2.45,"H", at=0.81, cex=2, las=1)



#
##
###
######################FIGURE 2B################################
###
##

#and finally alpha
#ALPHA
#again, we'll make some switches
#
mhia<-dan_au[which(dan_au$SPM_F < 0.3),]
#4721 male-biased genes
fhia<-dan_au[which(dan_au$SPM_F > 0.7),]
#1248 female
bofa<-dan_au[which(dan_au$SPM_F < 0.7 & dan_au$SPM_F > 0.3),]
#7529
mhiza<-dan_za[which(dan_za$SPM_F < 0.3),]
#189 male-biased genes
fhiza<-dan_za[which(dan_za$SPM_F > 0.7),]
#21 female
bofza<-dan_za[which(dan_za$SPM_F < 0.7 & dan_za$SPM_F > 0.3),]
#391 


mhizn<-dan_zn[which(dan_zn$SPM_F < 0.3),]
#184 male-biased genes
fhizn<-dan_zn[which(dan_zn$SPM_F > 0.7),]
#41 female
bofzn<-dan_zn[which(dan_zn$SPM_F < 0.7 & dan_zn$SPM_F > 0.3),]
#243

#strict
mhizn<-dan_zn[which(dan_zn$SPM_F < 0.15),]
#84 male-biased genes
fhizn<-dan_zn[which(dan_zn$SPM_F > 0.85),]
#15 female yiiiiiiikessss
bofzn<-dan_zn[which(dan_zn$SPM_F < 0.85 & dan_zn$SPM_F > 0.15),]
#369

#strict
mhia<-dan_au[which(dan_au$SPM_F < 0.15),]
#2794 male-biased genes
fhia<-dan_au[which(dan_au$SPM_F > 0.85),]
#655 female
bofa<-dan_au[which(dan_au$SPM_F < 0.85 & dan_au$SPM_F > 0.15),]
#10049

#strict
mhiza<-dan_za[which(dan_za$SPM_F < 0.15),]
#189 male-biased genes
fhiza<-dan_za[which(dan_za$SPM_F > 0.85),]
#21 female
bofza<-dan_za[which(dan_za$SPM_F < 0.85 & dan_za$SPM_F > 0.15),]
#391 

maa<-1-NItgCalc(mhia$dN,mhia$dS,mhia$pN.M,mhia$pS.M)
faa <-1-NItgCalc(fhia$dN,fhia$dS,fhia$pN.M,fhia$pS.M)
baa <-1-NItgCalc(bofa$dN,bofa$dS,bofa$pN.M,bofa$pS.M)

mzaa<-1-NItgCalc(mhiza$dN,mhiza$dS,mhiza$pN.M,mhiza$pS.M)
fzaa <-1-NItgCalc(fhiza$dN,fhiza$dS,fhiza$pN.M,fhiza$pS.M)
bzaa <-1-NItgCalc(bofza$dN,bofza$dS,bofza$pN.M,bofza$pS.M)

mzna<-1-NItgCalc(mhizn$dN,mhizn$dS,mhizn$pN.M,mhizn$pS.M)
fzna <-1-NItgCalc(fhizn$dN,fhizn$dS,fhizn$pN.M,fhizn$pS.M)
bzna <-1-NItgCalc(bofzn$dN,bofzn$dS,bofzn$pN.M,bofzn$pS.M)
#

#let's visualize 
dumx1<-c(-0.5, 0.5,2,3.5)
dumx2<-c(-0.3, 0.7, 2.2, 3.7)
dumx3<-c(-0.1, 0.9, 2.4, 3.9)
dumy1<-c(aza,mzaa, bzaa, fzaa)
dumy2<-c(azn,mzna, bzna, fzna)
dumy3<-c(aaa,maa, baa, faa)
#main =  expression(paste("Faster-Z in ", italic("Manduca sexta")))
plot(dumx1,dumy1,ylab="",las=1,xaxt='n', ylim=c(-1.65,0.38), xlab="", cex.axis=1.1, cex.lab=2, pch=16, xlim=c(-0.8,4), main =  expression(paste("Adaptive evolution across the ", italic("Danaus plexippus"), " genome")), cex=2, col="gold3",cex.main = 1.5)
points(dumx2,dumy2, pch=17, cex=2,col="yellowgreen")
points(dumx3,dumy3, pch=15, cex=2, col="gray27")

mtext(side = 3, line = 1, "B", cex=2, at = -1.1)
abline(v=0.11,lty=3)
abline(v=-0.63,lty=3)
points(-0.75, tza,pch=17, cex=2.75,col="yellowgreen")
points(-0.75, tza,pch=16, cex=2.15,col="gold3")
mtext(side = 2, line = 3, expression(alpha), las =1, cex=2, at = -0.8)

#wholeZ CI (-0.1973246,-0.05952206)
segments(-0.75,-0.1973246, -0.75,  tza,col="gold3")
segments(-0.75,tza, -0.75,  -0.05952206, col="yellowgreen")
segments(-0.85, -0.1973246, -0.65,  -0.1973246,col="gold3")
segments(-0.85,-0.05952206, -0.65,  -0.05952206, col="yellowgreen")

#anc CI
#(-0.277029, -0.08532541)
segments(-0.5, -0.277029, -0.5,   -0.08532541,col="gold3")
segments(-0.4, -0.277029, -0.6,  -0.277029,col="gold3")
segments(-0.4,-0.08532541, -0.6,  -0.08532541, col="gold3")

##neo CI (-0.1707368, 0.01237971)
segments(-0.3, -0.1707368, -0.3, 0.01237971, col="yellowgreen")
segments(-0.2, -0.1707368, -0.4, -0.1707368, col="yellowgreen")
segments(-0.2, 0.01237971, -0.4, 0.01237971, col="yellowgreen")

#auto CIs (-0.3466682, -0.2965795)
segments(-0.1, -0.3466682, -0.1, -0.2965795, col="gray27")
segments(-0.2, -0.3466682, 0, -0.3466682, col="gray27")
segments(-0.2, -0.2965795, 0, -0.2965795, col="gray27")

#M Z A
#( -0.1862542, 0.09125801)
segments(0.5, -0.1862542, 0.5,   0.09125801,col="gold3")
segments(0.4, -0.1862542, 0.6,  -0.1862542,col="gold3")
segments(0.4,0.09125801, 0.6,   00.09125801, col="gold3")


#M Z N (-0.2579251, 0.01972085)
segments(0.7, -0.2579251, 0.7, 0.01972085, col="yellowgreen")
segments(0.6, -0.2579251, 0.8, -0.2579251, col="yellowgreen")
segments(0.6, 0.01972085, 0.8, 0.01972085, col="yellowgreen")

#M A (-0.2973143, -0.2128369)
segments(0.9, -0.2973143, 0.9, -0.2128369, col="gray27")
segments(0.8, -0.2973143, 1, -0.2973143, col="gray27")
segments(0.8, -0.2128369, 1, -0.2128369, col="gray27")


# B Z A (-0.3905101, -0.1743549,col="gold3")
segments(2, -0.3905101, 2,-0.1743549,col="gold3")
segments(1.9, -0.3905101, 2.1, -0.3905101,col="gold3")
segments(1.9, -0.1743549, 2.1, -0.1743549,col="gold3")

#B Z N (-0.2335471, 0.0500527)
segments(2.2, -0.2335471, 2.2,  0.0500527, col="yellowgreen")
segments(2.1,-0.2335471, 2.3,  -0.2335471, col="yellowgreen")
segments(2.1,0.0500527, 2.3, 0.0500527, col="yellowgreen")

#B A (-0.3880204, -0.3241477)
segments(2.4, -0.3880204, 2.4,   -0.3241477, col="gray27")
segments(2.3,-0.3880204, 2.5,  -0.3880204, col="gray27")
segments(2.3,  -0.3241477, 2.5,   -0.3241477, col="gray27")

#F Z A #(-0.9687395, -0.1006715)
segments(3.5, -0.9687395, 3.5,-0.1006715,col="gold3")
segments(3.4, -0.9687395, 3.6, -0.9687395,col="gold3")
segments(3.4, -0.1006715, 3.6, -0.1006715,col="gold3")

#F Z N (-0.09668151, 0.3056676)
segments(3.7, -0.09668151, 3.7, 0.3056676, col="yellowgreen")
segments(3.6, -0.09668151, 3.8,-0.09668151, col="yellowgreen")
segments(3.6, 0.3056676, 3.8,  0.3056676, col="yellowgreen")

#F A (-0.4368034, -0.2739353)
segments(3.9, -0.4368034, 3.9,  -0.2739353, col="gray27")
segments(3.8, -0.4368034, 4, -0.4368034, col="gray27")
segments(3.8,  -0.2739353, 4,  -0.2739353, col="gray27")

mtext(side = 1, line = 1, "Male-biased", at=0.52, cex=1.4)
mtext(side = 1, line = 1, "Unbiased", at=2.24, cex=1.4)
mtext(side = 1, line = 1, "Female-biased", at=3.5, cex=1.4)
mtext(side = 1, line = 3, "Sex-bias class", at=2, cex=1.9)
mtext(side = 1, line = 3, "Whole Genome", at=-0.5, cex=1.9)
mtext(side = 1, line = 1, "Whole Z", at=-0.82, cex=1.4)
mtext(side = 1, line = 1, "Split Z", at=-0.3, cex=1.4)
legend(-0.61, -01.0, c("Ancestral Z", "Neo Z","Autosomes"), pch=c(16,17,15), col=c("gold3","yellowgreen","gray27"))


#########################

#can we show sex-bias on the Zs?

splitz<-table(dantran3$ZorA,dantran3$Sex.Bias)
autovec<-c(1248, 4721, 7529)
allzvec<-c(85, 463, 521)
#for strict
#autovec<-c(655, 2794, 10049)
#allzvec<-c(36,273, 760)
combz<-rbind(autovec,allzvec)
row.names(combz)<-c("auto","Z")
colnames(combz)<-c("F","M","UB")

allchi<-chisq.test(combz)
allchi
round(allchi$residuals,3)
#yes there's an effect of sex chromosomes on sex bias
#X-squared = 30.043, df = 2, p-value = 2.994e-07
#               F      M     UB
#auto  0.365 -1.191  0.808
#Z     -1.296  4.233 -2.870
#the residuals suggest this is mainly an excess of male-biased things on the Z, secondarily a deficit of unbiased genes
#and with a strict cutoff:
#X-squared = 17.12, df = 2, p-value = 0.0001916
#          F      M     UB
#auto  0.581 -0.899  0.332
#Z    -2.066  3.195 -1.179
#same results


splitchi<-chisq.test(splitz)
splitchi
round(splitchi$residuals,3)
#X-squared = 35.902, df = 4, p-value = 3.032e-07
#splitting into anc and neo is still signficant, but not radically moreso
#and continues to be significant under strict cutoffs:
#X-squared = 46.679, df = 4, p-value = 1.779e-09

#           Female   Male   Unbiased
#  A          0.365 -1.191      0.808
#  Z_anc -1.483  4.453     -2.970
#  Z_neo -0.279  1.352     -0.972
#here we see that the highst 3 residuals come from the anc-Z

#same with the strict cutoffs:
#         Female   Male Unbiased
#  A      0.581  -0.899    0.332
#  Z_anc -1.406   5.553   -2.602
#  Z_neo -1.528  -1.464    1.166


#and the final question: what's the effect ratio of Zs to As?
z4x<-fread("/Users/Andrew/Documents/Manduca_Demography/Dplex_4x_all_z_thetas.thetas",header=T,stringsAsFactors=F)
zth<-exp(z4x$Watterson)
a4x<-fread("/Users/Andrew/Documents/Manduca_Demography/Dplex_4x_auto_thetas.thetas",header=T,stringsAsFactors=F)
ath<-exp(a4x$Watterson)
median(zth)/median(ath)
#0.4478786
mean(zth)/mean(ath)
#0.6561195
#:o well that's surprising
#what about parts of the Z
az4x<-fread("/Users/Andrew/Documents/Manduca_Demography/Dplex_4x_anc_thetas.thetas",header=T,stringsAsFactors=F)
azth<-exp(az4x$Watterson)
nz4x<-fread("/Users/Andrew/Documents/Manduca_Demography/Dplex_4x_neo_thetas.thetas",header=T,stringsAsFactors=F)
nzth<-exp(nz4x$Watterson)
median(azth)/median(nzth)
#0.4480511
mean(azth)/mean(nzth)
#0.601671

median(azth)/median(ath)
#0.3809617
mean(azth)/mean(ath)
#0.5810242

median(nzth)/median(ath)
#0.8502639
mean(nzth)/mean(ath)
#0.9656843
