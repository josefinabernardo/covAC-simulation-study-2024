# Josi 31/07/2023
#
# exact sim does not work with gee
#
rm(list=ls(all=TRUE))
library(MASS)
library(OpenMx)
library(geepack)   # replaces gee()  ... tests are chi2 not t.
#
# function to calculate power given squared wald test
powchi=function(alpha, df, ncp) {
crit=qchisq(alpha,df,lower=F)
if (abs(ncp)<.0001) { ncp=0 }
power_=pchisq(crit,df,ncp,lower=F)
power_
}
#
# ------------------------------------------------------------------------- settings
#
seed=NA  # if NA random seed, if !is.na(NA) a chosen seed
if (!is.na(seed)) {set.seed(seed)}
#
nrep=2000  # number of replications ....
#
alpha=.05 # for power
# cmethod='exchangeable'  # gee error cov structure
cmethod='independence'  # gee error cov structure
# cmethod='unstructured'  # gee error cov structure
#
stzPRS=F # TRUE # standardize the prs
#
nmz=4000  # sample size MZ pairs
ndz=4000  # sample size DZ pairs
#
par_as=p_as=sqrt(.4)  # A
par_es=p_es=sqrt(.3)  # E
par_cs=p_cs=sqrt(.3)  # C
#
#             sources of A-C covariance effect size based on phenotypic variance
par_gs=p_gs= sqrt(c(0, .0025, .01)) # ,.01,.025))   # parent genotype to twins phenotype
par_bs=p_bs=  sqrt(c(0, .0025, .01)) #,.01,.025))   # twin 1 (2) genotype to twin 2 (1) phenotype ................ source of misfit
par_xs=p_xs=sqrt(.0)  # twin 1 (2) phenotype to twin 2 (1) phenotype - "sib interaction model"
#
#
nset=
length(p_as)*length(p_es)*length(p_cs)*
length(p_gs)*length(p_bs)*length(p_xs)
#
setkeep=matrix(NA,nset,30)   # to keep settings
reskeep=matrix(NA,nset,30)   # to keep results each data set
mxkeep=matrix(NA,nset,30) # openmx results
colnames(setkeep) = c(c('nmz','ndz','a','c','e','g','b','x','prs','A'), rep(NA, 20))
#
ii=0  # count sets
for (par_a in par_as) {
for (par_c in par_cs) {
for (par_e in par_es) {
#
for (par_g in par_gs) {
for (par_b in par_bs) {
for (par_x in par_xs) {
ii=ii+1 # count sets in factorial design
#
print(c(ii))
#
ng=100  # number of diallelic loci
ngp=10  # number of loci comprising polygenic score pgs (0 <= npg <= ng).
p_pgs=ngp/ng  # percentage of genetic variance explained by pgs
p_A=1-p_pgs # not explained = A without pgs effect
# e.g., if par_as^2 = .4, then this A variance is due to ng genes
#                         of .4*(npg/ng) is due to the PRS
#                         Given var(PH) = 1 (assuming no covAC), the PRS explained {.4*(npg/ng)}/1 of the phenotypic variance
#
setkeep[ii,1:10]=c(nmz,ndz, par_a, par_c, par_e, par_g, par_b, par_x, p_pgs, p_A)
# colnames(setkeep) = c('nmz','ndz','a','c','e','g','b','x','prs','A')
#
VA1=p_A; VP=p_pgs;VC=1; VE=1 # .... VA1+VP = par_as^2
##
# simulate data exactly:
# the sample MZ and DZ phenotypic covariance = the population matrices
#                 m    m   v    v   t1    t1    t2    t2
#                 1   2    3    4   5     6     7      8   9  10   11
SLmz=SLdz=diag(c(VA1, VP, VA1, VP, VA1/2, VP/2,VA1/2, VP/2,VC, VE, VE))
#
#SLmz[5,7]=SLmz[7,5]=VA1/2   #
#SLmz[6,8]=SLmz[8,6]=VP/2    #
#
# simulate the latent variables exactly
#
Ldz=mvrnorm(ndz, rep(0,11), Sigma=SLdz, emp=T)  # emp=T means exact data simulation, cov(Ldz) = SLdz
Lmz=mvrnorm(nmz, rep(0,11), Sigma=SLmz, emp=T)  # emp=T means exact data simulation, cov(Lmz) = SLmz  #
# build the exact simulated data DZ
Am_=Ldz[,1]
Pm=Ldz[,2] # prs
Af_=Ldz[,3]
Pf=Ldz[,4] # prs
At1r=Ldz[,5]
Pt1r=Ldz[,6]
At2r=Ldz[,7]
Pt2r=Ldz[,8]
#
At1_=.5*Am_+.5*Af_+At1r # A res t1
Pt1=.5*Pm+.5*Pf+Pt1r # Prs t1
At2_=.5*Am_+.5*Af_+At2r #t2
Pt2=.5*Pm+.5*Pf+Pt2r# t2
C=Ldz[,9]
E1=Ldz[,10]
E2=Ldz[,11]
#
Am=Am_+Pm
Af=Af_+Pf
At1=At1_+Pt1
At2=At2_+Pt2
#
Pht1=par_a*At1 + par_e*E1 + par_c*C +  par_g*(Am+Af) + par_b*(At2)
Pht2=par_a*At2 + par_e*E2 + par_c*C +  par_g*(Am+Af) + par_b*(At1)
Pht1=Pht1 + par_x*Pht2
Pht2=Pht2 + par_x*Pht1
# standardize prs
if (stzPRS) {
Pm=scale(Pm)
Pf=scale(Pf)
Pt1=scale(Pt1)
Pt2=scale(Pt2)
}
#
exdatdz=as.data.frame(cbind(Pm, Pf, Pt1, Pt2, Pht1, Pht2))  # exact data sim
exdatdzA=as.data.frame(cbind(Am, Af, At1, At2, Pht1, Pht2))
#
#
# build the exact simulated data MZ
Am_=Lmz[,1]
Pm=Lmz[,2]
Af_=Lmz[,3]
Pf=Lmz[,4]
At1r=Lmz[,5]
Pt1r=Lmz[,6]
At2r=Lmz[,5] # MZ MZ MZ MZ MZ Lmz[,7]
Pt2r=Lmz[,6] # MZ MZ MZ MZ MZ Lmz[,8]
#
At1_=.5*Am_+.5*Af_+At1r # A res
Pt1=.5*Pm+.5*Pf+Pt1r # Prs
At2_=At1_ #.5*Am_+.5*Af_+At2r
Pt2 = Pt1 # MZMZMZMZMZ .5*Pm+.5*Pf+Pt2r #
C=Lmz[,9]
E1=Lmz[,10]
E2=Lmz[,11]
#
Am=Am_+Pm
Af=Af_+Pf
At1=At1_+Pt1
At2=At1 # At2_+Pt2
#
Pht1=par_a*At1 + par_e*E1 + par_c*C +  par_g*(Am+Af) + par_b*(At2)
Pht2=par_a*At2 + par_e*E2 + par_c*C +  par_g*(Am+Af) + par_b*(At1)
Pht1=Pht1 + par_x*Pht2
Pht2=Pht2 + par_x*Pht1
#
## standardize prs
#
if (stzPRS) {
Pm=scale(Pm)
Pf=scale(Pf)
Pt1=scale(Pt1)
Pt2=scale(Pt2)
}
#
exdatmz=as.data.frame(cbind(Pm, Pf, Pt1, Pt2, Pht1, Pht2))  #
exdatmzA=as.data.frame(cbind(Am, Af, At1, At2, Pht1, Pht2))
#
colnames(exdatdz) =colnames(exdatmz) = c('pgsm','pgsf','pgst1','pgst2','pht1','pht2')
# keep following
round(cor(exdatmz),4) -> phrmz_
round(cor(exdatdz),4)      -> phrdz_
round(apply(exdatmz,2,var),4)[6] -> phvarmz_
round(apply(exdatdz,2,var),4)[6] -> phvardz_
#
# --------------------------------- end exact data simualtion
# exdatdz and exdatmz are dataframes, exact data simulation
# when analysing these data, we should "get out" what we "put in"
#
# exact simulated data in the wide (horizontal) format
#
phdatmz_e = as.data.frame(exdatmz)
phdatdz_e = as.data.frame(exdatdz)
colnames(phdatmz_e)=colnames(phdatdz_e) =c('pgsm','pgsf','pgst1','pgst2','pht1','pht2')
#
# add sum of mother and father prs and add mean of twins prs ... we need these additional variables to
#                                                                detect cov(AC)
#
addmz_e=cbind(phdatmz_e$pgsm+phdatmz_e$pgsf, (phdatmz_e$pgst1+phdatmz_e$pgst2)/2)
colnames(addmz_e) = c('pgsmf','mpgst')
phdatmz_e = cbind(phdatmz_e, addmz_e)
adddz_e=cbind(phdatdz_e$pgsm+phdatdz_e$pgsf, (phdatdz_e$pgst1+phdatdz_e$pgst2)/2)
colnames(adddz_e) = c('pgsmf','mpgst')
phdatdz_e = cbind(phdatdz_e, adddz_e)
# add sum and mean
 c(1,2,3,4,5,7,9,10) -> i1
 c(5,6,1,2,3,4,7,8) -> i2
#
 apply(phdatmz_e,2,var)[i2]
round(cor(phdatmz_e),3)[i2,i2]
#
apply(phdatdz_e,2,var)[i2]
round(cor(phdatdz_e),3)[i2,i2]
#
#
#
# [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
#
# Organize data in long format simulated data long format or vertical format
#
#
# Long format exact simulated data
phdatmzL_e = matrix(1,nmz*2,8)
phdatmzL_e[,2]=nmz + c(c(1:nmz),c(1:nmz))
phdatmzL_e[,3]=c(phdatmz_e$pht1,phdatmz_e$pht2)
phdatmzL_e[,4]=c(phdatmz_e$pgsm,phdatmz_e$pgsm)
phdatmzL_e[,5]=c(phdatmz_e$pgsf,phdatmz_e$pgsf)
phdatmzL_e[,6]=c(phdatmz_e$pgst1,phdatmz_e$pgst2)
phdatmzL_e[,7]=c(phdatmz_e$pgsmf,phdatmz_e$pgsmf)
phdatmzL_e[,8]=c(phdatmz_e$mpgst,phdatmz_e$mpgst)
#
ix_ = sort.int(phdatmzL_e[,2], index.return=T)
phdatmzL_e = phdatmzL_e[ix_$ix,]
colnames(phdatmzL_e)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
phdatmzL_e=as.data.frame(phdatmzL_e)
#
#
# Long format exact simulated data
#
phdatdzL_e = matrix(1,ndz*2,8)
phdatdzL_e[,2]=c(c(1:ndz),c(1:ndz))
phdatdzL_e[,3]=c(phdatdz_e$pht1,phdatdz_e$pht2)
phdatdzL_e[,4]=c(phdatdz_e$pgsm,phdatdz_e$pgsm)
phdatdzL_e[,5]=c(phdatdz_e$pgsf,phdatdz_e$pgsf)
phdatdzL_e[,6]=c(phdatdz_e$pgst1,phdatdz_e$pgst2)
phdatdzL_e[,7]=c(phdatdz_e$pgsmf,phdatdz_e$pgsmf)
phdatdzL_e[,8]=c(phdatdz_e$mpgst,phdatdz_e$mpgst)
#
ix_ = sort.int(phdatdzL_e[,2], index.return=T)
phdatdzL_e = phdatdzL_e[ix_$ix,]
colnames(phdatdzL_e)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
phdatdzL_e=as.data.frame(phdatdzL_e)
#
phdatL_e=rbind(phdatmzL_e, phdatdzL_e)
#
#
# simulated stochastically
#                wide        wide     long      long      long mz+dz
# data sets are  phdatdz and phdatdz, phdatdzL, phdatmzL, phdatL
#                wide          wide        long       long        long mz+dz
# simulated exactly
# data sets are  phdatdz_e and phdatdz_e, phdatdzL_e, phdatmzL_e, phdatL_e
#              pheno t1 pheno t2 mother    father   twin 1  nt twin1  twin2   nt twin2  m+f prs  mean twin prs
# colnames [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
#
# --------------------------------------------------- end data sim
# start the analyses
#
# regression analyses. based on simulated data (not exact) ...
#
#
# regression analyses. based on simulated data exact ... with power
# DZ twin 1 only ...
#
eM0dz=lm(pht1~pgst1, data=phdatdz_e) #)$coefficients 			#   just regression of pheno on prs
eM1dz=lm(pht1~pgsmf+pgst1, data=phdatdz_e) #)$coefficients 		#   with pgsmf ... this is equivalent to the transmitted / non-transmitted design
eM2dz=lm(pht1~mpgst+pgst1, data=phdatdz_e) #)$coefficients 		#   with mpgst ...
eM3dz=lm(pht1~pgsmf+mpgst+pgst1, data=phdatdz_e) # )$coefficients 	#   both detects two sources of cov(AC): twin interaction and cult transmission
#round(summary(eM3dz)$coefficients,5)
#
# DZ twin 1 and twin 2 ... switch to geeglm
#
egeeM0dzL=geeglm(ph~pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients #
egeeM1dzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients #
egeeM2dzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients  #
egeeM3dzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients  #
#
# exact twins mz + dz gee
# DZ and MZ twins
egeeM0mzdzL=geeglm(ph~pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients #
egeeM1mzdzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients #
egeeM2mzdzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  #
egeeM3mzdzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  #
#
# using exact data analyses
# get power exact power pgsmf, mpgst, pgsmf+mpgst: 3 tests ...  dz pairs,   dz + mz pairs
# get power exact power pgsmf: 1 tests ...  dz singles
#
# dz 1 pgsmf test ... test of mpgst does not apply given 1 dz
# model M1dz=lm(pht1~pgsmf+pgst1, data=phdatdz) ... 1dz test of pgsmf
tmp=summary(eM1dz)$coefficients
  eM1dzest=tmp[2,1]; eM1dzse=tmp[2,2]; eM1dztest=tmp[2,3]**2; eM1dzpw=powchi(alpha,1,eM1dztest)  # ....**2 t to chi2
reskeep[ii,1:3]=c(eM1dzest,eM1dzse,eM1dzpw) # est st and power
#
# dz 1+2 test power
# model geeM1dzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod, data=phdatdzL) #  wald test is t^2 already
# pgsmf
tmp=summary(egeeM1dzL)$coefficients
  egeeM1dzest=tmp[2,1]; egeeM1dzse=tmp[2,2]; egeeM1dztest=tmp[2,3]**2; egeeM1dzp=tmp[2,4]; egeeM1dzpw=powchi(alpha,1,egeeM1dztest)
  reskeep[ii,4:6]=c(egeeM1dzest,egeeM1dzse,egeeM1dzpw) # est st and power
# dz 1+2 test power
#geeM2dzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL)#)$coefficients  #
# mpgst
tmp=summary(egeeM2dzL)$coefficients
  egeeM2dzest=tmp[2,1]; egeeM2dzse=tmp[2,2]; egeeM2dztest=tmp[2,3]**2; egeeM2dzp=tmp[2,4]; egeeM2dzpw=powchi(alpha,1,egeeM2dztest)
  reskeep[ii,7:9]=c(egeeM2dzest,egeeM2dzse,egeeM2dzpw) # est st and power
#
# geeM3dzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL)#)$coefficients  #
# pgsmf in presence of mfpgs
tmp=summary(egeeM3dzL)$coefficients
egeeM3dzest1=tmp[2,1]; egeeM3dzse1=tmp[2,2]; egeeM3dztest1=tmp[2,3]**2; egeeM3dzp1=tmp[2,4]; egeeM3dzpw1=powchi(alpha,1,egeeM3dztest1)
reskeep[ii,10:12]=c(egeeM3dzest1,egeeM3dzse1,egeeM3dzpw1) # est st and power
# mpgst
egeeM3dzest2=tmp[3,1]; egeeM3dzse2=tmp[3,2]; egeeM3dztest2=tmp[3,3]**2; egeeM3dzp2=tmp[3,4]; egeeM3dzpw2=powchi(alpha,1,egeeM3dztest2)
reskeep[ii,13:15]=c(egeeM3dzest2,egeeM3dzse2,egeeM3dzpw2) # est st and power
#
# mz dz 1+2 test power
#egeeM1mzdzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients #
# pgsmf
tmp=summary(egeeM1mzdzL)$coefficients
  egeeM1mzdzest=tmp[2,1]; egeeM1mzdzse=tmp[2,2]; egeeM1mzdztest=tmp[2,3]**2; egeeM1mzdzp=tmp[2,4]; egeeM1mzdzpw=powchi(alpha,1,egeeM1mzdztest)
  reskeep[ii,16:18]=c(egeeM1mzdzest,egeeM1mzdzse,egeeM1mzdzpw) # est st and power
# dz 1+2 test power
# egeeM2mzdzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  #
# mpgst
tmp=summary(egeeM2mzdzL)$coefficients
  egeeM2mzdzest=tmp[2,1]; egeeM2mzdzse=tmp[2,2]; egeeM2mzdztest=tmp[2,3]**2; egeeM2mzdzp=tmp[2,4]; egeeM2mzdzpw=powchi(alpha,1,egeeM2mzdztest)
  reskeep[ii,19:21]=c(egeeM2mzdzest,egeeM2mzdzse,egeeM2mzdzpw) # est st and power
#
# egeeM3mzdzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  # #
# pgsmf in presence of mfpgs
tmp=summary(egeeM3mzdzL)$coefficients
egeeM3mzdzest1=tmp[2,1]; egeeM3mzdzse1=tmp[2,2]; egeeM3mzdztest1=tmp[2,3]**2; egeeM3mzdzp1=tmp[2,4]; egeeM3mzdzpw1=powchi(alpha,1,egeeM3mzdztest1)
reskeep[ii,22:24]=c(egeeM3mzdzest1,egeeM3mzdzse1,egeeM3mzdzpw1) # est st and power
# mpgst
egeeM3mzdzest2=tmp[3,1]; egeeM3mzdzse2=tmp[3,2]; egeeM3mzdztest2=tmp[3,3]**2; egeeM3mzdzp2=tmp[3,4]; egeeM3mzdzpw2=powchi(alpha,1,egeeM3mzdztest2)
reskeep[ii,25:27]=c(egeeM3mzdzest2,egeeM3mzdzse2,egeeM3mzdzpw2) # est st and power
#
#}} }} }}
#
ipow=c(seq(3,27,3))
colnames(reskeep) = c('dz1mfpgs_e','dz1mfpgs_s','dz1mfpgs1_p',
                       'dz2mfpgs_e','dz2mfpgs_s','dz2mfpgs_p',
                       'dz2mpgsm_e','dz2mpgsm_s','dz2mpgsm_p',
                       'dz2mfpgs_m_e','dz2mfpgs_m_s','dz2mfpgs_m_p',
                       'dz2mpgst_mf_e','dz2mpgst_mf_s','dz2mpgst_mf_p',
                       'mzdzmfpgs_e','mzdzmfpgs_s','mzdzmfpgs_p',
                       'mzdzmpgsm_e','mzdzmpgsm_s','mzdzmpgsm_p',
                       'mzdzmfpgs_m_e','mzdzmfpgs_m_s','mzdzmfpgs_m_p',
                       'mzdzmpgst_mf_e','mzdzmpgst_mf_s','mzdzmpgst_mf_p',NA,NA,NA)
#reskeep[,ipow]
#
# switch to OpenMx covariance structure modeling
#
#             1      2      3 1    4 2      5 3      6 4
vnamesdz=c('pht1','pht2','pgst1','pgst2','pgsm','pgsf')
vnamesmz=c('pht1','pht2','pgst1','pgsm','pgsf')
#
##         1       2    3       4      5      6      7       8  .... 7 and 8 are correlated 1 in MZs
w_mzdat=phdatmz_e[,vnamesmz]
w_dzdat=phdatdz_e[,vnamesdz]
# ==============================================================================================
#
# model covariance matrix:
#  GDFD'G' + Y |  GDFD'
#  FD'G'D      |  DFD'
#
# DFD' is the covariance matrix 2x2 of the parental prs
# D 2x2 diagonal contains the stds of the parental PRS
# F is the 2x2 correlation matrix of the parental PRS.
#
# G contains the regressions of offspring phenotype on parental and offspring prs
# Y is the residual phenotypic covariance matrix modeled using ACE model
#
mz1out=4 # remove 8th variable prs of mz twin 2.
Filter=diag(6)[-mz1out,]
#
RAdz=matrix(c(.5),4,4)
diag(RAdz)=1
RAdz[3,4]=RAdz[4,3]=0  # m f assortative mating
RAmz=RAdz
RAmz[1,2]=RAmz[2,1]=1 # MZ twins
RAfree=matrix(F,4,4)
RAlabels=matrix(NA,4,4)
#
M1 = mxModel("M1",
   mxMatrix(type='Full', nrow=5, ncol=6, free=F, values=Filter, labels=c(NA), name='Filter'),
# mean
   mxMatrix(type="Full", nrow=1, ncol=6,
#
            free=c(T,T,T,T,T,T),
            values=c(0,0,0,0,0,0),
            label=c("mph","mph","mprs","mprs","mprs","mprs"),
            name="expMeandz"),
            mxAlgebra(expression=(expMeandz%*%t(Filter)), name="expMeanmz"),
# prs stdevs matrix D
   mxMatrix(type="Diag", nrow=4, ncol=4,
            free=c(T,T,T,T),
            values=c(.7,.7,.7,.7),
            labels=c("sdp","sdp","sdp","sdp"), lbound=.01,
            name="D"),
# prs correlation matris F
    mxMatrix(type="Symm", nrow=4, ncol=4,
                       free=RAfree, values=RAdz, labels=RAlabels, name="Fdz"),
    mxMatrix(type="Symm", nrow=4, ncol=4,
                       free=RAfree, values=RAmz, labels=RAlabels, name="Fmz"),
# A matrix for P = A + C + E
    mxMatrix(type="Symm", nrow=2, ncol=2,
             free=c(F), values=RAdz[1:2,1:2], labels=c(NA), name="Adz"),
    mxMatrix(type="Symm", nrow=2, ncol=2,
             free=c(F), values=RAmz[1:2,1:2], labels=c(NA), name="Amz"),
# C matrix for P = A + C + E
    mxMatrix(type="Symm", nrow=2, ncol=2,
             free=c(F), values=c(1), labels=c(NA), name="C"),
# E
    mxMatrix(type="Iden", nrow=2, ncol=2, name="E"),
#
# A matrix for P = A + C + E
    mxMatrix(type="Symm", nrow=1, ncol=1,
             free=c(T), values=c(.33), labels=c("sa2"), name="Av"),
# C matrix for P = A + C + E
    mxMatrix(type="Symm", nrow=1, ncol=1,
             free=c(T), values=c(.33), labels=c("sc2"), name="Cv"),
# E matrix for P = A + C + E
    mxMatrix(type="Symm", nrow=1, ncol=1,
             free=c(T), values=c(.33), labels=c("se2"), name="Ev"),
#
    mxMatrix(type="Full", nrow=2, ncol=4,
             free=matrix(c(
              T,T,T,T,
              T,T,T,T),2,4,byrow=T),
             label=matrix(c(
             'a1','b1','g1','g1',
             'b1','a1','g1','g1'),2,4,byrow=T),
             values=matrix(c(
             .1, 0, .01,.01,
              0,.1, .01,.01),2,4,byrow=T), name='G'),
#
    mxAlgebra(expression=Av%x%Amz + Cv%x%C + Ev%x%E, name="Ymz"),
    mxAlgebra(expression=Av%x%Adz + Cv%x%C + Ev%x%E, name="Ydz"),
    mxAlgebra(expression=D%*%Fdz%*%t(D), name="Pdz"),
    mxAlgebra(expression=D%*%Fmz%*%t(D), name="Pmz"),
    mxAlgebra(expression=G%*%Pdz%*%t(G)+Ydz, name="Sdz1"),
    mxAlgebra(expression=G%*%Pmz%*%t(G)+Ymz, name="Smz1"),
#
 #   mxAlgebra(expression=G%*%D%*%Pmz, name="GDPmz"),
 #   mxAlgebra(expression=G%*%D%*%Pdz, name="GDPdz"),
    mxAlgebra(expression=G%*%Pmz, name="GDPmz"),
    mxAlgebra(expression=G%*%Pdz, name="GDPdz"),
    mxAlgebra(expression=Filter%*%rbind(cbind(Smz1,GDPmz),cbind(t(GDPmz),Pmz))%*%t(Filter), name="Smz"),
    mxAlgebra(expression=rbind(cbind(Sdz1,GDPdz),cbind(t(GDPdz),Pdz)), name="Sdz"),
#
mxCI(c('b1','a1','g1'))
#
 )

DZ =  mxModel('DZ',
      mxData( observed=w_dzdat, type="raw"),
   mxExpectationNormal( covariance="M1.Sdz", means="M1.expMeandz",
   dimnames=vnamesdz),  # the fit function
   mxFitFunctionML()
                                                                  )
MZ =  mxModel('MZ',
      mxData( observed=w_mzdat, type="raw"),
   mxExpectationNormal( covariance="M1.Smz", means="M1.expMeanmz",
   dimnames=vnamesmz),  # the fit function
   mxFitFunctionML()
                                                                  )

#
Model_1 <-  mxModel("model2", M1, MZ, DZ,
                  mxFitFunctionMultigroup( c("MZ","DZ"))
                                               )
#
#
# fit the model
Model_1out <- mxRun(Model_1, intervals=F)
#
summary(Model_1out)
sat_1out = mxRefModels(Model_1out, run=T)
mxCompare(sat_1out, Model_1out)
#
#
Model_1b = omxSetParameters(Model_1out, labels=c('b1'), free=F, values=c(0))
Model_1b_out <- mxRun(Model_1b, intervals=T)
mxCompare(Model_1out, Model_1b_out)

Model_1g = omxSetParameters(Model_1out, labels=c('g1'), free=F, values=c(0))
Model_1g_out <- mxRun(Model_1g, intervals=T)
mxCompare(Model_1out, Model_1g_out)

Model_1bg = omxSetParameters(Model_1out, labels=c('g1','b1'), free=F, values=c(0))
Model_1bg_out <- mxRun(Model_1bg, intervals=T)
mxCompare(Model_1out, Model_1bg_out)
#
mxRefModels(Model_1out, run=T) -> sat_1out
mxCompare(sat_1out, Model_1out)
#
 #
# -----------------------------------------------   vars
varnames=c('pht1','pht2')#
#
# the model to calculate expected summary statistics
# this is the twin model phenotypic
# [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
# [1] "pgsm"  "pgsf"  "pgst1" "pgst2" "pht1"  "pht2"  "pgsmf" "mpgst"
nphen1=1
nphen2=2
ACE1  = mxModel("ACE",
#
  # Matrices a, c, and e to store the a, c, and e path coefficients
   mxMatrix(type="Lower", nrow=nphen1, ncol=nphen1,
                  free=c(T), values=t(.3),
                  label=c("a11"),name="a"),
   mxMatrix(type="Lower", nrow=nphen1, ncol=nphen1,
                  free=c(T), values=c(.2),
                  label=c("c11"),name="c"),
   mxMatrix(type="Lower", nrow=nphen1, ncol=nphen1,
                  free=c(T), values=c(.7),
                  label=c("e11"),name="e"),
#
# Matrixes A, C, and E to compute A, C, and E variance components
#
#   mxAlgebra( expression=a %*% t(a), name="A"),    # a^2
#   mxAlgebra( expression=c %*% t(c), name="C"),    # c^2
#   mxAlgebra( expression=e %*% t(e), name="E"),    # e^2
   mxAlgebra( expression=a, name="A"),    # a^2
   mxAlgebra( expression=c, name="C"),    # c^2
   mxAlgebra( expression=e, name="E"),    # e^2
#
#
# Matrix expCovMZ for expected covariance matrix for MZ twins
   mxAlgebra( expression=
           rbind( cbind(A+C+E, A+C),
                  cbind(A+C, A+C+E)),name="expCovMZ"),
# Matrix expCovMZ for expected covariance matrix for DZ twins
   mxAlgebra( expression=
              rbind( cbind(A+C+E, .5%x%A+C),
              cbind(.5%x%A+C, A+C+E)),name="expCovDZ")
                                                         )
# a model the data, the fit function (MZ)
MZmodel=mxModel("MZ",
#
  # Matrix expMean for expected mean vector for MZ and DZ twins
   mxMatrix(type="Full", nrow=1, ncol=4, free=F, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=c(T,T,T,T), values=c(0,0,0,0),
        labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=T, values=c(0,0,0,0),
        labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
   mxMatrix(type="Full", nrow=1, ncol=1,
            free=c(T),values=c(0),label=c("b0"),
            name="Int"),
   mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
mxData(observed=phdatmz_e, type="raw"),
   mxExpectationNormal(covariance="ACE.expCovMZ",
                        means = "expMean", varnames),
    mxFitFunctionML()
                                        )
# a model the data, the fit function (DZ)
DZmodel=mxModel("DZ",
#
  # Matrix expMean for expected mean vector for MZ and DZ twins
   mxMatrix(type="Full", nrow=1, ncol=4, free=F, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=c(T,T,T,T), values=c(0,0,0,0),
        labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=T, values=c(0,0,0,0),
        labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
   mxMatrix(type="Full", nrow=1, ncol=1,
            free=c(T),values=c(0),label=c("b0"),
            name="Int"),
   mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
   mxData(observed=phdatdz_e, type="raw"),
   mxExpectationNormal(covariance="ACE.expCovDZ",
                        means = "expMean", varnames),
    mxFitFunctionML()
                                    )
Model_2 <-  mxModel("twinACE", ACE1, MZmodel, DZmodel, mxFitFunctionMultigroup( c("MZ","DZ") )
  #               mxAlgebra(MZ.objective + DZ.objective, name="minus2loglikelihood"),
  #               mxFitFunctionAlgebra("minus2loglikelihood")
                                               )
# fit the model
Model_2out <- mxRun(Model_2)
summary(Model_2out)
#
#
# sat
#
# -----------------------------------------------   vars
varnames=c('pht1','pht2')#
#
# the model to calculate expected summary statistics
# this is the twin model phenotypic
# [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
# [1] "pgsm"  "pgsf"  "pgst1" "pgst2" "pht1"  "pht2"  "pgsmf" "mpgst"
nphen1=1
nphen2=2
SAT1  = mxModel("SAT",
#
  # Matrices a, c, and e to store the a, c, and e path coefficients
   mxMatrix(type="Stand", nrow=nphen2, ncol=nphen2,
                  free=c(T), values=c(.5),
                  label=c("rmz"),name="Rmz"),
   mxMatrix(type="Stand", nrow=nphen2, ncol=nphen2,
                  free=c(T), values=c(.25),
                  label=c("rdz"),name="Rdz"),
   mxMatrix(type="Diag", nrow=nphen2, ncol=nphen2,
                  free=c(T), values=c(.7),
                  label=c("sd","sd"),name="SD"),
#
#
# Matrix expCovMZ for expected covariance matrix for MZ twins
   mxAlgebra( expression=
               SD%*%Rmz%*%SD,name="expCovMZ"),
# Matrix expCovMZ for expected covariance matrix for DZ twins
   mxAlgebra( expression=
              SD%*%Rdz%*%SD,,name="expCovDZ")
                                                         )
# a model the data, the fit function (MZ)
MZmodel=mxModel("MZ",
#
  # Matrix expMean for expected mean vector for MZ and DZ twins
   mxMatrix(type="Full", nrow=1, ncol=4, free=F, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=c(T,T,T,T), values=c(0,0,0,0),
        labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=T, values=c(0,0,0,0),
        labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
   mxMatrix(type="Full", nrow=1, ncol=1,
            free=c(T),values=c(0),label=c("b0"),
            name="Int"),
   mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
mxData(observed=phdatmz_e, type="raw"),
   mxExpectationNormal(covariance="SAT.expCovMZ",
                        means = "expMean", varnames),
    mxFitFunctionML()
                                        )
# a model the data, the fit function (DZ)
DZmodel=mxModel("DZ",
#
  # Matrix expMean for expected mean vector for MZ and DZ twins
   mxMatrix(type="Full", nrow=1, ncol=4, free=F, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=c(T,T,T,T), values=c(0,0,0,0),
        labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
   mxMatrix(type="Full", nrow=1, ncol=4, free=T, values=c(0,0,0,0),
        labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
   mxMatrix(type="Full", nrow=1, ncol=1,
            free=c(T),values=c(0),label=c("b0"),
            name="Int"),
   mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
   mxData(observed=phdatdz_e, type="raw"),
   mxExpectationNormal(covariance="SAT.expCovDZ",
                        means = "expMean", varnames),
    mxFitFunctionML()
                                    )
Model_3 <-  mxModel("twinSAT", SAT1, MZmodel, DZmodel, mxFitFunctionMultigroup( c("MZ","DZ") )
  #               mxAlgebra(MZ.objective + DZ.objective, name="minus2loglikelihood"),
  #               mxFitFunctionAlgebra("minus2loglikelihood")
                                               )
# fit the model
Model_3out <- mxRun(Model_3)
summary(Model_3out)
#
#
#
sat_2out = mxRefModels(Model_2out, run=T)
mxCompare(sat_2out, Model_2out)
#
#
Model_2g=omxSetParameters(Model_2out, labels='bpgsg', value=0, free=F)
Model_2g_out = mxRun(Model_2g)
#
Model_2b=omxSetParameters(Model_2out, labels='bpgsb', value=0, free=F)
Model_2b_out = mxRun(Model_2b)
# "bpgsb","bpgsg"
Model_2bg = omxSetParameters(Model_2out, labels=c('bpgsb','bpgsg'), free=F, values=c(0))
Model_2bg_out <- mxRun(Model_2bg)
#
ncp1bg=mxCompare(Model_1out, Model_1bg_out)[2,7]
pw1bg=powchi(alpha, 2, ncp1bg)
mxkeep[ii,1:2]=c(ncp1bg,pw1bg) # openmx results
#
ncp1g=mxCompare(Model_1out, Model_1g_out)[2,7]
pw1g=powchi(alpha,1,ncp1g)
mxkeep[ii,3:4]=c(ncp1g,pw1g) # openmx results
#
ncp1b=mxCompare(Model_1out, Model_1b_out)[2,7]
pw1b=powchi(alpha,1,ncp1b)
mxkeep[ii,5:6]=c(ncp1b,pw1b) # openmx results
#
ncp2bg=mxCompare(Model_2out, Model_2bg_out)[2,7]
pw2bg=powchi(alpha, 2, ncp2bg)
mxkeep[ii,7:8]=c(ncp2bg,pw2bg) # openmx results
#
ncp2g=mxCompare(Model_2out, Model_2g_out)[2,7]
pw2g=powchi(alpha,1,ncp2g)
mxkeep[ii,9:10]=c(ncp2g,pw2g) # openmx results
#
ncp2b=mxCompare(Model_2out, Model_2b_out)[2,7]
pw2b=powchi(alpha,1,ncp2b)
mxkeep[ii,11:12]=c(ncp2b,pw2b) # openmx results
#
}}}
}}}
#summary(egeeM3mzdzL)$coefficients
#anova(egeeM3mzdzL, egeeM1mzdzL)
#anova(egeeM3mzdzL, egeeM2mzdzL)
#
round(reskeep[,ipow],3)
round(mxkeep[,seq(2,12,by=2)],3)
