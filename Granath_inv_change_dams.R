##########################################################
# Disentangling drivers of temporal 
# change of macroinvertebrate diversity in urban ponds
#
# Granath et al. 
#
# Contact: Gustaf.Granath@gmail.com
#
##########################################################

# load some general packages
library(gsheet)
library(cowplot)
library(ggplot2)
library(raster)
library(rgdal)

# Extract land cover data from locations####
library(terra)
library(sf)

# read raster layers 
grass.rast <- rast("GRAC_1518_020m_E47N40_03035_v010.tif")
forest.rast <- rast("TCCM_1518_020m_E47N40_03035_v010.tif")
imp.rast <- rast("IMCC_1518_020m_E47N40_03035_v010.tif")

# load pond locations and transform to same projection
url <- 'https://docs.google.com/spreadsheets/d/1FqxnVY_nmfik9BAO45-RvYVXzM4CtFSK/edit#gid=1833472143'
cent.raw <- as.data.frame(gsheet2tbl(url))
cent = st_as_sf(cent.raw, coords = c("Longitude","Latitude"), remove = FALSE,
                crs= 4326)
cent = st_transform(cent, crs=3035)
# set buffer zone
cent <- st_buffer(cent, 250)             

# check plots
plot(grass.rast)
plot(st_geometry(cent), cex=3,add=TRUE)
plot(forest.rast)
plot(st_geometry(cent), cex=3,add=TRUE)
plot(imp.rast)
plot(st_geometry(cent), cex=3,add=TRUE)

# Grass: get % per class inside buffer
( classes <- sort(terra::unique(grass.rast)[,1]) )
v.grass <- terra::extract(grass.rast, vect(cent))
d.grass <- as.data.frame(do.call(rbind,tapply(v.grass[,2], v.grass$ID, function(x) { 
  prop.table(table(factor(x, levels=classes)))})))
names(d.grass) <- paste0("class", names(d.grass)) 

# Forest: get % per class inside buffer
( classes <- sort(terra::unique(forest.rast)[,1]) )
v.forest <- terra::extract(forest.rast, vect(cent))
d.forest <- as.data.frame(do.call(rbind,tapply(v.forest[,2], v.forest$ID, function(x) { 
  prop.table(table(factor(x, levels=classes)))})))
names(d.forest) <- paste0("class", names(d.forest)) 

# Imp: get % per class inside buffer
( classes <- sort(terra::unique(imp.rast)[,1]) )
v.imp <- terra::extract(imp.rast, vect(cent))
d.imp <- as.data.frame(do.call(rbind,tapply(v.imp[,2], v.imp$ID, function(x) { 
  prop.table(table(factor(x, levels=classes)))})))
names(d.imp) <- paste0("class", names(d.imp)) 

# make data frame
landuse = data.frame(site = cent$Name, grass.2015 = d.grass$class10 + d.grass$class22,
                     grass.2018 = d.grass$class10 + d.grass$class11,
                     forest.2015 = d.forest$class10 + d.forest$class2,
                     fores.2018 = d.forest$class10 + d.forest$class1,
                     imp.2015 = d.imp$class10 + d.imp$class2,
                     imp.2018 = d.imp$class10 + d.imp$class1)

# Load pond environmental data ###
url <- 'https://docs.google.com/spreadsheets/d/1FqxnVY_nmfik9BAO45-RvYVXzM4CtFSK/edit#gid=484274396'
env14 <- as.data.frame(gsheet2tbl(url))
url <- 'https://docs.google.com/spreadsheets/d/1FqxnVY_nmfik9BAO45-RvYVXzM4CtFSK/edit#gid=68965734'
env19 <- as.data.frame(gsheet2tbl(url))
not_samp <- which(env19[,1] %in% "Ängsholmsdammen") 
NArows <- c(which(apply(env19,1, function (x) anyNA(x))), not_samp)

env14 <- env14[-NArows,]
env19 <- env19[-NArows,]
pairs(env14[,2:NCOL(env14)]) # no large collinearity
pairs(env19[,2:NCOL(env19)]) # no large collinearity
# C,N,P some high values.

landuse.sub <- landuse[-NArows,]


# Species richness analyses ####
url <- 'https://docs.google.com/spreadsheets/d/1FqxnVY_nmfik9BAO45-RvYVXzM4CtFSK/edit#gid=2141698541'
sp14 <- read.csv(construct_download_url(url), skip=1, nrows = 35)
sp19 <- read.csv(construct_download_url(url), skip=39)
# lestes sp recorded 2014 but other species increases 2019 and 
# lestes sp is zero in 2019. So we merge to one species.
Lestes.sp.  <- apply(cbind(sp14$Lestes.dryas, sp14$Lestes.sp.,sp14$Lestes.sponsa),1, sum)
Lestes.sp.  <- ifelse(Lestes.sp.> 1,1, Lestes.sp.)
sp14$Lestes.sp. <- Lestes.sp.
Lestes.sp.  <- apply(cbind(sp19$Lestes.dryas, sp19$Lestes.sp.,sp19$Lestes.sponsa),1, sum)
Lestes.sp.  <- ifelse(Lestes.sp.> 1,1, Lestes.sp.)
sp19$Lestes.sp. <- Lestes.sp.
lestes.remove <- which(colnames(sp14) %in% c("Lestes.dryas", "Lestes.sponsa"))
sp14 <- sp14[,-lestes.remove]
sp19 <- sp19[,-lestes.remove]

sp14 <- sp14[-NArows,]
sp19 <- sp19[-NArows,]
no.species <- which(colSums(rbind(sp14[,-1],sp19[,-1]))==0)
sp14 <- sp14[,-(no.species+1)] # +1 to add first column with pond names
sp19 <- sp19[,-(no.species+1)] # +1 to add first column with pond names

sr.14 <- rowSums(sp14[,-1])
sr.19 <- rowSums(sp19[,-1])
sr.dat <- data.frame(sr = c(sr.14,sr.19), id=factor(rep(sp14$Name,2)), 
                     time =factor(rep(c("1","2"),each=length(sr.14))))
mean(sr.14); sd(sr.14) # 10.3 species, sd 4.8
mean(sr.19); sd(sr.19) # 11.0 species, sd=5.4

#__test differences and variation in species richness####
library(lme4)
library(lmerTest)
# Is there a change in species richness over time?
var.mod = lmerTest::lmer(sr ~ time + (1|id), sr.dat)
summary(var.mod) # No. P=0.47

# Is spatial variation different between the two years
var.test(sr ~ time, sr.dat, alternative = "two.sided")
# No. P=0.53

# check variance components
var.mod = lmer(sr ~ 1 + (1|id), sr.dat, 
               control=lmerControl(check.nobs.vs.nRE  = "ignore",
                                   check.nobs.vs.nlev = "ignore"))
summary(var.mod)
VarCorr(var.mod)$id[1]/(VarCorr(var.mod)$id[1] + attributes(VarCorr(var.mod))$sc^2)
attributes(VarCorr(var.mod))$sc^2/(VarCorr(var.mod)$id[1] + attributes(VarCorr(var.mod))$sc^2)


#__species 2014 and 2019####
sum(apply(sp14[,-1], 2, sum)>0) # 92 species 2014
sum(apply(sp19[,-1], 2, sum)>0) # 82 species 2019
#___overlap
sp.names.14 <- colnames(sp14[,-1])[apply(sp14[,-1], 2, sum)>0]
sp.names.19 <- colnames(sp19[,-1])[apply(sp19[,-1], 2, sum)>0]
sum(!(sp.names.14 %in% sp.names.19)) # 37 species unique in 2014
sum(!(sp.names.19 %in% sp.names.14)) # 27 species unique in 2019
sum(!(sp.names.14 %in% sp.names.19)) / sum(apply(sp14[,-1], 2, sum)>0) 
sum(!(sp.names.19 %in% sp.names.14)) / sum(apply(sp19[,-1], 2, sum)>0) 
# 40% of the species in 2014 were only found in 2014
# 33% of the species in 2019 were only found in 2019

# How many determined to species level?
all.names = c(sp.names.14[which(!(sp.names.14 %in% sp.names.19))], sp.names.19)
no.sps = all.names[-grep("sp", all.names)]
no.sps[grep(".", no.sps, fixed=TRUE)] # minus puella.pulchellum

# Look at the most common species
ncol(sp14[,-1]) #number of species
sum((colSums(sp14[,-1])>0) + (colSums(sp19[,-1])>0)==2) # species found both years

    # first total abundance
comm.sp <- sp14[,-1]+sp19[,-1]
tot.sum = apply(comm.sp, 2, sum)
tot.sum <- data.frame(records=tot.sum, species=
                                       names(tot.sum))
tot.sum$species = with(tot.sum, reorder(species, records, mean,
                                              decreasing = T))
sp.names.sums = as.character(tot.sum[rev(order(tot.sum$records)),
                                               "species"])

ggplot(tot.sum[rev(order(tot.sum$records))[1:20],], 
       aes(y=records, x=species, fill="records")) +
  geom_col(color="gray", fill="gray") +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=unit(c(1,1,1,2),'lines')
  ) +
  scale_x_discrete(labels = c(sp.names.sums)) +
  labs(x ="", y="Number of records")

# here repeatedly observed species (at same pond both years)
sites.found.twice <- apply(comm.sp, 2, function (x) sum(x==2) )
plot.sites.found.twice <- data.frame(sites=sites.found.twice, species=
                                       names(sites.found.twice))
plot.sites.found.twice$species = with(plot.sites.found.twice, 
                                      reorder(species, sites, mean,
                                      decreasing = T))
sp.names = as.character(plot.sites.found.twice[rev(order(plot.sites.found.twice$sites))[1:10],
                                  "species"])

#__plot Figure 3####
ggplot(plot.sites.found.twice[rev(order(plot.sites.found.twice$sites))[1:10],], 
       aes(y=sites, x=species, fill="sites")) +
  geom_col(color="gray", fill="gray") +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=unit(c(1,1,0,2),'lines')
        ) +
  scale_x_discrete(labels = c(sp.names)) +
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  labs(x ="", y="Number of sites")

# Pond area data ####
url <- 'https://docs.google.com/spreadsheets/d/1FqxnVY_nmfik9BAO45-RvYVXzM4CtFSK/edit#gid=1881497525'
area <- as.data.frame(gsheet2tbl(url))
area <- area[-NArows,1:2]
summary(lm(sr.14 ~ area[,2])) # no effect of dam area on species richness

# Temporal beta diversity (TBI) analyses####
#___TBI - env change ####
library(adespatial)

# scale and fix env variables
env.dat.fix <- function (mat1,mat2, non.neg=TRUE) {
n12 <- dim(mat1)[1]
tmp <- scale(rbind(mat1,mat2))
if(non.neg) tmp <- tmp - min(tmp)
mat1 <- tmp[1:n12,]
mat2 <- tmp[(n12+1):(2*n12),]
list(mat1=mat1, mat2=mat2)
}
pond.env.change = env.dat.fix(env14[,-1], env19[,-1])
tbi.env = TBI(pond.env.change[[1]], pond.env.change[[2]], method="euclidean")

#___TBI - land change ####
# make a 2015 matrix an a 2018 matrix to get TBI for land change
# forest no change so not included
landuse.sub <- landuse[-NArows,]
lc.2015 = landuse.sub[,c(2,6)] 
lc.2018 = landuse.sub[,c(3,7)]

lc.change = env.dat.fix(as.matrix(lc.2015), as.matrix(lc.2018))
tbi.lc = TBI(lc.change[[1]], lc.change[[2]], method="euclidean")
#___Plot TBI indices####
plot(x=tbi.lc[[1]], y=tbi.env[[1]], ylab="Env change", xlab="land change")
cor.test(tbi.lc[[1]], tbi.env[[1]])

plot(landuse.sub$grass.2015- landuse.sub$grass.2018, tbi.env[[1]])
plot(landuse.sub$imp.2015- landuse.sub$imp.2018, tbi.env[[1]])

#___TBI - beta div ####
pond.sp.change = env.dat.fix(sp14[,-1], sp19[,-1])
tbi.sp = TBI(sp14[,-1], sp19[,-1], method="jaccard", nperm=9999, test.t.perm=TRUE)
tbi.sp # no difference in species loss or gain over time (41% losses, 44% gain)
tbi.sp$p.TBI[tbi.sp$p.TBI<0.05] # barely below 0.05
tbi.sp$TBI[]
hist(tbi.sp$TBI)
# no sites are statistically very different

# Temporal versus spatial beta diversity####
# turnover, nested and total temporal beta diversity
library(betapart)
beta.temp <- beta.temp(sp14[,-1], sp19[,-1], index.family="jaccard")
mean(beta.temp[,1], na.rm=T) #turnover
mean(beta.temp[,2], na.rm=T) #nested
mean(beta.temp[,3]) # total

beta2014 <- beta.multi(sp14[,-1], index.family="jaccard")

# spatial beta diversity
jac.2014 <- beta.pair(sp14[,-1], index.family="jac")
jac.2014.turn <- as.matrix(jac.2014$beta.jtu)
jac.2014.nest <- as.matrix(jac.2014$beta.jne)
jac.2014.tot <- as.matrix(jac.2014$beta.jac)
# turnover
dis.site.turn =vector()
for (i in 1:NCOL(jac.2014.turn)) {
  dis.site.turn[i] <- mean(c(jac.2014.turn[i,], jac.2014.turn[,i]), na.rm=TRUE)
}
mean(dis.site.turn);sd(dis.site.turn)
mean(beta.temp[,1]);sd(beta.temp[,1])

t.test(beta.temp[,1], dis.site.turn)
# nested
dis.site.nest =vector()
for (i in 1:NCOL(jac.2014.nest)) {
  dis.site.nest[i] <- mean(c(jac.2014.nest[i,], jac.2014.nest[,i]), na.rm=TRUE)
}
mean(dis.site.nest);sd(dis.site.nest)
mean(beta.temp[,2]);sd(beta.temp[,2])

t.test(beta.temp[,2], dis.site.nest)
# total
dis.site.tot =vector()
for (i in 1:NCOL(jac.2014.tot)) {
  dis.site.tot[i] <- mean(c(jac.2014.tot[i,], jac.2014.tot[,i]), na.rm=TRUE)
}
mean(dis.site.tot);sd(dis.site.tot)
mean(beta.temp[,3]);sd(beta.temp[,3])

t.test(beta.temp[,3], dis.site.tot)

#___Fig 4 to compare spatial-temporal beta diversity####
spa.temp <- data.frame(turnover = c(dis.site.turn,beta.temp[,1]),
           nestedness = c(dis.site.nest,beta.temp[,2]),
           total = c(dis.site.tot,beta.temp[,3]),
           Type = rep(c("spatial", "temporal"), each = 30))
spa.temp <- data.frame(index = c(dis.site.turn,beta.temp[,1],
                                    dis.site.nest,beta.temp[,2],
                                    dis.site.tot,beta.temp[,3]),
                       component = rep(c("turnover", "nestedness",
                                         "total"), each = 60),
                       Type = rep(c("spatial", "temporal"),3, each = 30))
spa.temp$component <- relevel(factor(spa.temp$component), ref = "turnover")
 
ggplot(spa.temp, aes(y=index, x=component, fill=Type)) +
  geom_boxplot() +
  guides(colour=guide_legend(title=NULL), 
         fill=guide_legend(title=NULL), 
         shape="none") +
  theme_cowplot(12) +
  theme(#legend.position = c(0.35,0.9), 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.key.width = unit(1.2,"cm")) +
  labs(x="", y="Dissimilarity") +
  scale_fill_manual(values=c("white", "gray66"))

# comparing beta diversity between years
beta2014 <- beta.multi(sp14[,-1], index.family="jaccard")
beta2019 <- beta.multi(sp19[,-1], index.family="jaccard")
cbind(beta2014,beta2019)

# species lost 2014-2019
betapart.core(rbind(sp14[,-1],sp19[,-1]))


# Spatial autocorrelation ####
cent.sf.sub <- cent[-NArows,]
cent.sf.sub$tbi <- tbi.sp$TBI
ggplot(cent.sf.sub, aes(col = tbi, size = tbi)) +
  geom_sf() +
  scale_color_gradient2()

# look at ponds with most gains or most loss of species
gain.loss = tbi.sp$BCD.mat[,1]/tbi.sp$BCD.mat[,2]
cent.sf.sub$tbi.gain.loss <- ifelse(gain.loss>1,"gain", "loss")
ggplot(cent.sf.sub, aes(col = tbi.gain.loss, size = tbi)) +
  geom_sf() #

# make distance matrix
spa.dists <- as.matrix(dist(cbind(cent.sf.sub$Longitude, cent.sf.sub$Latitude)))
spa.dists.inv <- 1/spa.dists
diag(spa.dists.inv) <- 0

#__test for spat autocorr####
library(ape)
# check beta div and species gains/losses
Moran.I(cent.sf.sub$tbi, spa.dists.inv)
Moran.I(tbi.sp$BCD.mat[,1], spa.dists.inv)
Moran.I(tbi.sp$BCD.mat[,2], spa.dists.inv)
Moran.I(sp.change, spa.dists.inv)
# No evidence.

# check residuals
Moran.I(resid(lm(tbi.sp ~ tbi.env + tbi.lc + log(area), sem.dat)), 
        spa.dists.inv)
# some evidence beta diversity
Moran.I(resid(lm(sp.change ~ sr.14 + tbi.env + tbi.lc + log(area), sem.dat)), 
        spa.dists.inv)
# Not for richness

# SEMs####
#___species richness####
library(piecewiseSEM)
sp.change <- sr.14 - sr.19
sem.dat <- data.frame(sp.change = sp.change, sr.14 = sr.14, sr.19 = sr.19,
                      tbi.sp = tbi.sp$TBI, area=area[,2],
                      tbi.lc = tbi.lc$TBI, tbi.env = tbi.env$TBI,
                      imp = landuse.sub$imp.2015- landuse.sub$imp.2018,
                      grass = landuse.sub$grass.2015 - landuse.sub$grass.2018,
                      tbi.b = tbi.sp$BCD.mat[,1], tbi.c = tbi.sp$BCD.mat[,2])
model.sr <- psem(
  lm( tbi.env ~ tbi.lc + log(area), sem.dat), 
  lm(sr.19 ~ sr.14 + tbi.env + tbi.lc+ log(area), sem.dat))
summary(model.sr, .progressBar = F)
write.table(summary(model.sr, .progressBar = F)$coefficients, "sem_sr.txt")

plot(model.sr, show="std",add_edge_label_spaces=TRUE)  
coefs(model.sr, standardize = "scale")

#___temporal beta####
# total
model.tempbeta <- psem(
  lm( tbi.env ~ tbi.lc + log(area), sem.dat), 
  lm(tbi.sp ~  tbi.env + tbi.lc + log(area), sem.dat))
summary(model.tempbeta, .progressBar = F)
write.table(summary(model.tempbeta, .progressBar = F)$coefficients, "sem_tbeta_tot.txt")

# species gain
model.tempbeta.gain <- psem(
  lm( tbi.env ~ tbi.lc + log(area), sem.dat), 
  lm(tbi.c ~  tbi.env + tbi.lc + log(area), sem.dat))
summary(model.tempbeta.gain, .progressBar = F)
write.table(summary(model.tempbeta.gain, .progressBar = F)$coefficients, "sem_tbeta_gain.txt")

#species loss
model.tempbeta.loss <- psem(
  lm( tbi.env ~ tbi.lc + log(area), sem.dat), 
  lm(tbi.b ~  tbi.env + tbi.lc + log(area), sem.dat))
summary(model.tempbeta.loss, .progressBar = F)
write.table(summary(model.tempbeta.loss, .progressBar = F)$coefficients, "sem_tbeta_loss.txt")

# compare gains(c) and loss(b)
plot(tbi.c ~ tbi.b, sem.dat)
abline(a=0,b=1)

plot(model.tempbeta, show="std",add_edge_label_spaces=TRUE)  
coefs(model.tempbeta, standardize = "scale")


#Rarefaction curves####
# Rarefaction is done with the iNEXT package 
sp.mat.dat <- rbind(sp14,sp19) 
time =factor(rep(c("2014","2019"),each=nrow(sp14)))
sp.mer <- sp14[,-1] + sp19[,-1]
sp.mer[sp.mer>1] <- 1

#____Rarefaction curves. Fig 2####
# A curve is done for each year, and both
# For sample based rarefaction we make pres/absence data
sp.list <- list()
for (i in 1:length(levels(time))){
  sp.list[[i]] <- t(sp.mat.dat[time == unique(time)[i],-1])
}
sp.list[[3]] <- t(sp.mat.dat[,-1])
sp.list[[4]] <- t(sp.mer)
names(sp.list) <- c(as.character(unique(time)),"all", "merged")

# Run sample-based rarefaction up to 70 ponds
library(iNEXT)
out.raw <- iNEXT(sp.list, q = 0, datatype="incidence_raw", endpoint=70)

ggiNEXT(out.raw, type=1, color.var="Assemblage") +
  theme_bw(base_size = 18) #+
  theme(legend.position="None") 

# Fix data frame and add labels etc
plot.dat <- fortify(out.raw, type=1)
plot.dat$Time <- c(rep(c("2014", "2019"),
                        each=40),rep("Both", 30))
data.sub <- plot.dat[which(plot.dat$Method=="Observed"),]

plot.dat$Method[plot.dat$Method=="Observed"]="Rarefaction"
plot.dat$lty <- factor(plot.dat$Method, levels = c("Rarefaction", "Extrapolation"))

# Make plot
g <- ggplot(plot.dat, aes_string(x="x", y="y", colour="Time")) + 
  geom_point(aes_string(shape="Time"), size=5, data=data.sub) +
  scale_colour_manual(values=c("red", "black", "blue"))

g <- g + geom_line(aes_string(linetype="lty"), lwd=1.5) +
  guides(linetype="none",
         colour=guide_legend(title=NULL), 
         fill=guide_legend(title=NULL), 
         shape="none") +
  theme_cowplot(12) +
  theme(legend.position = c(0.35,0.9), 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.key.width = unit(1.2,"cm")) +
  labs(x="Number of sampled ponds", y="Species richness") +
  geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", 
                         fill="factor(Time)", colour="NULL"), alpha=0.2) +
  scale_fill_manual(values=c("lightpink2", "gray66","blue")) +
  scale_x_continuous(breaks = seq(5,70,5)) +
  scale_y_continuous(breaks = seq(10,170,10)) #+
g

#Multi species analyses####
library(mvabund)
library(permute)

#___mvabund####
sp.mat.dat <- rbind(sp14,sp19) 

# We focus on more common species here so
# we remove species that occurred only 5 times (or less) across years 
rare <- which(apply(sp.mat.dat[,-1], 2, sum) < 6)
sp.mat.dat <- sp.mat.dat[,-c(rare+1)]

sub.env <- data.frame(time=factor(rep(c("2014","2019"),each=nrow(sp14))),
                      #tbi.env = tbi.env$TBI, tbi.lc = tbi.lc$TBI,
                      #tbi.env = tbi.env$TBI, 
                      imp = c(landuse.sub$imp.2015,landuse.sub$imp.2018),
                      grass = c(landuse.sub$grass.2015,landuse.sub$grass.2018),
                      id = factor(rep(1:nrow(sp14),2)))

#__land covef changes####
mean(landuse.sub$forest.2015 - landuse.sub$fores.2018)
range(landuse.sub$forest.2015 - landuse.sub$fores.2018)
imp.c = landuse.sub$imp.2015 - landuse.sub$imp.2018
grass.c = landuse.sub$grass.2015 - landuse.sub$grass.2018
mean(imp.c); range(imp.c)
mean(grass.c); range(grass.c)
#____test time effect with mvabund####
sp.mat <- as.data.frame(sp.mat.dat[,-1])
abu <- mvabund(sp.mat)

# Control for between-year dependence in permutations
permID=shuffleSet(n=nrow(abu),nset=2000,control=how(block=sub.env$id))

mod_time <- manyglm(abu ~ factor(time), 
                    data = sub.env, family = binomial("cloglog"))
# Overall test
anova(mod_time, test = "score", bootID = permID,cor.type = "shrink")

#___specific species change####
summary.manyglm(mod_time, cor.type = "shrink", test = "score",
                bootID = permID,p.uni = "unadjusted")
summary.manyglm(mod_time, cor.type = "shrink", test = "score",
                bootID = permID,p.uni = "adjusted")

#Species colon-extinc prob####
#__VGAM####
library(VGAM)
sp.mat.dat.14 <- sp14[, -1]
sp.mat.dat.19 <- sp19[,-1]

rare <- which(apply(rbind(sp.mat.dat.14,sp.mat.dat.19), 2, sum) < 6)
sp.mat.dat.14 <- sp.mat.dat.14[,-c(rare)]
sp.mat.dat.19 <- sp.mat.dat.19[, -c(rare)]

sp.names = colnames(sp.mat.dat.14)
res =as.data.frame(matrix(nrow=length(sp.names),ncol=9))
colnames(res) <- c("names","ps", "odd", "col.ebcom", "ext.ebcom",
                "00","10","01","11") # row-column

for (i in 1:length(sp.names)) {
  spname = sp.names[i]
  ebcom = vglm(cbind(sp.mat.dat.14[,spname], sp.mat.dat.19[,spname]) ~ 1,
                    family=binom2.or(lmu="clogloglink", exchangeable=T),
                    trace=TRUE)
  
  ps=clogloglink(coef(ebcom, mat=T)[1,"clogloglink(mu1)"],
                 inv=TRUE)

  odd = loglink(coef(ebcom, mat=T)[1,"loglink(oratio)"],
              inv=TRUE)

  #colonization prob, Equ 13 in the 2009 paper
  col.ebcom = probColonization(ps, odd)
    #(sqrt( 1+4*(odd-1)*ps*(1-ps))-1) / (2*(odd-1)*(1-ps))
  
  #extinction prob, Equ 14 in the 2009 paper
  ext.ebcom = probExtinction(ps, odd)
    #(sqrt( 1+4*(odd-1)*ps*(1-ps))-1) / (2*(odd-1)*ps)

  
  res[i,1] <- spname
  res[i,2:9] <-  c(round(c(ps, odd, col.ebcom, ext.ebcom
                ),2),
  c(table(factor(sp.mat.dat.14[,spname], levels = 0:1), 
          factor(sp.mat.dat.19[,spname], levels = 0:1) )))#/ nrow(sp.mat.dat.19) 
  print(i)

}
res

# Check species
res[which(res$ext.ebcom==1),]
res[which(res$'01'==0),]
res[which(res$odd>1.01),]
res[which(res$ext.ebcom==1),]
plot(res$'00', res$odd)

# 3 species with low extinction but still rare
res[which(res$ext.ebcom<0.5 & res$ps < 0.3),]

colo = ggplot(res, aes(x=ps, y=col.ebcom, color=names))+
  geom_abline(slope=1, linetype = "dashed") +
  geom_point(legend=FALSE) +
  #lims(y=c(0,1), x=c(0,1)) +
  labs(y="Local colonization probability")+
  theme_cowplot() +
  theme(plot.margin=unit(c(1,1,0,2),'lines')) +
  guides(colour="none") + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,1),
                     breaks=seq(0,1,0.1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1),
                     breaks=seq(0,1,0.1))
ext = ggplot(res, aes(x=ps, y=ext.ebcom, color=names))+
  geom_point() +
  geom_abline(slope=-1, intercept=1, linetype="dashed")+
  labs(y="Local extinction probability")+
  theme_cowplot() +
  theme(plot.margin=unit(c(1,1,0,2),'lines')) +
  guides(colour="none") + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,1),
                     breaks=seq(0,1,0.1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1),
                     breaks=seq(0,1,0.1))

plot_grid(colo, ext, ncol=2,labels=c("A", "B"))


#___Functions for colonization-extinction calculations####
probColonization = function(prob, oratio=1) {
  LLL = max(length(prob), length(oratio))
  prob = rep(prob, len=LLL)
  ans = oratio = rep(oratio, len=LLL)
  index1 = (abs(oratio - 1) < 0.0001)
  ans[index1] = prob[index1]
  ans[!index1] = (sqrt(1 + 4*(oratio[!index1] - 1) * prob[!index1] *
                         (1-prob[!index1])) - 1) *
    (1 / (2 * (oratio[!index1] - 1) * (1-prob[!index1])))
  ans
}

probExtinction = function(prob, oratio=1) {
  LLL = max(length(prob), length(oratio))
  prob = rep(prob, len=LLL)
  ans = oratio = rep(oratio, len=LLL)
  index1 = (abs(oratio - 1) < 0.0001)
  ans[index1] = 1 - prob[index1]
  ans[!index1] = (sqrt(1 + 4*(oratio[!index1] - 1) * prob[!index1] *
                         (1-prob[!index1])) - 1) *
    (1 / (2 * (oratio[!index1] - 1) * prob[!index1]))
  ans
}
