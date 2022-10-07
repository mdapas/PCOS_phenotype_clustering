################################################
# Notes 09.29.2019                             #
################################################

# This script is designed to read in a tab-delimited text file with 
# trait values and assay method codes for the following traits by default:
# Age, BMI, Testosterone, SHBG, Insulin, Glucose, BDHEAS, LH, and FSH
#
# The default traits, column names, and other variables are defined in 
# the "Variables" section of this script.
# These variable can be modified as desired. 
#
# This clustering method log-normalizes trait distributions, then adjusts by age and assay, 
# then applies an inverse normal transformation (INT), in which the trait distributions are fit onto a 
# normal distribution. Clusters are therefore primarily driven by multicollinearity of input variables.
#
# The manuscript detailing these methods and applications has been published in PLOS Medicine: https://doi.org/10.1371/journal.pmed.1003132.
# Direct questions to Matthew Dapas at mdapas@uchicago.edu

################################################
# Load Libraries                               #
################################################

#install.packages('fpc')
#install.packages('stats')
#install.packages('gplots')
#install.packages('clusterSim')
#install.packages('factoextra')
#install.packages('FactoMineR')

library(fpc)
library(stats)
library(gplots)
library(clusterSim)
library(factoextra)
library(FactoMineR)


################################################
# Functions                                    #
################################################

rntransform <- function(formula,data,family=gaussian) {

    if ( is(try(formula,silent=TRUE),"try-error") ) {
      if ( is(data,"gwaa.data") ) data1 <- phdata(data)
      else if ( is(data,"data.frame") ) data1 <- data
      else stop("'data' must have 'gwaa.data' or 'data.frame' class")
      formula <- data1[[as(match.call()[["formula"]],"character")]]
    }

    var <- ztransform(formula,data,family)
    out <- rank(var) - 0.5
    out[is.na(var)] <- NA
    mP <- .5/max(out,na.rm=T)
    out <- out/(max(out,na.rm=T)+.5)
    out <- qnorm(out)
    out
  }

ztransform <- function(formula,data,family=gaussian) {
    if (missing(data)) {
      if(is(formula,"formula"))
        data <- environment(formula)
      else
        data <- environment()
      #		wasdata <- 0
    } else {
      if (is(data,"gwaa.data")) {
        data <- data@phdata
      }
      else if (!is(data,"data.frame")) {
        stop("data argument should be of gwaa.data or data.frame class")
      }
      #		attach(data,pos=2,warn.conflicts=FALSE)
      #		wasdata <- 1
    }

    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }

    if ( is(try(formula,silent=TRUE),"try-error") ) {
      formula <- data[[as(match.call()[["formula"]],"character")]]
    }

    if (is(formula,"formula")) {
      #		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
      mf <- model.frame(formula,data,na.action=na.pass,drop.unused.levels=TRUE)
      mids <- complete.cases(mf)
      mf <- mf[mids,]
      y <- model.response(mf)
      desmat <- model.matrix(formula,mf)
      lmf <- glm.fit(desmat,y,family=family)
      #		if (wasdata)
      #			mids <- rownames(data) %in% rownames(mf)
      #		else
      resid <- lmf$resid
      #		print(formula)
    } else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
      y <- formula
      mids <- (!is.na(y))
      y <- y[mids]
      resid <- y
      if (length(unique(resid))==1) stop("trait is monomorphic")
      if (length(unique(resid))==2) stop("trait is binary")
    } else {
      stop("formula argument must be a formula or one of (numeric, integer, double)")
    }
    y <- (resid-mean(resid))/sd(resid)
    #	if (wasdata==1) detach(data)
    tmeas <- as.logical(mids)
    out <- rep(NA,length(mids))
    out[tmeas] <- y
    out
  }


################################################
# Variables                                    #
################################################

input_file <- "~/input_file.txt" # File path to input table

## Set variables and labels
sample <- 'sample_id'; age <- 'age'; bmi <- 'bmi'; t <- 'T'; dheas <- 'dheas'
i0 <- 'i0'; g0 <- 'g0'; shbg <- 'shbg'; lh <- 'lh'; fsh <- 'fsh';
assay <- 'assay_method'
z <- 'z'; d <- '.'; u <- '_'

t.assay <- paste(t, assay, sep=u); dheas.assay <- paste(dheas, assay, sep=u); i0.assay <- paste(i0, assay, sep=u)
g0.assay <- paste(g0, assay, sep=u); shbg.assay <- paste(shbg, assay, sep=u); lh.assay <- paste(lh, assay, sep=u)
fsh.assay <- paste(fsh, assay, sep=u)

t.z <- paste(t, z, sep=d); lh.z <- paste(lh, z, sep=d); fsh.z <- paste(fsh, z, sep=d); i0.z <- paste(i0, z, sep=d)
g0.z <- paste(g0, z, sep=d); dheas.z <- paste(dheas, z, sep=d); shbg.z <- paste(shbg, z, sep=d); 

variables <- c(t, dheas, i0, g0, shbg, lh, fsh)
var_labels <- c('BMI','T', 'DHEAS', 'Ins0', 'Glu0', 'SHBG', 'LH', 'FSH')

cluster_input_cols <- c(sample, age, bmi, t, dheas, i0, g0, shbg, lh, fsh, t.assay, dheas.assay, i0.assay,
                g0.assay, shbg.assay, lh.assay, fsh.assay)
cluster_cols <- c(sample, age, bmi, t.z, dheas.z, i0.z, g0.z, shbg.z, lh.z, fsh.z)

dist_metric <- 'manhattan'
clust_method <- 'ward.D'


################################################
# Main                                         #
################################################

## Read in input file (must have columns according to cluster_input_cols)
m.df <- read.delim(input_file,na.strings=c("","#N/A", "missing"))

## Remove duplicate entries
m.df <- m.df[!m.df$sample_id %in% m.df[duplicated(m.df[,cluster_input_cols[2:length(cluster_cols)-1]], fromLast=T), sample],]

## Only keep samples with complete, non-zero quant. trait data
m.df <- m.df[complete.cases(m.df[,variables]), ]
m.df <- m.df[!rowSums(m.df[,variables]==0) >0, ]

## Remove data excluded from GWAS
#m.df <- m.df[m.df$Final_Exclusion_Filter != 1, ]

## Correct BMI for age
m.df$bmi <- log(m.df$bmi)
model <- glm(m.df$bmi ~ m.df$age, family=gaussian())
m.df$bmi <- m.df$bmi - fitted.values(model)

## Iterate through other quantitative traits
for (var in variables) {
  method <- paste(var, assay, sep=u)
  z <- paste(var,'.z',sep='')
  m.df[,method] <- as.character(m.df[,method])
  m.df[, z] <- log(m.df[, var])  # Log-normalize trait
  if (nlevels(as.factor(m.df[,method])) > 1) {
    model <- glm(m.df[,z] ~ m.df$age + m.df[,method], family=gaussian()) # Model trait, adjust for age + method
  } else {
    model <- glm(m.df[,z] ~ m.df$age, family=gaussian()) # Model trait, adjust for age + method
  }
  m.df[which(!is.na(m.df[,z])), z] <-  m.df[which(!is.na(m.df[,z])), z] - fitted.values(model)  # Get residuals
}

## Apply Inverse Normal Transformation to residuals
for (col in cluster_cols[2:length(cluster_cols)]){
  #print(c(col,shapiro.test(m.df[,col])[[2]] ))
  m.df[which(!is.na(m.df[,col])), col]  <- rntransform(as.numeric(m.df[which(!is.na(m.df[,col])), col]))
}

## Create clustering matrix
m.df <- m.df[,cluster_cols]
names(m.df) <- c(sample, age, var_labels)
cluster_matrix <- as.matrix(na.omit(m.df[, var_labels]))

## Perform hierarchical clustering 
hc <- hclust(dist((cluster_matrix), method=dist_metric), method=clust_method)
row_hc <- hclust(dist(t(df.12), method='manhattan'), method='ward.D')

## Define clusters
mycl <- cutree(hc, k=3)

## Get clustering stats
stats <- cluster.stats(dist((cluster_matrix), method=dist_metric),mycl)

## Get a color palette equal to the number of clusters
cluster_colors <- c("#5480C4FF",'grey',"#CF4D40D9")  # May want to re-order colors based on subtype

## Create vector of colors for side bar, declare colors to be used in heat map
myClusterSideBar <- cluster_colors[mycl]
colors = c(seq(-4.401,-1.201,length=40),seq(-1.2,1.2,length=30),seq(1.201,4.401,length=40))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 109)

## Draw the heat map
heatmap.2(t(cluster_matrix), main=paste("[TITLE]\nn=",dim(cluster_matrix)[1],sep=''), revC=T,
          Rowv=as.dendrogram(row_hc), Colv=as.dendrogram(hc), dendrogram="both", trace="none",col=my_palette,
          ColSideColors=myClusterSideBar, key=TRUE, key.title='', key.xlab='', margins=c(2,6), 
          labCol='',symm=F,symkey=F,symbreaks=F,breaks=colors)

#   Run clusterboot() with hclust to determine cluster stability
cboot.hclust <- clusterboot(dist((cluster_matrix), method=dist_metric), B=1000, clustermethod=disthclustCBI,
                            method=clust_method, k=3)
cboot.hclust$bootmean  # Jaccard scores for clusters

## Cluster numbers re-assigned in order to match cluster labels with number
#mean(cluster_matrix[mycl==1, 'BMI'])
#mean(cluster_matrix[mycl==2, 'BMI'])
#mean(cluster_matrix[mycl==3, 'BMI'])
#mycl[mycl==3] <- 4
#mycl[mycl==2] <- 3
#mycl[mycl==4] <- 2

## Add the cluster ID# to your data
clusters <- cbind(cluster_matrix, cluster_id=mycl)

## Merge the sample_ids
clusters <- merge(clusters, m.df, by=var_labels)

## Re-order the data frame
col_idx <- grep("cluster_id", names(clusters))
clusters <- clusters[, c(col_idx, (1:ncol(clusters))[-col_idx])]
col_idx <- grep(sample, names(clusters))
clusters <- clusters[, c(col_idx, (1:ncol(clusters))[-col_idx])]

## Write output table
#write.table(clusters, '~/cluster_output_table.txt', sep='\t', quote=F, row.names=F)


### PLOT: Reproductive vs. Metabolic box plots

## May want to re-order labels according to difference
sorted_vars <- var_labels

box_data <- list(clusters[clusters$cluster_id==1,sorted_vars[1]], clusters[clusters$cluster_id==2,sorted_vars[1]],
                 clusters[clusters$cluster_id==1,sorted_vars[2]], clusters[clusters$cluster_id==2,sorted_vars[2]],
                 clusters[clusters$cluster_id==1,sorted_vars[3]], clusters[clusters$cluster_id==2,sorted_vars[3]],
                 clusters[clusters$cluster_id==1,sorted_vars[4]], clusters[clusters$cluster_id==2,sorted_vars[4]],
                 clusters[clusters$cluster_id==1,sorted_vars[5]], clusters[clusters$cluster_id==2,sorted_vars[5]],
                 clusters[clusters$cluster_id==1,sorted_vars[6]], clusters[clusters$cluster_id==2,sorted_vars[6]],
                 clusters[clusters$cluster_id==1,sorted_vars[7]], clusters[clusters$cluster_id==2,sorted_vars[7]],
                 clusters[clusters$cluster_id==1,sorted_vars[8]], clusters[clusters$cluster_id==2,sorted_vars[8]])

fill_colors = c('#DFB6B1', '#7B96B6')  # light red, light blue
b_cols = c(rgb(0.71,.3,0.25,1), '#446689')  # dark red, dark blue
intervals = c(1.65,2.35, 3.65,4.35, 5.65,6.35, 7.65,8.35, 9.65,10.35, 11.65,12.35, 13.65,14.35, 15.65,16.35)
op <- par(mar = c(6,4,4,2) + 0.1, mgp=c(3,2,0))
boxplot(box_data, at=intervals, ylim=c(-2,2.3), ylab='Z', xaxt='n', yaxt='n', boxwex=.55,
        font.main='serif', cex.main=1.3, outline=F,staplelty=0, lwd = 2, boxlwd = 3, whisklty=0, range=1)
grid(); par(new=TRUE)
abline(h = 0, col = "indianred",lwd=1.75, lty=3); par(new=TRUE)
boxplot(box_data, at=intervals, ylim=c(-2,2.3), xaxt='n', border=b_cols, col=fill_colors, yaxt='n', boxwex=.55,
        outline=F, axes=F, whisklty =0, staplelty=0, lwd = 2, boxlwd = 3, range=1)
axis(side=2,at=seq(-2,2,0.5),labels=seq(-2,2,0.5), cex.axis=1, mgp=c(1,1,0))
lablist.x <- as.vector(sorted_vars)
text(c(2,4,6,8,10,12,14,16), par("usr")[3] - .03,labels = lablist.x, pos = 1, xpd = TRUE, font=1, cex=1.1)
legend(11,2.3,c('Metabolic','Reproductive'), fill=fill_colors, cex=.8, bg='gray98', ncol=2, text.width=1.6, x.intersp=0.5,y.intersp=.1)

#stripchart(box_data, at=intervals,vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = b_cols, cex=1.2)


### PCA PLOT ###
clusters$cluster_id <- as.factor(clusters$cluster_id)
clusters$group <- ifelse(clusters$cluster_id==1, "Metabolic", ifelse(clusters$cluster_id==2, "Reproductive", "Indeterminate"))

pca <- PCA(clusters[3:10], graph = T)
fviz_pca_biplot(pca, col.ind = clusters$group, palette = c('grey', rgb(0.81,0.3,0.25,.85), rgb(.33,.5,.77,1)),
                addEllipses = TRUE, label = "var", col.var = "black", repel = TRUE,
                legend.title = "Subtype", pointsize=2,#pointshape=19,
                xlab = paste("PC1 (", round(pca$eig[1,2],2), '%)', sep=''),
                ylab = paste("PC2 (", round(pca$eig[2,2],2), '%)', sep=''),
                title = "")
#dev.off()
