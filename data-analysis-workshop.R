
#############################################################
##### Script by Felix J. Hartmann, University of Zurich #####
#############################################################






############################
##### install packages #####
############################



source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") # upgrading Bioconductor (if you have problems with installing something)
biocLite("flowCore")    # for handling of fcs files in R
biocLite("ggplot2")     # for advanced data plotting
biocLite("Rtsne")       # for calculating tSNE
biocLite("RColorBrewer")        # additional color schemes
biocLite("reshape2")    # reorganizing data
biocLite("FlowSOM")     # FlowSOM clustering
biocLite("plyr")        # for data handling
biocLite("xlsx")
biocLite("flowWorkspace")
biocLite("ggcyto")      #plotting flowSets
biocLite("devtools")
devtools::install_github("JinmiaoChenLab/Rphenograph")






##########################
##### open libraries #####
##########################


    
library(ggplot2)
library(Rtsne)
library(RColorBrewer)
library(reshape2)
library(FlowSOM)
library(plyr)
library(xlsx)
library(flowWorkspace)





########################################
##### Importing a Workspace into R #####
########################################



# nice introductions
# http://bioconductor.org/packages/release/bioc/vignettes/flowCore/inst/doc/HowTo-flowCore.pdf
# http://bioconductor.org/packages/release/bioc/vignettes/openCyto/inst/doc/openCytoVignette.html
# https://www.bioconductor.org/help/course-materials/2014/BioC2014/OpenCytoPracticalComponent.html



# open the packages
library(flowCore)  
library(flowWorkspace)
library(ggcyto)



# specify the file locations
dataDir <- "files/ungated/fcs"
wsfile <- "files/ungated/workspace9.xml"



# read in the Workspace
ws <- openWorkspace(wsfile)
gs <- parseWorkspace(ws, path = dataDir,  sampNloc = "sampleNode", #transform = F,
                     compensation = NULL, name = 1)



# explore the gating-set 
plot(gs) #plot the hierarchy tree
getNodes(gs) #show all the cell populations(/nodes)
getPopStats(gs, format = "wide", statistic = "freq") #show the population statistics



# plotting gates
autoplot(gs[[1]], bins = 128) + theme(aspect.ratio = 1)
autoplot(gs, "cd45-live", bins = 128) + theme(aspect.ratio = 1)



# extracting a flowSet
fs <- getData(gs, "cd45-live")



# compensation for FACS data
#compMat <- getCompensationMatrices(gs[[1]])
#fs <- compensate(fs, compMat)



# exploring the data further
ggcyto(fs, aes(x = CD3, y = CD8a)) + geom_hex(bins = 128)
ggcyto(fs, aes(x = Time, y = DNA1)) + geom_hex(bins = 128)
ggcyto(fs, aes(x = CD3)) + geom_density()







#################################
##### exploring the dataset #####
#################################



# alternatively, read in fcs files directly
fs  <- read.flowSet(path = "files/gated/", pattern = ".fcs")



# exploring flowSets
fs
sampleNames(fs)
colnames(fs)
fs[[2]]
ff2 <- fs[[2]]
str(ff2)



# renaming channels
desc <- as.vector(ff2@parameters$desc)
desc
colnames(fs) <- desc
colnames(fs)
ff2 <- fs[[2]]
fs





#############################
##### data manipulation #####
#############################



# accessing the data of an fcs file
data2 <- exprs(ff2)
dim(data2)
head(data2)
class(data2)



# combine data into a matrix
data <- fsApply(fs, exprs)
head(data)
dim(data)



# make a gate.source vector
fs_rows <- fsApply(fs, nrow)[,1]
fs_rows
gate.source <- as.vector(x = NULL)
for(i in 1:length(fs_rows)) {
        temp.source <- rep(i, fs_rows[i])
        gate.source <- c(gate.source, temp.source)
}
tail(gate.source)
table(gate.source)



# add the gate.source vector to the file
data <- cbind(data, gate.source)
tail(data)



# plotting data
ggplot(data.frame(data[1:5000,]), aes(x = CD19, y = MHCII)) + 
  geom_point() + 
  theme(aspect.ratio = 1)



# transform data
asinh_scale <- 5
data.trans <- asinh(data / asinh_scale)
data.trans[,"gate.source"] <- data[,"gate.source"]



# plotting the transformed data
ggplot(data.frame(data.trans[1:5000,]), aes(x = CD19, y = MHCII)) + 
  geom_point() + 
  theme(aspect.ratio = 1)






####################################
##### percentile normalization #####
####################################

 

# percentile normalization for a random sequence of numbers
v <- 1:1000
v
quantile_value <- 0.999
quantile(v, quantile_value)



# calculating percentile for 1 vector
percentile.vector <- apply(data.trans, 2, function(x) quantile(x, quantile_value, names = F))
percentile.vector
data.trans <- t(t(data.trans) / as.numeric(percentile.vector))
data.trans[,"gate.source"] <- data[,"gate.source"]
head(data.trans)



# plotting the transformed data & normalized data
ggplot(data.frame(data.trans[1:5000,]), aes(x = CD19, y = MHCII)) + 
  geom_point() + 
  theme(aspect.ratio = 1)





#########################
##### tSNE analysis #####
#########################



# select clustering channels
colnames(data.trans)
clustering.channels <- colnames(data)[c(4:18, 20:27)]
clustering.channels



# subsample data for t-SNE
set.seed(123)
ix <- sample(1:nrow(data), 20000)
data_rtsne <- data.matrix(data.trans[ix,clustering.channels])
gate.source <- gate.source[ix]
dim(data_rtsne)
head(data_rtsne)



# run bh SNE
set.seed(123)
out_rtsne <- Rtsne(data_rtsne, dims = 2, max_iter = 1000, verbose = T, pca = F)



# prepare the tSNE data
tsne <- as.data.frame(out_rtsne$Y)
colnames(tsne) <- c("tSNE1", "tSNE2")



# plot tSNE 
ggplot(tsne, aes(x = tSNE1, y = tSNE2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1) +
  ggtitle("tSNE map") 



# prepare the expression data
data_rtsne <- cbind(data_rtsne, tsne, gate.source)
head(data_rtsne)
data_melt <- melt(data_rtsne, id.vars = c("tSNE1", "tSNE2", "gate.source"))



# define a color scheme
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                                 "yellow", "#FF7F00", "red", "#7F0000"))



# plot tSNEs with expression overlayed
ggplot(data_melt, aes(x = tSNE1, y = tSNE2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = jet.colors(100), limits = c(0,1)) +
  facet_wrap(~ variable, ncol = 8, scales = "free") +
  ggtitle("tSNE map with expression values") 






##############################
##### FlowSOM Clustering #####
##############################



# make a flowframe for FlowSOM
ff_new <- flowFrame(exprs = data.matrix(data_rtsne), desc = list(FIL = 1))
ff_new



# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clustering.channels)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]



# plot pie charts
out_fSOM <- UpdateNodeSize(out_fSOM, reset = TRUE)
FlowSOM::PlotStars(out_fSOM, view = "grid", markers = clustering.channels)
out_fSOM <- UpdateNodeSize(out_fSOM)
FlowSOM::PlotStars(out_fSOM, view = "MST", markers = clustering.channels)



# try the suggested automatic metaclustering method for a hint for k
auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 20)
max(auto_meta)



# select a k.value
choosen_k <- 12



# do a manual metaclustering
set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = choosen_k)
meta_results <- out_meta[labels]



# prepare the metaclustering data
data_meta <- cbind(tsne, meta_results, gate.source)
data_melt <- melt(data_meta, id.vars = c("tSNE1", "tSNE2", "gate.source"))



# plot tSNEs with clustering overlayed
ggplot(data_melt, aes(x = tSNE1, y = tSNE2, color = as.factor(value))) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggtitle("tSNE map with overlaid metaclustering") 






####################################
##### Other clustering methods #####
####################################



# phenograph
library(Rphenograph)
set.seed(123)
pheno_out <- Rphenograph(data_rtsne, k = 100)
p_labels <- factor(membership(pheno_out[[2]]))



# prepare the metaclustering data
data_meta_p <- cbind(tsne, p_labels, gate.source)
data_melt_p <- melt(data_meta_p, id.vars = c("tSNE1", "tSNE2", "gate.source"))



# plot tSNEs with clustering overlayed
ggplot(data_melt_p, aes(x = tSNE1, y = tSNE2, color = as.factor(value))) +
        geom_point(size = 0.05) +
        coord_fixed(ratio = 1) +
        #scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))) +
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        ggtitle("tSNE map with overlaid metaclustering") 





#####################################
##### analyse different samples #####
#####################################



# plot tSNEs for different groups
ggplot(data_melt, aes(x = tSNE1, y = tSNE2, color = as.factor(value))) +
  geom_point(size = 0.2) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))) +
  facet_wrap(~ gate.source, ncol = 3, scales = "free") + 
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggtitle("tSNE for different samples")





#########################
##### make heatmaps #####
#########################



# make combined expression matrix
ce <- cbind(data_rtsne, meta_results)



# go through all clusters and calculate mean for every channel
heat_mat <- matrix(, nrow = choosen_k, ncol = length(clustering.channels))
for(i in 1:choosen_k) {
    temp_mat <- ce[ce[,"meta_results"] == i, clustering.channels]
    heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}



# rename
rownames(heat_mat) <- paste("cluster", 1:choosen_k, sep = "")
colnames(heat_mat) <- clustering.channels  
heat_mat



# make a cluster heatmap
heatmap(heat_mat, scale = "none", Colv = T, Rowv = T,
          col = colorRampPalette(c("black", "gold"))(n = 9),
          breaks = seq(0, 1, by = 0.1111111))

 




###################################
##### show sample composition #####
###################################



# calculate the frequency of every cluster in each sample
fm <- ddply(data_meta, .(gate.source), function(x) {
        df <- data.frame(table(x$meta_results))
        prop <- (df$Freq/sum(df$Freq))*100})
colnames(fm) <- c("gate.source", paste("cluster", 1:choosen_k, sep = "."))
fm



# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- sampleNames(fs)
df_plot



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:(choosen_k+1)])



# make stacked bargraphs
qplot(samples, data = df_melt, geom = "bar", weight = value, fill = variable) 





###########################################################
##### statistical analysis of gated data/excel sheets #####
###########################################################



# read in excel sheet with data (eg exported from FlowJo)
freqs <- read.xlsx("files/frequencies.xlsx", sheetIndex = 1)
freqs
str(freqs)



# read in meta data
md <- read.xlsx("files/meta-data.xlsx", sheetIndex = 1)
md
str(md)



# make a combined dataset
dataset <- merge(freqs, md, by = "sample")
dataset



# calculate t-test
t.test(dataset$CD4 ~ dataset$genotype)
t.test(dataset$CD4[1:6], dataset$CD4[7:12])



# calculate t.tests for all frequencies in excel sheet
dataset_melt <- melt(dataset, id.vars = colnames(md))
dataset_melt
stat_out <- ddply(dataset_melt, .(variable), 
                  function(x) {c(p_value = t.test(x$value ~ x$genotype)$p.value)}
                  )
stat_out



# adjusting for multiple comparisons
stat_out$adjusted <- p.adjust(stat_out$p_value, method = "BH")
stat_out



# plot 1
ggplot(dataset_melt, aes(x = genotype, y = value)) +
        geom_boxplot() +
        facet_grid(~ variable, scales = "free")



# plot 2
ggplot(dataset_melt, aes(x = genotype, y = value)) +
        geom_boxplot() +
        facet_grid(organ ~ variable, scales = "free")



# more elaborate t-test
stat_out_2 <- ddply(dataset_melt, .(variable, organ), 
                  function(x) {c(p_value = t.test(x$value ~ x$genotype)$p.value)}
)
stat_out_2






##########################
##### Power analysis #####
##########################



# specify the things that you know or estimate
n <- 5          # sample size (per group)
delta <- 5      # difference between groups
sd <- 3         # standard deviation of the groups (same for both)
power <- NULL
        

   
# you need to specify three out of the four variables 
# and the last one will be returned
power.t.test(n = n,
             delta = delta,
             sd = sd,
             sig.level = 0.05,
             power = power,
             type = "two.sample",
             alternative = "two.sided")



# alternatively, check out the pwr package for power calculations
# for anovas and so on
# biocLite("pwr")

