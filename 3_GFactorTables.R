##############################################
#### Raphiel Murden
#### Create tables of G-factor by node in Power 264 system
#### 24 June 2020: 
##############################################

library(R.matlab)
library(lme4) 
library(lmerTest)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(matlab)
library(fields)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(ggpubr)
library(arsenal)
library(xtable)

# Change working directory depending on whether running on cluster or locally:
setwd('~/../../Dropbox/Multiband Study/Results/')

GFactors = readMat('Gfactor_power264_9p__9p_tf.mat')
allmb.gf.node = GFactors$allmb.gf.node
allmb.tf.gf.node = GFactors$allmb.tf.gf.node
mblist = sapply(GFactors$mblist, function(x) sub("_AP", "", sub("iso_", "", sub("rsfMRI_", "", x))))

dim(allmb.gf.node)
#############################################################################################
################################  No Temporal Filtering
#############################################################################################
median.mb.node = apply(allmb.gf.node, c(1,2), median, na.rm = T)
q1.mb.node = apply(allmb.gf.node, c(1,2), function(x) quantile(x, 0.25, na.rm = T))
q2.mb.node = apply(allmb.gf.node, c(1,2), function(x) quantile(x, 0.75, na.rm = T))

node.gf.table = as.data.frame(matrix(paste(round(median.mb.node, 2), " (", round(q1.mb.node, 2), ",", round(q2.mb.node, 2), ")", 
                                           sep = ""), ncol = 9), row.names = paste("Node", seq(264)))
colnames(node.gf.table) = mblist
node.gf.table = cbind(paste("Node", seq(264)), node.gf.table)
colnames(node.gf.table)[1] = "Node"

write.table(node.gf.table, file = "../Documents/GF Table Node-by-mb.csv", sep = ",", row.names = F, col.names = T)

median.mb.all = apply(median.mb.node, 2, median, na.rm = T)
q1.mb.all = apply(allmb.gf.node, 2, function(x) quantile(x, 0.25, na.rm = T))
q2.mb.all = apply(allmb.gf.node, 2, function(x) quantile(x, 0.75, na.rm = T))

overall.gf.table = as.data.frame(matrix(paste(round(median.mb.all, 2), " (", round(q1.mb.all, 2), ",", round(q2.mb.all, 2), ")", 
                                           sep = ""), ncol = 9))

#############################################################################################
######################################  WITH Temporal Filtering
#############################################################################################
median.mb.tf.node = apply(allmb.tf.gf.node, c(1,2), median, na.rm = T)
q1.mb.tf.node = apply(allmb.tf.gf.node, c(1,2), function(x) quantile(x, 0.25, na.rm = T))
q2.mb.tf.node = apply(allmb.tf.gf.node, c(1,2), function(x) quantile(x, 0.75, na.rm = T))

node.gf.tf.table = as.data.frame(matrix(paste(round(median.mb.tf.node, 2), " (", round(q1.mb.tf.node, 2), ",", round(q2.mb.tf.node, 2), ")", 
                                           sep = ""), ncol = 9), row.names = paste("Node", seq(264)))
colnames(node.gf.tf.table) = mblist
node.gf.tf.table = cbind(paste("Node", seq(264)), node.gf.tf.table)
colnames(node.gf.table)[1] = "Node"

write.table(node.gf.tf.table, file = "../Documents/GF Table Node-by-mb.tf.csv", sep = ",", row.names = F, col.names = T)

median.mb.tf.all = apply(median.mb.tf.node, 2, median, na.rm = T)
q1.mb.tf.all = apply(allmb.tf.gf.node, 2, function(x) quantile(x, 0.25, na.rm = T))
q2.mb.tf.all = apply(allmb.tf.gf.node, 2, function(x) quantile(x, 0.75, na.rm = T))

overall.gf.tf.table = as.data.frame(matrix(paste(round(median.mb.tf.all, 2), " (", round(q1.mb.tf.all, 2), ",", round(q2.mb.tf.all, 2), ")", 
                                              sep = ""), ncol = 9))
overall.table = rbind(overall.gf.table, overall.gf.tf.table)
rownames(overall.table) = c("9p", "9p + TFs")
colnames(overall.table) = mblist

xtable(overall.table)
                