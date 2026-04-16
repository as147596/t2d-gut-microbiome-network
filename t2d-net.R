library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(MicrobiotaProcess)
library(tidytree)
library(ggstar)
library(forcats)
library(curatedMetagenomicData)
library(NetCoMi)
library(tidyverse)
library(tidyr)
library(reshape2)
library(ggraph)
library(igraph)
meta_d<-sampleMetadata |>
  #filter(age >= 18) |>
  filter(body_site == "stool") |>
  filter(age_category=="adult")|>
  filter(study_name=="QinJ_2012")|>
  filter(antibiotics_current_use=="no")|>
  filter(country=="CHN")|>
  filter(disease == "T2D")

studies<-unique(meta_d$study_name)
meta_h<-sampleMetadata |>
  filter(study_name %in%studies) |>
  filter(age_category=="adult")|>
  filter(study_name=="QinJ_2012")|>
  filter(antibiotics_current_use=="no")|>
  filter(body_site == "stool") |>
  filter(disease == "healthy")

meta_all<-rbind(meta_d,meta_h)
write.csv(meta_all,"meta_t2d.csv")
T2D_data<-meta_all |> returnSamples("relative_abundance", rownames = "short")

mpse<-as.MPSE(T2D_data)
mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e5
mpse %<>% mp_rrarefy()
saveRDS(mpse,"T2D.rds")



mpse %<>% 
  mp_cal_alpha(.abundance=RareAbundance)

p1<-mpse %>% 
  mp_plot_alpha(
    .group=disease, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
  ) +
  scale_fill_manual(values=c("#9FC3E2", "#F9D5B2"), guide="none") +
  scale_color_manual(values=c("#9FC3E2", "#F9D5B2"), guide="none")
print(p1)
ggsave("alpha.pdf",plot = p1,width = 8,height = 3)
mpse %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=disease
  )


p3<-mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=disease,
    taxa.class = phylum,
    topn = 20,
    plot.group = TRUE
  )
ggsave("abund_bar_phylum.pdf",plot = p3)

p4<-mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=disease,
    taxa.class = genus,
    topn = 20,
    plot.group = TRUE
  )
ggsave("abund_bar_genus.pdf",width = 10,height = 6,plot = p4)


p5<-mpse %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = disease,
    taxa.class = class,
    relative = TRUE,
    topn = 20,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )


mpse %<>% 
  mp_decostand(.abundance=Abundance)
mpse %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")


p6<-mpse %>% mp_plot_dist(.distmethod = bray,.group = disease)

p7<-p6 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("orange", "deepskyblue"), 
    guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
    )
  ), 
  aes_var = time # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 3),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )


mpse %<>%
  mp_adonis(.abundance=hellinger, .formula=~disease, distmethod="bray", permutations=9999, action="add")
adonis<-mpse %>% mp_extract_internal_attr(name=adonis)

mpse %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")

adonis <-paste0("adonis R2: ",round(adonis$R2,2),"; P-value: ", adonis$`Pr(>F)`)

p8 <- mpse %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = disease, 
    .color = disease, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#9FC3E2", "#F9D5B2"), 
    #guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#9FC3E2", "#F9D5B2"),
    #guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    #guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )+labs(subtitle = adonis)+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
ggsave("beta.pdf",plot = p8,width = 4.5,height = 4.5)




mpse %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = disease,
    first.test.alpha = 0.01
  )

taxa.tree <- mpse %>% 
  mp_extract_tree(type="taxatree")


p1 <- ggtree(
  taxa.tree,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  )
# display the high light of phylum clade.
p2 <- p1 +
  geom_hilight(
    data = td_filter(nodeClass == "phylum"),
    mapping = aes(node = node, fill = label)
  )
# display the relative abundance of features(OTU)
p3 <- p2 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, disease, .fun=min),
      size = RelRareAbundanceBySample,
      fill = disease,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))
# display the tip labels of taxa tree
p4 <- p3 + geom_tiplab(size=2, offset=7.2)
# display the LDA of significant OTU.
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_disease,
      subset = !is.na(LDAmean)
    ),
    orientation = "y",
    offset = 0.3,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  )

# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p4 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_disease)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_disease,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))

p7<-p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)
ggsave("T2D/cir_plot.pdf",width = 10,height = 10,plot=p7)


p10<-mpse %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.1
  ) +
  scale_fill_manual(values=c("deepskyblue", "orange")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("deepskyblue", "orange")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999"
    )
  )
ggsave("T2D/tree_diff_group.pdf",plot=p10,width = 10,height = 10)


p11<-mpse %>%
  mp_plot_diff_cladogram(
    label.size = 2.5,
    hilight.alpha = .3,
    bg.tree.size = .5,
    bg.point.size = 2,
    bg.point.stroke = .25
  ) +
  scale_fill_diff_cladogram( # set the color of different group.
    values = c('deepskyblue', 'orange')
  ) +
  scale_size_continuous(range = c(1, 4))+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))
ggsave("tree_diff.pdf",plot = p11)

p12<-mpse %>%
  mp_plot_diff_boxplot(
    .group = disease,
  ) %>%
  set_diff_boxplot_color(
    values = c("#9FC3E2", "#F9D5B2"),
    guide = guide_legend(title=NULL,size=20)
  )
ggsave("diff_bar.pdf",plot = p12)

library(microbiomeMarker)
phy<-as.phyloseq(mpse)
colnames(phy@tax_table)<-c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
mm_lefse <- run_lefse(
  phy,
  wilcoxon_cutoff = 0.05,
  group = "disease",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2,
  taxa_rank = "Genus"
)
mm_lefse@marker_table$lda<-ifelse(mm_lefse@marker_table$enrich_group=="healthy",
                                  -mm_lefse@marker_table$ef_lda,
                                  mm_lefse@marker_table$ef_lda)
marker_table<-mm_lefse@marker_table
lefse_bar<-ggplot(marker_table,aes(x=lda,y=reorder(feature,lda),fill=enrich_group))+
  geom_col()+
  theme_test()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        axis.title = element_text(size=14,face = "bold"),
        axis.text = element_text(size=10),
        plot.margin = margin(t=5,b=15,l=40,r=20),
        legend.position = "right",
        plot.title = element_text(hjust = -0.14,vjust = 0),
        title = element_text(size = 16,face="bold"))+
  geom_text(data = marker_table[marker_table$enrich_group == "healthy",], aes(y = feature, x = 0.1, label = feature),
            hjust = 0, size = 3) +
  geom_text(data = marker_table[marker_table$enrich_group == "T2D",], aes(y = feature, x = -0.1, label = feature),
            hjust = 1, size = 3)+expand_limits(x=c(-6,5))+
  labs(x="LDA score",y="")+
  scale_fill_manual(values=c(T2D="#F9D5B2",healthy="#9FC3E2"))
lefse_bar
ggsave("genus_dif.pdf",width = 6,height = 5.5)

########network

T2D_data<-meta_all |> returnSamples("relative_abundance")

ab<-T2D_data@assays@data@listData[["relative_abundance"]]

taxa<-sapply(rownames(ab), function(x){
  tmp<-strsplit(x,"\\|")[[1]]
  tmp[grep("g__",tmp)]
})
ab<-as.data.frame(ab)
ab$taxa<-taxa
ab<-aggregate(.~taxa,data=ab,sum)
rownames(ab)<-ab$taxa
ab_genus<-ab[,-1]

keep<-apply(ab_genus, 1, function(x){
  mean(x!=0)>0.25
})
sum(keep)
ab_genus<-ab_genus[keep,]
rownames(ab_genus)<-gsub("g__","",rownames(ab_genus))
net <- netConstruct(data = t(ab_genus)[meta_all$disease=="T2D",],data2 = t(ab_genus)[meta_all$disease=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta_all$disease=="healthy")),
                    measure = "sparcc",
                    #measurePar = list(nlambda=10, 
                    #                  rep.num=10,
                    #                  Rmethod = "approx"),
                    normMethod = "clr", 
                    zeroMethod = "alrEM",
                    sparsMethod = "threshold",
                    thresh = 0.05,
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)
thresh=(sum(abs(net[["assoMat1"]])+abs(net[["assoMat2"]]))-2*nrow(net[["assoMat2"]]))/(2*nrow(net[["assoMat2"]])*(nrow(net[["assoMat2"]])-1))
net <- netConstruct(data = t(ab_genus)[meta_all$disease=="T2D",],data2 = t(ab_genus)[meta_all$disease=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0.1),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta_all$disease=="healthy")),
                    measure = "sparcc",
                    #measurePar = list(nlambda=10, 
                    #                  rep.num=10,
                    #                  Rmethod = "approx"),
                    normMethod = "clr", 
                    zeroMethod = "alrEM",
                    sparsMethod = "threshold",
                    thresh = 0.05,
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)

props <- netAnalyze(net, 
                    centrLCC = FALSE,
                    avDissIgnoreInf = TRUE,
                    sPathNorm = FALSE,
                    clustMethod = "cluster_fast_greedy",
                    hubPar = c("degree", "eigenvector"),
                    hubQuant = 0.9,
                    lnormFit = TRUE,
                    normDeg = FALSE,
                    normBetw = FALSE,
                    normClose = FALSE,
                    normEigen = FALSE)
pdf("network_compare.pdf",width=12,height = 8)
plot(props, 
     sameLayout = TRUE, 
     layoutGroup = 1,
     rmSingles = "inboth", 
     nodeSize = "mclr", 
     labelScale = FALSE,
     cexNodes = 0.5, 
     cexLabels = 0.7,
     cexHubLabels = 0.8,
     cexTitle = 1,
     # mar = c(20,10,10,20),
     groupNames = c("T2D","healthy"),
     hubBorderCol  = "gray40")

dev.off()
legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)

comp <- netCompare(props, 
                   permTest = FALSE, 
                   verbose = FALSE,
                   seed = 123456)

summ_res<-summary(comp, 
                  groupNames = c("T2D","healthy"),pAdjust=T,
                  #showCentr = c("degree", "between", "closeness"), 
                  showCentr="all",
                  numbNodes = 200)
degree_dif<-summ_res$topProps$topDeg[1:20,1:2]
degree_dif<-reshape2::melt(as.matrix(degree_dif))
degree_dif$value[degree_dif$Var2=="healthy"]<- -degree_dif$value[degree_dif$Var2=="healthy"]
degree_dif$Var1<-factor(degree_dif$Var1,levels = rev(levels(degree_dif$Var1)))
ggplot(degree_dif,aes(value,y=Var1,fill=Var2))+
  geom_col()+
  annotate("text",x=24,y=levels(degree_dif$Var1),
           label=rev(summ_res$topProps$topDeg$abs.diff.[1:20]),
           color="red")+
  annotate("text",x=24,y=21.5,label="abs.dif",fontface = "bold",size=4)+
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.1)))+
  geom_text(aes(label=abs(value)))+
  ggthemes::theme_base()+
  scale_fill_manual(name="",values = c(T2D="#F9D5B2",healthy="#9FC3E2"))+
  xlim(c(-20,28))+
  theme(plot.background = element_blank(),
        legend.text = element_text(size=12))+
  scale_x_continuous(limits = c(-20,28),breaks = seq(-20,25,10),labels = abs(seq(-20,25,10)))+
  labs(x="",y="")
ggsave("degree.pdf")


diff <- diffnet(net,discordThresh=0.8,alpha = 0.05,lfdrThresh = 0.2,
                diffMethod = "permute", nPerm = 1000,cores = 15,
                storeCountsPerm = TRUE, 
                fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                storeAssoPerm = TRUE,
                fileStoreAssoPerm = "assoPerm",
                adjust = "none",seed = 123456)
#saveRDS(diff,"specie_network_res.rds")
diff<-readRDS("specie_network_res.rds")

plot(diff,
     cexNodes = 2, 
     cexLegend = 0.8,
     cexTitle = 2,
     mar = c(10,10,10,20),layout=NULL,
     legendGroupnames = c("T2D", "healthy"),
     legendPos = c(1.1,0.8))

props_pears <- netAnalyze(net, 
                          clustMethod = "cluster_fast_greedy",
                          weightDeg = TRUE,
                          normDeg = FALSE,
                          gcmHeat = FALSE)


diffmat_sums <- rowSums(diff$diffAdjustMat)
diff_asso_names <- names(diffmat_sums[diffmat_sums > 0.25])
pdf("dif_net.pdf",width = 10,height = 5.5)
plot(props_pears, 
     nodeFilter = "names",
     nodeFilterPar = diff_asso_names,
     nodeColor = "gray",
     highlightHubs = FALSE,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "clr",
     edgeTranspHigh = 20,
     labelScale = FALSE,
     cexNodes = 1, 
     cexLabels = 0.8,
     cexTitle = 0.8,
     groupNames = c("T2D", "healthy"),
     hubBorderCol  = "gray40")
legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.1, cex = 1.5, lty = 1, lwd = 1, 
       bty = "n", horiz = TRUE)
dev.off()

n1<-net$edgelist1
n1<-n1[n1$v1%in%diff_asso_names|n1$v2%in%diff_asso_names,]
n2<-net$edgelist2
n2<-n2[n2$v1%in%diff_asso_names|n2$v2%in%diff_asso_names,]
g_1<-graph_from_data_frame(n1)
g_2<-graph_from_data_frame(n2)

igraph::edge_density(g_1)
igraph::edge_density(g_2)
mean_distance(g_1,weights = E(g_1)$diss)
mean_distance(g_2,weights = E(g_2)$diss)

diffedge<-diff$assoMat1-diff$assoMat2
diffpadj<-diff[["pAdjustMat"]]
diffpadj[upper.tri(diffpadj)]<-1
#diffedge<-diffedge[diff_asso_names,diff_asso_names]
#diffpadj<-diffpadj[diff_asso_names,diff_asso_names]
diffpadj[diffpadj>0.05]<-1
diffedge<-reshape2::melt(diffedge)
colnames(diffedge)<-c("m1","m2","diff")
diffpadj<-reshape2::melt(diffpadj)
colnames(diffpadj)<-c("m1","m2","p")
diffpadj<-diffpadj[diffpadj$p<0.05&diffpadj$m1!=diffpadj$m2,]
diffedge<-diffedge[paste(diffedge$m1,diffedge$m2,sep = ",")%in%paste(diffpadj$m1,diffpadj$m2,sep = ","),]
difnetwork<-data.frame(Var1=diffedge$m1,Var2=diffedge$m2,cor=diffedge$diff)
#keep<-c(intersect(grep("m$",difnetwork$m1),grep("\\[e\\]",difnetwork$m2)),intersect(grep("\\[e\\]",difnetwork$m1),grep("m$",difnetwork$m2)))
#difnetwork<-difnetwork[keep,]
nodes<-unique(c(difnetwork$Var1,difnetwork$Var2))


nodeatr<-data.frame(nodes = nodes,degree=summ_res[["topProps"]][["topDeg"]][nodes,3],
                    Between=summ_res[["topProps"]][["topBetw"]][nodes,3],nodetype="difnetwork")
nodeatr<-nodeatr[nodeatr$nodes%in%unique(c(difnetwork$Var1,difnetwork$Var2)),]
dif_graph<-graph_from_data_frame(d=difnetwork,vertices = nodeatr)

cluster<-cluster_walktrap(dif_graph)
plot(cluster,dif_graph,vertex.label=NA,
     edge.arrow.mode='-',    
     vertex.size = 5 )
cluster<-cluster$membership
V(dif_graph)$cluster<-factor(paste0("guild_",cluster),levels = paste0("guild_",1:11))

g <- sample_islands(length(table(cluster)), max(table(cluster)), 0.5, length(table(cluster)))
g <- igraph::simplify(g)
set.seed(123)
bb <- graphlayouts::layout_as_backbone(g, keep = 0.2)
pos_tmp<-bb$xy

pos<-data.frame()
nodelist<-V(dif_graph)$cluster%>%levels()
for(i in 1:length(table(cluster))){
  set.seed(123)
  tmp<-pos_tmp[((max(table(cluster))*(i-1))+1):(max(table(cluster))*i),]
  tmp<-tmp[sample(1:max(table(cluster)),sum(paste0("guild_",cluster)==nodelist[i])),,drop=F]
  pos<-rbind(pos,tmp)
}
rownames(pos)<-names(V(dif_graph))[order(V(dif_graph)$cluster)]
pos<-pos[names(V(dif_graph)),]
difgraph<-ggraph(dif_graph,
                 layout = "manual",
                 x = pos[, 1],
                 y = pos[, 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  ggforce::geom_mark_hull( 
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("guild_1"="orange3","guild_2"="orchid4","guild_3"="red","guild_4"="orangered4",
                                "guild_5"="lightsteelblue4","guild_6"="linen","guild_7"="orange","guild_8"="royalblue4",
                                "guild_9"="skyblue4","guild_10"="thistle","guild_11"="turquoise4","guild_12"="palegreen4","guild_13"="mediumseagreen"),
                     labels=c(paste("guild",1:2,sep="_"),"guild_3",paste("guild",4:13,sep="_"))) +
  scale_fill_manual(values = c("guild_1"="orange3","guild_2"="orchid4","guild_3"="red","guild_4"="orangered4",
                               "guild_5"="lightsteelblue4","guild_6"="linen","guild_7"="orange","guild_8"="royalblue4",
                               "guild_9"="skyblue4","guild_10"="thistle","guild_11"="turquoise4","guild_12"="palegreen4","guild_13"="mediumseagreen"),
                    labels=c(paste("guild",1:2,sep="_"),"guild_3",paste("guild",4:13,sep="_"))) +
  scale_edge_colour_gradientn(name = "difcor",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, 0.28),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "none")+
  guides(fill = guide_legend(ncol = 1))+
  theme_graph(base_family = "sans")+theme(legend.box = 'horizontal',
                                          legend.box.just = 'top')+
  theme(plot.title = element_text(hjust = -0.02,vjust = 4.5),title = element_text(size = 16,face="bold"))
difgraph
ggsave("dif_cluster_net.pdf")

library(WGCNA)
ab_difnet<-ab_genus[names(V(dif_graph)),]
MEList <- WGCNA::moduleEigengenes(t(ab_difnet), colors = V(dif_graph)$cluster)
MES0<-MEList$eigengenes
colnames(MES0)<-gsub("ME","",colnames(MES0))
design <- model.matrix(~0+meta_all$disease)  # 构建模型矩阵
colnames(design)= c("healthy", "T2D")
moduleTraitCor <- cor(MES0,design,use = "p")  # 计算模块特征向量与表型的相关系数矩阵
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,ncol(ab_difnet))  # 计算相关系数矩阵的p值
textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")  # 构建绘图时用的文本矩阵
dim(textMatrix)=dim(moduleTraitCor)  # 修改文本矩阵的维度，与相关系数矩阵相同
rownames(moduleTraitCor)<-gsub("ME","",rownames(moduleTraitCor))
pdf("MT_relationship.pdf",height = 4,width = 4)
#par(mar=c(6, 8.5, 3, 3))  # 设置绘图边距
labeledHeatmap(Matrix = moduleTraitCor,  # 绘制带标签的热图
               xLabels = colnames(design),  # x轴标签
               yLabels = names(MES0),  # y轴标签
               ySymbols = names(MES0),  # y轴符号
               colorLabels = FALSE,  # 不显示颜色标签
               colors = blueWhiteRed(50),  # 颜色范围
               textMatrix = textMatrix,  # 显示文本矩阵
               setStdMargins = FALSE,  # 不设置标准边距
               cex.text = 0.5,  # 文本大小
               zlim = c(-1,1),  # 颜色映射范围
               main = paste("Guild-trait relationships"))
dev.off()


ab_net_d<-ab_difnet[,meta_all$disease=="T2D"]
ab_net_h<-ab_difnet[,meta_all$disease!="T2D"]
h_net<-diff$assoMat1|>as.data.frame()
h_net[diff$pAdjustMat>0.05]<-0
h_net<-h_net[rownames(ab_net_h),rownames(ab_net_h)]
d_net<-diff$assoMat2|>as.data.frame()
d_net[diff$pAdjustMat>0.05]<-0
d_net<-d_net[rownames(ab_net_d),rownames(ab_net_d)]

ab_net_d1<-as.matrix(d_net)%*%as.matrix(ab_net_d)
ab_net_h1<-as.matrix(h_net)%*%as.matrix(ab_net_h)
ab_net_all<-cbind(ab_net_d1,ab_net_h1)|>as.data.frame()
ab_net_all<-ab_net_all[,meta_all$sample_id]
MEList <- WGCNA::moduleEigengenes(t(ab_net_all), colors = V(dif_graph)$cluster)
MES0<-MEList$eigengenes
colnames(MES0)<-gsub("ME","",colnames(MES0))
design <- model.matrix(~0+meta_all$disease)  # 构建模型矩阵
colnames(design)= c("healthy", "T2D")
moduleTraitCor1 <- cor(MES0,design,use = "p")  # 计算模块特征向量与表型的相关系数矩阵
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor1,ncol(ab_difnet))  # 计算相关系数矩阵的p值
textMatrix1 = paste(signif(moduleTraitCor1,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")  # 构建绘图时用的文本矩阵
dim(textMatrix1)=dim(moduleTraitCor1)  # 修改文本矩阵的维度，与相关系数矩阵相同
rownames(moduleTraitCor1)<-gsub("ME","",rownames(moduleTraitCor1))
pdf("MT_relationship1.pdf",height = 4,width = 4)
#par(mar=c(6, 8.5, 3, 3))  # 设置绘图边距
labeledHeatmap(Matrix = moduleTraitCor1,  # 绘制带标签的热图
               xLabels = colnames(design),  # x轴标签
               yLabels = names(MES0),  # y轴标签
               ySymbols = names(MES0),  # y轴符号
               colorLabels = FALSE,  # 不显示颜色标签
               colors = blueWhiteRed(50),  # 颜色范围
               textMatrix = textMatrix1,  # 显示文本矩阵
               setStdMargins = FALSE,  # 不设置标准边距
               cex.text = 0.5,  # 文本大小
               zlim = c(-1,1),  # 颜色映射范围
               main = paste("Guild-trait relationships"))
dev.off()


guild4<-names(V(dif_graph))[V(dif_graph)$cluster=="guild_4"]

meta_d<-sampleMetadata |>
  #filter(age >= 18) |>
  filter(body_site == "stool") |>
  filter(age_category=="adult")|>
  filter(study_name=="YuJ_2015")|>
  filter(antibiotics_current_use=="no")|>
  filter(country=="CHN")|>
  filter(disease == "T2D")

studies<-unique(meta_d$study_name)
meta_h<-sampleMetadata |>
  filter(study_name %in%studies) |>
  filter(age_category=="adult")|>
  filter(antibiotics_current_use=="no")|>
  filter(body_site == "stool") |>
  filter(disease == "healthy")

meta_all<-rbind(meta_d,meta_h)
T2D_data_test<-meta_all |> returnSamples("relative_abundance")

ab<-T2D_data_test@assays@data@listData[["relative_abundance"]]
taxa<-sapply(rownames(ab), function(x){
  tmp<-strsplit(x,"\\|")[[1]]
  tmp[grep("g__",tmp)]
})
ab<-as.data.frame(ab)
ab$taxa<-taxa
ab<-aggregate(.~taxa,data=ab,sum)
rownames(ab)<-ab$taxa
ab_genus<-ab[,-1]
rownames(ab_genus)<-gsub("g__","",rownames(ab_genus))
ab_4<-ab_genus[guild4,]
net <- netConstruct(data = t(ab_4)[meta_all$disease=="T2D",],data2 = t(ab_4)[meta_all$disease=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0.1),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta_all$disease=="healthy")),
                    measure = "sparcc",
                    #measurePar = list(nlambda=10, 
                    #                  rep.num=10,
                    #                  Rmethod = "approx"),
                    normMethod = "clr", 
                    zeroMethod = "alrEM",
                    sparsMethod = "threshold",
                    thresh = 0.05,
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)

ab_difnet<-ab_4
MEList <- WGCNA::moduleEigengenes(t(ab_difnet), colors = rep("guild_4",15))
MES0<-MEList$eigengenes
colnames(MES0)<-gsub("ME","",colnames(MES0))
design <- model.matrix(~0+meta_all$disease)  # 构建模型矩阵
colnames(design)= c("healthy", "T2D")
moduleTraitCor3 <- cor(MES0,design,use = "p")  # 计算模块特征向量与表型的相关系数矩阵
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor3,ncol(ab_difnet))  # 计算相关系数矩阵的p值
textMatrix3 = paste(signif(moduleTraitCor3,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")  # 构建绘图时用的文本矩阵
dim(textMatrix3)=dim(moduleTraitCor3)  # 修改文本矩阵的维度，与相关系数矩阵相同


ab_net_d<-ab_4[,meta_all$disease=="T2D"]
ab_net_h<-ab_4[,meta_all$disease!="T2D"]
h_net<-net$adjaMat1|>as.data.frame()
h_net<-h_net[rownames(ab_net_h),rownames(ab_net_h)]
d_net<-net$assoMat2|>as.data.frame()
d_net<-d_net[rownames(ab_net_d),rownames(ab_net_d)]
ab_net_d1<-as.matrix(d_net)%*%as.matrix(ab_net_d)
ab_net_h1<-as.matrix(h_net)%*%as.matrix(ab_net_h)
ab_net_all<-cbind(ab_net_d1,ab_net_h1)|>as.data.frame()
ab_net_all<-ab_net_all[,meta_all$sample_id]
MEList <- WGCNA::moduleEigengenes(t(ab_net_all), colors = rep("guild_4",15))
MES0<-MEList$eigengenes
colnames(MES0)<-gsub("ME","",colnames(MES0))
design <- model.matrix(~0+meta_all$disease)  # 构建模型矩阵
colnames(design)= c("healthy", "T2D")
moduleTraitCor4 <- cor(MES0,design,use = "p")  # 计算模块特征向量与表型的相关系数矩阵
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor4,ncol(ab_difnet))  # 计算相关系数矩阵的p值
textMatrix4 = paste(signif(moduleTraitCor4,2),"\n(",
                    signif(moduleTraitPvalue,1),")",sep = "")
rownames(moduleTraitCor3)<-rownames(moduleTraitCor4)<-"guild_4 "

heat_all<-rbind(cbind(moduleTraitCor,moduleTraitCor1),cbind(moduleTraitCor3,moduleTraitCor4))
text_mat<-rbind(cbind(textMatrix,textMatrix1),cbind(matrix(textMatrix3,ncol = 2),matrix(textMatrix4,ncol = 2)))
colnames(heat_all)[3:4]<-paste0(" ",colnames(heat_all))[3:4]
anno<-data.frame(row.names = colnames(heat_all),group=c("Pre-Convolution","Pre-Convolution","Post-Convolution","Post-Convolution"))
anno_row<-data.frame(row.names = c(paste0("guild_",1:11),"guild_4 "),cohort=c(rep("discovery",11),"validation"))
pdf("MT_relationship_all.pdf",height = 5,width = 7)
ComplexHeatmap::pheatmap(heat_all,cluster_rows = F,color = colorRampPalette(c("#67a9cf","#f7f7f7", "#ef8a62"))( 100 ),
                         cluster_cols = F,gaps_col = 2,display_numbers = text_mat,gaps_row = 11,
                         fontsize_row = 14,fontsize = 13,angle_col = "45",
                         annotation_row = anno_row,
                         annotation_col = anno,
                         annotation_colors = list(group=c("Pre-Convolution"='deepskyblue',
                                                                                "Post-Convolution"='orange'),
                                                  cohort=c(discovery="#b9d9cf",validation="#a2608f")),
                         name = " "
)
dev.off()
