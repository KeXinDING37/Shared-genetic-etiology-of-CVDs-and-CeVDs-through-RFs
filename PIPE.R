# PIPE
# load packages
library(PIPE)
library(tidyverse)
library(patchwork)

# specify built-in data placeholder
placeholder <- "http://www.comptransmed.pro/bigdata_pipe"

# create the output folder 'pipe_showcase'
outdir <- 'pipe_showcase'
if(!dir.exists(outdir)) dir.create(outdir)

data.file <- file.path(placeholder, "GWAS.cross.SNP.txt")
input_data <- read_delim(data.file, delim='\t') %>% filter(pubmedid==26974007) %>% distinct(snp_id_current,pvalue)

input_data <- data.frame(matrix(NA,21,2))
input_data[,1] <- c('rs17035646',
'rs880315',
'rs4151702',
'rs4135240',
'rs2107595',
'rs2891168',
'rs7859727',
'rs635634',
'rs2519093',
'rs532436',
'rs360153',
'rs3184504',
'rs1535791',
'rs9534454',
'rs4942561',
'rs4932373',
'rs35346340',
'rs76774446',
'rs73015024',
'rs55997232',
'rs112374545')
input_data[,2] <- 5e-10
colnames(input_data) <- c('snp_id_current','pvalue')

# LD.customised
LD.customised <- oRDS('GWAS_LD.EUR', placeholder=placeholder)

# All common SNPs (primarily sourced from dbSNP)
options(timeout=120) # please increase the timeout if failed to download
GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)

# All known genes (primarily sourced from UCSC)
GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)

# Functional interaction network (primarily sourced from STRING)
# restricted to high-quality interactions ("experiments" and "databases")
network.customised <- oRDS("STRING_fin", placeholder=placeholder)

# include.RGB
include.RGB <- c("PCHiC_PMID27863249_Activated_total_CD4_T_cells","PCHiC_PMID27863249_Macrophages_M0","PCHiC_PMID27863249_Macrophages_M1","PCHiC_PMID27863249_Macrophages_M2","PCHiC_PMID27863249_Megakaryocytes","PCHiC_PMID27863249_Monocytes","PCHiC_PMID27863249_Naive_B_cells","PCHiC_PMID27863249_Naive_CD4_T_cells","PCHiC_PMID27863249_Naive_CD8_T_cells","PCHiC_PMID27863249_Neutrophils","PCHiC_PMID27863249_Nonactivated_total_CD4_T_cells","PCHiC_PMID27863249_Total_B_cells","PCHiC_PMID27863249_Total_CD4_T_cells","PCHiC_PMID27863249_Total_CD8_T_cells","PCHiC_PMID31501517_DorsolateralPrefrontalCortex","PCHiC_PMID31501517_H1NeuralProgenitorCell","PCHiC_PMID31501517_Hippocampus","PCHiC_PMID31501517_LCL", "PCHiC_PMID31367015_astrocytes","PCHiC_PMID31367015_excitatory","PCHiC_PMID31367015_hippocampal","PCHiC_PMID31367015_motor", "PCHiC_PMID25938943_CD34","PCHiC_PMID25938943_GM12878")

# include.QTL
include.QTL <- c("eQTL_eQTLGen", "eQTL_eQTLCatalogue_Alasoo_2018_macrophage_IFNg","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_IFNgSalmonella","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_Salmonella","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_naive","eQTL_eQTLCatalogue_BLUEPRINT_Tcell","eQTL_eQTLCatalogue_BLUEPRINT_monocyte","eQTL_eQTLCatalogue_BLUEPRINT_neutrophil","eQTL_eQTLCatalogue_BrainSeq_brain","eQTL_eQTLCatalogue_CEDAR_microarray_Bcell_CD19","eQTL_eQTLCatalogue_CEDAR_microarray_Tcell_CD4","eQTL_eQTLCatalogue_CEDAR_microarray_Tcell_CD8","eQTL_eQTLCatalogue_CEDAR_microarray_monocyte_CD14","eQTL_eQTLCatalogue_CEDAR_microarray_neutrophil_CD15","eQTL_eQTLCatalogue_CEDAR_microarray_platelet","eQTL_eQTLCatalogue_Fairfax_2012_microarray_Bcell_CD19","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_IFN24","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_LPS2","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_LPS24","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_naive","eQTL_eQTLCatalogue_GENCORD_LCL","eQTL_eQTLCatalogue_GENCORD_Tcell","eQTL_eQTLCatalogue_GEUVADIS_LCL","eQTL_eQTLCatalogue_GTEx_V8_Brain_Amygdala","eQTL_eQTLCatalogue_GTEx_V8_Brain_Anterior_cingulate_cortex_BA24","eQTL_eQTLCatalogue_GTEx_V8_Brain_Caudate_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cerebellar_Hemisphere","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cerebellum","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cortex","eQTL_eQTLCatalogue_GTEx_V8_Brain_Frontal_Cortex_BA9","eQTL_eQTLCatalogue_GTEx_V8_Brain_Hippocampus","eQTL_eQTLCatalogue_GTEx_V8_Brain_Hypothalamus","eQTL_eQTLCatalogue_GTEx_V8_Brain_Nucleus_accumbens_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Putamen_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Spinal_cord_cervical_c1","eQTL_eQTLCatalogue_GTEx_V8_Brain_Substantia_nigra","eQTL_eQTLCatalogue_GTEx_V8_Nerve_Tibial","eQTL_eQTLCatalogue_GTEx_V8_Pituitary","eQTL_eQTLCatalogue_GTEx_V8_Whole_Blood","eQTL_eQTLCatalogue_Kasela_2017_microarray_Tcell_CD4","eQTL_eQTLCatalogue_Kasela_2017_microarray_Tcell_CD8","eQTL_eQTLCatalogue_Lepik_2017_blood","eQTL_eQTLCatalogue_Naranbhai_2015_microarray_neutrophil_CD16","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_Listeria","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_Salmonella","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_naive","eQTL_eQTLCatalogue_Quach_2016_monocyte_IAV","eQTL_eQTLCatalogue_Quach_2016_monocyte_LPS","eQTL_eQTLCatalogue_Quach_2016_monocyte_Pam3CSK4","eQTL_eQTLCatalogue_Quach_2016_monocyte_R848","eQTL_eQTLCatalogue_Quach_2016_monocyte_naive","eQTL_eQTLCatalogue_ROSMAP_brain_naive","eQTL_eQTLCatalogue_Schmiedel_2018_Bcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_CD4_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_Schmiedel_2018_CD4_Tcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_CD8_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_Schmiedel_2018_CD8_Tcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_NKcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_Tfh_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th117_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th17_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th1_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th2_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Treg_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Treg_naive","eQTL_eQTLCatalogue_Schmiedel_2018_monocyte_CD16_naive","eQTL_eQTLCatalogue_Schmiedel_2018_monocyte_naive","eQTL_eQTLCatalogue_Schwartzentruber_2018_sensory_neuron")

# prepare predictors
ls_pNode_genomic <- oPierSNPsAdv(data=input_data, LD.customised=LD.customised, significance.threshold=5e-8, GR.SNP=GR.SNP, GR.Gene=GR.Gene, distance.max=20000, decay.kernel="constant", include.QTL=include.QTL, include.RGB=include.RGB, network.customised=network.customised, placeholder=placeholder, verbose.details=F)

# do prioritisation（the initial round）
ls_pNode <- Filter(Negate(is.null), ls_pNode_genomic)
dTarget0 <- oPierMatrix(ls_pNode, displayBy="pvalue", aggregateBy="fishers", rangeMax=10, GR.Gene=GR.Gene, placeholder=placeholder)

# GSP (phase 2 or above drug targets for inflammatory disorders)
ig.EF.ChEMBL <- oRDS('ig.EF.ChEMBL', placeholder=placeholder)
codes_name <- c("psoriasis","ankylosing spondylitis","Crohn's disease","ulcerative colitis","sclerosing cholangitis")
codes_name <- c("cardiovascular disease","coronary artery disease","myocardial infarction",
                "heart failure","atrial fibrillation","Ischemic stroke","stroke",
                "hypertension","genetic hypertension",                
                "Disorder of lipid metabolism","hyperlipidemia",
                "Genetic obesity","obesity",
                "diabetes mellitus",
                "internal carotid artery stenosis")
ind <- match(codes_name, V(ig.EF.ChEMBL)$name)
GSP <- base::Reduce(union, V(ig.EF.ChEMBL)$GSP[ind])

# GSN
druggable_targets <- do.call(rbind, V(ig.EF.ChEMBL)$phase) %>% pull(target) %>% unique()
g1 <- oRDS('org.Hs.string_medium', placeholder=placeholder)
g2 <- oRDS('org.Hs.PCommons_UN', placeholder=placeholder)
g <- oCombineNet(list(g1,g2), combineBy='union', attrBy="intersect")
sGS <- oGSsimulator(GSP=GSP, population=druggable_targets, network.customised=g, neighbor.order=1, verbose=F)
GSN <- sGS$GSN
message(sprintf("GSP (%d), GSN (%d)", length(GSP), length(GSN)), appendLF=T)

# predictor matrix
df_predictor <- oPierMatrix(ls_pNode, displayBy="score", aggregateBy="none", GR.Gene=GR.Gene, placeholder=placeholder)

# informative predictors
sClass <- oClassifyRF(df_predictor, GSP=GSP, GSN=GSN, nfold=3, nrepeat=10, mtry=NULL, ntree=1000, cv.aggregateBy="fishers")
df_imp <- sClass$importance %>% as_tibble(rownames='predictor')
df_perf <- sClass$performance %>% as_tibble(rownames='predictor')
imp_perf <- df_imp %>% inner_join(df_perf, by='predictor') %>% mutate(group=str_replace_all(predictor,'_.*','')) %>% arrange(group,-MeanDecreaseAccuracy)
xintercept <- imp_perf %>% filter(group=='nGene') %>% pull(MeanDecreaseAccuracy)
yintercept <- imp_perf %>% filter(group=='nGene') %>% pull(MeanDecreaseGini)
zintercept <- imp_perf %>% filter(group=='nGene') %>% pull(auroc)
imp_perf2 <- imp_perf %>% filter(MeanDecreaseAccuracy>=0, MeanDecreaseGini>=0, direction=='+') %>% mutate(x=MeanDecreaseAccuracy>=xintercept, y=MeanDecreaseGini>=yintercept, z=auroc>=zintercept) %>% arrange(-MeanDecreaseAccuracy)

# predictors selected based the importance no worse than nGene
vec_predictor_selected <- imp_perf2 %>% filter(x|y|z) %>% pull(predictor)
vec_predictor_selected

# do prioritisation（based on informative predictors）
dTarget <- oPierMatrix(ls_pNode[vec_predictor_selected], displayBy="pvalue", aggregateBy="fishers", rangeMax=10, GR.Gene=GR.Gene, placeholder=placeholder, verbose=F)

# write into a file 'CRM_IND_priority.txt.gz' under the folder 'pipe_showcase'
dTarget$priority %>% select(name,rating,rank,seed,nGene,eGene,cGene,description) %>% write_delim(file.path(outdir,'CRM_IND_priority.txt.gz'), delim='\t')

# save into a file 'dTarget_IND.CRM.RData' under the folder 'pipe_showcase'
save(list='dTarget', file=file.path(outdir,'dTarget_IND.CRM.RData'))

# df_priority
df_priority <- dTarget$priority %>% column_to_rownames("name") %>% top_frac(1,rating) %>% transmute(priority=rating, rank=rank)

write.table(V(ig.EF.ChEMBL)$name, "name.csv", quote = F, row.names = F)

# prepare SET
codes_name <- c("cardiovascular disease","coronary artery disease","myocardial infarction",
"heart failure","atrial fibrillation","Ischemic stroke","stroke",
"hypertension","genetic hypertension",                
"Disorder of lipid metabolism","hyperlipidemia",
"Genetic obesity","obesity",
"diabetes mellitus",
"internal carotid artery stenosis")
ind <- match(codes_name, V(ig.EF.ChEMBL)$name)
codes <- c('CVD','CAD','MI','HF','AF','AS','AIS','HYP','gHYP','LM','HLIP',
           'GOS','OS','DM','IMT')

ls_gsp <- V(ig.EF.ChEMBL)$GSP[ind]
names(ls_gsp) <- codes

info <- ls_gsp %>% enframe() %>% transmute(id=name, member=value)

info <- info %>% mutate(name=codes_name, namespace='IND', distance=NA) %>% mutate(n=map_dbl(member,length)) %>% select(name,id,namespace,distance,member,n)

set <- list(info=info, stamp=as.Date(Sys.time()))
class(set) <- "SET"

# do leading prioritisation analysis (LPA)
set.seed(1020)

str(df_priority)
df_cvd <- data.frame(matrix(NA,12,2))
df_cvd[,1] <- 10
df_cvd[,2] <- 1
rownames(df_cvd) <- c('CASZ1', 'CDKN1A', 'TWIST1', 'CDKN2B', 'ABO', 
'SWAP70', 'SH2B3', 'LRCH1', 'FES', 'GOSR2', 'RPRML', 'LDLR')
colnames(df_cvd) <- c('priority','rank')

df <- df_priority[100:120, ]

eLPA <- oPierGSEA(df_priority, customised.genesets=set, size.range=c(10,5000), nperm=500000)

df_summary <- eLPA$df_summary %>% mutate(flag=ifelse(adjp<0.05,'Y','N')) %>% mutate(frac=nLead/nAnno) %>% filter(frac<=1) %>% arrange(tolower(setID))
df_leading <- eLPA$leading[df_summary %>% pull(setID)] %>% enframe() %>% transmute(setID=name,members=value) %>% mutate(rank_members=map_chr(members, function(x) str_c(names(x),' (',x,')', collapse=', ')))
df_summary_leading <- df_summary %>% inner_join(df_leading, by='setID')

# save into a file 'CRM_IND_lpa.txt' under the folder 'pipe_showcase'
df_summary_leading %>% select(-members) %>% write_delim(file.path(outdir,'CRM_IND_lpa.txt'), delim='\t')


pict <- df_summary_leading %>% mutate(pict=log10(nes * frac / adjp)) %>% select(setID, nes, adjp, frac, nAnno, nLead, pict) %>% arrange(tolower(setID))
gp_pict <- oPolarBar(pict %>% transmute(name=setID,value=pict), colormap='sci_locuszoom', size.name=9, size.value=4)

# saved into a file 'CRM_IND_pict.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_pict.pdf'), gp_pict, width=6.5, height=6.5)

leading.query <- df_priority %>% as_tibble(rownames='target') %>% top_frac(1,priority) %>% pull(target)
gp_leading <- oGSEAdotplot(eLPA, top='ulcerative colitis', x.scale="normal", leading=T, leading.edge.only=T, leading.query.only=T, leading.query=leading.query, colormap="jet.top", zlim=c(0,10), peak.color='black', leading.size=2, clab='Credit\nscore', leading.label.direction="left") + theme(plot.title=element_text(hjust=0.5,size=10), plot.subtitle=element_text(hjust=0.5,size=8))

# saved into a file 'CRM_IND_leading_uc.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_leading_uc.pdf'), gp_leading, width=6, height=8)

# Member genes for each disease defined as clinical proof-of-concept targets recovered at the leading prioritisation
data <- df_summary_leading %>% transmute(name=setID, members=rank_members)

# remove (rank)
data0 <- data %>% separate_rows(members, sep=", ") %>% mutate(members=str_replace_all(members, ' \\(.*', '')) %>% nest(data=-name) %>% mutate(members=map_chr(data, function(x) str_c(unlist(x), collapse=', '))) %>% select(-data)

# bigraph
big <- data %>% separate_rows(members, sep=", ") %>% mutate(flag=1) %>% select(members,name,flag) %>% pivot_wider(names_from=members, values_from=flag, values_fill=list(flag=0)) %>% column_to_rownames("name") %>% oBicreate() %>% oBigraph()
big <- big %>% oLayout("graphlayouts.layout_with_stress")
igraph::degree(big) -> V(big)$degree
if(1){
  V(big)$rating <- str_replace_all(V(big)$name, '.*\\(|\\).*', '') %>% as.numeric()
  ifelse(V(big)$type=='xnode', V(big)$name, "") -> V(big)$label
  ifelse(V(big)$type=='ynode' & V(big)$rating<=nrow(df_priority) * 0.02 | V(big)$type=='xnode', V(big)$name, "") -> V(big)$label
}
gp_big <- big %>% oGGnetwork(node.label="label", node.label.size=2, node.label.force=0.02, node.xcoord="xcoord", node.ycoord="ycoord", colormap='steelblue-steelblue', node.shape="type", node.size="degree", node.size.range=c(0.5,6),edge.color="grey80", edge.color.alpha=0.3, edge.arrow.gap=0, edge.curve=0.05)
## saved into a file 'CRM_IND_big.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_big.pdf'), gp_big, width=6, height=6)

# identify the minimum spanning tree (mst)
tig <- data %>% oTIG(method='empirical', empirical.cutoff=0.8)
tig <- tig %>% oLayout("gplot.layout.fruchtermanreingold")
gp_mst <- oGGnetwork(tig, node.label='name', node.label.size=2.5, node.label.color='steelblue', node.label.force=1, node.xcoord='xcoord', node.ycoord='ycoord', colormap="steelblue-steelblue", edge.size='weight', node.shape=18, node.size='n_members', node.size.title="Num of\ngenes", slim=c(0,25), node.size.range=c(2,5), edge.color='gray80', edge.color.alpha=0.3, edge.curve=0,edge.arrow.gap=0.01)
## saved into a file 'CRM_IND_mst.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_mst.pdf'), gp_mst, width=5, height=5)

library(ComplexHeatmap)

# df_priority_leading
data <- df_summary_leading %>% transmute(name=setID,members) %>% mutate(target=map(members, function(x) names(x))) %>% mutate(rating=map(members, function(x) x)) %>% select(name,target,rating) %>% unnest(c(target,rating))
df_x <- data %>% pivot_wider(names_from=name, values_from=rating)
df_y <- df_priority %>% as_tibble(rownames='target')
df_priority_leading <- df_y %>% inner_join(df_x, by='target') %>% arrange(-priority)

# write into a file 'CRM_IND_priority_leading.txt' under the folder 'pipe_showcase'
df_priority_leading %>% write_delim(file.path(outdir,'CRM_IND_priority_leading.txt'), delim='\t')

# ht_main
D <- df_priority_leading %>% mutate(rname=str_c(target,' (',rank,')')) %>% select(rname, priority) %>% column_to_rownames('rname')
col_fun <- circlize::colorRamp2(seq(0,6,length=64), oColormap('spectral')(64))
#vec_group <- df_priority_leading %>% pull(rank)
#vec_group <- as.numeric(cut(df_priority_leading$rank, breaks=seq(0,10000,100)))
vec_group <- as.numeric(cut(df_priority_leading$rank, breaks=ceiling(nrow(df_priority) * seq(0,0.5,0.05))))
ht_main <- Heatmap(D, name="Credit\nscore", col=col_fun, border=T, 
                   row_split=vec_group, cluster_rows=F, show_row_dend=T, row_title=NULL, show_row_names=T, row_names_gp=gpar(fontsize=5), row_names_rot=0, row_names_side="left", clustering_distance_rows="euclidean", clustering_method_rows="average",
                   cluster_columns=F, column_order=1:ncol(D), show_column_names=F, column_names_side="top", column_names_gp=gpar(fontsize=8), column_names_rot=0
)

# row_anno: ls_ha_mark
colormap <- "sci_locuszoom"
ncolor <- ncol(df_priority_leading) - 3
ls_ha_mark <- lapply(1:ncolor, function(j){
  node.color <- oColormap(colormap)(ncolor)[j]
  vec_at <- df_priority_leading %>% pull(j+3)
  vec_at[!is.na(vec_at)] <- 1
  vec_at[is.na(vec_at)] <- 0
  at <- which(vec_at==1)
  ha_mark_at <- rowAnnotation(
    D=anno_simple(vec_at, col=c("0"="grey95", "1"=node.color), width=unit(1.5,"mm")),
    mark=anno_mark(at=at, labels=df_priority_leading$target[at], labels_gp=gpar(fontsize=5,col="black"), lines_gp=gpar(col=node.color))
  )
  # D1, D2, D3, ...
  ha_mark_at@anno_list$D@label <- str_c('D',j)
  ha_mark_at
})
names(ls_ha_mark) <- colnames(df_priority_leading)[4:(3+ncolor)]
names(ls_ha_mark)
#ht <- ht_main + ls_ha_mark[[1]] + ls_ha_mark[[2]] + ls_ha_mark[[3]] + ls_ha_mark[[4]] + ls_ha_mark[[5]]
#draw(ht, heatmap_legend_side="left", annotation_legend_side="left")

# row_anno: ha_mark_approved
codes_name <- c("ankylosing spondylitis","Crohn's disease","psoriasis","sclerosing cholangitis","ulcerative colitis")
ChEMBL <- oRDS("ChEMBL_v33", placeholder=placeholder)
ChEMBL_subset <- ChEMBL %>% filter(efo_term %in% codes_name)
### df_target_label
df_target_label <- ChEMBL_subset %>% filter(phase==4,target_number<=5) %>% select(Symbol,pref_name_drug,efo_term) %>% nest(data=pref_name_drug) %>% mutate(drugs=map_chr(data, function(x) {
  y <- x %>% pull(pref_name_drug)
  if(length(y)>1){
    y <- c(y %>% head(1), '...')
    str_c(y,collapse=" | ")
  }else{
    str_c(y %>% head(1),collapse=" | ")
  }
})) %>% mutate(label=str_c(efo_term,' (',drugs,')')) %>% select(Symbol,label) %>% group_by(Symbol) %>% summarise(label=str_c(label,collapse='\n    ')) %>% ungroup() %>% mutate(label=str_c(Symbol,':\n    ',label))
### df_label
df_label <- df_priority_leading %>% select(target) %>% left_join(df_target_label, by=c('target'='Symbol'))
vec_at <- df_label$label
at <- which(!is.na(vec_at))
ha_mark_approved <- rowAnnotation(
  A=anno_simple(0+!is.na(vec_at), col=c("0"="grey90", "1"="black"), width=unit(3,"mm")),
  mark=anno_mark(at=at, labels=df_label$label[at], labels_gp=gpar(fontsize=5,col='black'), lines_gp=gpar(col='black'))
)
#ht <- ht_main + ha_mark_approved
#draw(ht, heatmap_legend_side="left", annotation_legend_side="left")

# do visualisation
## saved into a file 'CRM_IND_heatmap_therapeutics.pdf' under the folder 'pipe_showcase'
ht <- ht_main + ls_ha_mark[[1]] + ls_ha_mark[[2]] + ls_ha_mark[[3]] + ls_ha_mark[[4]] + ls_ha_mark[[5]] + ha_mark_approved
pdf(file.path(outdir, 'CRM_IND_heatmap_therapeutics.pdf'), width=6.7, height=9.5)
draw(ht, heatmap_legend_side="left", annotation_legend_side="left")
dev.off()