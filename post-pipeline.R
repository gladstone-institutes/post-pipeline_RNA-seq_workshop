rm(list=ls()) # clear workspace, just in case

# if need be, use setwd() to set your working directory to the folder that contains the data and code
setwd("~/git.repos/ilab.repos/post-pipeline_RNA-seq_workshop/")

# load libraries/packages we'll need
library(magrittr) # for programming pipes
library(edgeR) # for differential expression analysis
library(org.Mm.eg.db) # for mouse gene annotation
library(ggplot2) # to visualize parts of the linear model
library(reshape2) # to reformat data for plotting

# Load experimental metadata
targets <- read.delim("targets.txt", stringsAsFactors=FALSE)
targets

# Make a grouping variable that contains all possible combinations of experimental variables
group <- paste(targets$CellType, targets$Status, sep=".") %>% factor
group # Notice that by default the factor levels are alphabetically ordered
table(group) # 2 reps per group or experimental condition combination

# Load raw gene-wise counts
GenewiseCounts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz",
                             row.names="EntrezGeneID")

colnames(GenewiseCounts) <- colnames(GenewiseCounts) %>% substring(.,1,7)
colnames(GenewiseCounts) %<>% substring(.,1,7) # trim column names to match those in `targets`
dim(GenewiseCounts)
GenewiseCounts %>% head

# Create a DGEList object, which is the main object edgeR uses to store data and results
y <- DGEList(GenewiseCounts[,-1], # -1 removes the first column, which is not a column of read counts
             group=group, 
             samples=targets,
             genes=GenewiseCounts[,1,drop=FALSE]) # drop=FALSE is keeping it as a data.frame instead of simplifying to a vector
# NOTE: you don't necessarily need to specify group, sample, and genes at this stage. Often it's only group that's specified.

# Add gene symbols (i.e. abbreviated gene names)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENTREZID", column="SYMBOL")
# NOTE: the biomaRt package is also very useful for gene annotation

head(y$genes)
dim(y$genes) # 27179 rows (genes) & 2 columns
which(is.na(y$genes$Symbol)) %>% length
is.na(y$genes$Symbol) %>% sum

# DGEList objects can be subset by rows (genes) and columns (samples) just like a matrix or data frame
y <- y[!is.na(y$genes$Symbol), ] # remove genes without a symbol (i.e. abbreviated name)
y <- y[is.na(y$genes$Symbol) %>% not, ]

dim(y) # 26221 genes left

# Remove genes with very low counts
keep <- filterByExpr(y) # basically does this: keep <- rowSums(cpm(y) > 0.5) >= 2. See filterByExpr's documentation and Chen & Smyth (2016) for details.
table(keep) # keeping ~15K genes

# DGEList objects can be subset by rows (genes) and columns (samples) just like a matrix or data frame
y <- y[keep, , keep.lib.sizes=FALSE] # keep.lib.sizes=FALSE causes the library sizes to be recomputed after filtering, which is generally recommended

# Normalization ---------
y <- calcNormFactors(y)
str(y)
y %>% str
y$samples

# Use plotMDS for QC check and for data exploration ----------
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

# Explain experimental setup and the questions we'll be answering with the specific tests => **make this a slide**

# Specify the design matrix ------------
design <- model.matrix( ~ 0 + group )
colnames(design) # these are the default R design matrix column labels
colnames(design) <- levels(group) # rename columns to match the group labels we made above
design

# Estimate Dispersion (variability) ----------
y <- estimateDisp(y, design, robust=TRUE) # using robust=TRUE is recommended
# edgeR shares information across genes and conditions to better estimate dispersion. More details in the User's Guide and manuscripts.

# Test different hypotheses/comparisons using a joint model fit --------
# Fit the above model (as specified in the design matrix) to the data
fit <- glmQLFit(y, design, robust=TRUE) # another option is to use glmFit
head(fit$coefficients %>% exp)

fit %>% names
fit$fitted.values %>% head
fit$AveLogCPM %>% str

fit$design %>% colnames
edgeR::cpm(y, log = T) %>% head


grp.model.tests <- lapply(1:6, FUN = function(i) glmQLFTest(fit, coef = i) )

lapply(grp.model.tests, FUN = function(coef.test) topTags(coef.test) )

lapply(grp.model.tests, FUN = function(coef.test) decideTestsDGE(coef.test) %>% summary ) %>% 
  do.call(cbind, .)

grp.model.tests %>% str

ggplot(data = fit,
       aes())

# Test DE between the basal pregnant and lactating groups
# Construct the contrast for the comparison of interest
con <- makeContrasts(B.pregnant - B.lactating, levels = design) # Note the sign of the comparison
qlf <- glmQLFTest(fit, contrast = con) # another option is to use glmLRT (after first using glmFit)
# The coef option to glmQLFTest & glmLRT can be used to test if indiv. coefficients are ≠ 0

topTags(qlf) # view/extract top DE genes

summary(decideTests(qlf)) # Summary of DE genes

# This is an equivalent numerical way of specifying the same contrast: 
qlf <- glmQLFTest(fit, contrast = c(-1,1,0,0,0,0) ) # ...since we're contrasting the 2nd coefficient against the 1st
# NOTE: the contrasts *must* sum to 0

# The following is an example of how to use the F test to obtain a joint p-value
# for several contrasts/comparisons, in ANOVA fashion:
con <- makeContrasts(
  L.PvsL = L.pregnant - L.lactating,
  L.VvsL = L.virgin - L.lactating,
  L.VvsP = L.virgin - L.pregnant, levels=design)

anov <- glmQLFTest(fit, contrast=con)
topTags(anov)

con <- makeContrasts(
  pre = L.pregnant - B.pregnant,
  lac = L.lactating - B.lactating,
  vir = L.virgin - B.virgin, levels=design)

anov <- glmQLFTest(fit, contrast=con)
topTags(anov)

# EXERCISE:
# How would you setup the contrast to find genes that are DE between all luminal (L) cells versus all basal (B) cells?

LvB <- glmQLFTest(fit, contrast = c(-1/3,-1/3,-1/3,1/3,1/3,1/3) )
# LvB <- glmQLFTest(fit, 
#                   contrast = makeContrasts(contrasts = ("L.lactating+L.pregnant+L.virgin")/3 - 
#                                              ("B.lactating+B.pregnant+B.virgin"
#                                               )/3) ,
#                   levels = design))
summary(decideTests(LvB))

# Let's try a more typical model with an intercept, and main and interaction effects ------------

# Turn CellType and Status into factor variables with a specific level order,
# i.e. with an appropriate reference/baseline level
targets$CellType %<>% factor
targets$CellType %>% levels # The default alphanumeric sorting is the same as the order we want in this case

targets$Status %<>% factor(levels = c("virgin", "pregnant", "lactating"))
targets$Status <- targets$Status %>% factor(levels = c("virgin", "pregnant", "lactating"))
targets$Status %>% levels

design <- model.matrix( ~ 1 + CellType + Status + CellType:Status, data = targets)
design <- model.matrix( ~ 1 + CellType*Status, data = targets)
colnames(design)
design
# The intercept represents the combination of the reference levels for both predictor variables: B (Basal) & virgin

# Estimate Dispersion ----------
y <- estimateDisp(y, design, robust=TRUE)

fit <- glmQLFit(y, design, robust=TRUE) # using robust=TRUE is recommended
head(fit$coefficients %>% exp)
fit$AveLogCPM %>% head

# Let's build a quick summary of DE results of testing all the coefficients (except for the intercept)
int.model.tests <- lapply(2:6, FUN = function(i) glmQLFTest(fit, coef = i) )
int.model.tests %>% length
int.model.tests[[1]] %>% topTags
int.model.tests %>% names
lapply(int.model.tests, FUN = function(coef.test) decideTestsDGE(coef.test) %>% summary ) %>% 
  do.call(cbind, .)

# As you can see, this result agrees with the MDS plot above and makes sense biologically.

# However, it is often the case that main effects capture most of the variance in the data 
# (i.e. explain it well for most genes).
# In those situations the interaction terms have few significantly DE genes.

# If there's time, during Q & A we'll discuss R model formulae to analyze experimental designs 
# from your own research

#----------
# Additional Resources / Further Reading:
# This workshop was partly based on the edgeR User Guide and this publication:
# https://f1000research.com/articles/5-1438/v2
