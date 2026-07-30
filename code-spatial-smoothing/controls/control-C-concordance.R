## Control C (R2.8): concordance of emergent L-R domains with Chen 2022 (MOSTA)
## anatomical annotations for the SAME section (E16.5_E1S3), aligned to our grid.
## Positive control: if the domains from clustering the L-R hyperfield recover known
## anatomy, that is evidence they capture real tissue architecture.
suppressMessages(library(data.table))
setwd("/var/home/aoz/Dropbox/genetics/a-s-omics"); ctrl <- "data-outputs/mouse-embryo-raw/controls"
dom <- as.data.table(readRDS(file.path(ctrl,"observed-domain-labels.rds")))   # x,y,domain (20)
ann <- fread(file.path(ctrl,"e16-fov-annotation.csv"))                          # x,y,annotation
d <- merge(dom, ann, by=c("x","y")); cat("merged bins:", nrow(d), "\n")
a <- d$domain; b <- d$annotation

ARI <- function(a,b){ t<-table(a,b); n<-sum(t)
  ci<-function(x) sum(choose(x,2))
  idx<-ci(t); ai<-ci(rowSums(t)); bj<-ci(colSums(t))
  exp<-ai*bj/choose(n,2); mx<-(ai+bj)/2; (idx-exp)/(mx-exp) }
NMI <- function(a,b){ t<-table(a,b); n<-sum(t); pa<-rowSums(t)/n; pb<-colSums(t)/n; p<-t/n
  H<-function(pr) -sum(pr[pr>0]*log(pr[pr>0]))
  mi<-sum(ifelse(p>0, p*log(p/outer(pa,pb)), 0)); mi/sqrt(H(pa)*H(pb)) }

ari <- ARI(a,b); nmi <- NMI(a,b)
cat(sprintf("\nObserved: ARI = %.3f | NMI = %.3f  (0 = chance, 1 = perfect)\n", ari, nmi))

## per-domain dominant anatomy + purity
tab <- d[, .N, by=.(domain, annotation)]
pur <- tab[, .(dominant=annotation[which.max(N)], purity=max(N)/sum(N), n=sum(N)), by=domain][order(domain)]
cat("\nPer-domain dominant anatomy (purity = fraction of domain in its top tissue):\n")
print(pur)
cat(sprintf("\nmean domain purity (size-weighted) = %.1f%%\n", 100*sum(pur$purity*pur$n)/sum(pur$n)))

## permutation null: shuffle domain labels, recompute ARI
set.seed(1); nperm<-499; null<-numeric(nperm)
for(i in seq_len(nperm)) null[i]<-ARI(sample(a), b)
p <- (1+sum(null>=ari))/(1+nperm)
cat(sprintf("\nPermutation null ARI: mean=%.4f max=%.4f ; observed=%.3f ; p=%.4f\n",
            mean(null), max(null), ari, p))
saveRDS(list(ari=ari,nmi=nmi,purity=pur,p=p,null=null), file.path(ctrl,"control-C-result.rds"))
