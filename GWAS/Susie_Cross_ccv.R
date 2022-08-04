### Susie_Cross_ccv.R
library(methods) 
args <- commandArgs(trailingOnly=T)
chr <- as.numeric(args[1])
pos <- as.numeric(args[2])


summary_chr <- read.table(paste("summary_EURASN_chr",chr, sep=""), sep="\t", head=F)
names(summary_chr) <- c("SNP","BETA","P","SE","TEST","OTHER","TEST_FREQ","N","CHR","POS","Z")
summary_pos <- summary_chr[summary_chr$POS > pos-500000 & summary_chr$POS < pos+500000 ,]
summary_pos <- summary_pos[order(summary_pos$POS), ]
summary_pos$ld_name <-paste("X", gsub(":",".",summary_pos$SNP), "_",summary_pos$TEST,sep="")
rownames(summary_pos) <- summary_pos$ld_name

ld_EUR <- read.table(paste("EUR_BOTH_chr",chr,"_",pos,".raw",sep=""), sep="\t", head=T)
ld_EUR2 <-ld_EUR[,-(1:6)]
sd1 <-apply(ld_EUR2 ,2,sd)
ld_EUR2 <- ld_EUR2[,which(sd1!=0)]

for (i in 1:ncol(ld_EUR2)){
	ld_EUR2[,i][is.na(ld_EUR2[,i])]<-mean(ld_EUR2[,i],na.rm=TRUE)
}
ld_EUR3 <- ld_EUR2[,names(ld_EUR2) %in% as.character(summary_pos$ld_name)]

ld_ASN <- read.table(paste("ASN_BOTH_chr",chr,"_",pos,".raw",sep=""), sep="\t", head=T)
ld_ASN2 <-ld_ASN[,-(1:6)]
sd2 <-apply(ld_ASN2 ,2,sd)
ld_ASN2 <- ld_ASN2[,which(sd2!=0)]

for (i in 1:ncol(ld_ASN2)){
	ld_ASN2[,i][is.na(ld_ASN2[,i])]<-mean(ld_ASN2[,i],na.rm=TRUE)
}
ld_ASN3 <- ld_ASN2[,names(ld_ASN2) %in% as.character(summary_pos$ld_name)]


summary_pos_match <- summary_pos[rownames(summary_pos) %in% names(ld_EUR3) & rownames(summary_pos) %in% names(ld_ASN3),]
summary_pos_match$order <- 1:nrow(summary_pos_match)

ld_EUR4 <- ld_EUR3[,rownames(summary_pos_match)]
ld_ASN4 <- ld_ASN3[,rownames(summary_pos_match)]


R1 <- cor(ld_EUR4)
R2 <- cor(ld_ASN4)
R <- (247173*R1 + 139523*R2)/386696
N=386696

library(susieR)
set.seed(1)
fitted_rss <- susie_rss(z=summary_pos_match$Z, R=R, n=N, L = 10)
	#plot<-susie_plot(fitted_rss, y="PIP")
	#png("pip.png", width=1080, height=720)
	#susie_plot(fitted_rss, y="PIP")
	#dev.off()
#cs <- summary(fitted_rss)$cs

cs<-susie_get_cs(fitted_rss ,as.matrix(R),min_abs_corr = 0.2)
cs
length(cs$cs$L1)
summary_pos_match[summary_pos_match$POS==pos,]
write.table(summary_pos_match[summary_pos_match$order %in% cs$cs$L1,],paste("CROSS_chr",chr,"_",pos,"_CCV",sep=""),row.names=F,quote=F,col.names=F,sep="\t")

if(nrow(cs)==1){
	array <- strsplit(cs$variable[1],",")
	write.table(summary_pos_match[summary_pos_match$order %in% as.numeric(array[[1]]),], paste("CROSS_chr",chr,"_",pos,"_CCV",sep=""),row.names=F,quote=F,col.names=F,sep="\t")
}

if(nrow(cs)>1){
	for (i in 1:nrow(cs)){
		array <- strsplit(cs$variable[i],",")
		if(pos %in% summary_pos_match[summary_pos_match$order %in% as.numeric(array[[1]]),]$POS) {
			write.table(summary_pos_match[summary_pos_match$order %in% as.numeric(array[[1]]),], paste("CROSS_chr",chr,"_",pos,"_CCV",sep=""),row.names=F,quote=F,col.names=F,sep="\t")
		}		
	}
}

cs <-susie_get_cs(fitted_rss ,as.matrix(cor(ld_raw3)),min_abs_corr = 0.2)

