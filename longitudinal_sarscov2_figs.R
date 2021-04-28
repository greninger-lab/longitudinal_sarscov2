library(ggplot2)
library(ggExtra)
library(dplyr)
library(ggthemes)
library(tidyr)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(gtable)
library(grid)
library(ggmsa)
library(viridis)
library(ggridges)

setwd("/Users/gerbix/Documents/Michelle/longitudinal_sars_cov2/curated_set/variants_files")

# Pavi's variants table
old_variants_df <- read.table("All_allelefreqs_annotated_2020-09-26.csv",sep=",",header=TRUE) %>% mutate(af_diff =abs(af_diff))

# Variant table from 2021 Taylor reruns, and 2nd replicate table
variants_df <- read.csv("all_taylor_variants.csv", header = TRUE)
rep2_variants_df <- read.csv("rep2_variants.csv", header = TRUE)
# Convert columns to numerics
variants_df$AAFREQ <- as.numeric(as.character(variants_df$AAFREQ))
variants_df$TCOV <- as.numeric(as.character(variants_df$TCOV))
variants_df$NTPOS <- as.numeric(as.character(variants_df$NTPOS))
# Add pt and sample variable column by str splitting file name
variants_df <- variants_df %>% rowwise() %>% mutate(patient = strsplit(SAMPLE, split="_")[[1]][1])
variants_df<- variants_df %>% rowwise() %>% mutate(day = as.numeric(strsplit(strsplit(SAMPLE, split="_")[[1]][2],split="d"[[1]][1]))) %>% drop_na(TCOV)
rep2_variants_df <- rep2_variants_df %>% rowwise() %>% mutate(patient = strsplit(SAMPLE, split="_")[[1]][1])
rep2_variants_df<- rep2_variants_df %>% rowwise() %>% mutate(day = as.numeric(strsplit(strsplit(SAMPLE, split="_")[[1]][2],split="d"[[1]][1]))) %>% drop_na(TCOV)


# Find maximum allele differences within-host
#filt_variants_df <- filt_variants_df %>% group_by(patient,snpid) %>% add_tally() %>% filter(n>=2)
filt_variants_df <- variants_df %>% group_by(patient,snpid) %>% 
                    mutate(min_af = min(AAFREQ), max_af = max(AAFREQ), 
                           af_diff = case_when((min_af == max_af & min_af >= 0.9) ~ 0, 
                                               (min_af == max_af & min_af <= 0.9) ~ max_af,
                                               min_af != max_af ~ max_af - min_af)) %>% ungroup()

# Filter out positions with sequencing artefacts
filt_variants_df <- filt_variants_df %>% filter(NTPOS > 100) %>% filter(NTPOS < 29840) %>% filter(!grepl("dup",snpid)) %>%
                           filter(NTPOS != 13402) %>% filter(NTPOS!=24389) %>% 
                            filter(NTPOS != 24390) %>% filter(NTPOS!=21530) %>% filter(NTPOS != 21534) %>% filter(NTPOS!=21535) %>% filter(NTPOS != 21537) %>%
                            filter(NTPOS != 2464) %>% filter(NTPOS!=2465)  %>% filter(NTPOS!=2471) %>% filter(NTPOS != 23076) %>% filter(NTPOS!=23077)  %>% filter(NTPOS!=23078) %>% 
                            filter(NTPOS!=4655)  %>% filter(NTPOS != 5570) %>% filter(NTPOS!=21536) %>%
                            filter(NTPOS!=5765)  %>% filter(NTPOS != 5766) %>% filter(NTPOS!=14335) %>%
                            filter(!(NTPOS %in% c(632, 633,635,636,4968,19814,25438, 6218, 6222,6225,6227,11992,11994,22524,22773,1566,4274,25337,4656,4658,4661,4662,12582,14914,29199,29443,29445,29447,28676, 16256,16258,18142,18146,18148,518,519,520,4031,4027,7481,7487,7488,7490,7491,7492,7493,7499,7500,8883,8885,8886,9401,15071,20308,20309,20310,21209,21212,24161,26017))) %>% 
                            # homopolymeric regions?
                            filter (NTPOS != 11081) %>% filter(NTPOS != 11082) %>% filter(NTPOS != 11083) %>%  filter(NTPOS != 6700) %>% filter(NTPOS != 19989) %>% filter(NTPOS != 29056)

# Filter for replicate 2
rep2_filtvariants <- rep2_variants_df %>% filter(NTPOS > 100) %>% filter(NTPOS < 29840) %>% filter(!grepl("dup",snpid)) %>%
  filter(NTPOS != 13402) %>% filter(NTPOS!=24389) %>% 
  filter(NTPOS != 24390) %>% filter(NTPOS!=21530) %>% filter(NTPOS != 21534) %>% filter(NTPOS!=21535) %>% filter(NTPOS != 21537) %>%
  filter(NTPOS != 2464) %>% filter(NTPOS!=2465)  %>% filter(NTPOS!=2471) %>% filter(NTPOS != 23076) %>% filter(NTPOS!=23077)  %>% filter(NTPOS!=23078) %>% 
  filter(NTPOS!=4655)  %>% filter(NTPOS != 5570) %>% filter(NTPOS!=21536) %>%
  filter(NTPOS!=5765)  %>% filter(NTPOS != 5766) %>% filter(NTPOS!=14335) %>%
  filter(!(NTPOS %in% c(632, 633,635,636,4968,19814,25438, 6218, 6222,6225,6227,11992,11994,22524,22773,1566,4274,25337,4656,4658,4661,4662,12582,14914,29199,29443,29445,29447,28676, 16256,16258,18142,18146,18148,518,519,520,4031,4027,7481,7487,7488,7490,7491,7492,7493,7499,7500,8883,8885,8886,9401,15071,20308,20309,20310,21209,21212,24161,26017))) %>% 
  # homopolymeric regions?
  filter (NTPOS != 11081) %>% filter(NTPOS != 11082) %>% filter(NTPOS != 11083) %>%  filter(NTPOS != 6700) %>% filter(NTPOS != 19989) %>% filter(NTPOS != 29056)


# get rid of indels for stats and fig2, they are too messy
#filt_variants_df <- filt_variants_df %>% filter(AASUB != "-") %>% filter(AASUB != "fs")

# Filter at various allele frequencies and depths
filt_variants_df_50depth <- filt_variants_df %>% filter(AAFREQ >= 0.01) %>% filter(TCOV >=50)
filt_variants_df_5af_100depth <- filt_variants_df %>% filter(AAFREQ >= 0.05) %>% filter(TCOV>=100) 
filt_variants_df_1af_100depth <- filt_variants_df %>% filter(AAFREQ >= 0.01) %>% filter(TCOV>=100) 
filt_variants_df_95 <- filt_variants_df_5af_100depth %>% filter(AAFREQ < 0.95)


## Calculating some stats

# How many unique variants between 5-95%?
length(unique(filt_variants_df_95$snpid))
# Gene  distribution?
length(which(filt_variants_df_95$gene == "ORF1ab_polyprotein"))
length(which(filt_variants_df_95$gene == "ORF1a_polyprotein"))
length(unique(which(filt_variants_df_95$gene == "surface_glycoprotein")))

nsp_summary <- filt_variants_df_95 %>% group_by(gene, nsp) %>% summarise(num_variants = n_distinct(snpid))
sum(nsp_summary$num_variants)

# Stats for 1-5%?
low_freq_variants <- filt_variants_df %>% filter(AAFREQ >= 0.01) %>% filter(AAFREQ<=0.05) %>% filter(VCOV >=10)
length(unique(low_freq_variants$snpid))

# Stats for num variants in samples
low_freq_variants_counts <- low_freq_variants %>% group_by(SAMPLE) %>% summarize(num_variants = n_distinct(snpid))
mean(low_freq_variants_counts$num_variants)
filt_variants_df_95_counts <- filt_variants_df_95 %>% group_by(SAMPLE) %>% summarize(num_variants = n_distinct(snpid))
median(filt_variants_df_95_counts$num_variants)

### FIGURE GENERATION ###
## Figure 1, adapted from Pavi's code
ct_values <- read_csv("longitudinal_cts.csv")
#Sampling and Ct
p1<-ggplot(subset(ct_values,!ptid%in%c('P004','P006','P012')),
           aes(x=t_elap_seq,y=as.numeric(ct),group=ptid))+
  geom_line(colour='grey')+
  geom_point()+
  scale_y_reverse(limits=c(42,10))+
  ylab('Ct')+
  xlab('Time since first sequenced sample (d)')+
  theme_clean()+
  theme(plot.background=element_blank())

#Reads and genome completeness
plot_data<-ct_values %>% filter(!SpID%in%c('11151','4107'))
plot_data<-plot_data %>% 
  mutate(seq_qual=ifelse(is.na(conNsperc),'NA',
                         ifelse(conNsperc<5,'>95%',
                                ifelse(conNsperc<10,'90-95%','<90%'))),
         ct_bin=case_when(ct<20 ~ 'under 20',
                          ct>=20&ct<25 ~ '20-25',
                          ct>=25&ct<30 ~ '25-30',
                          ct>=30 ~ 'over 30',
                          TRUE ~ 'NA')) %>% 
  as.data.frame()
plot_data$ct_bin<-factor(plot_data$ct_bin,
                         levels=c('under 20','20-25','25-30','over 30'))
plot_data$seq_qual<-factor(plot_data$seq_qual,
                           levels=c('>95%','90-95%','<90%','NA'))

# moved p2 to supplementary figs
p2<-ggplot(data=subset(plot_data,seq_qual!='NA'),aes(x=ct_bin,y=100-conNsperc))+
  geom_jitter(aes(colour=raw_reads/1e6),width=0.1,height=0,size=2)+
  geom_jitter(data=subset(plot_data,seq_qual=='NA'),aes(y=79),colour='darkgrey',shape=1,width=0.2,height=0.2,size=2)+
  scale_color_binned('raw reads\n(millions)',n.breaks=4)+
  ylim(78,100)+
  theme_clean()+
  ylab('Genome completeness (%)')+xlab('Ct')+
  theme(plot.background=element_blank(),
        legend.position='top',
        legend.background=element_blank())
plot_grid(p1, p2, labels = c('A','B'), ncol=2,label_size = 11)
ggsave("Fig1_ct.pdf",p1,height=3,width=4,units=c("in"))

## Supplementary fig: compare afs across replicates
rep1_filtvariants <- filt_variants_df %>% filter(TCOV >= 50) %>% filter(VCOV >= 10) %>% filter(AASUB != "fs")
rep2_filtvariants <- rep2_filtvariants %>% filter(TCOV >= 50) %>% filter(VCOV >= 10)%>% filter(AASUB != "fs")
rep1_filtvariants <- rep1_filtvariants %>% unite("pat_day",c("patient","day"),remove=FALSE)
rep2_filtvariants <- rep2_filtvariants %>% unite("pat_day",c("patient","day"),remove=FALSE)

#rep1_filtvariants <- rep1_filtvariants %>% filter(pat_day %in% c(unique(rep2_filtvariants$pat_day)))
#rep2_filtvariants <- rep2_filtvariants %>% filter(pat_day %in% c(unique(rep1_filtvariants$pat_day)))

compare_replicate_af <- merge(x=rep1_filtvariants, y = rep2_filtvariants, #all=T,
                              by=c("gene","AAPOS","AAREF","AASUB","NTPOS","snpid","nsp","NSPPOS","NSPREF","NSPSUB",
                                   "patient","day"))
#compare_replicate_af <- compare_replicate_af %>% filter(!patient %in% c("P020","P002","P004","P013"))
#compare_replicate_af[is.na(compare_replicate_af)] <- 0

replicates_plot <- ggplot(compare_replicate_af, aes(x=AAFREQ.x, y=AAFREQ.y)) + 
  geom_smooth(method = "lm", se = FALSE, color="mediumpurple1", fill="aquamarine",alpha=0.1,size=1.5)   + 
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color="black", linetype = "dotted", size=0.5, alpha = 1) + 
  geom_jitter(width=0.02,height=0.02, alpha=0.7) + 
    theme_clean() + labs(x = "Replicate 1 Allele Frequency", y = "Replicate 2 Allele Frequency") + 
  theme(plot.background = element_blank())
  
sup_fig <- plot_grid(p2, replicates_plot, labels=c("A","B"))#,align="h",axis="t")
ggsave("SFig_replicate_af_plot.pdf", sup_fig, height=4,width=8.5,units="in")


## Figure 2 : Combined figure with frequencies and marginal histogram, adapted from Pavi's code
fig2<-ggplot(filt_variants_df_1af_100depth,aes(x=(NTPOS),y=AAFREQ))+
  geom_point(aes(colour=abs(af_diff),size=log10(TCOV)),alpha=0.75)+
  scale_color_viridis_c(expression(paste('Maximum Change in AF')),
                        option="C",end=0.9,begin=0, direction=-1, breaks=c(0,0.25,0.5,0.75), labels=c("0.00","0.25","0.50","0.75"),limits=c(0,1)) + 
  #scale_color_viridis_c(expression(paste('abs(','\u0394','af)')),
  #                      option='C',direction=-1,end=.9,
  #                      breaks=c(0,0.25,0.5,0.75))+
  scale_size_binned('log10(depth)',range=c(0.1,4),limits=c(2,4.5),breaks=c(2,3,4))+
  xlab('Nt pos')+ylab('Variant allele frequency')+
  theme_clean()+
  scale_x_continuous(breaks=c(0,5000,10000,15000,20000,25000,30000))+
  geom_rug(data=filt_variants_df %>% filter(TCOV>=100&AAFREQ>=0.95),sides='top',
           alpha=.2,position='jitter',length=unit(0.05,'npc'),colour='black') +
  coord_cartesian(clip = "off") +
  theme(plot.margin=margin(1, 0, 0, 0, "cm"),
        plot.background=element_blank(),
        legend.box='vertical',
        legend.spacing=unit(0,'cm'),
        legend.background=element_rect(colour=NA), 
        legend.title = element_text(face="plain"))
fig2 <- ggMarginal(fig2,type='histogram',binwidth=500,size=2.5,margins='x')
ggsave("Fig2_plot.pdf", fig2, height=6,width=8,units="in")

## Looking at common variants and evolving variants
# Find most common variants in our dataset
filt_variants_df_1af_100depth$AAFREQ <- as.numeric(as.character(filt_variants_df_1af_100depth$AAFREQ))
common_variants <- filt_variants_df_1af_100depth %>% group_by(snpid,gene, AAREF, AAPOS, AASUB, nsp, NSPREF, NSPPOS, NSPSUB, .drop=FALSE) %>% 
      summarise(min_AF = round(min(AAFREQ),2), max_AF = round(max(AAFREQ),2), num_patients = n_distinct(patient), num_samples = n_distinct(SAMPLE)) %>%
      filter(num_samples >=15) %>% filter(max_AF >= 0.1)
common_variants <- common_variants[order(common_variants$num_samples),] %>% map_df(rev)
common_variants[] <- lapply(common_variants, gsub, pattern = "_", replacement = " ", fixed=TRUE)

# Formatting
common_variants <- common_variants %>% unite("AAChange", AAREF:AASUB, remove=TRUE, sep="") %>%
                    unite("AA Change", gene:AAChange, remove=TRUE,sep = " ") %>%
                    unite("nspchange", NSPREF:NSPSUB, remove=TRUE, sep="") %>%
                    unite("NSP Change", nsp:nspchange, remove=TRUE, sep=" ") %>%
                    unite("AF Range", min_AF:max_AF, remove=TRUE, sep = " - ")

write.table(common_variants, file = "common_variants.csv", sep = ",", row.names=FALSE, quote=FALSE)

# Find evolved variants (those with allele frequency change of >20%)
evolved_variants <- filt_variants_df_1af_100depth %>% filter(af_diff >= 0.15) %>% group_by(patient, snpid,gene, AAREF, AAPOS, AASUB, .drop=FALSE) %>%
                    # Pt6 and 12 don't have longitudinal data
                    filter(VCOV >= 10) %>% #filter(max_af >= 0.25) %>%
                    filter(patient != "P012") %>% filter(patient != "P006") %>%
                    filter(af_diff!=1) %>% unite("AAChange", c("AAREF","AAPOS","AASUB"), remove=TRUE, sep="") %>%
                    unite("AA Change", gene:AAChange, remove=TRUE,sep = " ") %>%
                    unite("nspchange", c("NSPREF","NSPPOS","NSPSUB"), remove=TRUE, sep="") %>%
                    unite("NSP Change", nsp:nspchange, remove=TRUE, sep=" ")

#write.table(evolved_variants, file = "evolved_variants2.csv", sep = ",", row.names=FALSE, quote=FALSE)

# evolved_variants$snpid <- factor(evolved_variants$snpid, levels=unique(evolved_variants$snpid[order(evolved_variants$patient,evolved_variants$NTPOS)]))
# 
# evolved_variants_plot <- ggplot(evolved_variants, aes(x=day, y=AAFREQ, color=snpid, group = snpid)) +
#   geom_point(alpha = 0.75) + geom_line(alpha=0.75) + facet_wrap(~ patient, ncol=5) +
#   theme_clean() + labs(x = "Days Since First Sample", y = "Allele Frequency") +
#   theme(legend.position = "None")
# 
# other_variants <- filt_variants_df_1af_100depth %>% filter(AAFREQ <= 0.75) %>% filter(AAFREQ >= 0.2) %>%
#   filter(af_diff<=0.2) %>%
#   filter(patient != "P012") %>% filter(patient != "P006") %>%
#   filter(af_diff!=1) %>% unite("AAChange", c("AAREF","AAPOS","AASUB"), remove=TRUE, sep="") %>%
#   unite("AA Change", gene:AAChange, remove=TRUE,sep = " ") %>%
#   unite("nspchange", c("NSPREF","NSPPOS","NSPSUB"), remove=TRUE, sep="") %>%
#   unite("NSP Change", nsp:nspchange, remove=TRUE, sep=" ")
# other_variants <- sample_n(other_variants,28) 
# 
# write.table(other_variants, file = "other_variants.csv", sep = ",", row.names=FALSE, quote=FALSE)

# Check against GISAID database here
full_gisaid_variants <- read_tsv("GISAID_variant_surveillance.tsv") # GISAID table downloaded ~April 1 2021
# evolved_variants table, slightly manually cleaned up
gisaid_variants <- read_csv("evolved_variants3.csv")

# Find relative frequency (988472 is total # of GISAID consensuses at time of download)
gisaid_variants$percent_GISAID <- gisaid_variants$num_gisaid_consensuses / 988472 * 100
gisaid_variants_plot <- ggplot(gisaid_variants,aes(x=NTPOS,y=percent_GISAID, size=af_diff)) + geom_point(alpha = 0.75, color="darkslateblue") +
                        theme_clean() +
                        scale_x_continuous(breaks=c(0,5000,10000,15000,20000,25000,30000))+
                        scale_y_continuous(limits=c(0,0.25)) +
                        scale_size_continuous(limits=c(0.2,0.8),breaks=c(0.8,0.6,0.4,0.2)) + 
                        #scale_color_viridis(discrete=FALSE, limits = c(0.2,0.8),breaks=c(0.2,0.4,0.6,0.8)) +
                        theme(legend.background=element_rect(colour=NA,fill=NA), legend.title = element_text(face="plain"), 
                              legend.spacing.y = unit(0,"cm"), plot.background=element_blank()) +
                        labs(size="Maximum Change in\nAllele Frequency\n",#, color = "",
                            y = "Relative Frequency of Variant\nin GISAID Consensuses (%)", x="Nucleotide Position of Variant") + 
                        guides(size=guide_legend(order=1))

# For a subset of our GISAID_variants (>20% AF difference), let's pull the consensus sequence matches
gisaid_variants <- gisaid_variants %>% filter(af_diff >= 0.4)
gisaid_consensuses <- data.frame()
for(variant in unique(gisaid_variants$gisaid_mut)) {
  print(paste("Looking for GISAID consensuses for",variant,"...",sep=" "))
  subset <- filter(full_gisaid_variants, grepl(variant,`AA Substitutions`))
  if(variant=="Spike_V143del") {
    subset <- filter(subset, grepl("Spike_Y144del",`AA Substitutions`))
  }
  subset$query_variant <- variant
  print(nrow(subset))
  subset$num_gisaid_samps <- nrow(subset)
  gisaid_consensuses <- rbind(gisaid_consensuses,subset)
}

#fresh_gisaid_consensuses <- gisaid_consensuses
#gisaid_consensuses <- fresh_gisaid_consensuses
gisaid_consensuses$`Collection date` <- as.POSIXct.Date(as.Date(gisaid_consensuses$`Collection date`))
gisaid_consensuses <- gisaid_consensuses[!is.na(gisaid_consensuses$`Collection date`),]

gisaid_consensuses <- separate(gisaid_consensuses, Location, into=c("continent","region","subregion"), sep=" / ") %>%
  mutate(region = recode(region,"United Kingdom" = "UK"), subregion = recode(subregion,"England" = "Great Britain"))

gisaid_consensuses <- gisaid_consensuses %>% group_by(query_variant, `Collection date`) %>% add_tally(name = "Num_Date")
gisaid_consensuses <- separate(gisaid_consensuses, `Collection date`, c("Year","Month","Day"), remove = FALSE)
gisaid_consensuses <- gisaid_consensuses %>% group_by(Year, Month,query_variant) %>% mutate(Num_Month = n_distinct(`Accession ID`))
gisaid_consensuses <- unite(gisaid_consensuses, Month_Date, c("Month","Year"), sep="/")

gisaid_consensuses$query_variant <- factor(gisaid_consensuses$query_variant, levels = c("NSP3_T749I","NSP3_T1036I","NSP12_M601I","NSP12_A660S","NSP12_G852D","NS3_Y211N","E_S68F","NS8_H40Y"))
gisaid_consensuses$query_variant <- gsub("_", ": ", gisaid_consensuses$query_variant)
gisaid_consensuses$percent_GISAID <- gisaid_consensuses$num_gisaid_samps / 988472 * 100

# ridges plot
# gisaid_consensuses_plot <- ggplot(data=gisaid_consensuses,aes(x=`Collection date`,y=query_variant, 
#                                                               height = stat(density), fill = stat(density))) + 
#   geom_density_ridges_gradient(alpha = 0.75, stat="density") +
#   scale_fill_viridis_c(name = "# Samples", option="C") + labs(y="Variant") + theme_clean() + theme(legend.position = "None")
# ggsave("GISAID_consensuses_ridgeplots.pdf",gisaid_consensuses_plot, height=8, width=6, units="in")
# 
# # categorical scatter? nope
# gisaid_consensuses_plot <- ggplot(data=gisaid_consensuses,aes(x=`Collection date`,y=query_variant, color = query_variant)) + 
#   geom_jitter(width=0.1,height=0.3, alpha = 0.5) + scale_color_cyclical(values = c("#4040B0", "#9090F0")) + 
#   theme_clean() + theme(legend.position = "None")
# 
# # area chart
# gisaid_consensuses_plot <- ggplot(data=gisaid_consensuses,aes(x=`Collection date`,fill=query_variant, color=query_variant)) + 
#                             geom_histogram(alpha=0.7, bins = 13) + theme_clean() + theme(legend.position = "right") + 
#                             scale_fill_brewer(palette="Dark2") + scale_color_brewer(palette="Dark2")



gisaid_consensuses <- gisaid_consensuses %>% filter(Month_Date != "04/2021") %>% filter(Month_Date != "NA")
date_list <- c("03/2020","04/2020","05/2020","06/2020","07/2020","08/2020","09/2020","10/2020","11/2020","12/2020","01/2021","02/2021","03/2021")
gisaid_consensuses$Month_Date <- factor(gisaid_consensuses$Month_Date, levels=date_list)
gisaid_consensuses$query_variant <- factor(gisaid_consensuses$query_variant, levels = c("NSP3: T749I","NSP3: T1036I","NSP12: M601I","NSP12: A660S","NSP12: G852D","NS3: Y211N","E: S68F","NS8: H40Y"))

gisaid_consensuses <- gisaid_consensuses %>% group_by(query_variant) %>% mutate(num_continents = n_distinct(continent))
gisaid_consensuses <- gisaid_consensuses %>% group_by(query_variant) %>% mutate(num_subregions = n_distinct(subregion))
gisaid_consensuses <- gisaid_consensuses %>% group_by(query_variant,Month_Date) %>% mutate(num_continents_month = n_distinct(continent))


gisaid_consensuses_plot <- ggplot(data=gisaid_consensuses,aes(x=Month_Date, fill=num_continents_month)) + geom_bar(alpha=0.8) + 
                            facet_wrap(~ query_variant,ncol=4) + theme_clean() + scale_y_continuous(labels = scales::comma) + 
                            scale_x_discrete(breaks = c("03/2020","06/2020","09/2020","12/2020","03/2021")) + 
                            scale_fill_viridis(limits = c(1,7), breaks = c(1,3,5,7)) + 
                            theme(strip.text.y.right = element_text(angle = 0), panel.spacing = unit(1,"lines"), plot.background=element_blank(),
                                  legend.background = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) + 
                            labs(fill = "Total #\nContinents", x = "Month", y = "# of GISAID Consensuses")

fig3_gisaid <- plot_grid(gisaid_variants_plot, gisaid_consensuses_plot, labels=c("A","B"),ncol=1, rel_heights = c(0.8, 1))
ggsave("GISAID_fig3.pdf",fig3_gisaid,height=9,width=8.5,units="in")

# gisaid_consensuses_plot <- ggplot(data=gisaid_consensuses,aes(x=`Collection date`)) + 
#                             facet_grid(rows =vars(query_variant)) + 
#                             geom_histogram(alpha=0.7,color="white",bins=13) + 
#                             theme_clean() + theme(legend.position = "None") + 
#                             theme(strip.text.y.right = element_text(angle=0)) + 
#                             scale_x_datetime(breaks = date_breaks("3 months"), labels = date_format("%Y-%b"))


median(unique(gisaid_consensuses$num_gisaid_samps))

#spike v143del, y144del

ggsave("GISAID_fig_plot.pdf", gisaid_variants_plot, height=5,width=6,units="in")

## num_variants
num_variants_1 <- filt_variants_df %>% group_by(patient,day) %>% summarise(num_variants = n_distinct(snpid))
num_variants_5 <- filt_variants_df %>% filter(AAFREQ>=0.05)  %>% group_by(patient,day) %>% summarise(num_variants = n_distinct(snpid))
num_variants_1$freq <- ">1%"
num_variants_5$freq <- ">5%"
num_variants <- rbind(num_variants_1, num_variants_5)
ggplot(num_variants, aes(x=day,y=num_variants, color=freq)) + geom_point() + geom_smooth()

## the deletion variant that evolved >20%, spike in pt16
spike_deletions <- Biostrings::readDNAStringSet("pt16_spike_deletion_variants.fasta")
variants_order <- c("NC_045512.2", "P016_0_999_97.50","P016_0_110_1.87","P016_0_12_0.2","P016_0_9_0.15","P016_0_9_0.14","P016_0_5_0.09",
                    "P016_0_3_0.05","P016_3_999_78.1","P016_3_930_20.7","P016_3_52_1.2","P016_6_999_99.8")

spike_deletions_plot <- ggmsa(spike_deletions, none_bg = TRUE, char_width = 0.7, order = c(variants_order), seq_name = TRUE, consensus_views = TRUE, ref="NC_045512.2") + 
  theme(axis.text.y = element_blank()) #+ ylim(0,9.)

spike_dels <- separate(data.frame(variants_order), variants_order, c("Patient","Day","Count","Rf"), sep = "_") %>% filter(Patient != "NC")
spike_dels$y_axis <- 11-1:nrow(spike_dels) + 0.5

spike_dels$Rf <- as.numeric(as.character(spike_dels$Rf))
spike_dels_bar <- ggplot(spike_dels, aes(x=y_axis, y = Rf, fill=Day)) + geom_col() + coord_flip() + theme_nothing() + 
    scale_x_continuous(limits=c(-1,12)) + scale_y_sqrt(limits=c(0,100)) + 
    scale_fill_brewer(palette="Dark2") + 
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    geom_text(aes(label=Rf),position = position_dodge(width = 0.9), hjust = -0.20)

#plot_grid(spike_deletions_plot, spike_dels_bar, align = "h", ncol = 2)

ggsave("spike_deletions2.pdf",spike_deletions_plot, units = "in", width = 6, height = 4)
ggsave("spike_dels_bar2.pdf",spike_dels_bar, units = "in", width = 2, height=2)


## Let's look at evolving variants in spike
## Only take pts with min spike coverage > 100 (no 4, 5, 7,)
spike_variants <- filt_variants_df_50depth %>% filter(gene=="surface_glycoprotein")
spike_variants$AAPOS <- as.numeric(as.character(spike_variants$AAPOS))
spike_variants <- (spike_variants[,!names(spike_variants) %in% c("min_af", "max_af","af_diff")])
spike_variants <- distinct(spike_variants)

# Fill empty days with 0
fill_with_zero <- function(df, day_list, patient_name) {
  df2 = df
  for(aa_change in unique(df$snpid)) {
    subset <- df %>% filter(snpid==aa_change)
    nuc_pos <- subset$NTPOS[[1]]
    sample_name <- subset$SAMPLE[[1]]
    aapos <- subset$AAPOS[[1]]
    aaref <- subset$AAREF[[1]]
    aasub <- subset$AASUB[[1]]
    for(day in day_list) {
      if(!(day %in% subset$day)) {
        to_add <- data.frame("SAMPLE" = patient_name, "gene" = "surface_glycoprotein", "AAPOS" = aapos, "AAREF" = aaref, "AASUB" = aasub, "TCOV" = 0, "VCOV" = 0, "AAFREQ" = 0, "NTPOS" = nuc_pos,
                             "snpid" = aa_change, "nsp" = "", "NSPPOS" = "", "NSPREF" = "", "NSPSUB" = "", "patient" = patient_name, day = day)
        df2 <- rbind(df2,to_add)
      }
    }
  }
  return(df2)
}

for(patient_name in unique(variants_df$patient)) { 
  print(patient_name)
  day_list = unique(variants_df$day[variants_df$patient == patient_name])
  print(day_list)
  sample_subset <- spike_variants %>% filter(patient == patient_name)
  spike_variants <- rbind(spike_variants, fill_with_zero(sample_subset, day_list, patient_name))
}

spike_variants <- distinct(spike_variants)
a <- data.frame()

# Grab correct change in AF
for(patient_name in unique(variants_df$patient)) {
  print(patient_name)
  day_list = unique(variants_df$day[variants_df$patient == patient_name])
  min_day <- (min(day_list))
  max_day <- (max(day_list))
  subset <- spike_variants %>% filter(patient == patient_name)
  for(unique_variant in unique(subset$snpid)) {
    subset2 <- subset %>% filter(snpid == unique_variant)
    af_diff <- subset2$AAFREQ[subset2$day==max_day] - subset2$AAFREQ[subset2$day==min_day]
    subset2$af_diff <- af_diff
    a <- rbind(a,subset2)
  }
}

spike_variants <- a

#spike_variants$plasma[spike_variants$patient=="P001"] <- "TRUE"
#spike_variants$plasma[spike_variants$patient=="P018"] <- "TRUE"

# P012 and P006 don't have longitudinal data; P011 later timepoints no bueno coverage
spike_variants <- spike_variants %>% filter(patient!="P012") %>% filter(patient != "P006") %>% filter(patient != "P004")
                  #filter(patient != "P011") %>% 

spike_variants_filt <- spike_variants %>% filter(AAFREQ >= 0.05)
spike_variants_filt$TCOV <- as.numeric(as.character(spike_variants_filt$TCOV))
spike_variants_filt <- spike_variants_filt %>% filter(TCOV >= 10)
spike_variants_filt <- spike_variants_filt %>% group_by(snpid,patient) %>% add_tally() 

change_in_af <- ggplot(spike_variants, aes(x=AAPOS, y=AAFREQ, color=af_diff, group = snpid)) + 
  #scale_colour_brewer(type = "div", palette = "RdBu", direction=-1) + 
  scale_color_gradient2_tableau(palette="Classic Area Red-Green",limits=c(-1,1)) + 
  geom_point() + geom_line() + facet_wrap(~ patient, ncol=5) +
  theme_clean() + labs(color = "Net Change\n in AF\n", x = "Amino Acid Residue", y = "Allele Frequency") + 
  scale_x_continuous(breaks=c(0,500,1000,1500))+
  theme(plot.margin=margin(1, 0, 0, 0, "cm"),
        plot.background=element_blank(),
        legend.box='vertical',
        #legend.spacing=unit(0,'cm'),
        legend.background=element_rect(colour=NA), 
        legend.title = element_text(face="plain", size = 10),
        legend.title.align = 0.5,
        legend.text = element_text(size = 10))

ggsave("spike_change_in_af.pdf",change_in_af, units = "in", width = 8, height = 6)

spike_variants_over20 <- spike_variants %>% group_by(patient,snpid) %>% 
                      mutate(min_af = min(AAFREQ), max_af = max(AAFREQ), 
                      af_diff = case_when((min_af == max_af & min_af >= 0.9) ~ 0, 
                               (min_af == max_af & min_af <= 0.9) ~ max_af,
                                min_af != max_af ~ max_af - min_af)) %>% ungroup() %>%
                      filter(af_diff >= 0.2)

# manually reviewed for accuracy
spike_variants_over20 <- spike_variants_over20 %>% filter(!(AAPOS %in% c("1259","321","404","140","320","398","399","804","867", "143")))
spike_variants_over20 <- spike_variants_over20 %>% unite("AA_Change",c("AAREF","AAPOS","AASUB"),remove=TRUE, sep="")

spike_variants_over20$AA_Change[spike_variants_over20$AA_Change=="I144-"] <- "143-145del"

spike_variants_over20$AA_Change <- factor(spike_variants_over20$AA_Change, levels=unique(spike_variants_over20$AA_Change[order(spike_variants_over20$patient)]))

spike_mutations_plot <- ggplot(spike_variants_over20, aes(x=day, y=AAFREQ, color=AA_Change, group = AA_Change)) + 
  #scale_color_manual(values=c("black","red"), aesthetics="color") +
  geom_point(alpha = 0.75) + geom_line(alpha=0.75) + facet_wrap(~ patient, ncol=6) +
  theme_clean() + labs(x = "Days Since First Sample", y = "Allele Frequency") + 
  theme(legend.position = "right", legend.background=element_rect(colour=NA), legend.title = element_blank())

ggsave("spike_mutations_changes.pdf", spike_mutations_plot, units = "in", width = 8, height = 2)

