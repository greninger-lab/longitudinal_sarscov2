library(tidyverse)
library(ggthemes) 
library(cowplot)



rm(list=ls())
library(tidyverse); library(ggthemes) ; library(cowplot)
options(dplyr.width=Inf)

setwd('/Users/vikas/Downloads/scratch/9_1_combined/')

long_metadata<-read.csv('/Users/vikas/Downloads/scratch/9_1_combined/longitudinal_seqd_filtered_2020-09-22_withpid.csv',stringsAsFactors=F) %>% 
  mutate(t_int=as.character(round(t_elapsed))) %>% 
  filter(SpID!='4097')
all_rpms<-read.csv('/Users/vikas/Downloads/scratch/9_1_combined/all_rpms.csv',stringsAsFactors=F)
clin_rel<-read.csv('/Users/vikas/Downloads/scratch/9_1_combined/9_1_combined_clinically_relevant.csv',stringsAsFactors=F) %>% 
  pivot_longer(cols = starts_with("P0"),
               names_to='pat_clomp',
               values_to='rpm') %>% 
  mutate(ptid=str_split(pat_clomp,'_',simplify=T)[,1],
         t=gsub('d','',str_split(pat_clomp,'_',simplify=T)[,2]))

read_breakdown<-read_csv('/Users/vikas/Downloads/scratch/9_1_combined/read_comparison.csv') %>% 
  # mutate(`Other classified`=`Classified reads`-(
  #   `Fungal reads`+`Bacterial reads`+`Viral reads`+`Protozoan reads`+`Artificial reads`)) %>% 
   select(-c(`Classified reads`,`Microbial reads`, `Protozoan reads`, `Artificial reads`)) %>% 
  pivot_longer(cols=-Name,
               names_to='class',
               values_to='reads')%>% 
  mutate(ptid=str_split(Name,'_',simplify=T)[,1],
         t=gsub('d','',str_split(Name,'_',simplify=T)[,2]))

clomp_merged<-left_join(long_metadata,read_breakdown,by=c('ptid','t_int'='t')) %>% 
  group_by(ptid,t_int)

clin_rel$t_seq<-''
for(i in 1:nrow(clin_rel)){
  clin_rel$t_seq[i]<-ifelse(grepl('P006_0d_a|P012_0d_a',clin_rel$pat_clomp[i]),1,
                            ifelse(grepl('P006_0d_b|P012_0d_b',clin_rel$pat_clomp[i]),2,
                                   long_metadata$t_seq[long_metadata$ptid==clin_rel$ptid[i]&
                                                         long_metadata$t_int==clin_rel$t[i]]))
}
read_breakdown$t_seq<-''
for(i in 1:nrow(read_breakdown)){
  read_breakdown$t_seq[i]<-ifelse(grepl('P006_0d_a|P012_0d_a',read_breakdown$Name[i]),1,
                                  ifelse(grepl('P006_0d_b|P012_0d_b',read_breakdown$Name[i]),2,
                                         long_metadata$t_seq[long_metadata$ptid==read_breakdown$ptid[i]&
                                                               long_metadata$t_int==read_breakdown$t[i]]))
}



# Panel A 
total_reads<-ggplot(read_breakdown %>% filter(class=='Number of raw reads'),aes(x=t_seq,y=reads))+
  geom_bar(stat='identity',fill='darkred')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  ylab('Number of total reads')+
  facet_wrap(~ptid,scales='free_x',nrow=1)+
  theme_clean()+
  theme(plot.background=element_blank(),
        strip.text = element_text(size=5),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank()) +
  scale_y_log10(expand = c(0, 0), limits = c(1,1e8)) 
total_reads
dev.off()

# Panel B 
plot_data<-read_breakdown %>% 
  group_by(Name) %>% 
  mutate(perc_raw=reads/reads[class=='Number of raw reads']) %>% 
  filter(class!='Number of raw reads')
read_breakdown_plot<-ggplot(plot_data ,aes(x=t_seq,y=reads))+
  geom_bar(aes(fill=class),stat='identity',position='fill')+
  ylab('% of raw reads')+
  facet_wrap(~ptid,scales='free_x',nrow=1)+
  theme_clean()+
  scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.background=element_blank(),
        legend.position='bottom',
        legend.key.size=unit(5,'pt'),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.text = element_text(size=5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7))
read_breakdown_plot
dev.off()

# Panel C 
clin_rel %>% 
  group_by(name) %>% 
  summarize(total_rpm_all=sum(rpm),
            n_pts=length(unique(ptid[rpm>10])),
            n_samps=length(unique(pat_clomp[rpm>10]))) %>% 
  arrange(desc(total_rpm_all)) %>% 
  as.data.frame()


ptids<-c('P004','P005','P006','P015')
plot_data<-clin_rel %>% 
  filter(name%in%c('Staphylococcus aureus','Moraxella catarrhalis')) %>% 
  filter(ptid%in%ptids)
plot_data$rpm[plot_data$rpm<1]<-0
#png('Clomp_longitudinal.png',width=4,height=2,res=300,units='in')
clinically_relevant_plot<-ggplot(plot_data,aes(x=as.factor(t_seq),y=rpm,group=name))+
  geom_point(aes(colour=name,shape=name),size=2)+
  geom_line(aes(colour=name))+xlab('Sample number')+
  scale_y_log10(limits=c(1,1e4))+
  theme_clean()+
  facet_wrap(~ptid,drop=T,scales='free_x',nrow=2)+
  scale_color_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ylab('RPM') + 
  theme(legend.position='bottom',
        plot.background=element_blank(),
        legend.text=element_text(size=9),
        legend.title=element_blank(),
        legend.background=element_rect(color=NA),
        legend.box.spacing=unit(0,'pt'))
clinically_relevant_plot
dev.off()



# Combining plots and saving

panel<-plot_grid(read_breakdown_plot, total_reads,ncol = 1, nrow = 2, labels = c('A','B'))
panel

AB<-align_plots(read_breakdown_plot, total_reads, align = 'v', axis = 'l')
bottom_row<-plot_grid(AB[[2]], AB[[1]], ncol = 1, nrow = 2, labels = c('A','B'))
bottom_row

new_panel<-plot_grid(bottom_row, clinically_relevant_plot, rel_widths = c(3,1), labels = c('', 'C'))
new_panel


#VIKAS: latest revision (8/3/20)
final_panel<-plot_grid(read_breakdown_plot,clinically_relevant_plot,ncol = 2, nrow = 1, labels = c('A','B'),rel_widths = c(3,1))
final_panel
ggsave(plot = final_panel, filename = 'draft_panel_4.pdf', height = 5, width = 9)


supp_plot<-total_reads
ggsave(plot = supp_plot, filename = 'draft_supp_1.pdf', height = 5, width = 9)


