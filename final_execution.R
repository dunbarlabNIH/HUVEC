setwd("~/Documents/NIH/HUVEC/final")
source("final_functions.R")

load_new_space=FALSE
if(load_new_space){
  rm(list=ls())
  source("final_functions.R")
  data.zj25.lib7=read.combined("ZJ25_data/ZJ25_combined_lib7_final.txt")
  data.zj25.lib10=read.combined("ZJ25_data/ZJ25_combined_lib10_final.txt")
  data.zj25.lib9=read.combined("ZJ25_data/ZJ25_combined_lib9_final.txt")
  
  data.zj25.lib7=barcodetrackR::threshold(data.zj25.lib7)
  data.zj25.lib10=barcodetrackR::threshold(data.zj25.lib10)
  data.zj25.lib9=barcodetrackR::threshold(data.zj25.lib9)
  
  readme.zj25.lib7=read.delim("ZJ25_data/ZJ25_combined_lib7_final.txt_README.txt",skip=3,header=TRUE,row.names=1)
  readme.zj25.lib10=read.delim("ZJ25_data/ZJ25_combined_lib10_final.txt_README.txt",skip=3,header=TRUE,row.names=1)
  readme.zj25.lib9=read.delim("ZJ25_data/ZJ25_combined_lib9_final.txt_README.txt",skip=3,header=TRUE,row.names=1)
  
  zj25.key=read.delim("ZJ25_data/ZJ25_key_final.txt",row.names=NULL,header=TRUE)
  zj25.samples=as.vector(read.delim("ZJ25_data/ZJ25_core_samples.txt", header=FALSE)$V1)
  zj25.samples.filename=filename(zj25.samples,zj25.key)
  zj25.t=zj25.samples.filename[1:7];zj25.b=zj25.samples.filename[8:14];zj25.mono=zj25.samples.filename[15:21]
  zj25.nk=zj25.samples.filename[22:28];zj25.grans=zj25.samples.filename[29:35]
  zj25.time=c(1,2,4.5,5,20,44,51)
  
  data.m11.lib7=read.combined("M11021072_data/M11021072_combined_lib7_final.txt")
  data.m11.lib11=read.combined("M11021072_data/M11021072_combined_lib11_final.txt")
  
  m11.key=read.delim("M11021072_data/M11021072_key_final.txt",row.names=NULL,header=T)
  m11.samples=as.vector(read.delim("M11021072_data/M11021072_core_samples.txt", header=FALSE)$V1)
  m11.samples.filename=filename(m11.samples,m11.key)
  m11.t=m11.samples.filename[1:14];m11.b=m11.samples.filename[15:28];m11.mono=m11.samples.filename[29:42]
  m11.nk=m11.samples.filename[43:56];m11.grans=m11.samples.filename[57:70]
  m11.time=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
  
  data.m11.lib7=barcodetrackR::threshold(data.m11.lib7[,c(m11.samples.filename,"M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq","M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq","M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq","M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq","M11021072_15m__Grpos__BMR_CD175_sampled_I527_S22_L001_R1_001.fastq","M11021072_15m__Grpos__BML_CD175_sampled_I528_S23_L001_R1_001.fastq","M11021072_4m__Grans_BM_LT19_sampled_LT19_i519_S39_L005_R1_001.fastq")])
  data.m11.lib11=barcodetrackR::threshold(data.m11.lib11[,c(m11.samples.filename,"M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq","M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq","M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq","M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq","M11021072_15m__Grpos__BMR_CD175_sampled_I527_S22_L001_R1_001.fastq","M11021072_15m__Grpos__BML_CD175_sampled_I528_S23_L001_R1_001.fastq","M11021072_4m__Grans_BM_LT19_sampled_LT19_i519_S39_L005_R1_001.fastq")])
  
  readme.m11.lib7=read.delim("M11021072_data/M11021072_combined_lib7_final.txt_README.txt",skip=3,header=TRUE,row.names=1)
  readme.m11.lib11=read.delim("M11021072_data/M11021072_combined_lib11_final.txt_README.txt",skip=3,header=TRUE,row.names=1)
  
  for(readme in ls(pattern="^readme")){
    temp=get(readme)
    
    temp=temp[!is.na(rownames(temp)),]
    temp=temp[-nrow(temp),]
    
    assign(readme,temp)
  }
  rm(readme)
  
  save.image("final_data.Rdata")

}else{
  load("final_data.Rdata")
}

#############################################################################################
#Figure 3A#
#############################################################################################

pdf("figures/f3a_zj25_lib9.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.zj25.lib9[,zj25.samples.filename],names=zj25.samples,label_size=10,n_clones=20,grid=FALSE,cellnote_size = 0)
dev.off()
pdf("figures/f3a_zj25_lib7.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.zj25.lib7[,zj25.samples.filename],names=zj25.samples,label_size=10,n_clones = 20,grid=FALSE,cellnote_size = 0)
dev.off()

pdf("figures/f3a_m11_lib11.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.m11.lib11[,m11.samples.filename],names=m11.samples,label_size=10,n_clones=20,grid=FALSE,cellnote_size = 0)
dev.off()
pdf("figures/f3a_m11_lib7.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.m11.lib7[,m11.samples.filename],names=m11.samples,label_size=10,n_clones = 20,grid=FALSE,cellnote_size = 0)
dev.off()

#############################################################################################
#Figure 3B#
#############################################################################################

zj25.combined_data_list=list(data.zj25.lib7,data.zj25.lib9)
zj25.combined_readme_list=list(readme.zj25.lib7,readme.zj25.lib9)
pdf("figures/f3b_zj25.pdf",height=20,width=20)
barcode_Xlib(combined_data_list=zj25.combined_data_list,combined_readme_list=zj25.combined_readme_list,
             samples=zj25.samples,samples.filenames=zj25.samples.filename)
dev.off()

m11.combined_data_list=list(data.m11.lib7,data.m11.lib11)
m11.combined_readme_list=list(readme.m11.lib7,readme.m11.lib11)
pdf("figures/f3b_m11.pdf",height=20,width=20)
barcode_Xlib(combined_data_list=m11.combined_data_list,combined_readme_list=m11.combined_readme_list,
             samples=m11.samples,samples.filenames=m11.samples.filename)
dev.off()

#############################################################################################
#Figure 3C#
#############################################################################################
pdf("figures/f3c_zj25_lib7.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.zj25.lib7[,c("SP9_R28_ZJ25_4_5m_CD34pos_leftBM.fastq",
                                                   "SP9_R43_ZJ25_4_5m_grans_leftBM.fastq",
                                                  "SP9_R29_ZJ25_4_5m_CD34pos_rightBM.fastq",
                                                  "SP9_R44_ZJ25_4_5m_grans_rightBM.fastq",
                                                  "ZJ25_51m__CD34pos__BML_CD175_sampled_I507_S7_L001_R1_001.fastq",
                                                  "ZJ25_51m__Grpos__BML_CD175_sampled_I514_S11_L001_R1_001.fastq",
                                                  "ZJ25_51m__CD34pos__BMR_CD175_sampled_I506_S6_L001_R1_001.fastq",
                                                  "ZJ25_51m__Grpos__BMR_CD175_sampled_I512_S10_L001_R1_001.fastq")],
                                 n_clones=10,grid = FALSE,cellnote_size = 0,label_size = 10,
                                 names=c("4.5m BML","4.5m BML Grans","4.5m BMR","4.5m BMR Grans","51m BML","51m BML Grans","51m BMR","51m BMR Grans"))
dev.off()

pdf("figures/f3c_zj25_lib9.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.zj25.lib9[,c("SP9_R28_ZJ25_4_5m_CD34pos_leftBM.fastq",
                                                   "SP9_R43_ZJ25_4_5m_grans_leftBM.fastq",
                                                   "SP9_R29_ZJ25_4_5m_CD34pos_rightBM.fastq",
                                                   "SP9_R44_ZJ25_4_5m_grans_rightBM.fastq",
                                                   "ZJ25_51m__CD34pos__BML_CD175_sampled_I507_S7_L001_R1_001.fastq",
                                                   "ZJ25_51m__Grpos__BML_CD175_sampled_I514_S11_L001_R1_001.fastq",
                                                   "ZJ25_51m__CD34pos__BMR_CD175_sampled_I506_S6_L001_R1_001.fastq",
                                                   "ZJ25_51m__Grpos__BMR_CD175_sampled_I512_S10_L001_R1_001.fastq")],
                                 n_clones=10,grid = FALSE,cellnote_size = 0,label_size = 10,
                                 names=c("4.5m BML","4.5m BML Grans","4.5m BMR","4.5m BMR Grans","51m BML","51m BML Grans","51m BMR","51m BMR Grans"))
dev.off()

pdf("figures/f3c_m11_lib7.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.m11.lib7[,c("M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq",
                                                   "M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq",
                                                   "M11021072_4m__Grans_BM_LT19_sampled_LT19_i519_S39_L005_R1_001.fastq",
                                                   "M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq",
                                                  "M11021072_15m__Grpos__BML_CD175_sampled_I528_S23_L001_R1_001.fastq",
                                                   "M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq",
                                                   "M11021072_15m__Grpos__BMR_CD175_sampled_I527_S22_L001_R1_001.fastq"
                                                  )],
                                 n_clones=10,grid = FALSE,cellnote_size = 0,label_size = 10,
                                 names=c("4m BML","4m BMR","4m BM Grans","15m BML","15m BML Grans","15m BMR","15m BMR Grans"))
dev.off()

pdf("figures/f3c_m11_lib11.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.m11.lib11[,c("M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq",
                                                  "M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq",
                                                  "M11021072_4m__Grans_BM_LT19_sampled_LT19_i519_S39_L005_R1_001.fastq",
                                                  "M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq",
                                                  "M11021072_15m__Grpos__BML_CD175_sampled_I528_S23_L001_R1_001.fastq",
                                                  "M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq",
                                                  "M11021072_15m__Grpos__BMR_CD175_sampled_I527_S22_L001_R1_001.fastq")],
                              n_clones=10,grid = FALSE,cellnote_size = 0,label_size = 10,
                              names=c("4m BML","4m BMR","4m BM Grans","15m BML","15m BML Grans","15m BMR","15m BMR Grans"))
dev.off()

#############################################################################################
#Figure 4A#
#############################################################################################
pdf("figures/f4a_zj25_t.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,zj25.t],"libX"=data.zj25.lib9[,zj25.t]),zj25.time)
dev.off()

pdf("figures/f4a_zj25_b.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,zj25.b],"libX"=data.zj25.lib9[,zj25.b]),zj25.time)
dev.off()

pdf("figures/f4a_zj25_mono.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,zj25.mono],"libX"=data.zj25.lib9[,zj25.mono]),zj25.time)
dev.off()

pdf("figures/f4a_m11_t.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,m11.t],"libX"=data.m11.lib11[,m11.t]),m11.time)
dev.off()

pdf("figures/f4a_m11_b.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,m11.b],"libX"=data.m11.lib11[,m11.b]),m11.time)
dev.off()

pdf("figures/f4a_m11_mono.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,m11.mono],"libX"=data.m11.lib11[,m11.mono]),m11.time)
dev.off()
#############################################################################################
#Figure 4B#
#############################################################################################

pdf("figures/f4b_zj25_nk.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,zj25.nk],"libX"=data.zj25.lib9[,zj25.nk]),zj25.time)
dev.off()

pdf("figures/f4b_zj25_grans.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,zj25.grans],"libX"=data.zj25.lib9[,zj25.grans]),zj25.time)
dev.off()

pdf("figures/f4b_m11_nk.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,m11.nk],"libX"=data.m11.lib11[,m11.nk]),m11.time)
dev.off()

pdf("figures/f4a_m11_grans.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,m11.grans],"libX"=data.m11.lib11[,m11.grans]),m11.time)
dev.off()

pdf("figures/f4a_zj25_4_cd34L.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,rep("SP9_R28_ZJ25_4_5m_CD34pos_leftBM.fastq",2)],"libX"=data.zj25.lib9[,rep("SP9_R28_ZJ25_4_5m_CD34pos_leftBM.fastq",2)]),c(1,max(zj25.time)))+xlab("Time=4.5")
dev.off()

pdf("figures/f4a_zj25_4_cd34R.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,rep("SP9_R29_ZJ25_4_5m_CD34pos_rightBM.fastq",2)],"libX"=data.zj25.lib9[,rep("SP9_R29_ZJ25_4_5m_CD34pos_rightBM.fastq",2)]),c(1,max(zj25.time)))+xlab("Time=4.5")
dev.off()

pdf("figures/f4a_zj25_51_cd34L.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,rep("ZJ25_51m__CD34pos__BML_CD175_sampled_I507_S7_L001_R1_001.fastq",2)],"libX"=data.zj25.lib9[,rep("ZJ25_51m__CD34pos__BML_CD175_sampled_I507_S7_L001_R1_001.fastq",2)]),c(1,max(zj25.time)))+xlab("Time=51")
dev.off()

pdf("figures/f4a_zj25_51_cd34R.pdf",height=20,width=20)
stacked_2(list("lib7"=data.zj25.lib7[,rep("ZJ25_51m__CD34pos__BMR_CD175_sampled_I506_S6_L001_R1_001.fastq",2)],"libX"=data.zj25.lib9[,rep("ZJ25_51m__CD34pos__BMR_CD175_sampled_I506_S6_L001_R1_001.fastq",2)]),c(1,max(zj25.time)))+xlab("Time=51")
dev.off()

pdf("figures/f4a_m11_4_cd34L.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,rep("M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq",2)],"libX"=data.m11.lib11[,rep("M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq",2)]),c(1,max(m11.time)))+xlab("Time=4")
dev.off()

pdf("figures/f4a_m11_4_cd34R.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,rep("M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq",2)],"libX"=data.m11.lib11[,rep("M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq",2)]),c(1,max(m11.time)))+xlab("Time=4")
dev.off()

pdf("figures/f4a_m11_15_cd34L.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,rep("M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq",2)],"libX"=data.m11.lib11[,rep("M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq",2)]),c(1,max(m11.time)))+xlab("Time=15")
dev.off()

pdf("figures/f4a_m11_15_cd34R.pdf",height=20,width=20)
stacked_2(list("lib7"=data.m11.lib7[,rep("M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq",2)],"libX"=data.m11.lib11[,rep("M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq",2)]),c(1,max(m11.time)))+xlab("Time=15")
dev.off()

#############################################################################################
#Figure 4C#
#############################################################################################

pdf("figures/f4c_m11_cum.pdf",height=5,width=5)
unicum(list("lib7"=data.m11.lib7,"lib11"=data.m11.lib11),samples=c(m11.grans,m11.t,m11.b,m11.mono,m11.nk),
       type="Cumulative",time=m11.time)
dev.off()

pdf("figures/f4c_m11_un.pdf",height=5,width=5)
unicum(list("lib7"=data.m11.lib7,"lib11"=data.m11.lib11),samples=c(m11.grans,m11.t,m11.b,m11.mono,m11.nk),
       type="Unique",time=m11.time)
dev.off()

pdf("figures/f4c_zj25_cum.pdf",height=5,width=5)
unicum(list("lib7"=data.zj25.lib7,"lib9"=data.zj25.lib9),samples=c(zj25.grans,zj25.t,zj25.b,zj25.mono,zj25.nk),
       type="Cumulative",time=zj25.time)
dev.off()

pdf("figures/f4c_zj25_un.pdf",height=5,width=5)
unicum(list("lib7"=data.zj25.lib7,"lib9"=data.zj25.lib9),samples=c(zj25.grans,zj25.t,zj25.b,zj25.mono,zj25.nk),
       type="Unique",time=zj25.time)
dev.off()

library(stringr)

#############################################################################################
#Supplemental Table#
#############################################################################################
pdf("figures/supp_zj25_t.pdf",height=20,width=20)
stacked_3(list("lib7"=data.zj25.lib7[,zj25.t],"libX"=data.zj25.lib9[,zj25.t],"Other"=data.zj25.lib10[,zj25.t]
               ),zj25.time)
dev.off()

pdf("figures/supp_zj25_b.pdf",height=20,width=20)
stacked_3(list("lib7"=data.zj25.lib7[,zj25.b],"libX"=data.zj25.lib9[,zj25.b],"Other"=data.zj25.lib10[,zj25.b]
               ),zj25.time)
dev.off()

pdf("figures/supp_zj25_mono.pdf",height=20,width=20)
stacked_3(list("lib7"=data.zj25.lib7[,zj25.mono],"libX"=data.zj25.lib9[,zj25.mono],"Other"=data.zj25.lib10[,zj25.mono]
               ),zj25.time)
dev.off()

pdf("figures/supp_zj25_nk.pdf",height=20,width=20)
stacked_3(list("lib7"=data.zj25.lib7[,zj25.nk],"libX"=data.zj25.lib9[,zj25.nk],"Other"=data.zj25.lib10[,zj25.nk]
               ),zj25.time)
dev.off()

pdf("figures/supp_zj25_grans.pdf",height=20,width=20)
stacked_3(list("lib7"=data.zj25.lib7[,zj25.grans],"libX"=data.zj25.lib9[,zj25.grans],"Other"=data.zj25.lib10[,zj25.grans]
               ),zj25.time)
dev.off()

pdf("figures/supp_zj25_cum.pdf",height=5,width=5)
unicum(list("lib7"=data.zj25.lib7,"lib9"=data.zj25.lib9,"Other"=data.zj25.lib10),samples=c(zj25.grans,zj25.t,zj25.b,zj25.mono,zj25.nk),
       type="Cumulative",time=zj25.time)
dev.off()

pdf("figures/supp_zj25_un.pdf",height=5,width=5)
unicum(list("lib7"=data.zj25.lib7,"lib9"=data.zj25.lib9,"Other"=data.zj25.lib10),samples=c(zj25.grans,zj25.t,zj25.b,zj25.mono,zj25.nk),
       type="Unique",time=zj25.time)
dev.off()

zj25.combined_data_list=list(data.zj25.lib7,data.zj25.lib9,data.zj25.lib10)
zj25.combined_readme_list=list(readme.zj25.lib7,readme.zj25.lib9,readme.zj25.lib10)
pdf("figures/supp_zj25.pdf",height=20,width=20)
barcode_Xlib(combined_data_list=zj25.combined_data_list,combined_readme_list=zj25.combined_readme_list,
             samples=zj25.samples,samples.filenames=zj25.samples.filename)
dev.off()


barcodetrackR::cor_plot(data.zj25.lib7[,zj25.samples.filename],names=NULL)

x=NULL
for(i in 1:7){
  x=c(x,seq(i,35,by=7))
}

temp.zj25.samples=zj25.samples.filename[x]

names=NULL
names$time=rep(zj25.time,each=5)
names$cell=rep(c("T","B","Mono","NK","Grans"),times=length(zj25.time))
names$full=paste(names$time,"m ",names$cell,sep='')

pdf("figures/supp_zj25_pearson_lib7.pdf",height=5,width=5)
barcodetrackR::cor_plot(data.zj25.lib7[,temp.zj25.samples],names=names$full,method_corr = "pearson")
dev.off()

pdf("figures/supp_zj25_pearson_lib9.pdf",height=5,width=5)
barcodetrackR::cor_plot(data.zj25.lib9[,temp.zj25.samples[-28]],names=names$full[-28],method_corr = "pearson",your_title = "Note: had to delete 44m mono")
dev.off()


pdf("figures/supp_zj25_lib10.pdf",height=20,width=20)
barcodetrackR::barcode_ggheatmap(data.zj25.lib10[,zj25.samples.filename],names=zj25.samples,label_size=10,n_clones = 20,grid=FALSE,cellnote_size = 0)
dev.off()


x=NULL
for(i in 1:14){
  x=c(x,seq(i,70,by=14))
}

# temp.m11.samples=m11.samples.filename[x]
# temp.m11.samples=temp.m11.samples[-(31:60)];temp.m11.samples=temp.m11.samples[-(11:15)]
# 
# names=NULL
# names$time=rep(c(1:2,4:6,13,14),each=5)
# names$cell=rep(c("T","B","Mono","NK","Grans"),times=7)
# names$full=paste(names$time,"m ",names$cell,sep='')
# 
# pdf("figures/supp_m11_pearson_lib7wo3.pdf",height=5,width=5)
# barcodetrackR::cor_plot(data.m11.lib7[,temp.m11.samples],names=names$full,method_corr = "pearson")
# dev.off()
# 
# pdf("figures/supp_m11_pearson_lib11wo3.pdf",height=5,width=5)
# barcodetrackR::cor_plot(data.m11.lib11[,temp.m11.samples],names=names$full,method_corr = "pearson")
# dev.off()

temp.m11.samples=m11.samples.filename[x]
temp.m11.samples=temp.m11.samples[-(31:60)];temp.m11.samples=temp.m11.samples[-(26:30)]

names=NULL
names$time=rep(c(1:5,13,14),each=5)
names$cell=rep(c("T","B","Mono","NK","Grans"),times=7)
names$full=paste(names$time,"m ",names$cell,sep='')

pdf("figures/supp_m11_pearson_lib7wo6.pdf",height=5,width=5)
barcodetrackR::cor_plot(data.m11.lib7[,temp.m11.samples],names=names$full,method_corr = "pearson")
dev.off()

pdf("figures/supp_m11_pearson_lib11wo6.pdf",height=5,width=5)
barcodetrackR::cor_plot(data.m11.lib11[,temp.m11.samples],names=names$full,method_corr = "pearson")
dev.off()

#############################################################################################
#NGS Table#
#############################################################################################

readme.zj25=list(readme.zj25.lib7,readme.zj25.lib9)

i=1
zj25.unique.lib7=NULL
zj25.unique.lib9=NULL
for(sample in filename){
  temp=data.zj25.lib7[,sample]
  temp=sum(temp>0)
  zj25.unique.lib7[i]=temp
  
  temp=data.zj25.lib9[,sample]
  temp=sum(temp>0)
  zj25.unique.lib9[i]=temp
  
  i=i+1
}

table.zj25=data.frame(Seq_total_reads=readme.zj25[[1]][filename,]$READS,
                      Percent_mapped_to_Lib7_Nonexpanded=readme.zj25[[1]][filename,]$MAP..,
                      Unique_barcodes_lib7=zj25.unique.lib7,
                      Percent_mapped_to_Lib9_HUVEC=readme.zj25[[2]][filename,]$MAP..,
                      Unique_barcodes_lib9=zj25.unique.lib9,
                      Assigned_threshold=readme.zj25[[1]][filename,]$MAP..)

zj25.table.names=zj25.samples

zj25.table.names=str_replace(zj25.table.names, "ZJ25 ", "")
zj25.table.names=str_replace(zj25.table.names, "LT(..)", "")
zj25.table.names=str_replace(zj25.table.names, "SP(..)", "")
zj25.table.names=str_replace(zj25.table.names, "SS(.)", "")
zj25.table.names=str_replace(zj25.table.names, "CD(...)", "")
zj25.table.names=str_replace(zj25.table.names, "cell", "")
zj25.table.names=gsub("Gr ","Grans",zj25.table.names)
zj25.table.names=gsub(" ","",zj25.table.names)
zj25.table.names=gsub("m","m ",zj25.table.names)

rownames(table.zj25)=zj25.table.names

readme.m11=list(readme.m11.lib7,readme.m11.lib11)

i=1
m11.unique.lib7=NULL
m11.unique.lib11=NULL
for(sample in m11.samples.filename){
  temp=data.m11.lib7[,sample]
  temp=sum(temp>0)
  m11.unique.lib7[i]=temp
  
  temp=data.m11.lib11[,sample]
  temp=sum(temp>0)
  m11.unique.lib11[i]=temp
  
  i=i+1
}

table.m11=data.frame(Seq_total_reads=readme.m11[[1]][m11.samples.filename,]$READS,
                     Percent_mapped_to_Lib7_Nonexpanded=readme.m11[[1]][m11.samples.filename,]$MAP..,
                     Unique_barcodes_lib7=m11.unique.lib7,
                     Percent_mapped_to_Lib9_HUVEC=readme.m11[[2]][m11.samples.filename,]$MAP..,
                     Unique_barcodes_lib11=m11.unique.lib11,
                     Assigned_threshold=readme.m11[[1]][m11.samples.filename,]$MAP..)

m11.table.names=m11.samples

m11.table.names=str_replace(m11.table.names, "M11021072 ", "")
m11.table.names=str_replace(m11.table.names, "(.)LT(..)", "")
m11.table.names=str_replace(m11.table.names, "(.)SP(..)", "")
m11.table.names=str_replace(m11.table.names, "(.)SS(.)", "")
m11.table.names=str_replace(m11.table.names, "(.)CD(...)", "")
m11.table.names=str_replace(m11.table.names, "cell", "")
m11.table.names=str_replace(m11.table.names, "CELL", "")
m11.table.names=str_replace(m11.table.names, "(.)IY(..)", "")
m11.table.names=gsub("Gr ","Grans",m11.table.names)
m11.table.names=gsub("hm",".5m",m11.table.names)
m11.table.names=gsub(" ","",m11.table.names)
m11.table.names=gsub("m","m ",m11.table.names)


rownames(table.m11)=m11.table.names

write.table(table.zj25,file="tables/NGS_zj25_by_sample.txt",col.names=NA,quote=FALSE,sep='\t')
write.table(table.m11,file="tables/NGS_m11_by_sample.txt",col.names=NA,quote=FALSE,sep='\t')

#############################################################################################
#Dynamics Table#
#############################################################################################
table2.zj25=as.data.frame(matrix(ncol=2,nrow=length(zj25.samples.filename)))
colnames(table2.zj25)=c("biased_lib7","biased_lib9")
rownames(table2.zj25)=zj25.table.names

table2.m11=as.data.frame(matrix(ncol=2,nrow=length(m11.samples.filename)))
colnames(table2.m11)=c("biased_lib7","biased_lib11")
rownames(table2.m11)=m11.table.names

create.table.bias=function(table,data,filename){
  for(library in names(data)){
    for(timepoint in 1:(length(filename)/5)){
      index=(0*length(filename)/5+(timepoint))
      temp.data=data[[library]]
      
      temp.one=temp.data[,index]
      temp.others=temp.data[,c((1*length(filename)/5+(timepoint)),(2*length(filename)/5+(timepoint)),
                               (3*length(filename)/5+(timepoint)),(4*length(filename)/5+(timepoint)))]
      temp.barcodes=rownames(temp.data)[temp.one>10*apply(temp.others,1,max)]
      
      temp.number=length(temp.barcodes)
      temp.contribution=100*sum(temp.data[temp.barcodes,index])
      
      table[index,paste("biased_",library,sep="")]=temp.contribution
      
      
      index=(1*length(filename)/5+(timepoint))
      temp.one=temp.data[,index]
      temp.others=temp.data[,c((0*length(filename)/5+(timepoint)),(2*length(filename)/5+(timepoint)),
                               (3*length(filename)/5+(timepoint)),(4*length(filename)/5+(timepoint)))]
      temp.barcodes=rownames(temp.data)[temp.one>10*apply(temp.others,1,max)]
      
      temp.number=length(temp.barcodes)
      temp.contribution=100*sum(temp.data[temp.barcodes,index])
      
      table[index,paste("biased_",library,sep="")]=temp.contribution
      
      index=(2*length(filename)/5+(timepoint))
      temp.one=temp.data[,index]
      temp.others=temp.data[,c((1*length(filename)/5+(timepoint)),(0*length(filename)/5+(timepoint)),
                               (3*length(filename)/5+(timepoint)),(4*length(filename)/5+(timepoint)))]
      temp.barcodes=rownames(temp.data)[temp.one>10*apply(temp.others,1,max)]
      
      temp.number=length(temp.barcodes)
      temp.contribution=100*sum(temp.data[temp.barcodes,index])
      
      table[index,paste("biased_",library,sep="")]=temp.contribution
      
      index=(3*length(filename)/5+(timepoint))
      temp.one=temp.data[,index]
      temp.others=temp.data[,c((1*length(filename)/5+(timepoint)),(2*length(filename)/5+(timepoint)),
                               (0*length(filename)/5+(timepoint)),(4*length(filename)/5+(timepoint)))]
      temp.barcodes=rownames(temp.data)[temp.one>10*apply(temp.others,1,max)]
      
      temp.number=length(temp.barcodes)
      temp.contribution=100*sum(temp.data[temp.barcodes,index])
      
      table[index,paste("biased_",library,sep="")]=temp.contribution
      
      index=(4*length(filename)/5+(timepoint))
      temp.one=temp.data[,index]
      temp.others=temp.data[,c((1*length(filename)/5+(timepoint)),(2*length(filename)/5+(timepoint)),
                               (3*length(filename)/5+(timepoint)),(0*length(filename)/5+(timepoint)))]
      temp.barcodes=rownames(temp.data)[temp.one>10*apply(temp.others,1,max)]
      
      temp.number=length(temp.barcodes)
      temp.contribution=100*sum(temp.data[temp.barcodes,index])
      
      table[index,paste("biased_",library,sep="")]=temp.contribution
    }
  }
  
  return(table)
}

table2.zj25=create.table.bias(table2.zj25,list("lib7"=barcodetrackR::barcode_ggheatmap(data.zj25.lib7[,zj25.samples.filename],n_clones=20,printtable=TRUE),
                                               "lib9"=barcodetrackR::barcode_ggheatmap(data.zj25.lib9[,zj25.samples.filename],n_clones=20,printtable=TRUE)),zj25.samples.filename)
table2.m11=create.table.bias(table2.m11,list("lib7"=barcodetrackR::barcode_ggheatmap(data.m11.lib7[,m11.samples.filename],n_clones=20,printtable=TRUE),
                                             "lib11"=barcodetrackR::barcode_ggheatmap(data.m11.lib11[,m11.samples.filename],n_clones=20,printtable=TRUE)),m11.samples.filename)

table3.zj25=as.data.frame(matrix(nrow=5,ncol=6));rownames(table3.zj25)=c("T","B","Mono","NK","Grans")
colnames(table3.zj25)=c("ave_lib7","ave_lib9","exp_lib7","exp_lib9","biased_lib7","biased_lib9")

table3.m11=as.data.frame(matrix(nrow=5,ncol=6));rownames(table3.m11)=c("T","B","Mono","NK","Grans")
colnames(table3.m11)=c("ave_lib7","ave_lib11","exp_lib7","exp_lib11","biased_lib7","biased_lib11")

mean0=function(vector){
  if(sum(vector)==0 & !is.na(sum(vector))){
    return(0)
  }else if(sum(is.na(vector))==length(vector)){
    return(NA)
  }else{
    vector=vector[vector>0 &!is.na(vector) & !is.nan(vector) & is.finite(vector)]
    return(mean(vector))
  }
}

create.table.bias.small=function(table.small,table.large,data,filename){
  for(library in names(data)){
    temp.data=data[[library]]
    
    i=0
    for(celltype in rownames(table.small)){
      
      temp.data.cell=temp.data[,((length(filename)/5)*i+1):((length(filename)/5)*(i+1))]
      
      temp.ave=apply(temp.data.cell,1,mean0)
      
      table.small[celltype,paste("biased_",library,sep='')]=
        paste(round(mean(table.large[((length(filename)/5)*i+1):((length(filename)/5)*(i+1)),paste("biased_",library,sep="")]),digits=2),"%",sep="")
      
      table.small[celltype,paste("ave_",library,sep='')]=paste(round(10000*mean(temp.ave),digits=2),"x10^4",sep="")
      
      temp.exp=matrix(NA,ncol=(ncol(temp.data.cell)-1),nrow=nrow(temp.data.cell))
      for(timepoint in 1:(ncol(temp.data.cell)-1)){
        
        temp.exp[(temp.data.cell[,(timepoint+1)]>0 & temp.data.cell[,(timepoint)]>0),timepoint]=(temp.data.cell[,(timepoint+1)]/temp.data.cell[,timepoint])[(temp.data.cell[,(timepoint+1)]>0 & temp.data.cell[,(timepoint)]>0)]
      }
      
      temp.exp=apply(temp.exp,1,mean0)
      table.small[celltype,paste("exp_",library,sep='')]=mean0(temp.exp)
      
      i=i+1
    }
  }
  
  return(table.small)
}

table3.zj25=create.table.bias.small(table3.zj25,table2.zj25,list("lib7"=barcodetrackR::barcode_ggheatmap(data.zj25.lib7[,zj25.samples.filename],n_clones=20,printtable=TRUE),
                                                                 "lib9"=barcodetrackR::barcode_ggheatmap(data.zj25.lib9[,zj25.samples.filename],n_clones=20,printtable=TRUE)),zj25.samples.filename)
table3.m11=create.table.bias.small(table3.m11,table2.m11,list("lib7"=barcodetrackR::barcode_ggheatmap(data.m11.lib7[,m11.samples.filename],n_clones=20,printtable=TRUE),
                                                              "lib11"=barcodetrackR::barcode_ggheatmap(data.m11.lib11[,m11.samples.filename],n_clones=20,printtable=TRUE)),m11.samples.filename)

write.table(table3.zj25,file="tables/dynamics_zj25_by_cell.txt",col.names=NA,quote=FALSE,sep='\t')
write.table(table3.m11,file="tables/dynamics_m11_by_cell.txt",col.names=NA,quote=FALSE,sep='\t')


#############################################################################################
#BM bias Table#
#############################################################################################
bm.bias=function(monkey,data,bm.samples,other.samples){
  
  if(monkey=="ZJ25"){
    final.table=data.frame(matrix(NA,ncol=6,nrow=4))
    colnames(final.table)=c("N_periphery_biased","N_neither_biased","N_BM_biased","Cont_BMB_in_grans","Cont_BMB_in_LBM","Cont_BMB_in_RBM")
    rownames(final.table)=c("4m_lib7","4m_lib9","51m_lib7","51m_lib9")
  }else if(monkey=="M11"){
    final.table=data.frame(matrix(NA,ncol=6,nrow=4))
    colnames(final.table)=c("N_periphery_biased","N_neither_biased","N_BM_biased","Cont_BMB_in_grans","Cont_BMB_in_LBM","Cont_BMB_in_RBM")
    rownames(final.table)=c("4m_lib7","4m_lib11","15m_lib7","15m_lib11")
  }
  library.index=1
  for(library in names(data)){
    temp.data=data[[library]]
    
    temp.data=barcodetrackR::barcode_ggheatmap(temp.data[,c(bm.samples,other.samples)],n_clones=20,printtable=TRUE)
    
    temp.data=apply(temp.data,2,function(x){x/sum(x)})
    i.index=1
    for(i in 0:(length(bm.samples)/2-1)){
      
      bm.data=temp.data[,bm.samples[c(i*2+1,i*2+2)]]
      other.data=temp.data[,other.samples[(i*length(other.samples)/5+1):(i*length(other.samples)/5+5)]]
      
      temp.barcode=vector(length=nrow(temp.data))
      
      temp.comparison=cbind(bm.data,apply(other.data,1,max))
      temp.barcode[apply(temp.comparison,1,function(x){sum(x[1:(length(x)-1)]>10*x[length(x)])>0})]="BM"
      final.table$N_BM_biased[seq(i.index,4,by=2)[library.index]]=sum(temp.barcode=="BM")
      bm.bar=rownames(bm.data)[temp.barcode=="BM"]
      final.table$Cont_BMB_in_LBM[seq(i.index,4,by=2)[library.index]]=paste(round(100*apply(bm.data[bm.bar,],2,sum)[1],digits=2),"%",sep="")
      final.table$Cont_BMB_in_RBM[seq(i.index,4,by=2)[library.index]]=paste(round(100*apply(bm.data[bm.bar,],2,sum)[2],digits=2),"%",sep="")
      final.table$Cont_BMB_in_grans[seq(i.index,4,by=2)[library.index]]=paste(round(100*sum(other.data[bm.bar,5]),digits=2),"%",sep="")
      
      temp.comparison=cbind(other.data,apply(bm.data,1,max))
      temp.barcode[apply(temp.comparison,1,function(x){sum(x[1:(length(x)-1)]>10*x[length(x)])>0})]="PERIPHERAL"
      final.table$N_periphery_biased[seq(i.index,4,by=2)[library.index]]=sum(temp.barcode=="PERIPHERAL")
      
      
      temp.barcode[temp.barcode != "BM" & temp.barcode != "PERIPHERAL"]="NEITHER"
      final.table$N_neither_biased[seq(i.index,4,by=2)[library.index]]=sum(temp.barcode=="NEITHER")
      
      
      i.index=i.index+1
    }
    
    library.index=library.index+1
  }
  return(final.table)
  
}
m11.bm=c("M11021072_4m__CD34p_BML_LT19_sampled_LT19_i521_S41_L005_R1_001.fastq","M11021072_4m__CD34p_BMR_LT19_sampled_LT19_i520_S40_L005_R1_001.fastq","M11021072_15m__CD34pos__BML_CD175_sampled_I523_S19_L001_R1_001.fastq","M11021072_15m__CD34pos__BMR_CD175_sampled_I522_S18_L001_R1_001.fastq")
temp.by.time=NA
for(i in 1:(length(m11.samples.filename)/5)){temp.by.time=c(temp.by.time,c(seq(i,length(m11.samples.filename),by=length(m11.samples.filename)/5)))};temp.by.time=temp.by.time[-1]
m11.other=m11.samples.filename[temp.by.time]
m11.other=m11.other[c((1+(5*(4-1))):(5+(5*(4-1))),(1+(5*(14-1))):(5+(5*(14-1))))]
m11.bm.bias=bm.bias("M11",list("lib7"=data.m11.lib7,"lib11"=data.m11.lib11),m11.bm,m11.other)

zj25.bm=c("SP9_R28_ZJ25_4_5m_CD34pos_leftBM.fastq","SP9_R29_ZJ25_4_5m_CD34pos_rightBM.fastq","ZJ25_51m__CD34pos__BML_CD175_sampled_I507_S7_L001_R1_001.fastq","ZJ25_51m__CD34pos__BMR_CD175_sampled_I506_S6_L001_R1_001.fastq")
temp.by.time=NA
for(i in 1:(length(zj25.samples.filename)/5)){temp.by.time=c(temp.by.time,c(seq(i,length(zj25.samples.filename),by=length(zj25.samples.filename)/5)))};temp.by.time=temp.by.time[-1]
zj25.other=zj25.samples.filename[temp.by.time]
zj25.other=zj25.other[c((1+(5*(3-1))):(5+(5*(3-1))),(1+(5*(7-1))):(5+(5*(7-1))))]
zj25.bm.bias=bm.bias("ZJ25",list("lib7"=data.zj25.lib7,"lib9"=data.zj25.lib9),zj25.bm,zj25.other)

write.table(zj25.bm.bias,file="tables/BMbias_zj25.txt",col.names=NA,quote=FALSE,sep='\t')
write.table(m11.bm.bias,file="tables/BMbias_m11.txt",col.names=NA,quote=FALSE,sep='\t')



