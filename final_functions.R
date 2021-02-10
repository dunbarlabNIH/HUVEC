#############################################################################################
#Utilities#
#############################################################################################

read.combined=function(filename){
  return(read.delim(filename,header=TRUE,sep='\t',row.names=1))
}

filename=function(samples,key){
  all.samples.filename=as.vector(unlist(key$FILENAME))
  all.samples.givenname=as.vector(unlist(key$GIVENNAME))
  
  key.samples=NULL
  
  for(i in 1:length(samples)){
    for(j in 1:length(all.samples.filename)){
      if(samples[i]==all.samples.givenname[j]){
        key.samples=c(key.samples,all.samples.filename[j])
      }
    }
  }
  return(key.samples)
}

custom_log = function(x, log_choice, vary_log_min){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  if(vary_log_min) {
    x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  } else {
    x[x == -Inf] <- (log(100/4000000, log_choice)-1)
  }
  return(x)
}

require(magrittr)
require(dplyr)

#############################################################################################
#Plots#
#############################################################################################


barcode_Xlib=function(combined_data_list,combined_readme_list,samples,samples.filenames,colors=c("red","blue","chartreuse4"),
                      n_clones=10,grid=FALSE,log_choice=exp(1)){
  
  library(ComplexHeatmap)
  library(circlize)
  
  for(check in 1:length(combined_data_list)){
    temp_data=combined_data_list[[check]]
    sample_barcode=rownames(temp_data)[1]
    libid_data=substring(sample_barcode,1,6)
    
    temp_readme=combined_readme_list[[check]]
    libid_readme=temp_readme$LIBID[1]
    
    if(libid_data==libid_readme){
      print(paste("List element",check,"matches"))
    }else{
      stop("The input order of combined files and readme files are not consistent")
    }
  }
  
  master_combined=cbind(combined_data_list[[1]][,samples.filenames],"color"=LETTERS[1])
  master_readme=data.frame(combined_readme_list[[1]]$'MAP..',row.names = rownames(combined_readme_list[[1]]))
  
  
  for(combined_library in 2:length(combined_data_list)){
    master_combined=rbind(master_combined,
                          cbind(combined_data_list[[combined_library]][,samples.filenames],
                                "color"=LETTERS[combined_library]))
    master_readme=cbind(master_readme,data.frame(combined_readme_list[[combined_library]]$'MAP..',row.names=rownames(combined_readme_list[[combined_library]])))
  }
  
  temp=rowSums(master_readme)
  mapping=master_readme/temp*100
  

  barcode_colors=data.frame("color"=master_combined$color,row.names = rownames(master_combined))

  heatmap_combined=barcodetrackR::barcode_ggheatmap(master_combined[,samples.filenames],n_clones=n_clones,printtable=TRUE)
  top_clones=rownames(heatmap_combined)
  colnames(heatmap_combined)=samples
  
  heatmap_combined=custom_log(heatmap_combined, log_choice, FALSE)
  actual_scale = c(log(100/4000000, log_choice) - 1, log(100/4000000, log_choice), log(0.001, log_choice), log(0.01,log_choice), log(0.1, log_choice), 0)
  your_scale = scales::rescale(actual_scale, to = c(0,1))
  your_labels = c("Defined 0", "0.0025% ", "0.1%", "1%", "10%", "100%")
  your_limits = c(log(100/4000000, log_choice)-1,0)
  
  heatmap_colors=data.frame(color=barcode_colors[top_clones,],row.names=top_clones)

  if(length(colors)==2){
    barcode_source = HeatmapAnnotation( df = heatmap_colors, col = list(color = c("A"=colors[1],"B"=colors[2])), which = 'row')
  }
  if(length(colors)==3){
    barcode_source = HeatmapAnnotation( df = heatmap_colors, col = list(color = c("A"=colors[1],"B"=colors[2],"C"=colors[3])), which = 'row')
  }

  
  mapping=mapping[samples.filenames,]
  rownames(mapping)=samples
  
  column_ha = HeatmapAnnotation(mapping = anno_barplot(as.matrix(mapping), 
                                                       axis = TRUE,
                                                       axis_param=(axis_side = "right"),
                                                       gp = gpar( fill = colors[1:length(combined_data_list)] , 
                                                                  width = unit(2, "cm")),
                                                       show_legend = TRUE),show_annotation_name = TRUE,
                                annotation_name_offset = unit(1, "cm"),which = "column")
  
  if(!grid){grid_choice=NA}else{grid_choice="black"}

  plot=Heatmap(as.matrix(heatmap_combined),show_row_dend = FALSE,show_column_dend = FALSE,
          show_row_names = FALSE,row_order = top_clones, cluster_columns = FALSE,
          col = colorRamp2(actual_scale, c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")),
          rect_gp = gpar(col = grid_choice), top_annotation = column_ha
  )+ barcode_source

  return(plot)
  
}

stacked=function(list_data,time_points,colors,clone_selection_option=c("last","all"),n_top_clones=5,create_legends=FALSE){
  nlib=length(list_data)
  
  top_barcodes=NULL;list_top=NULL;list_palette=NULL;total_lib=NULL;legend=NULL
  total_lib=matrix(data=0,nrow=1,ncol=length(time_points))
  if(clone_selection_option=="last"){
    for(unilib in names(list_data)){
      unilib_data=list_data[[unilib]]
      
      ntime=ncol(unilib_data)
      temp_data=unilib_data[,ntime]
      
      top_barcodes[[unilib]]=rownames(unilib_data)[rank(-temp_data)<=n_top_clones]
      
      top_data=unilib_data[top_barcodes[[unilib]],]
      other_data=unilib_data[!(rownames(unilib_data) %in% top_barcodes[[unilib]]),]
      other_data_comb=colSums(other_data)
      
      list_top[[unilib]]=rbind(top_data,other_data_comb)
      rownames(list_top[[unilib]])=c(paste(unilib,top_barcodes[[unilib]]),paste(unilib,"other"))
      
      n_barcodes=nrow(list_top[[unilib]])-1
      
      if(RColorBrewer::brewer.pal.info[colors[[unilib]],]$maxcolors<n_barcodes){
        temp_n=RColorBrewer::brewer.pal.info[colors[[unilib]],]$maxcolors
      }else{
        temp_n=n_barcodes
      }
      
      getPalette = colorRampPalette(RColorBrewer::brewer.pal(temp_n,colors[[unilib]]))
      temp_colors=sapply(n_barcodes,getPalette)
      list_palette[[unilib]]=data.frame(c(sample(temp_colors,n_barcodes),"NA"),row.names=c(paste(unilib,top_barcodes[[unilib]]),paste(unilib,"other")))
      
      if(create_legends){
        legend[[unilib]]=barplot(table(rep(1,times=100)+runif(100),rep(1,times=100)),col=as.vector(sapply(100,getPalette)),border=NA)
      }
      
      total_lib=rbind(total_lib,total_lib[nrow(total_lib),])
      total_lib=rbind(total_lib,colSums(unilib_data)+total_lib[nrow(total_lib),])
    }
  }
  if(clone_selection_option=="all"){
    for(unilib in names(list_data)){
      unilib_data=list_data[[unilib]]
      
      ntime=ncol(unilib_data)
      
      top_barcodes[[unilib]]=rownames(barcodetrackR::barcode_ggheatmap(unilib_data,n_clones=n_top_clones,printtable=TRUE))
      
      top_data=unilib_data[top_barcodes[[unilib]],]
      other_data=unilib_data[!(rownames(unilib_data) %in% top_barcodes[[unilib]]),]
      other_data_comb=colSums(other_data)
      
      list_top[[unilib]]=rbind(top_data,other_data_comb)
      rownames(list_top[[unilib]])=c(paste(unilib,top_barcodes[[unilib]]),paste(unilib,"other"))
      
      n_barcodes=nrow(list_top[[unilib]])-1
      
      
      if(RColorBrewer::brewer.pal.info[colors[[unilib]],]$maxcolors<n_barcodes){
        temp_n=RColorBrewer::brewer.pal.info[colors[[unilib]],]$maxcolors
      }else{
        temp_n=n_barcodes
      }
      
      getPalette = colorRampPalette(RColorBrewer::brewer.pal(temp_n,colors[[unilib]]))
      temp_colors=sapply(n_barcodes,getPalette)
      list_palette[[unilib]]=data.frame(c(sample(temp_colors,n_barcodes),"NA"),row.names=c(paste(unilib,top_barcodes[[unilib]]),paste(unilib,"other")))
      
      if(create_legends){
        legend[[unilib]]=barplot(table(rep(1,times=100)+runif(100),rep(1,times=100)),col=as.vector(sapply(100,getPalette)),border=NA)
      }
      
      total_lib=rbind(total_lib,total_lib[nrow(total_lib),])
      total_lib=rbind(total_lib,colSums(unilib_data)+total_lib[nrow(total_lib),])
    }
  }
  
  total_lib=total_lib[-1,]
  rownames(total_lib)=paste(rep(names(list_data),each=2),rep(c("start","end"),times=nlib))
  
  plotting_data=list_top[[names(list_top)[1]]]
  palette_data=list_palette[[names(list_top)[1]]]
  for(unilib in names(list_top)[2:nlib]){
    plotting_data=rbind(plotting_data,list_top[[unilib]])
    palette_data=rbind(palette_data,list_palette[[unilib]])
  }
  colnames(palette_data)="values"
  plotting_data=apply(plotting_data,2,function(x){x/sum(x)})
  total_lib=apply(total_lib,2,function(x){x/max(x)})
  
  colnames(plotting_data)=time_points;colnames(total_lib)=time_points
  
  temp_list=lapply(1:ncol(plotting_data),function(col){data.frame(prop=plotting_data[,col],month=colnames(plotting_data)[col],barcode=row.names(plotting_data))})
  full=NULL;for(row in 1:length(temp_list)){full=rbind(full,temp_list[[row]])};rownames(full)=NULL
  
  temp_list=lapply(1:ncol(total_lib),function(col){data.frame(prop=total_lib[,col],month=colnames(total_lib)[col],barcode=row.names(total_lib))})
  full_lib=NULL;for(row in 1:length(temp_list)){full_lib=rbind(full_lib,temp_list[[row]])};rownames(full_lib)=NULL
  full_lib$lib=rep(names(list_data),each=2)
  full_lib$color=paste("dark",rep(rep(gsub('.{1}$', '',tolower(as.vector(unlist(colors)))),each=2),times=length(time_points)));
  full_lib$level=c("start","end")
  full_lib=full_lib%>%dplyr::arrange(level);full_lib[full_lib$level=="start",]=full_lib[full_lib$level=="start",]%>%dplyr::arrange(desc(month))
  
  if(all.equal(rownames(palette_data),as.vector(full$barcode[1:(nrow(full)/ntime)]))){
    library(ggplot2)
    
    full$barcode_color=palette_data$values
    full$barcode=paste(seq(from=0,to=1,length.out=nrow(plotting_data)),full$barcode)
    
    temp_lib=NULL;i=1
    for(unilib in names(list_data)){
      temp_lib=c(temp_lib,rep(unilib,1+length(top_barcodes[[unilib]])))
      legend[[paste(unilib,"total")]]=barplot(table(100,1),col=unique(full_lib$color)[i],border=NA)
      i=i+1
    }
    
    full$lib=temp_lib
    
    y=function(m,x,b){return(m*x+b)};x=function(m,y,b){return((y-b)/m)}
    xlim=c(min(time_points)+.04,max(time_points)-.04);ylim=c(0.5,99.7);nline=40;m=(ylim[2]-ylim[1])/(xlim[2]-xlim[1])
    
    pattern_df=data.frame(x0=xlim[1],y0=ylim[1],xend=xlim[2],yend=ylim[2])
    pattern_df=rbind(pattern_df,data.frame(x0=pattern_df$x0[1]+(xlim[2]-xlim[1])/nline,y0=ylim[1],xend=xlim[2],yend=pattern_df$yend[1]-(ylim[2]-ylim[1])/nline))
    pattern_df=rbind(pattern_df,data.frame(x0=xlim[1],y0=pattern_df$y0[1]+(ylim[2]-ylim[1])/nline,xend=pattern_df$xend[1]-(xlim[2]-xlim[1])/nline,yend=ylim[2]))
    for(i in seq(4,2*nline,by=2)){
      pattern_df=rbind(pattern_df,data.frame(x0=pattern_df$x0[i-2]+(xlim[2]-xlim[1])/nline,y0=ylim[1],xend=xlim[2],yend=pattern_df$yend[i-2]-(ylim[2]-ylim[1])/nline))
      pattern_df=rbind(pattern_df,data.frame(x0=xlim[1],y0=pattern_df$y0[i-1]+(ylim[2]-ylim[1])/nline,xend=pattern_df$xend[i-1]-(xlim[2]-xlim[1])/nline,yend=ylim[2]))
    }
    
    print_plot=ggplot(full) +
      
      geom_polygon(data=full_lib,aes(x=as.numeric(as.vector(month)),y=100*prop,group=lib),fill=full_lib$color)+
      geom_segment(data = pattern_df, aes(x = x0, y = y0, xend = xend, yend = yend),size=2,lineend = "round",color="grey48")+
      
      geom_area(data=full,aes(x=as.numeric(as.vector(month)), y=prop*100, group=barcode, fill=barcode),position=position_stack(reverse=T))+
      scale_fill_manual(values = as.vector(full$barcode_color))+
      theme(legend.position="none")+
      scale_x_continuous(limits=c(min(as.numeric(as.vector(full$month))),max(as.numeric(as.vector(full$month)))))+
      xlab("Time (months)")+
      geom_line(data=full,aes(x=as.numeric(as.vector(month)), y=prop*100, group=barcode),colour="black", size=.2, alpha=.4,position = position_stack(reverse = T))+
      
      geom_line(data=full_lib[(full_lib$prop != 0 & full_lib$prop != 1 & !duplicated(paste(full_lib$prop,full_lib$month))) ,],aes(x=as.numeric(as.vector(month)),y=100*prop,group=lib),size=2,color="black")+
      
      ylab("Percent Contribution \n of Barcodes")+ylim(0,100)
    
    if(create_legends){
      return(legend)
    }else{
      return(print_plot)
    }
  }else{print("Internal error")}
}
stacked_2=function(list_data,time_points){
  return(stacked(list_data,time_points,colors=list("lib7"="Reds","libX"="Blues"),clone_selection_option = "last",n_top_clones = 20,create_legends=FALSE))
}

stacked_3=function(list_data,time_points){
  return(stacked(list_data,time_points,colors=list("lib7"="Reds","libX"="Blues","Other"="Greens"),clone_selection_option = "last",n_top_clones = 20,create_legends=FALSE))
}

unicum=function(data,samples,type,time){
  nlib=length(data)
  celltypes=length(samples)/length(time)
  
  custom.rbind=function(a,b){
    shared=rownames(a)[rownames(a) %in% rownames(b)]
    comb.shared=as.data.frame(a[shared,]+b[shared,])
    comb.not=as.data.frame(rbind(a[!(rownames(a) %in% shared),],b[!(rownames(b) %in% shared),]))
    
    return(rbind(comb.shared,comb.not))
  }
  
  plot.data=data.frame(matrix(NA,ncol=length(time),nrow=length(names(data))));colnames(plot.data)=time;rownames(plot.data)=names(data)
  i=1
  for(lib in names(data)){
    temp.data.list=NULL
    temp.data=data[[lib]][,samples]
    for(cell in 1:celltypes){
      temp.data.list[[cell]]=temp.data[,((cell-1)*length(time)+1):((cell)*length(time))]
      colnames(temp.data.list[[cell]])=time
    }
    temp.comb.lib=temp.data.list[[1]]
    for(cell in 2:celltypes){
      temp.comb.lib=custom.rbind(temp.comb.lib,temp.data.list[[cell]])
    }
    
    if(type=="Cumulative"){
      plot.data[lib,]=as.vector(barcodetrackR::barcodecount(temp.comb.lib,months=time,count="cumulative"))
    }else if(type=="Unique"){
      plot.data[lib,]=as.vector(barcodetrackR::barcodecount(temp.comb.lib,months=time,count="unique"))
    }
    
    if(lib=="lib7"){
      color="red"
    }else if(lib=="Other"){
      color="chartreuse4"
    }else{
      color="blue"
    }
    if(i==1){
      plot(y=plot.data[1,],x=time,ylim=c(0,10000),type="l",col=color,xlab="Month",ylab=paste(type,"Barcodes"),lwd=5,cex=5)
    }else{
      points(y=plot.data[i,],x=time,type="l",col=color,lwd=5)
    }
    i=i+1
  }
  
}
