#Install necessary packages
package.list<-c("magick",
                "png",
                "this.path",
                "stringr",
                "ggplot2",
                "dplyr",
                "utils",
                "tcltk",
                "RColorBrewer",
                "e1071")
for(i in 1:length(package.list)){
  tryCatch(find.package(package.list[i]),
           error = function(e) install.packages(package.list[i],repos="http://lib.stat.cmu.edu/R/CRAN/"))
}

library(magick)
library(png)
library(this.path)
library(stringr)
library(ggplot2)
library(dplyr)
library(utils)
library(tcltk)
library(RColorBrewer)
library(e1071)

#Set seed
set.seed(123)

#Set image to be X by X px
img_size<-100

dir<-this.dir()

#Load parms
parms<-read.table(paste(dir,"/Parms.txt",sep=""),sep=":")

band_height<-as.numeric(str_trim(parms[which(parms$V1=="HeightBwBands"),]$V2))
band_width<-as.numeric(str_trim(parms[which(parms$V1=="WidthBwBands"),]$V2))
#save_obj<-as.character(str_trim(parms[which(parms$V1=="Objects"),]$V2))
#  keep_obj<-str_split(save_obj,pattern=",")[[1]]
alpha<-as.numeric(str_trim(parms[which(parms$V1=="Alpha"),]$V2))
use_mean<-as.logical(str_trim(parms[which(parms$V1=="FindWithMean"),]$V2))
run_id<-as.character(str_trim(parms[which(parms$V1=="RunID"),]$V2))
sample_size<-as.numeric(str_trim(parms[which(parms$V1=="Sample"),]$V2))
skew_status<-as.logical(str_trim(parms[which(parms$V1=="Skew"),]$V2))

save_obj<-tclVar("ALL")

#Create save location
run_loc<-paste(dir,"/RUNS/",as.Date(Sys.time()),"_",run_id,sep="")
dir.create(run_loc)

#Sets save location
plots_dir<-run_loc
debug_dir<-run_loc

#Get object metric data  
all_obj_metrics<-data.frame()

#-----------START IMAGE CYCLE----------------

#Select files
raw_files<-as.character(choose.files(multi = TRUE))

#Curate file names
all_files<-str_sub(raw_files,start=str_locate_all(raw_files,str_sub(as.Date(Sys.time()),0,4))[[1]][,1][[1]])
file<-all_files[str_detect(all_files,".tsv",negate=TRUE)]

#Get .tsv files
data_logs<-all_files[str_detect(all_files,".tsv",negate=FALSE)]

#Cycle counter increases after each run
cycle<-1
for(sel_file in file){  
  
  #------------------------------------------------
  #-----------------PROCESS IMAGES-------------------
  #------------------------------------------------
  
  x<-image_read(raw_files[cycle])
  
  #Convert image to png and resize to 100x100 px, ignore aspect ratio
  img_conv<-image_convert(x,format="png",colorspace="gray")
  img_res<-image_scale(img_conv,paste(img_size,"x",img_size,"!",sep=""))
  
  #Ensure image is 8 bits
  img_quant<-image_quantize(img_res,max=256)
  #img_quant<-image_modulate(img_quant,brightness=200)
  
  #Save as png and reload with png package
  setwd(paste(dir,"/ProcBands/",sep=""))
  image_write(img_quant,paste("PROCESSED_",sel_file,"_TEST.png",sep=""))
  img_png<-readPNG(paste(dir,"/ProcBands/PROCESSED_",sel_file,"_TEST.png",sep=""))
  
  #Number rows and col 1:100
  rownames(img_png)<-1:img_size
  colnames(img_png)<-1:img_size
  
  #Get vector of px brightness values with corresponding pixel position
  colvec<-c()
  for(i in 1:img_size){
    colvec<-c(colvec,as.vector(img_png[i,]))
  }
  
  #Create index for each pixel position
  colvec_data<-data.frame(Index=1:length(colvec),Signal=colvec)
  
  #Identify row for each index
  rowvec<-c()
  for(i in 1:img_size){
    tmpvec<-rep(i,times=img_size)
    rowvec<-c(rowvec,tmpvec)
  }
  colvec_data$Row=rowvec
  
  #Identify column for each index
  coln<-c(1:(img_size/10))
  colvec_data$Col=rep(1:img_size,times=img_size)
  
  #Only get objects for cycle 1
  if(cycle==1){
  
    #------------------------------------------------
    #------GENERATE DISTRIBUTION---------------------
    #------------------------------------------------
    
    #By central limit theorem, sampling distribution of sample means are ~normally distributed
    means<-c()
    for(i in 1:(5000)){
      x<-sample(colvec_data$Index,size=5000)
      sampled<-colvec_data[x,]
      sample_mean<-mean(sampled$Signal)
      means<-c(means,sample_mean)
    }
    
    #Get cumulative probability
    mean_frame<-as.data.frame(table(means))
    mean_frame$means<-as.numeric(as.character(mean_frame$means))
    mean_frame$Prob<-mean_frame$Freq/sum(mean_frame$Freq)
    mean_frame$CumProb<-NA
    
    for(i in 1:length(mean_frame$means)){
      x<-subset(mean_frame,means<=mean_frame[i,]$means)
      mean_frame[i,]$CumProb<-sum(x$Prob)
    }
    
    #Get upper lim, alpha=0.05
    upper_mean<-mean_frame[max(which(mean_frame$CumProb<=(1-alpha))),]$means
    real_alpha<-1-mean_frame[max(which(mean_frame$CumProb<=(1-alpha))),]$CumProb
    
    #Get lower lim, used for background, always a=~0.05
    lower_mean<-mean_frame[max(which(mean_frame$CumProb<=0.001)),]$means
    lower_alpha<-mean_frame[max(which(mean_frame$CumProb<=0.001)),]$CumProb
    
    #Check normal distribution
    setwd(debug_dir)
    qplot(means)+
      geom_histogram(col="black",fill="skyblue")+
      ylab("Frequency")+
      xlab("Sample Mean (j=5000, n=5000)")+
      ggtitle(paste("Sampling Distribution of Sample Means\nShapiro p=",round(shapiro.test(means)$p.value,digits=3),sep=""))+
      geom_vline(xintercept=upper_mean,col="red")+
      geom_text(aes(x=upper_mean,y=200,label=paste("Mean>",round(upper_mean,digits=4),"\np<",round(real_alpha,digits=4),sep="")),hjust=0,nudge_x=0.0001,col="red")+
      geom_vline(xintercept=lower_mean,col="blue")+
      geom_text(aes(x=lower_mean,y=200,label=paste("Mean<",round(lower_mean,digits=4),"\np<",round(lower_alpha,digits=4),sep="")),hjust=1,nudge_x=-0.0001,col="blue")+
      theme_bw()
    ggsave("1_DistributionCheck.png",width=7,height=5)
    
  #----Take 1000 samples of 100 px to cover entire image
    #Aggregate all samples w/ mean > upper_mean (p<0.05)
    #Remove px w/ signal below mean of aggregates samples
  
  band_sets<-data.frame()
  back_data<-data.frame()
  for(i in 1:((10*img_size^2)/sample_size)){
    x<-sample(colvec_data$Index,size=sample_size)
    #Get mean of x
    sampled<-colvec_data[x,]
    tmp_mean<-mean(sampled$Signal)
    if(tmp_mean>upper_mean){
      band_sets<-rbind(band_sets,sampled)
    } else if(tmp_mean<lower_mean){
      back_data<-rbind(back_data,sampled)
    }
  }

  og_kurt<-kurtosis(band_sets$Signal)
  
  #Ensure right-skew
    #Images taken during wash at Texp=3636 tend to have 10<kurtosis<15
  
  if(skew_status==TRUE){
    while(kurtosis(band_sets$Signal)<10){
      band_sets$Signal<-band_sets$Signal^1.1
    }
    
    while(kurtosis(band_sets$Signal)>20){
      band_sets$Signal<-band_sets$Signal^0.9
    }
  }
  
    if(use_mean==FALSE){
    band_signals<-as.numeric(levels(as.factor(band_sets$Signal)))
    sse_df<-data.frame()
    for(i in 2:length(band_signals)){
      x<-subset(band_sets,Signal<band_signals[i])
        x_sse<-sum((x$Signal-mean(x$Signal))^2)/sd(x$Signal)
      y<-subset(band_sets,Signal>=band_signals[i])
        y_sse<-sum((y$Signal-mean(y$Signal))^2)
      tot_sse<-sum(x_sse,y_sse)
      tmp<-data.frame(Split=band_signals[i],SSE=tot_sse)
      sse_df<-rbind(sse_df,tmp)
    }
    sse_df<-sse_df[order(sse_df$SSE,decreasing=FALSE),]
    split<-sse_df$Split[1]
  } else{
    split<-mean(band_sets$Signal)
  }
  
  #Plot mean cutoff
  setwd(debug_dir)
  qplot(band_sets$Signal)+
    geom_histogram(col="black",fill="skyblue")+
    geom_vline(xintercept=split,col="red")+
    xlab("Signal")+
    ylab("Frequency")+
    ggtitle(paste("Significant Sample Signal Distribution, Kurtosis=",round(kurtosis(band_sets$Signal),digits=2),", Og=",round(og_kurt,digits=2),sep=""))+
    geom_text(aes(x=split,y=10000,label=paste("Split=",round(split,digits=2),sep="")),nudge_x=0.01,hjust=0,col="red")+
    theme_bw()
  ggsave("2_BandDataDistributionCheck.png",width=7,height=5)
  
  band_sets_cur<-subset(band_sets,Signal>split)
  back_sets_cur<-subset(back_data,Signal<mean(back_data$Signal))
  
  #plot(band_sets_cur$Col,band_sets_cur$Row,xlim=c(1,100),ylim=c(100,1),cex=0.1)
  
  #Remove duplicates
  bands<-band_sets_cur[-which(duplicated(band_sets_cur$Index)==TRUE),]
  back<-back_sets_cur[-which(duplicated(back_sets_cur$Index)==TRUE),]
  
  #plot(bands$Col,bands$Row,xlim=c(1,100),ylim=c(100,1))
  
  #------------------------------------------------
  #-----------IDENTIFY OBJECTS---------------------
  #------------------------------------------------
  
  #Identify objects according to clustering
    #Object is group of px within 4 rows and cols of each other
  object_count<-1
  band_rows<-as.numeric(levels(as.factor(bands$Row)))
  obj_num<-c()
  for(i in 1:length(band_rows)){
    #Get previous row
    x<-band_rows[i-1]
    #Get current row
    y<-band_rows[i]
    
    #Get height
    h<-abs(y-x)
    
    if(h>band_height&&is.na(h)==FALSE&&length(h)>0){
      object_count<-object_count+1
    } 
    
    obj_num<-c(object_count,obj_num)
  }
  #Separate obj by row
  obj_rows<-data.frame(Row=band_rows,Obj=rev(obj_num))
  
  #Create object database
  obj_database<-data.frame()
  
  #Separate by col
  for(i in levels(as.factor(obj_rows$Obj))){
    x<-subset(obj_rows,Obj==i)
    obj_data<-subset(bands,Row%in%x$Row)
    
    object_count_2<-1
    obj_num_2<-c()
    band_cols<-as.numeric(levels(as.factor(obj_data$Col)))
    for(m in 1:length(band_cols)){
      a<-band_cols[m-1]
      b<-band_cols[m]
      
      width=abs(a-b)
      
      if(width>band_width&&is.na(width)==FALSE&&length(width)>0){
        object_count_2<-object_count_2+1
      } 
      
      obj_num_2<-c(obj_num_2,object_count_2)
    }
    tmp_obj<-data.frame(Col=band_cols,Obj=obj_num_2)
    tmp_obj$ObjID<-paste(levels(as.factor(x$Obj)),"-",tmp_obj$Obj,sep="")
    sep_tmp<-left_join(obj_data,tmp_obj,by="Col")
    obj_database<-rbind(obj_database,sep_tmp)
    
    #End row obj cycle
  }
  
  obj_list<-levels(as.factor(obj_database$ObjID))
  
 #Plot obj database
  gradient<-colorRampPalette(c("yellow3", "darkorange","orangered","mediumpurple2","royalblue","black"))(length(obj_list))
  setwd(plots_dir)
  png("3a_Objects.png")
  plot(x=NULL,y=NULL,xlim=c(1,100),ylim=c(100,-3),
       ylab="Row",xlab="Col")
  for(i in obj_list){
    x<-subset(obj_database,ObjID==i)
    points(x=x$Col,y=x$Row,cex=1.5,pch=19,col=gradient[which(obj_list%in%i)])
    text(x=mean(x$Col),y=min(x$Row)-4,i,col=gradient[which(obj_list%in%i)])
  }
  dev.off()
  
  setwd(plots_dir)
  png("3b_Background.png")
  #Plot obj database
  plot(x=back$Col,y=back$Row,xlim=c(1,100),ylim=c(100,1),
       ylab="Row",xlab="Col")
  dev.off()
  
  #Prompt user for object input
    
  #Open object image
  shell.exec(paste(plots_dir,"/3a_Objects.png",sep=""))
  
    #Input function
    enter_fun<-function(h){
      .GlobalEnv$save_obj_in<-tclvalue(save_obj)
      .GlobalEnv$keep_obj<-str_split(save_obj_in,pattern=",")[[1]]
      tkdestroy(win)
    }
    
    #Create window
    win<-tktoplevel()
    tkwm.geometry(win,"200x70")
    in_frame<-tkframe(win)
    label<-tklabel(in_frame,text="Enter Object IDs")
    tkpack(label,in_frame)
    input<-tkentry(in_frame,textvariable=save_obj)
    tkpack(input,in_frame)
    enter<-tkbutton(in_frame,text="Enter",command=enter_fun)
    tkpack(enter,in_frame)
    tkwait.window(win)
    
  #End if cycle=1
  }
  
  #------------------------------------------------
  #-----------PROCESS OBJECTS---------------------
  #------------------------------------------------
  
  if("ALL"%in%keep_obj){
    proc_data<-obj_database
    gradient_cur<-gradient
  } else{
    proc_data<-obj_database[which(obj_database$ObjID%in%keep_obj),]
    gradient_cur<-gradient[which(obj_list%in%keep_obj)]
  }
  
  metric_data<-left_join(colvec_data,proc_data,by="Index")
  metric_data_cur<-metric_data[which(!is.na(metric_data$ObjID)),]
  
  #Get area and signal
  obj_metrics<-data.frame()
  for(i in levels(as.factor(metric_data_cur$ObjID))){
    x<-subset(metric_data_cur,ObjID==i)
    signal<-sum(x$Signal.x)
    area<-length(x$Index)
    tmp<-data.frame(Object=i,Signal=signal,Area=area)
    obj_metrics<-rbind(obj_metrics,tmp)
  }
  obj_metrics$Cycle=cycle
  
  #------------------------------------------------
  #-----------PROCESS BACKGROUND---------------------
  #------------------------------------------------
  
  #Get back signal and area
  back_signal<-round(sum(colvec_data[which(colvec_data$Index%in%back$Index),]$Signal),digits=3)
  back_area<-length(back$Index)
  back_Imean<-back_signal/back_area
  
  #Add to df
  obj_metrics$BackgroundSignal<-back_signal
  obj_metrics$B<-round(obj_metrics$Area*back_Imean,digits=3)
  
  #Order output
  ordered_metrics<-data.frame(Cycle=obj_metrics$Cycle,
                              Object=obj_metrics$Object,
                              RawSignal=round(obj_metrics$Signal,digits=3),
                              BackgroundSignal=obj_metrics$BackgroundSignal,
                              Area=obj_metrics$Area,
                              B=obj_metrics$B)
  
  all_obj_metrics<-rbind(all_obj_metrics,ordered_metrics)
  
  cycle<-cycle+1
}

#Curate obj metrics
all_obj_metrics$Object<-paste("ID: ",all_obj_metrics$Object,sep="")

#Calculate background-adj signal
all_obj_metrics$Signal<-all_obj_metrics$RawSignal-all_obj_metrics$B

#Calculate SBR
all_obj_metrics$SBR<-round(all_obj_metrics$Signal/all_obj_metrics$B,digits=3)

#Add rel_sig if needed
all_obj_metrics$RelSig<-NA

#Reorder columns
all_obj_metrics<-all_obj_metrics[,c(1,2,7,9,3,4,5,6,8)]

if(length(file)>1){
  
  #Get relative signal
  for(i in obj_list){
    tmp_id<-paste("ID: ",i,sep="")
    x<-subset(all_obj_metrics,Object==tmp_id)
    for(m in x$Cycle){
      y<-subset(x,Cycle==m)
      rel_sig<-round(y$Signal/x$Signal[1],digits=2)
      all_obj_metrics[which(all_obj_metrics$Object==tmp_id&all_obj_metrics$Cycle==m),]$RelSig<-rel_sig
    }
  }
  
  
  setwd(plots_dir)
  ggplot(all_obj_metrics,aes(as.factor(Cycle),Signal,color=Object,group=Object))+
    geom_point()+
    geom_line()+
    geom_text(aes(as.factor(Cycle),Signal*1.05,label=round(Signal,digits=2)),show.legend=FALSE)+
    scale_y_continuous(limits=c(0,NA),n.breaks=10)+
    scale_color_manual(values=gradient_cur)+
    xlab("Cycle")+
    theme_bw()
  ggsave("4_Signal.png",width=7,height=5)
  
  setwd(plots_dir)
  ggplot(all_obj_metrics,aes(as.factor(Cycle),SBR,color=Object,group=Object))+
    geom_point()+
    geom_line()+
    geom_text(aes(as.factor(Cycle),SBR*1.05,label=round(SBR,digits=2)),show.legend=FALSE)+
    scale_y_continuous(limits=c(0,NA),n.breaks=10)+
    scale_color_manual(values=gradient_cur)+
    xlab("Cycle")+
    theme_bw()
  ggsave("5_SBR.png",width=7,height=5)
  
  setwd(plots_dir)
  ggplot(all_obj_metrics,aes(as.factor(Cycle),BackgroundSignal,color=Object,group=Object))+
    geom_point()+
    geom_line()+
    geom_text(aes(as.factor(Cycle),BackgroundSignal*1.05,label=round(BackgroundSignal,digits=2)),show.legend=FALSE)+
    scale_y_continuous(limits=c(0,NA),n.breaks=10)+
    scale_color_manual(values=gradient_cur)+
    xlab("Cycle")+
    theme_bw()
  ggsave("6_BackgroundSignal.png",width=7,height=5)
  
  setwd(plots_dir)
  ggplot(all_obj_metrics,aes(as.factor(Cycle),RelSig,color=Object,group=Object))+
    geom_point()+
    geom_line()+
    geom_text(aes(as.factor(Cycle),RelSig*1.05,label=round(RelSig,digits=2)),show.legend=FALSE)+
    scale_y_continuous(limits=c(0,1.1),n.breaks=10)+
    scale_color_manual(values=gradient_cur)+
    xlab("Cycle")+
    ylab("Relative Signal")+
    theme_bw()
  ggsave("7_RelativeSignal.png",width=7,height=5)
  
}

setwd(plots_dir)
write.csv(all_obj_metrics,"8_Metrics.csv",row.names = FALSE)
write.table(parms,"Parms.txt")

#-----------ANALYZE DATA LOGS------------------

if(length(data_logs)>0){
  
  run_data<-data.frame()
  for(i in 1:length(data_logs)){
    log<-data_logs[i]
    file_path<-paste(str_sub(raw_files[1],end=str_locate_all(raw_files[1],"2022")[[1]][,2][[1]]-4),log,sep="")
    tmp<-as.data.frame(read.table(file = file_path, sep = '\t', header = FALSE))
    if(length(which(tmp[1,]%in%"tec_temp"))==0){
      tmp$V12<-NA
    }
    tmp$Cycle=i
    run_data<-rbind(run_data,tmp)
  }
  names(run_data)<-run_data[1,]
  names(run_data)[length(names(run_data))]<-"Cycle"
  run_data<-run_data[-1,]
  
  run_data$T<-as.double(run_data$T)
  run_data$fluid<-as.factor(run_data$fluid)
  run_data$P0<-as.double(run_data$P0)
  run_data$P1<-as.double(run_data$P1)
  run_data$flow_temp_0<-as.double(run_data$flow_temp_0)
  run_data$flow_temp_1<-as.double(run_data$flow_temp_1)
  run_data$Q0<-as.double(run_data$Q0)
  run_data$Q1<-as.double(run_data$Q1)
  run_data$Cycle<-as.numeric(run_data$Cycle)
  
  run_data$T<-1:length(run_data$T)
  
  #dP, Q and flow temp
  run_data$dP<-as.numeric(abs(run_data$P1-run_data$P0))
  run_data$Q<-as.numeric(abs(run_data$Q1-run_data$Q0)/2)
  run_data$FlowTemp<-as.numeric(abs(run_data$flow_temp_1-run_data$flow_temp_0)/2)
  
  setwd(plots_dir)
  ggplot(run_data,aes(x=T,y=dP,color=fluid))+
    geom_point()+
    scale_color_discrete(name="Reagent")+
    scale_y_continuous(limits=c(0,4),n.breaks=10)+
    xlab("Time (s)")+
    ggtitle(run_id)+
    theme_bw()
  ggsave("9_RunStats.png",width=7,height=5)
  
  #Get medians per cycle
  vars<-c("dP","Q","FlowTemp")
  sum_stats<-data.frame()
  for(i in levels(as.factor(run_data$Cycle))){
    x<-subset(run_data,Cycle==i&fluid=="Ab2")
    for(m in vars){
      y<-x[,which(names(x)%in%m)]
      tmp<-data.frame(Cycle=i,Var=m,Value=mean(y))
      sum_stats<-rbind(sum_stats,tmp)
    }
  }
  
  #Get regression for each object
  reg_data<-data.frame()
  for(i in levels(as.factor(all_obj_metrics$Object))){
    x<-subset(all_obj_metrics,Object==i)
    for(m in vars){
      y<-subset(sum_stats,Var==m)
      fit<-lm(x$RelSig~y$Value)
      sum<-summary(fit)
      adj_rsqr<-sum$adj.r.squared
      if(is.null(sum$fstatistic)==FALSE){
        high_tail_p_value<-pf(sum$fstatistic[1],sum$fstatistic[2],sum$fstatistic[3],lower.tail = FALSE)[[1]]
        low_tail_p_value<-pf(sum$fstatistic[1],sum$fstatistic[2],sum$fstatistic[3],lower.tail = TRUE)[[1]]
        p_value<-min(c(high_tail_p_value,low_tail_p_value))
      } else{
        p_value<-1
      }
      tmp<-data.frame(Object=i,
                      Var=m,
                      Rsqr=round(adj_rsqr,digits=3),
                      P_value=round(p_value,digits=3))
      reg_data<-rbind(reg_data,tmp)
    }
  }
  
  setwd(plots_dir)
  write.csv(sum_stats,"10_SummaryStats.csv",row.names = FALSE)
  write.csv(reg_data,"11_RegressionOutput.csv",row.names = FALSE)
  
  if(length(which(reg_data$P_value<0.05))>0){
  
  reg_data$Signif<-NA
  reg_data[which(reg_data$P_value<0.05),]$Signif<-TRUE
  reg_data[which(reg_data$P_value>=0.05),]$Signif<-FALSE
  
  signif_data<-subset(reg_data,Signif==TRUE)
  
    for(i in signif_data$Object){
      x<-subset(signif_data,Object==i)
      for(m in x$Var){
        #Get signal data for object
        met<-subset(all_obj_metrics,Object==i)$RelSig
        
        #Get run data for var
        var_stats<-subset(sum_stats,Var==m)$Value
        
        #Get model
        lm_model<-lm(met~var_stats)
        lm_sum<-summary(lm_model)
        coef<-as.data.frame(lm_sum$coefficients)
        eqt<-paste("y= ",round(coef$Estimate[2],digits=3),"x+",round(coef$Estimate[1],digits=3),sep="")
        p<-round(pf(lm_sum$fstatistic[1],lm_sum$fstatistic[2],lm_sum$fstatistic[3],lower.tail = TRUE)[[1]],digits=3)
        rsqr<-round(lm_sum$r.squared,digits=3)
        
        tmp_obj<-paste(str_extract_all(i,"[:digit:]")[[1]],collapse="-")
        
        setwd(plots_dir)
        plot_data<-data.frame(RelSig=met,VarData=var_stats)
        ggplot(plot_data,aes(x=VarData,y=RelSig,group=1))+
          geom_point()+
          geom_smooth(method="lm",lwd=0.5)+
          scale_y_continuous(limits=c(0,NA),n.breaks=10)+
          geom_text(aes(x=min(VarData),y=0.1,label=paste(eqt,"\nRsqr= ",rsqr,"\np= ",p,sep="")),hjust=0,col="blue")+
          ylab("Relative Signal")+
          xlab(m)+
          theme_bw()
        ggsave(paste("12__",tmp_obj,"_",m,".png",sep=""),width=7,height=5)
      }
    }
  }
  
  for(i in vars){
    x<-subset(sum_stats,Var==i)

    setwd(plots_dir)
    ggplot(x,aes(x=as.factor(Cycle),y=Value,group=1))+
    geom_point()+
    geom_line()+
    xlab("Cycle")+
    geom_text(aes(x=as.factor(Cycle),y=Value*1.05,label=round(Value,digits=3)))+
    scale_y_continuous(limits=c(0,NA),n.breaks=10)+
    ylab(i)+
    theme_bw()
    ggsave(paste("13_",i,".png",sep=""),width=7,height=5)
}
  
  #END IF
}

#Open run
shell.exec(plots_dir)
