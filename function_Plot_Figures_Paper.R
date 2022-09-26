
## To plot time series of proba and ToE
fastplot_time_series<-function(array_prob_68, array_prob_95,p_i2, p_i1,
                               array_dotted_68, array_dotted_95, ylim_=c(-0.01,0.2), 
                               main_=NaN,
                               coord_baseline=NaN,coord_sub_ts=NaN,subplot_name=NaN,
                               length_sw=30){
  toe_red_68=col=rgb(1,0,0, alpha=0.2)
  tmp_toe_68=compute_ToE_matrix_from_array(array_prob_68,length_sw, "1850_2100", coord_baseline)
  tmp_toe_95=compute_ToE_matrix_from_array(array_prob_95,length_sw, "1850_2100", coord_baseline)
  toe_red_95=col=rgb(1,0,0, alpha=0.1)
  
  plot(array_prob_68[p_i2,p_i1,coord_sub_ts,1],type='l',ylim=ylim_,
       ylab="Proba.",  xaxt="n", xlab="Sliding windows", main=paste0(main_),
       cex.main=2, cex.lab=1.7, cex.axis=1.7)
  axis(1,xlab_pos,xlab_name, cex.axis=1.7)
  polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(array_prob_68[p_i2, p_i1,coord_sub_ts,2],rev(array_prob_68[p_i2,p_i1,coord_sub_ts,3])),col= toe_red_68, border = FALSE)
  polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(array_prob_95[p_i2, p_i1,coord_sub_ts,2],rev(array_prob_95[p_i2,p_i1,coord_sub_ts,3])),col= toe_red_95, border = FALSE)
  #abline(v=1,col=rgb(0,0,0, alpha=0.1),lwd=2)
  abline(h=min(array_prob_68[p_i2,p_i1,coord_baseline,2]), col=toe_red_68, lwd=2)
  abline(h=max(array_prob_68[p_i2,p_i1,coord_baseline,3]), col=toe_red_68, lwd=2)
  abline(v=tmp_toe_68$coord_TOE[p_i2,p_i1]-coord_sub_ts[1]+1, col=toe_red_68, lwd=2)
  
  abline(h=min(array_prob_95[p_i2,p_i1,coord_baseline,2]), col=toe_red_95, lwd=2, lty=2)
  abline(h=max(array_prob_95[p_i2,p_i1,coord_baseline,3]), col=toe_red_95, lwd=2, lty=2)
  abline(v=tmp_toe_95$coord_TOE[p_i2,p_i1]-coord_sub_ts[1]+1, col=toe_red_95, lwd=2, lty=2)
  
  # lines(array_dotted_68[p_i2,p_i1,coord_sub_ts,2], lty=2)
  # lines(array_dotted_68[p_i2,p_i1,coord_sub_ts,3], lty=2)
  # lines(array_dotted_95[p_i2,p_i1,coord_sub_ts,2], lty=2)
  # lines(array_dotted_95[p_i2,p_i1,coord_sub_ts,3], lty=2)
  if(!is.na(tmp_toe_68$TOE[p_i2,p_i1])){
    toe1=paste0(tmp_toe_68$TOE[p_i2,p_i1]," (", tmp_toe_68$TOE_text[p_i2,p_i1], ")")
  }else{
    toe1=paste0("n/a")
  }  
  if(!is.na(tmp_toe_95$TOE[p_i2,p_i1])){
    toe2=paste0(tmp_toe_95$TOE[p_i2,p_i1]," (", tmp_toe_95$TOE_text[p_i2,p_i1], ")")
  }else{
    toe2=paste0("n/a")
  }
  
  coord_toe_68=rep(NaN,length(coord_sub_ts))
  coord_toe_95=rep(NaN,length(coord_sub_ts))
#  k=0
#  for(i in coord_sub_ts){
#    k=k+1
#    if(array_prob_68[p_i2,p_i1,222,1]>=array_prob_68[p_i2,p_i1,i,3] |
#       array_prob_68[p_i2,p_i1,222,1]<=array_prob_68[p_i2,p_i1,i,2]){
#      coord_toe_68[k]=k
#    }
#    if(array_prob_95[p_i2,p_i1,222,1]>=array_prob_95[p_i2,p_i1,i,3] |
#       array_prob_95[p_i2,p_i1,222,1]<=array_prob_95[p_i2,p_i1,i,2]){
#      coord_toe_95[k]=k
#    }
#  }
  # abline(v=coord_toe, col=toe_red_68)
  leg.txt <- c(paste0("ToE 68%: ",  toe1),  paste0("ToE 95%: ",  toe2))
  legend("topleft",leg.txt,col=c(toe_red_68, toe_red_95),lty=c(1, 1),
         lwd=c(2,2), cex=1.5, bg = "white")
#  points(coord_toe_68,rep(ylim_[1]+0.01,length(coord_toe_68)),pch=3, col=toe_red_68)
#  points(coord_toe_95,rep(ylim_[1],length(coord_toe_95)),pch=3, col=toe_red_95)
  legend("topright", paste0(subplot_name),cex=2.5, bty="n")
}



### Fastplot ToE matrix
fastplot_toe_mat<-function(array_prob,p_i2, p_i1, main_=NaN,
                       coord_baseline=NaN,coord_sub_ts=NaN, subplot_name="NaN",
                       zlim_=c(2020,2086), length_sw=30){
  toe_red=col=rgb(1,0,0, alpha=0.2)
  tmp_toe=compute_ToE_matrix_from_array(array_prob,length_sw, "1850_2100",
                                        coord_baseline)
  image.plot(1:19, 1:19,postproc_image.plot(tmp_toe$TOE),zlim=zlim_,
               col=rdylbuPal, main=paste0(main_),
              cex.main=2, cex.lab=1.7, cex.axis=1.7, 
                            ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n",
                            yaxt="n",
                            axis.args=list(cex.axis=1.5))
  axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  bool_infTOE=which(postproc_image.plot(tmp_toe$TOE)<zlim_[1],arr.ind=TRUE)
  points(bool_infTOE, pch=20, col="red")
  legend("topleft", c(subplot_name), cex=2.5, bty="n")
 
} 

fastplot_difftoe_mat<-function(array_prob_margdep, array_prob_marg,
                               array_prob_dep, p_i2, p_i1, main_=NaN,
                       coord_baseline=NaN,coord_sub_ts=NaN){
  tmp_toe_margdep=compute_ToE_matrix_from_array(array_prob_margdep,30, "1850_2100",
                                        coord_baseline)
  tmp_toe_marg=compute_ToE_matrix_from_array(array_prob_marg,30, "1850_2100",
                                        coord_baseline)
  tmp_toe_dep=compute_ToE_matrix_from_array(array_prob_dep,30, "1850_2100",
                                        coord_baseline)
  diff_marg=tmp_toe_marg$TOE-tmp_toe_margdep$TOE
  diff_dep=tmp_toe_dep$TOE-tmp_toe_margdep$TOE
  image.plot(1:19, 1:19,postproc_image.plot(tmp_toe_margdep$TOE),zlim=c(2020,2086),
               col=rdylbuPal, main=paste0("Marg.-dep."), 
                            ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n",
                            yaxt="n",
                            cex.main=2, cex.lab=1.7, cex.axis=1.7,
                            axis.args=list(cex.axis=1.5))
  axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  bool_infTOE=which(postproc_image.plot(tmp_toe_margdep$TOE)<2020,arr.ind=TRUE)
  points(bool_infTOE, pch=20, col="red")
  legend("topleft", c("(a)"), cex=2.5, bty="n")
  image.plot(1:19,
             1:19,postproc_image.plot(diff_dep),zlim=c(-100,100),
               col=rbPal, main=paste0("Dep.-(Marg.-dep.)"), 
                            ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n",
                            yaxt="n", cex.main=2, cex.lab=1.7, cex.axis=1.7,
                            axis.args=list(cex.axis=1.5))
  axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  bool_infTOE=which(postproc_image.plot(diff_dep)>100,arr.ind=TRUE)
  points(bool_infTOE, pch=20, col="red")
  legend("topleft", c("(b)"), cex=2.5, bty="n")
 
  image.plot(1:19,
             1:19,postproc_image.plot(diff_marg),zlim=c(-100,100),
               col=rbPal, main=paste0(main_), 
                            ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n",
                            yaxt="n",
                            cex.main=2, cex.lab=1.7, cex.axis=1.7,
                            axis.args=list(cex.axis=1.5))
  axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),
       cex.axis=1.7)
  bool_infTOE=which(postproc_image.plot(diff_marg)>100,arr.ind=TRUE)
  points(bool_infTOE, pch=20, col="red")
  legend("topleft", c("(c)"), cex=2.5, bty="n")
 

}



## To compute_Contrib
compute_Contrib_matrix_from_array<-function(array_margdep,array_marg, array_dep, coord_Ref){
  tmp_margdep=array_margdep[,,,1]
  tmp_marg=array_marg[,,,1]
  tmp_dep=array_dep[,,,1]
  p00<-replicate(dim(array_margdep)[3], array_margdep[,,coord_Ref,1], simplify="array")
  effet_total=tmp_margdep-p00
  effet_marg=tmp_marg- p00
  effet_dep=tmp_dep- p00
  effet_int= effet_total-effet_dep-effet_marg
  FAR_margdep= effet_total/tmp_margdep
  FAR_marg= effet_marg/tmp_margdep
  FAR_dep= effet_dep/tmp_margdep
  FAR_int= effet_int/tmp_margdep
  Ecart_relat_margdep= effet_total/p00
  Ecart_relat_margdep[is.infinite(Ecart_relat_margdep)]<-NaN
  Ecart_relat_marg= effet_marg/p00
  Ecart_relat_dep= effet_dep/p00
  Ecart_relat_int= effet_int/p00
  Contrib_marg=effet_marg/effet_total*100
  Contrib_dep=effet_dep/effet_total*100
  Contrib_int=effet_int/effet_total*100
  return(list(Contrib_marg=Contrib_marg, 
              Contrib_dep=Contrib_dep, 
              Contrib_int=Contrib_int,
              Ecart_relat_margdep=Ecart_relat_margdep,
              Ecart_relat_marg= Ecart_relat_marg,
              Ecart_relat_dep=Ecart_relat_dep,
              Ecart_relat_int=Ecart_relat_int,
              FAR_margdep=FAR_margdep,
              FAR_marg=FAR_marg,
              FAR_dep=FAR_dep,
              FAR_int=FAR_int))
}

## To plot Contrib time series, matrix and barplot
fastplot_Contrib<-function(input,p_i2,p_i1, plot_ts=NaN, plot_matrix=NaN,
                           plot_barplot=NaN,coord_sub_ts=NaN,
                           subplot_name=c("(d)", "(e)", "(f)")){
  Med_C_marg<-apply(input$Contrib_marg,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
  Med_C_dep<-apply(input$Contrib_dep,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
  Med_C_int<-apply(input$Contrib_int,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
  tmp_Contrib_marg=input$Contrib_marg[p_i2, p_i1,coord_sub_ts]
  tmp_Contrib_dep=input$Contrib_dep[p_i2, p_i1,coord_sub_ts]
  tmp_Contrib_int=input$Contrib_int[p_i2, p_i1,coord_sub_ts]
  coord_Contrib_marg_high<- which(tmp_Contrib_marg>150)
  coord_Contrib_marg_low<- which(tmp_Contrib_marg<(-150))
  coord_Contrib_dep_high<- which(tmp_Contrib_dep>150)
  coord_Contrib_dep_low<- which(tmp_Contrib_dep<(-150))
  coord_Contrib_int_high<- which(tmp_Contrib_int>150)
  coord_Contrib_int_low<- which(tmp_Contrib_int<(-150))
  tmp_Contrib_marg[(tmp_Contrib_marg>150 | tmp_Contrib_marg<(-150))]<-NaN
  tmp_Contrib_dep[(tmp_Contrib_dep>150 | tmp_Contrib_dep<(-150))]<-NaN
  tmp_Contrib_int[(tmp_Contrib_int>150 | tmp_Contrib_int<(-150))]<-NaN
  #col_main=c("dodgerblue","darkorange", "chartreuse3", "indianred1")
  if(plot_ts==TRUE){
    plot(input$FAR_margdep[p_i2, p_i1,coord_sub_ts], xaxt="n", ylab="Bivar. FAR", xlab="Sliding windows",type='l',
         ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7)
    abline(h=0)
    axis(1,xlab_pos,xlab_name, cex.axis=1.7)
    lines(input$FAR_marg[p_i2, p_i1,coord_sub_ts], type="l",
          col=col_main[1])
    lines(input$FAR_dep[p_i2, p_i1,coord_sub_ts], type="l",
          col=col_main[2])
    lines(input$FAR_int[p_i2, p_i1,coord_sub_ts], type="l",
          col=col_main[3])
    legend("bottomright", c(expression(paste(Delta, "P"^"FAR", sep="")),
                            expression(paste(Delta, "M"^"FAR", sep="")),
                            expression(paste(Delta,"D"^"FAR", sep="")), 
                            expression(paste(Delta, "I"^"FAR", sep=""))),
           lty=c(1,1,1), col=c("black", col_main[1],
                               col_main[2], col_main[3]), cex=1.5, bg="white")
    legend("topright", c("(d)"), cex=2.5, bty="n")
    if(sum(is.nan(input$Ecart_relat_margdep[p_i2, p_i1,coord_sub_ts]))==length(input$Ecart_relat_margdep[p_i2, p_i1,coord_sub_ts])){
      plot.new()
    }else{
      plot(input$Ecart_relat_margdep[p_i2, p_i1,coord_sub_ts], xaxt="n", 
           ylab="Relat. diff.", xlab="Sliding windows",type='l', ylim=c(-0.5,
                                                                        2),
                                                                          cex.axis=1.7,
                                                                          cex.lab=1.7)
      abline(h=0)
      axis(1,xlab_pos,xlab_name, cex.axis=1.7)
      lines(input$Ecart_relat_marg[p_i2, p_i1,coord_sub_ts], type="l", col=col_main[1])
      lines(input$Ecart_relat_dep[p_i2, p_i1,coord_sub_ts], type="l", col=col_main[2])
      lines(input$Ecart_relat_int[p_i2, p_i1,coord_sub_ts], type="l", col=col_main[3])
      legend("topleft",  c(expression(paste(Delta, "P"^"r. diff", sep="")),
                            expression(paste(Delta, "M"^"r. diff", sep="")),
                            expression(paste(Delta,"D"^"r. diff", sep="")), 
                            expression(paste(Delta, "I"^"r. diff", sep=""))),
             lty=c(1,1,1), col=c("black", col_main[1], col_main[2],
                                 col_main[3]), cex=1.5, bg="white")
     legend("topright", c("(e)"), cex=2.5, bty="n")
     }
    plot(tmp_Contrib_marg,ylim=c(-150,150), col=col_main[1], type="l",
         ylab="Contribution (%)", xlab="Sliding windows", xaxt="n",
         cex.axis=1.7, cex.lab=1.7)
    axis(1,xlab_pos,xlab_name, cex.axis=1.7)
    lines(tmp_Contrib_dep,ylim=c(-100,100), col=col_main[2], type='l')
    lines(tmp_Contrib_int,ylim=c(-100,100), col=col_main[3], type='l')
    abline(h=0,lty=2)
    abline(h=-100,lty=2)
    abline(h=100,lty=2)
#    lines(199:201,rep(Med_C_marg[p_i2, p_i1],3), pch=15, col=col_main[1], lwd=3)
#    lines(199:201,rep(Med_C_dep[p_i2, p_i1],3), pch=15, col=col_main[2], lwd=3)
#    lines(199:201,rep(Med_C_int[p_i2, p_i1],3), col=col_main[3], lwd=3)
    abline(h=Med_C_marg[p_i2, p_i1], col=col_main[1], lwd=2, lty=2)
    abline(h=Med_C_dep[p_i2, p_i1],  col=col_main[2], lwd=2, lty=2)
    abline(h=Med_C_int[p_i2, p_i1], col=col_main[3], lwd=2, lty=2)
 
    points(coord_Contrib_marg_high,rep(150,length(coord_Contrib_marg_high)),pch=8,
           col=col_main[1])
    points(coord_Contrib_marg_low,rep(-150,length(coord_Contrib_marg_low)),pch=8,
           col=col_main[1])
    points(coord_Contrib_dep_high,rep(150,length(coord_Contrib_dep_high)),pch=8,
           col=col_main[2])
    points(coord_Contrib_dep_low,rep(-150,length(coord_Contrib_dep_low)),pch=8,
           col=col_main[2])
    points(coord_Contrib_int_high,rep(150,length(coord_Contrib_int_high)),pch=8,
           col=col_main[3])
    points(coord_Contrib_int_low,rep(-150,length(coord_Contrib_int_low)),pch=8,
           col=col_main[3])
    legend("bottomright", c(expression(paste("Contrib"[Delta]["M"], sep="")),
                            expression(paste("Contrib"[Delta]["D"], sep="")),
                            expression(paste("Contrib"[Delta]["I"], sep=""))),
           lty=c(1,1,1), col=c( col_main[1],
                                col_main[2], col_main[3]), cex=1.5, bg="white")
    legend("topright", c("(f)"), cex=2.5, bty="n")
  }
  if(plot_matrix==TRUE){
    image.plot(1:19,1:19,postproc_image.plot(Med_C_marg),xaxt="n",
               yaxt="n",zlim=c(-30, 130),col=brbgPal, xlab=x_axis_lab,
               ylab=y_axis_lab,
               main=expression(~bold(paste("Median Contrib"[Delta]["M"],
                                           sep=""))),
               cex.main=2, cex.lab=1.7, cex.axis=1.7,
               axis.args=list(cex.axis=1.5))
    #points(which(postproc_image.plot(Med_C_marg)>=50,arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_marg)>=postproc_image.plot(Med_C_dep),arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_marg)<=(-30),arr.ind=TRUE), pch=20,
           col=brbgPal[1])
    axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
         cex.axis=1.7)
    axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),
         cex.axis=1.7)
    legend("topleft", c("(d)"), cex=2.5, bty="n")
    image.plot(1:19,1:19,postproc_image.plot(Med_C_dep),xaxt="n",
               yaxt="n",zlim=c(-30, 130),col=brbgPal, xlab=x_axis_lab,
               ylab=y_axis_lab,
               main=expression(~bold(paste("Median Contrib"[Delta]["D"],
                                           sep=""))),
               cex.main=2, cex.lab=1.7, cex.axis=1.7,
               axis.args=list(cex.axis=1.5))
    #points(which(postproc_image.plot(Med_C_dep)>=50,arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_dep)>=postproc_image.plot(Med_C_marg),arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_dep)<=(-30),arr.ind=TRUE), pch=20,
           col=brbgPal[1])
    axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
         cex.axis=1.7)
    axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),cex.axis=1.7)
    legend("topleft", c("(e)"), cex=2.5, bty="n")
    image.plot(1:19,1:19,postproc_image.plot(Med_C_int),xaxt="n",
               yaxt="n",zlim=c(-30, 130),col=brbgPal, xlab=x_axis_lab,
               ylab=y_axis_lab,
               main=expression(~bold(paste("Median Contrib"[Delta]["I"],
                                           sep=""))),
               cex.main=2, cex.lab=1.7, cex.axis=1.7,
               axis.args=list(cex.axis=1.5))
    axis(1,seq(2,18,by=4),round(quant_Init_index1_CNRMCM6[seq(2,18,by=4)],1),
         cex.axis=1.7)
    axis(2,seq(2,18,by=4),round(quant_Init_index2_CNRMCM6[seq(2,18,by=4)],1),
         cex.axis=1.7)
    legend("topleft", c("(f)"), cex=2.5, bty="n")
  }
  if(plot_barplot==TRUE){
    mat_count=matrix(NaN,ncol=1, nrow=3)
    rownames(mat_count)<-c("Marg.", "Dep.", "Int.")
    colnames(mat_count)<-c("Model CNRMCM6")#paste0(Model ", Mod))
    mat_count[,1]=c(Med_C_marg[p_i2,p_i1], Med_C_dep[p_i2,p_i1],
                    Med_C_int[p_i2,p_i1])
    barplot(mat_count, col=c(col_main[1], col_main[2], col_main[3]), beside=TRUE,
            ylab="Median contribution (%)", ylim=c(-150,150), las=1)
    legend("topright",
           legend = c("Marg.", "Dep.", "Int."),
           fill = c(col_main[1], col_main[2], col_main[3]))
    abline(h=0)
    abline(h=100,lty=2)
    #abline(h=50, lty=2)
    abline(h=-150, lty=2)
  }
}

#To plot Contrib matrix fast
fastplot_Contrib_matrix<-function(Med_C_marg, Med_C_dep, Med_C_int,
                                  plot_matrix=TRUE,coord_sub_ts=NaN,
                                  is_TardiveFrost=FALSE, subplot_name=c("(a)",
                                                                        "(b)",
                                                                        "(c)")){
  if(is_TardiveFrost==FALSE){
    mat_x_axis=round(probs_to_eval_index1_CNRMCM6[seq(2,18,by=2)],2)
    mat_y_axis=round(probs_to_eval_index2_CNRMCM6[seq(2,18,by=2)],2)
  }else{
    mat_x_axis=round(quant_Init_index1_CNRMCM6[seq(2,18,by=2)],2)
    mat_y_axis=round(quant_Init_index2_CNRMCM6[seq(2,18,by=2)],2)
  }

 if(plot_matrix==TRUE){
    image.plot(1:19,1:19,postproc_image.plot(Med_C_marg),xaxt="n",
               yaxt="n",zlim=c(-30, 130),col=brbgPal, xlab=x_axis_lab,
               ylab=y_axis_lab, main=expression(~bold(paste("Median Contrib"[Delta]["M"],
                                           sep=""))),
               cex.main=2, cex.axis=1.7, cex.lab=1.7,
               axis.args=list(cex.axis=1.7)) 
   #points(which(postproc_image.plot(Med_C_marg)>=50,arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_marg)>=postproc_image.plot(Med_C_dep),arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_marg)<=(-30),arr.ind=TRUE), pch=20,
           col=brbgPal[1])
    axis(1,seq(2,18,by=2),mat_x_axis, cex.axis=1.7)
    axis(2,seq(2,18,by=2),mat_y_axis, cex.axis=1.7)
    legend("topleft", subplot_name[1], bty="n", cex=2.5)
    image.plot(1:19,1:19,postproc_image.plot(Med_C_dep),xaxt="n",
               yaxt="n",zlim=c(-30, 130),col=brbgPal, xlab=x_axis_lab, ylab=y_axis_lab
               , cex.main=2, cex.axis=1.7, cex.lab=1.7,
               axis.args=list(cex.axis=1.7)
               ,main=expression(~bold(paste("Median Contrib"[Delta]["D"],
                                           sep=""))))
    #points(which(postproc_image.plot(Med_C_dep)>=50,arr.ind=TRUE), pch=2)
    points(which(postproc_image.plot(Med_C_dep)<=(-30),arr.ind=TRUE), pch=20,
           col=brbgPal[1])
    points(which(postproc_image.plot(Med_C_dep)>=postproc_image.plot(Med_C_marg),arr.ind=TRUE), pch=2)
    axis(1,seq(2,18,by=2),mat_x_axis, cex.axis=1.7)
    axis(2,seq(2,18,by=2),mat_y_axis, cex.axis=1.7)
    legend("topleft", subplot_name[2], bty="n", cex=2.5)
    image.plot(1:19,1:19,postproc_image.plot(Med_C_int),xaxt="n",
               yaxt="n",zlim=c(-30, 130),col=brbgPal, xlab=x_axis_lab, ylab=y_axis_lab,
               cex.main=2, cex.axis=1.7, cex.lab=1.7,
               axis.args=list(cex.axis=1.7)
               ,main=expression(~bold(paste("Median Contrib"[Delta]["I"],
                                           sep=""))))
    axis(1,seq(2,18,by=2),mat_x_axis, cex.axis=1.7)
    axis(2,seq(2,18,by=2),mat_y_axis, cex.axis=1.7)
    legend("topleft", subplot_name[3], bty="n", cex=2.5)
 }
}

## To plot multilines proba, ToE matrix and Sum
fastplot_Indiv<-function(list_array_prob,p_i2, p_i1, array_dotted=NaN,
                         ylim_=c(0,0.2), 
                         main_=NaN, coord_baseline=NaN,coord_sub_ts=NaN, plot_proba=FALSE, plot_toe_mat=FALSE,
                         plot_sum=FALSE, plot_iqr=FALSE, label_68_95="NaN",
                         is_TardiveFrost=FALSE, subplot_name="NaN",
                         subplot_pos="bottomleft",
                         zlim_=c(2020,2086)){
  if(is_TardiveFrost==FALSE){
    mat_x_axis=round(probs_to_eval_index1_CNRMCM6[seq(2,18,by=2)],2)
    mat_y_axis=round(probs_to_eval_index2_CNRMCM6[seq(2,18,by=2)],2)
  }else{
    mat_x_axis=round(quant_Init_index1_CNRMCM6[seq(2,18,by=2)],2)
    mat_y_axis=round(quant_Init_index2_CNRMCM6[seq(2,18,by=2)],2)
  }
  tmp_mat_mod=matrix(NaN,ncol=length(names(list_array_prob)),
                     nrow=length(coord_sub_ts))
  toe_red_68=col=rgb(1,0,0, alpha=0.2)
  toe_red_95=col=rgb(1,0,0, alpha=0.05)
  k=0
  if(plot_proba==TRUE){
    plot(list_array_prob[[1]][p_i2,p_i1,coord_sub_ts,1],type='l',col="white",
         ylim=ylim_, ylab="Proba.",xlab="Sliding windows", xaxt="n",
         main=paste0(main_),
         cex.main=2, cex.axis=1.7, cex.lab=1.7)
    axis(1,xlab_pos,xlab_name, cex.axis=1.7)
    count_ToE=0
    indiv_ToE=c()
    for(i_Mod in names(list_array_prob)){
      k=k+1
      tmp_proba=list_array_prob[[i_Mod]][p_i2,p_i1,coord_sub_ts,1]
      tmp_mat_mod[,k]<-tmp_proba
      lines(tmp_proba, type='l', col=toe_red_68)
      tmp_ci_low_68=list_array_prob[[i_Mod]][p_i2,p_i1,coord_sub_ts,2]
      tmp_ci_high_68=list_array_prob[[i_Mod]][p_i2,p_i1,coord_sub_ts,3]
      polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(tmp_ci_low_68,rev(tmp_ci_high_68)),
              col= toe_red_95, border = FALSE)
      tmp_toe=compute_ToE_matrix_from_array(list_array_prob[[i_Mod]],30, "1850_2100", coord_baseline)
      if(!is.na(tmp_toe$coord_TOE[p_i2,p_i1])){
        indiv_ToE=c(indiv_ToE,tmp_toe$TOE[p_i2,p_i1])
        count_ToE=count_ToE+1
        print(paste0(count_ToE, ", Mod ", i_Mod, " has a ToE of: ", tmp_toe$TOE[p_i2,p_i1]))
      }
      abline(v=tmp_toe$coord_TOE[p_i2,p_i1]-coord_sub_ts[1]+1, col=toe_red_68,lty=1)
    }
    lines(apply(tmp_mat_mod, 1, mean), col="black", lty=2)
    #polygon(c(coord_baseline,rev(coord_baseline)),c(rep(-1000,length(coord_baseline)),rev(rep(6000,length(coord_baseline)))),
    #        col=rgb(0,0,0, alpha=0.1), border = FALSE)
    if(!is.null(median(indiv_ToE))){
      toe1=paste0(trunc(median(indiv_ToE))," (", trunc(median(indiv_ToE))-15,"-", trunc(median(indiv_ToE))+14, ")")
      abline(v=median(indiv_ToE)-1864-coord_sub_ts[1]+1, col="red",lty=1)
      leg.txt <- c(paste0("Median ToE ", label_68_95,"%: ",  toe1),paste0("Nb. models w. ToE ", label_68_95,"%: ",  count_ToE))
      legend("topleft",leg.txt,col=c("red",NaN),lty=c(1,NaN), cex=1.5,
             bg="white")
   }else{
      leg.txt <- c(paste0("Median ToE ", label_68_95,"%: n/a"),paste0("Nb. models w. ToE ", label_68_95,"%: ",  count_ToE))
      legend("topleft",leg.txt,col=c("red",NaN),lty=c(1,NaN), cex=1.5, bg="white")
    }
    legend(subplot_pos, subplot_name, cex=2.5, bty="n")
 }
  ### ToE-Indiv
  list_toe=array(NaN, dim=c(19,19,length(names(list_array_prob))))
  k=0
  for(i_Mod in names(list_array_prob)){
    k=k+1
    list_toe[,,k]=compute_ToE_matrix_from_array(list_array_prob[[i_Mod]],30, "1850_2100", coord_baseline)$TOE
  }
  res_median=apply(list_toe,c(1,2),median,na.rm=TRUE)
  if(plot_toe_mat==TRUE){
   image.plot(1:19, 1:19,postproc_image.plot(res_median),zlim=zlim_, col=rdylbuPal, main=paste0(main_), 
               ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n", yaxt="n",
               cex.main=2, cex.lab=1.7, cex.axis=1.7,
               axis.args=list(cex.axis=1.5))
    axis(1,seq(2,18,by=2),mat_x_axis, cex.axis=1.7)
    axis(2,seq(2,18,by=2),mat_y_axis, cex.axis=1.7)
    bool_infTOE=which(postproc_image.plot(res_median)<zlim_[1],arr.ind=TRUE)
    points(bool_infTOE, pch=20, col="red")
    legend("topleft", subplot_name, cex=2.5, bty="n")
  }
  if(plot_sum==TRUE){
    #lab.breaks=c(0,1,3,5,7,9,11, 13, 15)
    #breaks = c(0,0.9, 2.9, 4.9,6.9, 8.9,10.9, 12.9, 14.9)
    #lab.breaks=c(0,1,3,5,7,9,11, 13)
    #breaks = c(0,0.9 , 2.9, 4.9 ,6.9, 8.9 ,10.9, 12.9)
    lab.breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12, 13)
    breaks = c(0,0.9,1.9, 2.9,3.9, 4.9,5.9 ,6.9,7.9, 8.9,9.9 ,10.9,11.9, 12.9,
               13.9)

    res_sum=apply(!is.na(list_toe),c(1,2),sum,na.rm=TRUE)
    image.plot(1:19, 1:19,postproc_image.plot(res_sum),
               #lab.breaks=lab.breaks,
               breaks= breaks,
               #axis.args=list( at=c(0.5,seq(1,14,2)+1),labels=c(0, "1-2",
               #                                                    "3-4",
               #                                                    "5-6",
               #                                                    "7-8",
               #                                                    "9-10",
               #                                                    "11-12",
               #                                                    "13-14")),
               axis.args=list( at=c(0.5,seq(1,13,1)+0.4), labels=c(0:13),
                              cex.axis=1.5),
               col=c("white",magma13), zlim=c(0,13),
               ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n", yaxt="n",
               main=paste0(main_), cex.main=2, cex.lab=1.7, cex.axis=1.7)
    axis(1,seq(2,18,by=2),mat_x_axis, cex.axis=1.7)
    axis(2,seq(2,18,by=2),mat_y_axis, cex.axis=1.7)
    legend("topleft", subplot_name, cex=2.5, bty="n")
  }
  if(plot_iqr==TRUE){
    res_iqr=apply(list_toe,c(1,2),function(x){quantile(x,probs=0.75,na.rm=TRUE)-quantile(x,
                                                                                         probs=0.25,
                                                                                         na.rm=TRUE)})
    image.plot(1:19, 1:19,postproc_image.plot(res_iqr),
               col=rbPal[251:500], zlim=c(0,110),
               ylab=y_axis_lab, xlab=x_axis_lab, xaxt="n", yaxt="n",
               main=paste0(main_),
               cex.axis=1.7, cex.lab=1.7, cex.main=2,
               axis.args=list(cex.axis=1.5))
    axis(1,seq(2,18,by=2),mat_x_axis, cex.axis=1.7)
    axis(2,seq(2,18,by=2),mat_y_axis, cex.axis=1.7)
    legend("topleft", subplot_name, cex=2.5, bty="n")
  }
}

## To plot Full proba and ToE
fastplot_Full_time_series<-function(array_prob_68, p_i2, p_i1, array_dotted_68, 
                                    main_=NaN, ylim_=c(-0.01, 0.2),
                                    coord_baseline=NaN,coord_sub_ts=NaN,
                                    label_68_95="NaN", subplot_name="NaN",
                                    subplot_pos="bottomleft"){
  toe_red_68=col=rgb(1,0,0, alpha=0.2)
  tmp_toe_68=compute_ToE_matrix_from_array(array_prob_68,30, "1850_2100", coord_baseline)
  #tmp_toe_95=compute_ToE_matrix_from_array(array_prob_95,30, "1850_2100", coord_baseline)
  toe_red_95=col=rgb(1,0,0, alpha=0.1)
  
  plot(array_prob_68[p_i2,p_i1,coord_sub_ts,1],type='l',ylim=ylim_,
       ylab="Proba.",  xaxt="n", xlab="Sliding windows", main=paste0(main_),
       cex.main=2, cex.lab=1.7, cex.axis=1.7)
  axis(1,xlab_pos,xlab_name, cex.axis=1.7)
  polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(array_prob_68[p_i2, p_i1,coord_sub_ts,2],rev(array_prob_68[p_i2,p_i1,coord_sub_ts,3])),col= toe_red_68, border = FALSE)
  #polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(array_prob_95[p_i2, p_i1,coord_sub_ts,2],rev(array_prob_95[p_i2,p_i1,coord_sub_ts,3])),col= toe_red_95, border = FALSE)
#  abline(v=1,col=rgb(0,0,0, alpha=0.1),lwd=2)
  abline(h=min(array_prob_68[p_i2,p_i1,coord_baseline,2]), col=toe_red_68, lwd=2)
  abline(h=max(array_prob_68[p_i2,p_i1,coord_baseline,3]), col=toe_red_68, lwd=2)
  abline(v=tmp_toe_68$coord_TOE[p_i2,p_i1]-coord_sub_ts[1]+1, col=toe_red_68, lwd=2)
  
  #abline(h=min(array_prob_95[p_i2,p_i1,coord_baseline,2]), col=toe_red_95, lwd=2, lty=2)
  #abline(h=max(array_prob_95[p_i2,p_i1,coord_baseline,3]), col=toe_red_95, lwd=2, lty=2)
  #abline(v=tmp_toe_95$coord_TOE[p_i2,p_i1]-coord_sub_ts[1]+1, col=toe_red_95, lwd=2, lty=2)
  
  # lines(array_dotted_68[p_i2,p_i1,coord_sub_ts,2], lty=2)
  # lines(array_dotted_68[p_i2,p_i1,coord_sub_ts,3], lty=2)
  # lines(array_dotted_95[p_i2,p_i1,coord_sub_ts,2], lty=2)
  # lines(array_dotted_95[p_i2,p_i1,coord_sub_ts,3], lty=2)
  if(!is.na(tmp_toe_68$TOE[p_i2,p_i1])){
    toe1=paste0(tmp_toe_68$TOE[p_i2,p_i1]," (", tmp_toe_68$TOE_text[p_i2,p_i1], ")")
  }else{
    toe1=paste0("n/a")
  }  
  #if(!is.na(tmp_toe_95$TOE[p_i2,p_i1])){
  #  toe2=paste0(tmp_toe_95$TOE[p_i2,p_i1]," (", tmp_toe_95$TOE_text[p_i2,p_i1], ")")
  #}else{
  #  toe2=paste0("n/a")
  #}
  
  coord_toe_68=rep(NaN,length(coord_sub_ts))
  #coord_toe_95=rep(NaN,length(coord_sub_ts))
#  k=0
#  for(i in coord_sub_ts){
#    k=k+1
#    if(array_prob_68[p_i2,p_i1,222,1]>=array_prob_68[p_i2,p_i1,i,3] |
#       array_prob_68[p_i2,p_i1,222,1]<=array_prob_68[p_i2,p_i1,i,2]){
#      coord_toe_68[k]=k
#    }
#   # if(array_prob_95[p_i2,p_i1,222,1]>=array_prob_95[p_i2,p_i1,i,3] |
#   #    array_prob_95[p_i2,p_i1,222,1]<=array_prob_95[p_i2,p_i1,i,2]){
#   #   coord_toe_95[k]=k
#   # }
#  }
  # abline(v=coord_toe, col=toe_red_68)
  leg.txt <- c(paste0("ToE ", label_68_95, "%: ",  toe1))#,  paste0("ToE 95%: ",  toe2))
  legend("topleft",leg.txt,col=c(toe_red_68),lty=c(1), lwd=c(2), cex=1.5,
         bg="white")
#  points(coord_toe_68,rep(ylim_[1],length(coord_toe_68)),pch=3, col=toe_red_68)
  legend(subplot_pos, subplot_name, cex=2.5, bty="n")
  #points(coord_toe_95,rep(-0.01,length(coord_toe_95)),pch=3, col=toe_red_95)
  # text(1, 0.10, paste0(subplot_name),cex=2)
}


### With Ruban en option
fastplot_Contrib_Indiv<-function(list_Contrib, p_i2,p_i1, coord_sub_ts=NaN,
                                 subplot_name=c("(a)","(b)","(c)")){
  col_main=c("black","dodgerblue","darkorange", "gray",  "indianred1")
  #chartreuse3
  col_k=0
  eval(parse(text=paste0("plot(list_Contrib[['CNRMCM6']]$FAR_margdep[p_i2,
                         p_i1,coord_sub_ts], xaxt='n', ylab='Bivar. FAR',
xlab='Sliding windows',type='l',ylim=c(-1,1),col='white', cex.axis=1.7,
cex.lab=1.7, cex.main=2)")))
  axis(1,xlab_pos,xlab_name, cex.axis=1.7)
  abline(h=0)
  
  for(v in c("margdep", "marg", "dep","int")){
    col_k=col_k+1
    mat_far=matrix(NaN,ncol=length(names(list_Contrib)), nrow=length(coord_sub_ts))
    k=0
    for(i_Mod in names(list_Contrib)){
      k=k+1
      eval(parse(text=paste0("mat_far[,k]<-list_Contrib[[i_Mod]]$FAR_", v,"[p_i2, p_i1,coord_sub_ts]")))
    }
    lines(apply(mat_far,1,median),col=col_main[col_k])
    if(v=="int"){col_second=rgb(0.1, 0.8,0.3, alpha=0.1)}
    if(v=="marg"){col_second=rgb(0, 0,1, alpha=0.1)}
    if(v=="dep"){col_second=rgb(1, 0.3,0, alpha=0.1)}
    if(v=="margdep"){col_second=rgb(0, 0,0, alpha=0.1)}
    mat_far[is.infinite(mat_far) & mat_far<0]<-(-1000)
    mat_far[is.infinite(mat_far) & mat_far>0]<-(1000)
    # polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(apply(mat_far,1,min,na.rm=TRUE),rev(apply(mat_far,1,max,na.rm=TRUE))),
    #         col= col_second, border = FALSE)
  }
  legend("bottomright", c(expression(paste(Delta, "P"^"FAR", sep="")),
                            expression(paste(Delta, "M"^"FAR", sep="")),
                            expression(paste(Delta,"D"^"FAR", sep="")), 
                            expression(paste(Delta, "I"^"FAR", sep=""))),
         lty=c(1,1,1, 1), col=c("black", col_main[2], col_main[3],
                                col_main[4]), cex=1.5, bg="white")
  legend("topright", subplot_name[1], cex=2.5, bty="n")
  col_k=0
  eval(parse(text=paste0("plot(list_Contrib[['CNRMCM6']]$Ecart_relat_margdep[p_i2,
                         p_i1,coord_sub_ts], xaxt='n', ylab='Evol. relat.',
xlab='Sliding windows',type='l',ylim=c(-0.5,2),col='white', cex.axis=1.7,
cex.lab=1.7, cex.main=2)")))
  axis(1,xlab_pos,xlab_name, cex.axis=1.7)
  abline(h=0)
  for(v in c("margdep", "marg", "dep", "int")){
    col_k=col_k+1
    mat_far=matrix(NaN,ncol=length(names(list_Contrib)), nrow=length(coord_sub_ts))
    k=0
    for(i_Mod in names(list_Contrib)){
      k=k+1
      eval(parse(text=paste0("mat_far[,k]<-list_Contrib[[i_Mod]]$Ecart_relat_", v,"[p_i2, p_i1,coord_sub_ts]")))
    }
    lines(apply(mat_far,1,median),col=col_main[col_k])
    if(v=="int"){col_second=rgb(0.1, 0.8,0.3, alpha=0.1)}
    if(v=="marg"){col_second=rgb(0, 0,1, alpha=0.1)}
    if(v=="dep"){col_second=rgb(1, 0.3,0, alpha=0.1)}
    if(v=="margdep"){col_second=rgb(0, 0,0, alpha=0.1)}
    mat_far[is.infinite(mat_far) & mat_far<0]<-(-1000)
    mat_far[is.infinite(mat_far) & mat_far>0]<-(1000)
    # polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(apply(mat_far,1,min,na.rm=TRUE),rev(apply(mat_far,1,max,na.rm=TRUE))),
    #         col= col_second, border = FALSE)
  }
  legend("topleft",  c(expression(paste(Delta, "P"^"r. diff", sep="")),
                            expression(paste(Delta, "M"^"r. diff", sep="")),
                            expression(paste(Delta,"D"^"r. diff", sep="")), 
                            expression(paste(Delta, "I"^"r. diff", sep=""))),
         lty=c(1,1,1, 1), col=c("black",col_main[2], col_main[3], col_main[4]),
         cex=1.5, bg="white")
  legend("topright", subplot_name[2], cex=2.5, bty="n")
  tmp_Contrib<- eval(parse(text=paste0("list_Contrib[['CNRMCM6']]$Contrib_marg[p_i2, p_i1,coord_sub_ts]")))
  eval(parse(text=paste0("plot(tmp_Contrib, xaxt='n', ylab='Contrib.',
                         xlab='Sliding
                         windows',type='l',ylim=c(-150,150),col='white',
                         cex.axis=1.7, cex.lab=1.7, cex.main=2)")))
  axis(1,xlab_pos,xlab_name, cex.axis=1.7)
  col_k=1
  for(v in c("marg", "dep","int")){
    col_k=col_k+1
    tmp_Contrib<- eval(parse(text=paste0("list_Contrib[['CNRMCM6']]$Contrib_", v,"[p_i2, p_i1,coord_sub_ts]")))
    coord_Contrib_high<- which(tmp_Contrib>150)
    coord_Contrib_low<- which(tmp_Contrib<(-150))
    tmp_Contrib[(tmp_Contrib>150 | tmp_Contrib<(-150))]<-NaN
    mat_far=matrix(NaN,ncol=length(names(list_Contrib)), nrow=length(coord_sub_ts))
    k=0
    for(i_Mod in names(list_Contrib)){
      k=k+1
      tmp_Contrib_Mod=eval(parse(text=paste0("list_Contrib[['",i_Mod, "']]$Contrib_", v,"[p_i2, p_i1,coord_sub_ts]")))
      # coord_Contrib_high<- which(tmp_Contrib_Mod>150)
      # coord_Contrib_low<- which(tmp_Contrib_Mod<(-150))
      # tmp_Contrib_Mod[(tmp_Contrib_Mod>150 | tmp_Contrib_Mod<(-150))]<-NaN
      # eval(parse(text=paste0("lines(tmp_Contrib_Mod, ,type='l',col=col_main[",col_k,"])")))
      eval(parse(text=paste0("mat_far[,k]<-tmp_Contrib_Mod")))
    }
    lines(apply(mat_far,1,median,na.rm=TRUE), col=col_main[col_k])
    abline(h=median(apply(mat_far,1,median,na.rm=TRUE),na.rm=TRUE), col=col_main[col_k],lty=2,lwd=2)
    if(v=="marg"){col_second=rgb(0, 0,1, alpha=0.1)}
    if(v=="dep"){col_second=rgb(1, 0.3,0, alpha=0.1)}
    if(v=="int"){col_second=rgb(0.1, 0.8,0.3, alpha=0.1)}
    # polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(apply(mat_far,1,min,na.rm=TRUE),rev(apply(mat_far,1,max,na.rm=TRUE))),
    #         col= col_second, border = FALSE)
  }
  abline(h=-100,lty=2)
  abline(h=100,lty=2)
  abline(h=0,lty=2)
  legend("bottomright", c(expression(paste("Contrib"[Delta]["M"], sep="")),
                            expression(paste("Contrib"[Delta]["D"], sep="")),
                            expression(paste("Contrib"[Delta]["I"], sep=""))),
         lty=c(1,1, 1), col=c(col_main[2], col_main[3], col_main[4]), cex=1.5,
         bg="white")
  legend("topright", subplot_name[3], cex=2.5, bty="n")
}



fastplot_Boxplot_ToE<-function(list_array_prob_Indiv_margdep, 
                               list_array_prob_Indiv_marg,
                               list_array_prob_Indiv_dep,
                               list_array_prob_Full_margdep,
                               list_array_prob_Full_marg,
                               list_array_prob_Full_dep,
                               p_i2, p_i1, array_dotted=NaN,
                               ylim_=c(1980,2086),
                               main_=NaN, coord_baseline=NaN,coord_sub_ts=NaN,label_68_95="NaN",
                               is_TardiveFrost=FALSE, subplot_name="NaN"){
  nb_Mod_final_IndivEnsemble=length(names(get(paste0("list_array_prob_Indiv_margdep"))))
  for(v in c("margdep", "marg", "dep")){
    count_ToE=0
    indiv_ToE=c()
    for(i_Mod in names(get(paste0("list_array_prob_Indiv_", v)))){
      print(v)
      tmp_proba=get(paste0("list_array_prob_Indiv_", v))[[i_Mod]][p_i2,p_i1,coord_sub_ts,1]
      tmp_toe=compute_ToE_matrix_from_array(get(paste0("list_array_prob_Indiv_", v))[[i_Mod]],30, "1850_2100", coord_baseline)
      indiv_ToE=c(indiv_ToE,tmp_toe$TOE[p_i2,p_i1])
      if(!is.na(tmp_toe$coord_TOE[p_i2,p_i1])){
        count_ToE=count_ToE+1
      }
    }
    print(count_ToE)
    assign(paste0("tmp_TOE_Indiv_", v), indiv_ToE)
  }
  df <- data.frame( id = c(rep("Indiv",nb_Mod_final_IndivEnsemble))
                    ,
                    #        rep("Full", nb_Mod_final_FullEnsemble)),
                    Margdep = c(tmp_TOE_Indiv_margdep),
                    Marg = c(tmp_TOE_Indiv_marg),
                    Dep = c(tmp_TOE_Indiv_dep))
  proportion <- c(sum(!is.na(tmp_TOE_Indiv_margdep)),
                  sum(!is.na(tmp_TOE_Indiv_marg)),
                  sum(!is.na(tmp_TOE_Indiv_dep)))/(13)
  boxplot(df[,-1],boxwex=proportion, boxfill= rdylgnPal[350], 
          medcol=rdylgnPal[300],border=rdylgnPal[400], staplelty = 1,
                  lwd=1, lty=1, boxlwd = 1, ylab="ToE",
          outline=FALSE,  names=c("Marg.-dep.","Marg.", "Dep."), ylim=ylim_) 
  ### Full
  for(v in c("margdep", "marg", "dep")){
    full_ToE=c()
    for(i_Mod in names(get(paste0("list_array_prob_Full_", v)))){
      print(v)
      tmp_proba=get(paste0("list_array_prob_Full_", v))[[i_Mod]][p_i2,p_i1,coord_sub_ts,1]
      tmp_toe=compute_ToE_matrix_from_array(get(paste0("list_array_prob_Full_", v))[[i_Mod]],30, "1850_2100", coord_baseline)
      full_ToE=c(full_ToE,tmp_toe$TOE[p_i2,p_i1])
    }
    assign(paste0("tmp_TOE_Full_", v), full_ToE)
    print(full_ToE)
  }
  segments(0.85,tmp_TOE_Full_margdep, x1 = 1.15, lwd=3, col=rdylgnPal[500])
  #points(1.25,tmp_TOE_Full_margdep, pch=8, col=rdylgnPal[500],
  #       lwd=2)
#  points(1.75,tmp_TOE_Full_marg, pch=8, col=rdylgnPal[500], lwd=2)
#  points(3.25,tmp_TOE_Full_dep, pch=8, col=rdylgnPal[500], lwd=2)
  segments(1.85,tmp_TOE_Full_marg, x1 = 2.15, lwd=3, col=rdylgnPal[500])
  segments(2.85,tmp_TOE_Full_dep, x1 = 3.15, lwd=3, col=rdylgnPal[500])
  legend("bottomright",
         legend = c("Indiv-Ensemble", "Full-Ensemble"),
         col = c(rdylgnPal[350], rdylgnPal[500]),
         pch=c(15, NaN), lty=c(NaN,1), lwd=c(NaN,2))
  legend("topright", subplot_name, bty="n", cex=2)
}







