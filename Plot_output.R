# Script that creates basic figures from RINGS_V3 model output.
# A more comprehensive script related to the publication
# Friend et al 2022 is in the publication's supplementary materials.

#Authors: Andrew Friend, Annemarie Eckes-Shephard


res= 300
cex.legend = 0.5
cex.lab = 0.7
file_path <- ""
obs_path  <- file_path"
figs_path <- paste0(file_path,"Figures/")

# Parameters taken from the literature and calculated from observations
# used for conversions in this script
ttl = 50.1 # cell tangential Length [mu m]
tal = 2680. # Cell axial length [ mu m]
##############################################################################
#Functions used in data aggregation and plotting

##################################################################################
#bin_and_running_mean
#input: model_config_output numeric vector,
#input: bin_size numeric, desired size used for binning across the radial file
#input: maxlength_radial_file numeric, maximal radial file length that has come
# out of the simulations, it is used together with the bin size to derive the
# total number of bins for the aggregation.
#input: varname, string, whether binning is done over cell mass density ("den")
# or cell lengths ("L")
#output: a list of three numeric vectors with either cell length ('L') or
# cell mass density ('den'), together with cell distance from the cambium ('Di')
# and number of datapoints in the bin ('n')
bin_and_running_mean <- function(model_config_output,bin_size,maxlength_radial_file,varname){
  # author of original code: Andrew F, translated to R, Annemarie ES
  #binning:
  bs = bin_size
  mlen = maxlength_radial_file
  den_a = vector()
  pos_a = vector()
  ns_a   = vector()
  a = .0
  b = a + bs
  while (b < mlen+bs){
    n = length(t(model_config_output[varname]))
    i = 1
    j = .0
    dena_sum = .0
    while (i <=n){
      d = model_config_output$Di[i]
      # print(i)
      # print(d)
      if (d >= a & d < b){
        dena_sum = dena_sum + model_config_output[i,varname]
        j = j + 1.
      }
      i = i + 1
    }
    if (j > 0){
      den_a <- append(den_a,dena_sum/j)
      ns_a   <- append(ns_a,j)  # record how many data points were actually used to derive the number in this bin
      pos_a <- append(pos_a,a+(b-a)/2.)
    }

    a = a + bs
    b = b + bs
  }

  n = length(den_a)
  run = 3
  den_ar = vector()
  den_ar = append(den_ar,den_a[1])
  # running means across bins:
  for (i in 2:(n-1)){
    sden = .0
    for (j in (i-1):(i+1)){
      sden = sden + den_a [j]
    }
    den_ar = append(den_ar,(sden/run))
  }
  den_ar = append(den_ar,den_a[n])

  #prep output:
  if(varname=="L"){
    out <- list(den_ar,pos_a,ns_a)
    names(out)[1:3] <- c("L","Di","n")
    return(out)
  }
  if(varname=="den"){
    out <- list(den_ar,pos_a,ns_a)
    names(out)[1:3] <- c("den","Di","n")
    return(out)
  }

}
###################makeTransparent
#To show where overlap between predicted cells is greatest
#note: always pass alpha on the 0-255 scale
makeTransparent <- function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
################################ add_legend
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
###################add_cellNumberColumn
# quick and dirty addition of column to highlight each radial file. Not used here but maybe
# of use for people who want to compare against that output in their own studies.
# year-specific binning and plotting
# input: model_output RINGS_v2 model output from fort.97 table read in.
# output model output, with additional column
add_cellNumberColumn <- function(model_output){
  pos=1
  rfile = 1
  model_output$celln[1] <- pos
  model_output$row[1] <- rfile
  for(i in 2:dim(model_output)[1]){
    pos=pos+1
    if(model_output$Di[i-1] < model_output$Di[i]){
      model_output$celln[i] <- pos
      model_output$row[i] <- rfile
    }else{
      pos = 1
      rfile = rfile+1
      model_output$celln[i] <- pos
      model_output$row[i] <- rfile
    }
  }
  return(model_output)
}



###########FIGURE1:#############################################################
#FIGURE1:
# we normalis every simulated set of radial files for each model configuration run
bs = 0.02 # bin size (frac)
mlen = 1.0 # length to process (frac)


# Di = a tracheid's distance from the cambium [mu]
# L  = a tracheid's length ( radial diameter), this includes cell walls [mu m]
# M  = a tracheid's mass [g]
#load data for each model configuration run:

full_model <- read.table(paste0(file_path,'fort.97'))
names(full_model) <- c("Di", "L ", "M")
# add another cellnumber column, for analysis with tgram
full_model <- add_cellNumberColumn(full_model)


# some data manipulation and aggregation:

#turn radial files around
full_model$Di                 <-  1- full_model$Di


  full_model$Vc  <-  full_model$L * ttl * tal / 1.E12 # ml    eq 11 in manuscript
  full_model$den <- 1.E-6 *  full_model$M /  full_model$Vc # g/cm3
  # create binned trachs:
  collect_stdTrach_single <- bin_and_running_mean(model_config_output = full_model,
                                                       bin_size = bs,maxlength_radial_file = mlen,varname="den")

# plot:

jpeg(filename=paste0(figs_path,'annual_densprof_combined.jpeg'),units="mm",width=90, height=50, res=res)

par(mfrow=c(1,1),mar=c(0,2,0,0),oma=c(1,0,0.02,4),cex.axis=0.7, cex.lab=0.7, cex.main=1, cex.sub=1,ps=8,xpd=TRUE)

#Fig 1a)

# start with all simulated density values that were simulated across all 20 years:
plot(full_model$Di,full_model$den,pch=16,cex=0.12,ylab="",xlab="",yaxt="n",col=makeTransparent('grey',alpha=50),
     ylim=c(0,1.2),xaxt="n")

# add all the tracheidograms across all 20 years to plot
# grey lines are the (binned) dens - tracheidograms of the different years at that site.
# 100 radial files were modelled for each year between the simulation years 1976 to 1995.
# (note "row" is equivalent to an individually simulated radial file)
for(row in 1:20){
  n=100
  tmp <- bin_and_running_mean(model_config_output = full_model[which(full_model$row > row*n-99 & full_model$row <= row*n),],
                              bin_size = bs,maxlength_radial_file = mlen,varname="den")
  lines(tmp$Di, tmp$den,type="l",col="black",lwd=0.2)
}

# add binned tracheidogram derived from all years' simulations:
lines(collect_stdTrach_single$Di,collect_stdTrach_single$den,lwd=1)



## layout and legend:
axis(side = 2,lwd = 0, line = -0.8, las = 2)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, tck = .015, labels = NA)



axis(side = 1, tck = -.015, labels = NA)
axis(side = 1,lwd = 0, line = -1.3,cex=0.7)
mtext(side = 2 , expression(Density ~(g ~cm^-3)),line=-1 ,outer = TRUE,cex = cex.lab)

mtext(outer=FALSE,side=1,expression(Distance~from~inner~edge~(fraction~of~ring)),line=0.2,cex = cex.lab)




# Add legend between a&b, outside plot region
add_legend(0.62,0.25,xpd=TRUE,  legend = c("density/cell","binned \nprofiles"), title = "Model output",
           pch = c(16,NA), lty = c(NA,1), lwd= c(NA,1.5), col='grey',box.lwd = 0, cex = 0.7, bg = "transparent")



dev.off()
################################################################################



############APPENDIX FIGURE 1###################################################
#APPENDIX FIGURE 1 - plots all binned dens - radial files and the simulated densities
# from Figure 1b) from the model configurations that were not included in the graph
# plot all binned tracheidograms for the appendix, this time with a graph each
# year containing all model configurations' output.

pdf(file = paste0(figs_path,'annual_density_profiles.pdf'), width = 10, height = 6)

#par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(4,4,0.5,0.5))
par(mfrow = c(4,6), mar = c(0,0,0,0), oma = c(4,4,0.5,0.5))

years <- seq(1976,1995,1)

for(row in 1:20){
  c=1

      tmp_save <- bin_and_running_mean(model_config_output = full_model[which(full_model$row > row*n-99 & full_model$row <= row*n),], bin_size = bs,maxlength_radial_file = mlen,varname="den")

      plot(full_model[which(full_model$row > row*n-99 & full_model$row <= row*n),]$Di,full_model[which(full_model$row > row*n-99 & full_model$row <= row*n),]$den,pch=16,cex=0.5,ylab="",xlab="",yaxt="n",col="grey",
           ylim=c(0,1.2),xlim=c(0,1),axes = FALSE)
      lines(tmp_save$Di, tmp_save$den,type="l",col="black",lwd=1.2)


  #axes on the outer edges only

  #if(row==1|| row==5 || row==9 ||row==13||row==17){
  if(row==1|row==7|row==13|row==19){
    axis(side=2,las=2)
  }
  #if(row>=17){
  if(row>=15){
    axis(side=1,las=2)
  }

  box()
  par(new=TRUE)
  mtext(paste(years[row]),side=1,line=-1.7,adj=0.4)


}
# add legend
plot.new()
legend("center",c("Single cells","Binned profile"),
       lty=c(NA,1),pch=c(16,16),pt.cex=c(2,0),lwd=c(0,1),col=c("grey","black"),
       text.font=c(1,0.5),bty="n",cex=0.9)

mtext(outer=TRUE,side=1,"Distance from inner edge (fraction of ring)",line=2.4)
mtext(outer=TRUE,side=2,expression(Density ~(g ~cm^-3)),line=2.1)


dev.off()
################################################################################



############FIGURE2:##############################################


#load data for each model configuration run:

full_model <- read.table(paste0(file_path,'fort.94'))
names(full_model) <- c("Di", "L", "M")


bs = 40. # bin size (um)
mlen = 1850. # length to process (um)
maxlen=mlen

  full_model$Vc  <- full_model$L * ttl * tal / 1.E12 # ml    eq 11 in manuscript
  full_model$den <- 1.E-6 *  full_model$M /  full_model$Vc # g/cm3
  # create binned trachs:
  collect_stdTrach_single <- bin_and_running_mean(model_config_output = full_model,
                                                       bin_size = bs,maxlength_radial_file = mlen,varname="den")


#plot:

jpeg(filename=paste0(figs_path,'annual_densprof_single.jpeg'),units="mm",width=90, height=50, res=res)

par(mfrow=c(1,1),mar=c(0,2,0,0),oma=c(1,0,0.02,4),cex.axis=0.7, cex.lab=0.7, cex.main=1, cex.sub=1,ps=8,xpd=TRUE)

# start with all simulated density positions:
plot(full_model$Di,full_model$den,
     pch=16,cex=0.5,col=makeTransparent('grey',alpha=50),
     xlim=c(maxlen,0),ylim=c(0,1.01),ylab='',xlab='',xaxt='n',yaxt='n')

#add binned tracheidograms
lines( collect_stdTrach_single$Di,
       collect_stdTrach_single$den,type="l",col="black", lwd=2)


axis(side = 2,lwd = 0, line = -0.8, las = 2)
axis(side = 2, tck = -.015, labels = NA)
axis(1, at = seq(0,maxlen,by=200), labels = NA,tck= .015)

mtext(outer=FALSE,side=2,expression(Cell~mass~density),line=1.3,cex= cex.lab)
mtext(outer=FALSE,side=2,expression((g ~cm^-3)),line=0.8,cex= cex.lab,adj=0.5)


axis(1, at = seq(0,maxlen,by=200), labels = NA,tck= -.015, tick = TRUE)
axis(1, at = seq(0,maxlen,by=200),line = -1.3,lwd=0)


mtext(outer=FALSE,side=1,expression(Distance~from~cambium~(~mu*m)),line=0.2,cex= cex.lab)

# Add legend between a&b, outside plot region
add_legend(0.62,0.25,xpd=TRUE,  legend = c("density/cell","binned \nprofiles"),
           pch = c(16,NA), lty = c(NA,1), lwd= c(NA,1.5), col='grey',box.lwd = 0, cex = 0.7, bg = "transparent")



dev.off()
################################################################################


maxlen = mlen = 1910 #file length to process(um)
bs =100 #bin size (um)


#load data for each model configuration run:

full_model <- read.table(paste0(file_path,'fort.94'))
names(full_model) <- c("Di", "L", "M")


  full_model$Vc  <- full_model$L * ttl * tal / 1.E12 # ml    eq 11 in manuscript
  full_model$den <- 1.E-6 *  full_model$M /  full_model$Vc # g/cm3
  # create binned trachs:
  collect_stdTrach_single <- bin_and_running_mean(model_config_output = full_model,
                                                       bin_size = bs,maxlength_radial_file = mlen,varname="L")




#plot:
jpeg(filename=paste0(figs_path,'Cell_lengthprofile_single.jpeg'),units="mm",width=90, height=50, res=res)

par(mfrow=c(1,1),mar=c(0,2,0,0),oma=c(1,0,0.02,4),cex.axis=0.7, cex.lab=0.7, cex.main=1, cex.sub=1,ps=8,xpd=TRUE)

# add all simulated length positions:
plot(full_model$Di,full_model$L,pch=16,cex=0.8,col=makeTransparent('grey',alpha=50),xlim=c(maxlen,0),ylim=c(0,160),ylab='',xlab='',xaxt='n',yaxt='n')

# add binned tracheidograms
lines(collect_stdTrach_single$Di, collect_stdTrach_single$L,type="l",col="black", lwd=2)

axis(side = 2,lwd = 0, line = -0.8, las = 2)
axis(side = 2, tck = -.015, labels = NA)

mtext(outer=FALSE,side=2,expression(Cell~length ~(mu*m)),line=1,cex = cex.lab)

# Add legend between a&b, outside plot region
add_legend(0.62,0.25,xpd=TRUE,  legend = c("density/cell","binned \nprofiles"),
           pch = c(16,NA), lty = c(NA,1), lwd= c(NA,1.5), col='grey',box.lwd = 0, cex = 0.7, bg = "transparent")


axis(1, at = seq(0,maxlen,by=200), labels = NA,tck= -.015, tick = TRUE)
axis(1, at = seq(0,maxlen,by=200),line = -1.3,lwd=0)

mtext(outer=FALSE,side=1,expression(Distance~from~cambium~(~mu*m)),line=0.2,cex= cex.lab)

# note that grey cells  to the right of 0 are alive cambial cells for the next season.
dev.off()
################################################################################

####################IMPORTANT!!
# for the below figures to render, run the model with
#Tornetrask
#in driver.txt


############FIGURE maximum density simulations against obs###########################

jpeg(filename=paste0(figs_path,'MXD.jpeg'),units="mm",width=90, height=90, res=res)
par(mfrow=c(1,1),oma=c(1,1,0.1,0.1),mar=c(1,1,0,0),cex.axis=0.7, cex.lab=0.7,ps=8)


# calculate temperatures on the fly from forcing data:
out_temp= read.table(paste0(file_path,'torn_clm_crujrav2.2.txt'))
names(out_temp) <- c("kyr_clm", "kt", "Temp" )

jja = list() #mean jja temperature for each year (degC)
nyr = 2004 - 1901 + 1
nmo = 12
ndays = c(31,28,31,30,31,30,31,31,30,31,30,31)
m = 1
sum = 0
for (i in 1:nyr){ # for each year
  for (j in 1:nmo){ # for each month
    for (k in 1:ndays[j]){ # for days in the month
      for (l in 1:4){ #each day has four 6 hour timesteps, so loop over them, too
        # if the month is june july august, sum up temperatures
        if ( (j == 6) | (j == 7)| (j == 8)){
          sum = sum + out_temp$Temp[m] - 273.15
        }
        m = m + 1
      }

    }
  }
  jja = append( jja,(sum/ (4*(ndays[6]+ndays[7]+ndays[8]) ) ) )
  sum = 0
}


#extract years of relevance for plotting jja
jja_obs = vector()

l_obs = nyr
ib = 1901
ie = 2004
j = ib
for (i in 1:l_obs){
  if (j <= ie){
    jja_obs= append(jja_obs,jja[i])
  }
  j = j + 1
}

### load up Torneträsk width and density observation data
#Source:
#Grudd, H. Torneträsk tree-ring width and density ad 500–2004: a test of climatic
#sensitivity and a new 1500-year reconstruction of north Fennoscandian summers.
# Clim. Dyn. 31, 843–857, DOI:10.1007/s00382-007-0358-2 (2008).
#  ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/reconstructions/
#  europe/sweden/tornetrask-temperature2008.txt

out2 = read.table(paste0(obs_path,"torn_tree_data.txt"))
names(out2) = c("Year", "MXDrcs", "MXDlo", "MXDlowh", "MXDlowt",
                "TRWrcs", "TRWlow", "TRWlowb", "TRWlowt")
out2 <- out2[which(out2$Year >= 1901),]

obs = list()
pre =  list()

l_obs = length(out2$Year)
ib = 1901
ie = 2004
j = ib
for (i in  1:l_obs){
  if (out2$Year[i]==j & j <= ie){
    # rcs = regional curve standardisation  ( see Grudd 2008 above)
    #MXD = maximum density
    obs = append(obs, out2$MXDrcs[i] * 0.615)
  }

  j = j + 1
}

# add model output
out = read.table(paste0(file_path,"TRW.dat"))
names(out) = c("kyr", "trw", "mxd")

pre <- list()
l_pre = length(out$kyr)
ib = 1901
ie = 2004
j = ib
for (i in 1:l_pre){
  if (out$kyr[i] == j & j <= ie){

    pre = append(pre,out$mxd[i])
    j = j + 1
  }
}

plot(jja_obs,obs,col=makeTransparent('dark green',alpha=200),type='p',pch=16,
     ylim=c(0.4,0.8),ylab="",xlab="",yaxt="n",xaxt="n",cex=0.8)
lines(type='p',jja_obs,pre, col='grey',pch=16,cex=0.8)

mtext(side=1, 'Mean June-July-August temperature (degree C)',line=0.7,
      cex = cex.lab)
mtext(side=2, expression(paste("Maximum density (g cm"^"-3",")")),
      line=1.1,cex = cex.lab)
axis(side = 2,lwd = 0, line = -0.7, las = 2)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, tck = -.015, labels = NA)
axis(side = 1,lwd = 0, line = -1.05)



##add regression lines
#prepare for lm function:
data = data.frame(jja_obs=as.vector(unlist(jja_obs)),obs=as.vector(unlist(obs))
                  ,pre=as.vector(unlist(pre)))
abline(lm(obs~jja_obs,data=data),col=makeTransparent('dark green',alpha=180))
abline(lm(pre~jja_obs,data=data),col='grey')

add_legend(-0.81,1.05, legend=c("Observed","Simulated"),
           bty = "n", col = c(makeTransparent('dark green',alpha=200),
                              "grey"),
           lty=c(1,1),cex=cex.legend)



dev.off()

######TRW-plot

# prepare observations:

obs = list()
pre =  list()

l_obs = length(out2$Year)
ib = 1901
ie = 2004
j = ib
for (i in  1:l_obs){
  if (out2$Year[i]==j & j <= ie){
    # rcs = regional curve standardisation ( see Grudd 2008 above)
    # TRW = tree ring width
    obs = append(obs, out2$TRWrcs[i] )
  }

  j = j + 1
}


l_pre = length(out$kyr)
ib = 1901
ie = 2004
j = ib
for (i in 1:l_pre){
  if (out$kyr[i] == j & j <= ie){
    pre = append(pre,(out$trw[i]/1000))
    j = j + 1
  }
}

############FIGURE RING WIDTH sim against obs ####################

jpeg(filename=paste0(figs_path,'TRW.jpeg'),units="mm",width=90, height=90, res=res)
par(mfrow=c(1,1),oma=c(1,1,0.1,0.1),mar=c(1,1,0,0),cex.axis=0.7, cex.lab=0.7,ps=8)

plot(jja_obs,obs,col=makeTransparent('dark green',alpha=200),type='p',pch=16,
     ylim=c(0.5,2.1),ylab="",xlab="",yaxt="n",xaxt="n",cex=0.8)
lines(type='p',jja_obs,pre, col='grey',pch=16,cex=0.8)

mtext(side=1, 'Mean June-July-August temperature (degree C)',line=0.7,
      cex = cex.lab)
mtext(side=2, "Ring width (mm)",line=1.1,cex = cex.lab)
axis(side = 2,lwd = 0, line = -0.7, las = 2)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, tck = -.015, labels = NA)
axis(side = 1,lwd = 0, line = -1.05)



##add regression lines
#prepare for lm function:
data = data.frame(jja_obs=as.vector(unlist(jja_obs)),obs=as.vector(unlist(obs)),
                  pre=as.vector(unlist(pre)))
abline(lm(obs~jja_obs,data=data),col=makeTransparent('dark green',alpha=180))
abline(lm(pre~jja_obs,data=data),col='grey')

add_legend(-0.81,1.05, legend=c("Observed","Simulated"),
           bty = "n", col = c(makeTransparent('dark green',alpha=200),
                              "grey"),
           lty=c(1,1),cex=cex.legend)



dev.off()





