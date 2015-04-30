#  libgc.R
#  
#  Copyright 2014 Christian Diener <christian@giiku.com>
#  
#  MIT license. See LICENSE for more information.

OD.min = 0.3
OD.max = 1.0

read.gc = function(name, time_scale=24)
{
    ftype = tolower( strsplit(name, '\\.')[[1]][2] )
    print(ftype)
    if( ftype=="csv" ) data = read.csv(name, header=T)
    else if( ftype=="txt" ) data = read.table(name, header=1)
    else stop("Not a recognized file type. Must be txt or csv!")
	
	data$Time = data$Time*time_scale
	return(data)
}

subsample = function(gcs, n_points)
{
    idx = seq(1:nrow(gcs), length.out=n_points)
    
    return( gcs[idx,] )
}

get_factors = function(name)
{
	name = gsub("\\.\\d+", "", name)
	name = gsub("\\.", "-", name)
	splits = strsplit(name, "_")[[1]]
	
	if (length(splits)<2) splits = c("WT", "C", NA)
	else if (length(splits)<3) splits = c( splits, NA )
	else splits = c(splits[1], splits[2], 
					paste(splits[3:length(splits)], collapse=' '))
	
	return( splits )
}

names_conv = function(name_list)
{
	splits = sapply(name_list, get_factors)
	df_info = data.frame(strain=splits[1,], treatment=splits[2,], info=splits[3,])
	
	# Resort by appearance
	for( i in 1:ncol(df_info) ) 
		df_info[,i] = factor( df_info[,i], levels=unique(df_info[,i]) )
	
	return( df_info )
}

normalize = function(data)
{
	ODs = as.matrix(data[,-1])
	ODs = apply(ODs, 2, function(col) col-min(col)) 
	data[,-1] = ODs
	return(data)
}

join = function(data_set1, data_set2, ref_ids=2:3) 
{
	gcs1 = data_set1
	gcs2 = data_set2
  
	control1 = rowMeans(gcs1[,ref_ids])
	control2 = rowMeans(gcs2[,ref_ids])
	t1 = sum(control1<max(control1)/2)
	print(t1)
	t2 = sum(control2<max(control2)/2)
	shift = t1-t2
	
	if(shift<1) {
		data_set1 = data_set1[1:(nrow(data_set1)+shift),]
		data_set2 = data_set2[shift:-1,]
		return(cbind(data_set1,data_set2[,-1]))
	}
	else if(shift>1) {
		data_set1 = data_set1[-shift:-1,]
		data_set2 = data_set2[1:(nrow(data_set2)-shift),]
		return(cbind(data_set2,data_set1[,-1]))
	}
	else return(cbind(data_set1,data_set2[,-1]))
}


lm_stats = function(OD, Time)
{
	# Calculate the area under the curve
	fun = approxfun(Time, OD, method="linear")
	ABC = integrate(fun, min(Time), max(Time), rel.tol=1e-3)$value
	
	good = which(OD>OD.min & OD<OD.max)
	if(length(good)>4)
	{
		logOD = log(OD[good])
		Time = Time[good]
		model = lm(logOD ~ Time)
		stats = c( coef(model), anova(model)[1,5], length(good), ABC)
		
		# Negative growth rates are invalid in our interpretation
		if( stats[2]<0 ) stats[2]=0
	}
	else stats = c(0, 0, 1, length(good), ABC)
	
	return(stats)
}

lm.stats.multi = function(data)
{
	Time = data$Time
	if(max(Time)<3) Time = Time*24
	probes = names(data[,-1])
	stats = t( apply(as.matrix(data[,-1]), 2, function(od) lm_stats(od, Time)) )
	stats = as.data.frame(stats)
	names(stats) = c("intercept", "growth.rate", "pval", "n_reg", "area")
	
	info = names_conv(probes)
	stats = cbind(stats, info)
	rownames(stats) = NULL
	
	return(stats)
}
