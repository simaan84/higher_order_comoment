library(lubridate)
library(data.table)
library(plyr)
library(rvest)
library(parallel)
library(quantmod)

rm(list = ls())
gc()


###########################################
########## S&P 500 DATA ###################



#### FIND THE S&P 500 CONST. 
theurl <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies"

file <- read_html(theurl)
tables <- html_nodes(file, "table")

# try extract each into a data.frame
tables.list <- numeric()
for (i in 1:length(tables)) { 
  tables.list[i] <- list(try(html_table(tables[i], fill = TRUE),silent = T))
}

# the table we are looking for should have columns for volume and price: 
find.sym <- lapply(tables.list, function(x) grep("Symbol",x,ignore.case = T)   )
find.sec <- lapply(tables.list, function(x) grep("Security",x,ignore.case = T)   )
locate.table <- which.max(sapply(find.sym, length) + sapply(find.sec, length) )

ds <- tables.list[[locate.table]][[1]]

#### LOAD CRSP DATA ######
file.i <- "/home/simaan/Dropbox/Data/DATA/CRSP_1960_2019_d.csv"

select_var <- c("PERMCO","PERMNO","date","COMNAM","SHROUT","RET","PRC","TICKER")
DS <- fread(file.i,select = select_var)

# keep tickers that show up in the wiki page for being part of S&P
DS <- DS[DS$TICKER %in% ds$Symbol,]
gc()

DS$RET <- as.numeric(DS$RET)
DS$date <- ymd(DS$date)

# let's pivot the data
DS <- unique(DS)
DS_N <- DS[,.N, by = c("date","PERMNO")]
table(DS_N$N)

DS_sub <- unique(DS[,list(PERMNO,date,RET),])
DS_sub <- DS_sub[order(DS_sub$PERMNO,DS_sub$date),]
DS_sub <- na.omit(DS_sub)
DS_sub[,N :=  .N, by = c("date","PERMNO")]
table(DS_sub$N)

# control for duplicates
DS_sub <- DS_sub[,lapply(.SD, mean),by = c("PERMNO","date"),.SDcols = "RET"]
DS_sub[,N :=  .N, by = c("date","PERMNO")]
table(DS_sub$N)
DS_sub$N <- NULL
rm(DS); gc()

DS_list <- dlply(DS_sub,"PERMNO",data.frame)

pivot_data <- function(ds_i) {
  ret_xts <- ds_i$RET
  names(ret_xts) <- ds_i$date 
  ret_xts <- as.xts(ret_xts)
  names(ret_xts) <- paste("PERMNO",unique(ds_i$PERMNO),sep = "_")
  return(ret_xts)
}

xts_list <- mclapply(DS_list,pivot_data,mc.cores = 16)

track_dates <- lapply(xts_list,function(x) range(date(x)) )
keep_last <- sapply(track_dates, function(x) x[2] ==  "2019-12-31")
# in total there are companies that available at the end
sum(keep_last)

# consider those companies available before
keep_first <- sapply(track_dates, function(x) x[1] <=  "1990-01-01")
keep_stocks <- keep_first&keep_last
sum(keep_stocks)

xts_list2 <- xts_list[keep_stocks]
xts_list2 <- lapply(xts_list2, function(x) x[date(x) >= "1990-01-01",]  )
total_obs <- sapply(xts_list2,nrow)
xts_list2 <- xts_list2[total_obs == max(total_obs)]

xts_all <- Reduce(merge.xts,xts_list2)
xts_all <- na.omit(xts_all)




########## ADD FAMA-FRENCH PUBLIC DATA #######
FF_file <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/100_Portfolios_10x10_Daily_CSV.zip"
temp <- tempfile()
download.file(FF_file,temp)
unz_files <- unzip(temp)
ds <- read.csv(unz_files,skip = 18)
flag_obs <- grep("Average Equal Weighted Returns -- Daily",ds[,1],ignore.case = T)
ds <- ds[1:(flag_obs-1),]
names(ds)[1] <- "date"
ds$date <- ymd(ds$date)
ds <- ds[ds$date >= "1990-01-01",]
ds[,-1] <- data.frame(apply(ds[,-1], 2, as.numeric))

FF_data <- ds
rownames(FF_data) <- FF_data$date
FF_data$date <- NULL
FF_data <- as.xts(FF_data)
FF_data <- FF_data/100 # report returns in decimals
FF_data <- FF_data["1990/2019",]
rm(ds); gc()


# save both data sets into a list
all_data_list <- list(xts_all,FF_data)

###############################################################################################################
