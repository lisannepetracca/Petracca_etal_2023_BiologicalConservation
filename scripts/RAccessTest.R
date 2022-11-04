# Abby Bratt
# aebratt@uw.edu
# last updated September 2020
# Example code for initiating connected between R and MS access databases

# install.packages(“RODBC”)
library(RODBC)
# note that this requires Rcpp
# also had to update R

# steps to take before this (if file is accdb)
# installing the driver described in step 3 here:
  # https://www.r-bloggers.com/getting-access-data-into-r/?
# Ursus has the 2019 version of access but the 2016 driver worked for me. 
  # Downloaded here: https://www.microsoft.com/en-us/download/details.aspx?id=54920
# additional resources
  # https://db.rstudio.com/odbc/

# open connection
con <- odbcConnect("SHLAabund") # name of database in working directory
con # check it out
sqlTables(con) # list all tables
sqlTables(con, tableType = "TABLE")$TABLE_NAME # restrict to relevant tables
# load data
SHLA_detail <- sqlFetch(con, "Tbl:SHLA_detail")
SHLA_Event <- sqlFetch(con, "Tbl:SHLA_Event")
site <- sqlFetch(con, "tlu:site")

# if you wanted to process data do it here
# blahblahblah

# if you wanted to write data back to con do it here
test <- site
sqlSave(con, dat = test)

# it worked!
bar <- sqlFetch(con, "test")

# close connections
odbcCloseAll()
