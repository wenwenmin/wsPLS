library(xtable)
# print(xtable(avg_table_data1, digits = 3))

newTable = function(Res_data1){
  avg = Res_data1$table_avg
  rownames(avg) = paste(rownames(avg),".avg")
  sd  = Res_data1$table_sd
  rownames(sd) = paste(rownames(sd),".sd")
  table = round(rbind(avg,sd),3)
  table2 = table[c(1,6,2,7,3,8,4,9,5,10),]
  return(table)
}

print(xtable(newTable(Res_data3), digits = 3))