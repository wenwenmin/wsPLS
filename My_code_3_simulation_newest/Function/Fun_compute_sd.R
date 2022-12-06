Fun_comput_sd = function(Res_data1){
  table_list = Res_data1$table_list 
  n_rep = length(table_list)
  
  mean_mat = Res_data1$table_avg
  sd_mat   = 0*mean_mat
  
  n_row = dim(mean_mat)[1]
  n_col = dim(mean_mat)[2]

  for(i in 1:n_row){
    for(j in 1:n_col){
      for(k in 1: n_rep){
        a = table_list[[k]][i,j]
        a_avg = mean_mat[i,j]
        sd_mat[i,j] = sd_mat[i,j] + (a-a_avg)^2
      }
    }
  }
  
  sd_mat2 = sd_mat
  sd_mat2 = sd_mat2/n_rep
  sd_mat2 = sqrt(sd_mat2)
  return(sd_mat2)
}