prep_data_indv_pred <- function(data, id, timeVar){
  
  id_unique <- unique(data[, id])
  
  data_out <- data.frame()
  
  for(i in id_unique){
    
    data_i <- data[data[, id] == i, ]
    
    for(j in 1:nrow(data_i)){
      data_i_j <- data_i[1:j, ]
      data_i_j[, id] <- paste0(i, "_", data_i_j[nrow(data_i_j), timeVar]) %>% rep(nrow(data_i_j))
      data_out <- rbind(data_out, data_i_j)
    }
    
  }
  
  return(data_out)
  
}