
for (nam in names(cors.list)) {        
  write(nam, file = outfile, append = TRUE)        
  write.table(round(cors.list[[nam]],3), file = outfile, append = TRUE) 
  write("", file = outfile, append = TRUE)        
}          
     





