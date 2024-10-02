library("xml2")
library("dplyr")
#library("tidyverse")
library("ggplot2")

substituteColIfMissing=function(data,cnam){
  if(length(cnam)==0){
    return (data)
  }
  if(is.vector(cnam) && length(cnam)>1){
    for (i in 1:length(cnam)) { 
      data=substituteColIfMissing(data,as.character(cnam[i]))
    }
    return (data)
  }
  if(!(cnam %in% colnames(data)))
  {
    data <- data %>%
      mutate(!!cnam :=NA)
  }
  return(data) 
}
addColbySubscript=function(row,content){
  if(grepl("_", content, fixed=TRUE)){
    cnam=unlist(strsplit(content, "_"))[1]
    size=unlist(strsplit(content, "_"))[2]
    #print(paste(content,cnam,size))
    row <- row %>%
      mutate(!!cnam :=size)
    #print(paste(cnam,size))
  }
  return(row)
}
extractDataFromName=function(row,file){
  parts <- unlist(strsplit(basename(file), "\\."))
  model=parts[1]
  row <- row %>%
    mutate("model":= model)
  scenario=parts[2]
  if(model=="zeroconf_dl"){
    scenario=parts[3]
    row=addColbySubscript(row,parts[4])
  }
  if(model=="philosophers"||model=="pnueli-zuck"){
    scenario=""
    row$n=parts[2]
  }
  if(model=="mer"){
    scenario=""
    row=addColbySubscript(row,parts[2])
  }
  if(startsWith(model,"firewire")){
    scenario=unlist(strsplit(parts[2], "_"))[1]
    row=addColbySubscript(row,parts[2])
  }
  if(model=="csma"){
    row=addColbySubscript(row,parts[3])
    row=addColbySubscript(row,parts[4])
  }
  if(model=="consensus"){
    row=addColbySubscript(row,parts[3])
  }
  
  software=parts[length(parts)-1]
  row <- row %>%
    mutate("software":= software)
  row <- row %>%
    mutate("scenario":= scenario)
  return(row)
}

#Converts the xml data to a dataframe containing the run information
XML2Data=function(data){
  # Extract data from the "result" element
  result_nodes <- xml_find_all(data, "//result/run")
  dataFrame=data.frame()
  # Loop through the "run" nodes to extract data
  for (node in result_nodes) {
    # Extract data from subnodes named "column title"
    title_nodes <- xml_find_all(node, ".//column[@title]")
    row=data.frame(c(1))
    cpu=1 # CPU column numbers differe. we set them to 1 and 2
    #Convert "column" nodes to columns in a dataframe
    for (entry in title_nodes) {
      key=xml_attr(entry, "title")
      if(startsWith(key,"cputime-cpu")){
        key=paste("cputime-cpu",cpu,sep="")
        cpu=cpu+1
      }
      value=xml_attr(entry, "value")
      
      if(key=="probability"){
        value=as.numeric(value)#ifelse(is.na(as.numeric(value)),-1,as.numeric(value))
      }
      row <- row %>%
        mutate(!!key := value)
    }
    #if terminated this column does not exist in xml and has to be supplied
    if(!("terminationreason" %in% colnames(row)))
    {
      row <- row %>%
        mutate("terminationreason":= NA)
      row <- row %>%
        mutate("exitsignal":= NA)
    }
    #if not terminated this column does not exist in xml and has to be supplied
    if(!("returnvalue" %in% colnames(row)))
    {
      row <- row %>%
        mutate("returnvalue":= NA)
    }
    if(!("cputime-cpu2" %in% colnames(row)))
    {
      row <- row %>%
        mutate("cputime-cpu2":= NA)
    }
    if(!("distribution" %in% colnames(row)))
    {
      row <- row %>%
        mutate("distribution":= NA)
    }
    
    # Extract data from the "name" attribute of the "run" node
    attribute_name <- xml_attr(node, "name")
    row <- row %>%
      mutate("file" := attribute_name)
    #Extract model size from Filename
    parts <- unlist(strsplit(attribute_name, "\\."))
    
    #model=parts[length(parts)-3]
    #row <- row %>%
    #  mutate("model" := model)
    size=parts[length(parts)-2]
    if(grepl("_", size, fixed=TRUE)){
      size=unlist(strsplit(size, "_"))[2]
    }
    
    if(is.na(as.numeric(size))){
      print(paste(size,attribute_name,parts))
    }
    size=ifelse(is.na(as.numeric(size)),-1,as.numeric(size))
    row <- row %>%
      mutate("size" := size)
    row=extractDataFromName(row,attribute_name)
    #print("#################################")
    #print(colnames(row))
    #print("BREAK")
    #print(colnames(dataFrame))
    #print("#################################")
    
    row=substituteColIfMissing(row,setdiff(colnames(dataFrame),colnames(row)))
    dataFrame=substituteColIfMissing(dataFrame,setdiff(colnames(row),colnames(dataFrame)))
    dataFrame=rbind(dataFrame,row)
  }
  dataFrame$c.1.=NULL
  return(dataFrame)
}

#read all xml.bz2 in the working directory and add them to a common data frame
file_list <- list.files("../results",pattern = "\\.xml\\.bz2$")
data=data.frame()
for(file in file_list){
  path=paste("../results/",file,sep="")
  #CAVE currently only works for storm. Modest comming soon
  parts <- unlist(strsplit(basename(file), "[\\.]"))
  engine=parts[1]
  engine_scenario=parts[length(parts)-1]
  result=XML2Data(read_xml(path))
  result$engine=engine
  result$engine_scenario=engine_scenario
  if(nrow(data)==0){
    data=result
  }else{
    result=substituteColIfMissing(result,setdiff(colnames(data),colnames(result)))
    data=substituteColIfMissing(data,setdiff(colnames(result),colnames(data)))
    data=rbind(data,result)
  }
}
data$engine=as.factor(data$engine)
data$engine_scenario=as.factor(data$engine_scenario)
data$fullscenario=as.factor(paste(data$model,data$scenario))

for( fact in levels(data$fullscenario)){
  #plot probability over size
  # Create a line graph
  filtered=data%>%filter(fullscenario==fact)%>%filter(!is.na(probability)&!is.na(size))%>%filter(status=="done")
  
  plot=ggplot(filtered, aes(x = size, y = probability, color = engine, group = engine), fill="white") +
    geom_line() +
    labs(x = "Size", y = "Probability") +
    scale_y_continuous(limits = c(0, 1)) +
    ggtitle(fact)+
    theme_classic()#+ theme(plot.background = element_rect(fill = "white"))
  ggsave(paste(fact,".png",sep=""))
}

