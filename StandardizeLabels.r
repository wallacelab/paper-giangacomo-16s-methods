# Helper function to standardize sample and treatment labels (so don't have to recode in each script) 

library(phyloseq)

# Takes a phyloseq object and standardizes sample.type and treatment columns of the metadata (sample_data)
standardize_labels = function(mydata, type=c('extraction','amplification')){
    metadata = sample_data(mydata)
    
    if(type=='extraction'){
        metadata$sample.type = sapply(as.character(metadata$sample.type), switch,
                                      "Leaf-Arabidopsis"="Arabidopsis Leaf", 
                                      "Leaf-Corn"="Maize Leaf", 
                                      "Leaf-Soybean"="Soybean Leaf", 
                                      "Soil-Soil SS1"="Soil",
                                      NA) # NA catches anything that didn't match
        metadata$sample.type = factor(metadata$sample.type, levels=c("Arabidopsis Leaf", "Maize Leaf", "Soybean Leaf", "Soil"))
        
        metadata$treatment = sapply(as.character(metadata$treatment), switch,
                                        ExtracNAmp="ExtractNAmp", 
                                        MoBioPowerSoil="PowerSoil",
                                        QiagenDNeasyPlant="DNeasyPlant",
                                        ZymoEasyDna = "EasyDNA",
                                        KazuBuffer = "KazuBuffer",
                                        NA) # NA catches anything that didn't match
        metadata$treatment = factor(metadata$treatment)
    }
    else if(type=='amplification'){
        metadata$sample.type = sapply(as.character(metadata$sample.type), switch,
                                      "leaf-maize"="Maize Leaf", 
                                      "leaf-soybean"="Soybean Leaf", 
                                      "defined-community"="Defined Community", 
                                      "soil-clay"="Soil 1",
                                      "soil-flowerbed"="Soil 2",
                                      "blank"="blank",
                                      "water"="water",
                                      NA) # NA catches anything that didn't match
        
        metadata$treatment = sapply(as.character(metadata$treatment), switch,
                                    BlockingOligos_v3v4="BO_3/4", 
                                    BlockingOligos_v5v7="BO_5/7",
                                    BlockingOligos_v5v7_noLinkers="BO_5/7",
                                    Discriminating = "Discriminating",
                                    PNAs= "PNA", 
                                    Universal="Universal",
                                    NA) # NA catches anything that didn't match
        
        
    }else{
        stop("Illegal sample type provided to standardize_labels")
    }
    
    sample_data(mydata) = metadata
    return(mydata)
}
