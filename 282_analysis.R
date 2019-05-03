library(tidyverse)

# load in data to bind into a data frame 
train_y <- read_delim('~/Downloads/TCGA_mRNA_RNA/y_train.csv', ',', col_names = F)
test_y <- read_delim('~/Downloads/TCGA_mRNA_RNA/y_test.csv', ',', col_names = F)

sample_labels_train <- read_delim('~/Downloads/TCGA_mRNA_RNA/train_X_samples.csv', ',', col_names=F)
sample_labels_test <- read_delim('~/Downloads/TCGA_mRNA_RNA/test_X_samples.csv', ',', col_names=F)

train_X <- read_delim('~/Downloads/TCGA_mRNA_RNA/X_train.csv', ',', col_names = F)
test_X <- read_delim('~/Downloads/TCGA_mRNA_RNA/X_test.csv', ',', col_names = F)

X <- bind_rows(train_X, test_X)

train_X_encoded <- read_delim('~/Downloads/TCGA_mRNA_RNA/train_X_encoded.csv', ',', col_names=F)
test_X_encoded <- read_delim('~/Downloads/TCGA_mRNA_RNA/test_X_encoded.csv', ',', col_names=F)

train_X_encoded_umap <- read_delim('~/Downloads/TCGA_mRNA_RNA/embedding.csv', ',', col_names=F)
train_X_raw_umap <- read_delim('~/Downloads/TCGA_mRNA_RNA/raw_embedding.csv', ',', col_names=F)

# compute cancer labels 
train_cancer <- floor(train_y/2)
test_cancer <- floor(test_y/2)

cancers <- c('BRCA', 'COAD', 'KIRC', 'KIRP', 'LGG', 'LUAD', 'READ', 'UCEC')
assays <- c('RNASeq2GeneNorm', 'mRNAArray')

train_interpolated_data <- list()
for(i in seq(from=1, to=15, by=2)){
  for(j in 1:2){
    train_interpolated_data[[i+j-1]] <- read_delim(paste0('~/Downloads/TCGA_mRNA_RNA/train_X_', 
                                                          cancers[ceiling(i/2)], '_', assays[j], 
                                                          '_', 'interpolated.csv', collapse = ''), 
                                                   ',', col_names = F)
  }
}

test_interpolated_data <- vector('list', 16)
for(i in seq(from=1, to=15, by=2)){
  for(j in 1:2){
    test_interpolated_data[[i+j-1]] <- read_delim(paste0('~/Downloads/TCGA_mRNA_RNA/test_X_', 
                                                          cancers[ceiling(i/2)], '_', assays[j], 
                                                          '_', 'interpolated.csv', collapse = ''), 
                                                   ',', col_names = F)
  }
}

#bind these lists and delete them to free up memory 

train_interpolated_df <- bind_rows(train_interpolated_data)
rm(train_interpolated_data)

test_interpolated_df <- bind_rows(test_interpolated_data)
rm(test_interpolated_data)

interpolated_df <- bind_rows(train_interpolated_df, test_interpolated_df)

bind_dataframe <- function(data, sample_labels, y_labels){
  #check that dimensions match up since R will silently recycle 
  stopifnot(dim(data)[0]==dim(sample_labels)[0])
  stopifnot(dim(y_labels)[0]==dim(sample_labels)[0])
  
  participant <- substr(sample_labels$X1, 9, 12)
  # cancer samples are 01-10 normal samples are 11-19
  sample_type <- ifelse(substr(sample_labels$X1, 14, 15)<10, 1, 0)
  
  df <- bind_cols(sample=sample_labels, 
                  participant = participant, 
                  sample_type = sample_type, 
                    y=y_labels, 
                    cancer=floor(pull(y_labels)/2), 
                    data)
  colnames(df)[1:5] <- c('SampleID', 'Participant', 'CancerStatus', 'Dataset', 'Cancer')
  return(df)
}


# compute metadata for all samples 

dataset <- tibble(c(rep(1, 3943), rep(0,982)))
names(dataset) <- c('TrainSplit')

metadata <- bind_cols(bind_rows(sample_labels_train, sample_labels_test), bind_rows(train_y,test_y), dataset)

metadata <- bind_cols(metadata,Participant=substr(metadata$X1, 9, 12), CancerStatus=ifelse(substr(metadata$X1, 14, 15)<10, 1, 0) )
names(metadata)[1:2] <- c('ID', 'CancerAssay')


# find samples with paired data 

replicates <- metadata %>% filter(CancerStatus==1) %>% select(-c(ID, TrainSplit, CancerStatus)) %>% group_by(Participant) %>% 
         summarise(nassays=n(), replicates=nassays>1) %>% filter(replicates) %>% select(Participant)

replicate_metadata <- metadata[metadata$Participant %in% pull(replicates),]

# find rna-seq mrna paired data 

paired_replicates <- replicate_metadata %>% group_by(Participant) %>% summarise(NumAssays=n_distinct(CancerAssay), PairedData=NumAssays>1) %>%
                      filter(PairedData) %>% select(Participant)

paired_replicate_metadata <- replicate_metadata[replicate_metadata$Participant %in% pull(paired_replicates),]

#remove 3+ pairs for easier analysis 
binary_pairs <- paired_replicate_metadata %>% group_by(Participant) %>% summarise(n=n(), junk=n>2) %>% filter(!junk) %>% select(Participant)

paired_replicate_metadata <- paired_replicate_metadata[paired_replicate_metadata$Participant %in% pull(binary_pairs),]

#sort the paired_replicate_metadata 
sorted_paired_replicate_metadata <- paired_replicate_metadata %>% arrange(desc(Participant))

# calculate distance between interpolated samples and GT samples 
# as percentage of the GT norms 

#holds the cancer assay, trainsplit, and resulting distance
interpolate_tibble <- vector('list', dim(sorted_paired_replicate_metadata)[1])

calc_interpolate_distance <- function(input_sample, paired_sample){
  #find index of input and paired sample 
  input_index <- which(metadata$ID==input_sample$ID & metadata$CancerAssay==input_sample$CancerAssay &
                        metadata$TrainSplit==input_sample$TrainSplit & metadata$CancerStatus==input_sample$CancerStatus)
  paired_index <- which(metadata$ID==paired_sample$ID & metadata$CancerAssay==paired_sample$CancerAssay &
                            metadata$TrainSplit==paired_sample$TrainSplit & metadata$CancerStatus==paired_sample$CancerStatus)
  #input_sample has been interpolated into paired samples domain 
  input_vector <- interpolated_df[input_index,]
  # get the paired GT from original data 
  paired_vector <- X[paired_index,]
  
  relative_norm <- rel_norm(input_vector,paired_vector)
  
  return(tibble(Assay=input_sample$CancerAssay, Split=input_sample$TrainSplit, RelNorm=relative_norm))
}

# correlation study 

calc_correlation <- function(input_sample, paired_sample){
  #find index of input and paired sample 
  input_index <- which(metadata$ID==input_sample$ID & metadata$CancerAssay==input_sample$CancerAssay &
                         metadata$TrainSplit==input_sample$TrainSplit & metadata$CancerStatus==input_sample$CancerStatus)
  paired_index <- which(metadata$ID==paired_sample$ID & metadata$CancerAssay==paired_sample$CancerAssay &
                          metadata$TrainSplit==paired_sample$TrainSplit & metadata$CancerStatus==paired_sample$CancerStatus)
  #input_sample has been interpolated into paired samples domain 
  input_vector <- unlist(X[input_index,])
  # get the paired GT from original data 
  paired_vector <- unlist(X[paired_index,])
  
  correlation <- cor(input_vector,paired_vector, method='spearman')
  
  return(correlation)
}

cor_list <- c()
for(i in seq(from=1, to=(dim(sorted_paired_replicate_metadata)[1]), by=2)){
  print(i)
  cor_list <- c(cor_list,calc_correlation(sorted_paired_replicate_metadata[i,],
                                                       sorted_paired_replicate_metadata[ifelse(i%%2==0, i-1, i+1),]))
}

cor_df <- data.frame(Correlation=cor_list, Cancer=sorted_paired_replicate_metadata$CancerAssay[seq(1, 1696, 2)])
correlation_histogram <- ggplot(cor_df, aes(x=Correlation, Color=factor(Cancer)))+
                        geom_histogram(binwidth = .01)+
                        ggtitle("Spearman Correlation Between Paired Microarray and RNA-Seq Samples")+
                        ylab('Count')+theme_minimal()

rel_norm <- function(x, y){
  sum((x-y)**2)/sum((y)**2)
}

for(i in 242:(dim(sorted_paired_replicate_metadata)[1])){
  print(i)
  interpolate_tibble[[i]] <- calc_interpolate_distance(sorted_paired_replicate_metadata[i,],
                                                       sorted_paired_replicate_metadata[ifelse(i%%2==0, i-1, i+1),])
}

interpolate_tibble_df <- bind_rows(interpolate_tibble)

# plot the latent space

umap_df <- bind_dataframe(train_X_encoded_umap, sample_labels_train, train_y)

long_cancers <- c('Breast invasive carcinoma', 'Colon adenocarcinoma',
                  'Kidney renal clear cell carcinoma',
                  'Kidney renal papillary cell carcinoma',
                  'Brain Lower Grade Glioma', 
                  'Lung adenocarcinoma',
                  'Rectum adenocarcinoma',
                  'Uterine Corpus Endometrial Carcinoma')

latent_space_plot <- ggplot(umap_df, aes(x=X12, y=X2, color=factor(Cancer), shape=factor(Dataset%%2)))+
  geom_point(alpha=.3)+
  ggtitle('UMAP of Latent Space Encoding')+
  scale_color_discrete(name='Cancer Type', labels=cancers)+
  scale_shape_discrete(name='Assay', labels=c('RNA-seq', 'mRNA Array'))+
  xlab('UMAP 1')+
  ylab('UMAP 2')+
  theme_light()

# plot just the mRNA arrays in the latent space 
mrna_umap_df <- umap_df %>% filter(Dataset%%2==1)
mRNA_latent_space_plot <- ggplot(mrna_umap_df, aes(x=X12, y=X2, color=factor(Cancer), shape=factor(2)))+
  geom_point(alpha=.9)+
  ggtitle('UMAP of mRNA Arrays in Latent Space')+
  scale_color_discrete(name='Cancer Type', labels=cancers)+
  scale_shape_discrete(name='Assay', labels=c('mRNA Array'))+
  xlab('UMAP 1')+
  ylab('UMAP 2')+
  theme_light()

Rel# plot the raw data 
raw_umap_df <- bind_dataframe(train_X_raw_umap, sample_labels_train, train_y)

raw_latent_space_plot <- ggplot(raw_umap_df, aes(x=X12, y=X2, color=factor(Cancer), shape=factor(Dataset%%2)))+
  geom_point(alpha=.3)+
  ggtitle('UMAP of TCGA Data')+
  scale_color_discrete(name='Cancer Type', labels=cancers)+
  scale_shape_discrete(name='Assay', labels=c('RNA-seq', 'mRNA Array'))+
  xlab('UMAP 1')+
  ylab('UMAP 2')+
  theme_light()

# plot interpolation ability 

interpolate_into_mRNA_plot <- ggplot(interpolate_tibble_df %>% filter(Assay%%2==0), aes(factor(Assay), sqrt(RelNorm), fill=factor(Assay)))+
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2)+
  ggtitle('Relative Norms of Interpolations from RNA-seq to mRNA Array')+
  scale_fill_discrete(name='Cancer Type', labels=cancers)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab('Relative Error in Norm of Interpolation')+
  theme_light()

interpolate_into_Rnaseq_plot <- ggplot(interpolate_tibble_df %>% filter(Assay%%2==1), aes(factor(Assay), sqrt(RelNorm), fill=factor(Assay)))+
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2)+
  ggtitle('Relative Norms of Interpolations from mRNA Array to RNA-seq')+
  scale_fill_discrete(name='Cancer Type', labels=cancers)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab('Relative Error in Norm of Interpolation')+
  theme_light()






