#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(dplyr)

f <- read.delim(file=args[1], header=TRUE)
print(nrow(f))
t <- f %>% filter(!seqID %in% "unclassified")
t <- t %>% mutate(V9 = hitLength/queryLength) %>% arrange(desc(V9)) %>% filter(V9 > 0.50) %>% select(readID)
a <- anti_join(f,t) %>% select(readID) %>% distinct()
print(nrow(a))
write.table(a, file="not_contam.txt",col.names=FALSE ,row.names=FALSE,sep="\t", quote = FALSE)


