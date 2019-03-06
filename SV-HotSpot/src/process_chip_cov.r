#!/gapp/x64linux/opt/R3.1.2/bin/Rscript

# Find regions/peaks whose SVs altered expression of nearby genes
# Written by Abdallah Eteleeb & Ha Dang

require(data.table)

args = commandArgs(T)

chip.file = args[1]
window.size = as.numeric(args[2])
out.dir = args[3]

#################### function to compute the average chip over window
avgChipOverWindows <- function(ch, chr, w=NULL){
  s = window.size
  imin = min(ch$start)
  imax = max(ch$end)
  win = data.frame(chrom=chr, start=seq(imin,imax,s))
  win$stop = win$start + s - 1
  win$pos = (win$start + win$stop)/2
  win$mean.cov = 0
  print(paste0('Averaging chip data coverage over ', s,'-window for chromosome ',chr))
  for (i in 1:nrow(win)){
    #cat ('.')
    chi = ch[ch$pos < win$stop[i] & ch$pos >= win$start[i],]
    win$mean.cov[i] = mean(chi$cov)
  }
  cat('\n')
  return(win)
}

### read chip-seq data 
if (file.exists((chip.file))) {
    cat ('Reading chip-seq coverage file ...\n')
    chip.data <- fread(chip.file)
    chip.data$pos = (chip.data$start+chip.data$end)/2
    chrs <- unique(chip.data$chrom)
    avg.chip.cov = NULL

    #Initiate the bar
    pb <- txtProgressBar(min = 0, max = nrow(res), style = 3)
    for (c in chrs) {
        ch.chr.data = chip.data[chip.data$chrom==c, ]
        avg.chip.cov = rbind(avg.chip.cov, avgChipOverWindows(ch.chr.data, c))
        #Update the progress bar
        setTxtProgressBar(pb, i)
    }
    
} else {
     stop('chip-seq file was not found!')
}

### write resutls 
avg.chip.cov = avg.chip.cov[, c('chrom','start','stop', 'mean.cov')]
colnames(avg.chip.cov) = c('chrom','start','end', 'cov')
write.table(avg.chip.cov, file=paste0(out.dir,'/processed_chip_coverage.tsv'), sep="\t", row.names=F, quote = F)
